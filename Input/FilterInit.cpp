/*
 * FilterInit.cpp
 *
 *  Created on: Oct 17, 2026
 *      Author: Aly Masani
 */

/**
 ********************************************************************************
 * @file    FilterInit.cpp
 * @author  Aly Masani
 * @date    Oct 17, 2026
 * @brief   UKF initial conditions from the tare. Replaces HALO::initializeHALO
 *          without Eigen or the scenario database
 ********************************************************************************
 */

/************************************
 * INCLUDES
 ************************************/
#include "Inc/FilterInit.hpp"
#include "Inc/AxisAlignment.hpp"
#include "Inc/FilterProfile.hpp"
#include "../Estimate/Inc/ScenarioData.hpp"
#include <cmath>

/************************************
 * PRIVATE MACROS AND DEFINES
 ************************************/
namespace {
// Legacy HALO::initializeHALO used (from covarianceCalc.py sim residuals):
//   Q = {{74777.41, 4458.13, -2164.91}, {4458.13, 2413.02, -4.52}, {-2164.91, -4.52, 503.78}}
//   R = {{15438.09, 1528.51, 727.68}, {1528.51, 5005.79, -469.20}, {727.68, -469.20, 613.20}}
// Q(0,0) = 74777 m^2 per step is a ~273 m altitude std, so the legacy filter mostly ignored
// its model. Replaced by the diagonal starting values below; tune against replayed flights.

// Process noise as spectral densities (per second), the UKF multiplies by dt every step.
// Acceleration gets the most: the scenario model is least sure at burnout and drogue / main deploy
constexpr Mat3 DEFAULT_Q = {{
	{0.1f, 0.0f, 0.0f},		// m^2/s
	{0.0f, 1.0f, 0.0f},		// (m/s)^2/s
	{0.0f, 0.0f, 100.0f}	// (m/s^2)^2/s
}};

// Measurement noise floors for the Everest measurement (alt, vel, acc). The UKF raises the
// altitude term to the live dynamite sigma when that is larger
// Tuned with the PC replay (Tools/replay): NIS / NEES near 3 on the Stormfly sims
constexpr Mat3 DEFAULT_R = {{
	{0.05f, 0.0f, 0.0f},	// (0.22 m)^2
	{0.0f, 0.25f, 0.0f},	// (0.5 m/s)^2, the dynamite velocity is smoothed, not measured
	{0.0f, 0.0f, 0.5f}		// (0.7 m/s^2)^2, IMU vibration in flight
}};

/**
 * @brief Variance of the average of up to two sensors with stds std0 / std1
 */
float CombinedVariance(bool valid0, float std0, bool valid1, float std1)
{
	if (valid0 && valid1) {
		return (std0 * std0 + std1 * std1) / 4.0f;
	}
	if (valid0) {
		return std0 * std0;
	}
	if (valid1) {
		return std1 * std1;
	}
	return 0.0f;
}

float Floor(float value, float minimum)
{
	return (value > minimum) ? value : minimum;
}
}

/************************************
 * FUNCTION DEFINITIONS
 ************************************/
FilterInitState FilterInit::Build(const TareResult& tare)
{
	FilterInitState init{};

	// Pad altitude, averaged over the valid baros
	uint8_t baroCount = 0;
	for (uint8_t i = 0; i < 2; i++) {
		if (tare.baroValid[i]) {
			init.groundAltitude += tare.groundAltitude[i];
			baroCount++;
		}
	}
	init.valid = (baroCount > 0);
	if (baroCount > 0) {
		init.groundAltitude /= baroCount;
	}

	// Filter runs relative to the pad, at rest
	init.x0[0] = 0.0f;
	init.x0[1] = 0.0f;
	init.x0[2] = 0.0f;

	// Altitude noise: average of the two baros
	const float altVar = CombinedVariance(tare.baroValid[0], tare.altitudeStd[0],
										  tare.baroValid[1], tare.altitudeStd[1]);

	// Acceleration noise along body +Z (signed permutation, so |rotated std| is exact)
	float accStdBody[2][3];
	for (uint8_t i = 0; i < 2; i++) {
		AxisAlignment::ApplyToImu(i, tare.accelStd[i], accStdBody[i]);
	}
	const float accVar = CombinedVariance(tare.imuValid[0], std::fabs(accStdBody[0][2]),
										  tare.imuValid[1], std::fabs(accStdBody[1][2]));

	init.P0.m[0][0] = Floor(altVar, P0_MIN_ALT_VAR);
	init.P0.m[1][1] = P0_VEL_VAR;
	init.P0.m[2][2] = Floor(accVar, P0_MIN_ACC_VAR);

#if APF_COVARIANCE == APF_COVARIANCE_TABER	// Input/Inc/FilterProfile.hpp
	// Taber profiles: the legacy HALO::initializeHALO matrices from covarianceCalc.py.
	// Legacy added Q once per 1/3 s step; this UKF adds Q * dt, so Q is scaled by the legacy
	// rate (3 /s) to add the same amount per step
	init.Q = {{
		{74777.41f * 3.0f, 4458.13f * 3.0f, -2164.91f * 3.0f},
		{4458.13f * 3.0f, 2413.02f * 3.0f, -4.52f * 3.0f},
		{-2164.91f * 3.0f, -4.52f * 3.0f, 503.78f * 3.0f}
	}};
	init.R = {{
		{15438.09f, 1528.51f, 727.68f},
		{1528.51f, 5005.79f, -469.20f},
		{727.68f, -469.20f, 613.20f}
	}};
	init.Qafter = init.Q;
#elif APF_COVARIANCE == APF_COVARIANCE_SCENARIO
	// Q before / after apogee from how far the closest cases sit apart (Tools/preprocess_cases.py),
	// never below DEFAULT_Q on the diagonal
	for (uint8_t i = 0; i < 3; i++) {
		for (uint8_t j = 0; j < 3; j++) {
			init.Q.m[i][j] = SCENARIO_Q_BEFORE[i][j];
			init.Qafter.m[i][j] = SCENARIO_Q_AFTER[i][j];
		}
		init.Q.m[i][i] = Floor(init.Q.m[i][i], DEFAULT_Q.m[i][i]);
		init.Qafter.m[i][i] = Floor(init.Qafter.m[i][i], DEFAULT_Q.m[i][i]);
	}
	init.R = DEFAULT_R;
#else
	init.Q = DEFAULT_Q;
	init.Qafter = DEFAULT_Q;
	init.R = DEFAULT_R;
#endif
	for (uint8_t i = 0; i < 3; i++) {
		for (uint8_t j = 0; j < 3; j++) {
			init.Q.m[i][j] *= APF_Q_SCALE;	// Input/Inc/FilterProfile.hpp, 1 unless sweeping
			init.Qafter.m[i][j] *= APF_Q_SCALE;
			init.R.m[i][j] *= APF_R_SCALE;
		}
	}

	return init;
}
