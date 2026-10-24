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
#include <cmath>

/************************************
 * PRIVATE MACROS AND DEFINES
 ************************************/
namespace {
// Q and R copied from HALO::initializeHALO (legacy HALO.cpp), produced by covarianceCalc.py
// from sim residuals. REVIEW: Q(0,0) = 74777 m^2 is a ~273 m altitude std, which looks
// far too large; confirm with whoever generated these before flight
constexpr Mat3 LEGACY_Q = {{
	{74777.41f, 4458.13f, -2164.91f},
	{4458.13f,  2413.02f,    -4.52f},
	{-2164.91f,   -4.52f,   503.78f}
}};

constexpr Mat3 LEGACY_R = {{
	{15438.09f, 1528.51f,   727.68f},
	{1528.51f,  5005.79f,  -469.20f},
	{727.68f,   -469.20f,   613.20f}
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

	init.Q = LEGACY_Q;
	init.R = LEGACY_R;

	return init;
}
