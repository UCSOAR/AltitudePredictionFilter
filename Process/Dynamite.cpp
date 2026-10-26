/*
 * Dynamite.cpp
 *
 *  Created on: Oct 24, 2026
 *      Author: Aly Masani
 */

/**
 ********************************************************************************
 * @file    Dynamite.cpp
 * @author  Aly Masani
 * @date    Oct 24, 2026
 * @brief   PROCESS stage, modelled on the legacy EverestTask dynamite() / recalculateGain():
 *          every source gets a trust weight from how close it is to the predicted
 *          altitude, the weighted blend is the Everest measurement for the UKF and
 *          Launch Detection
 ********************************************************************************
 */

/************************************
 * INCLUDES
 ************************************/
#include "Inc/Dynamite.hpp"
#include "Barometric.hpp"
#include "PrintUtils.hpp"
#include "SystemDefines.hpp"
#include <cmath>

/************************************
 * PRIVATE MACROS AND DEFINES
 ************************************/
namespace {
float MaxF(float a, float b)
{
	return (a > b) ? a : b;
}

const char* const SOURCE_NAMES[DYN_SRC_COUNT] = {"BARO0", "BARO1", "IMU0", "IMU1"};
}

/************************************
 * FUNCTION DEFINITIONS
 ************************************/
void Dynamite::Init(const TareResult& tare)
{
	Reset();

	for (uint8_t i = 0; i < 2; i++) {
		if (!tare.imuValid[i]) {
			continue;
		}
		madgwick[i].Init(tare.gravityBody[i]);

		// Largest axis noise, accelStd is in sensor axes
		for (uint8_t axis = 0; axis < 3; axis++) {
			accelStd[i] = MaxF(accelStd[i], tare.accelStd[i][axis]);
		}
	}

	initialized = true;
}

void Dynamite::Reset()
{
	for (uint8_t i = 0; i < 2; i++) {
		madgwick[i] = Madgwick{};
		imuAltitude[i] = 0.0f;
		imuMovingTime[i] = 0.0f;
		imuPrevAccel[i] = 0.0f;
		accelStd[i] = 0.0f;
	}
	baroDistrustUntilMs = 0;
	initialized = false;
	hasFused = false;
	measurement = EverestMeasurement{};
	status = DynamiteStatus{};
}

const EverestMeasurement& Dynamite::Update(const SensorFrame& calibrated, const TareResult& tare, bool inFlight)
{
	if (!initialized) {
		Init(tare);
	}

	MeasureBaros(calibrated, tare);
	MeasureImus(calibrated, tare, inFlight);

	bool anyValid = false;
	for (uint8_t i = 0; i < DYN_SRC_COUNT; i++) {
		anyValid = anyValid || status.sourceValid[i];
	}
	measurement.valid = anyValid;
	if (anyValid) {
		Fuse(calibrated.dt);
	}

	return measurement;
}

/**
 * @brief Altitude above the pad from each baro, each against its own tare reference.
 *        Sigma is the altitude noise the tare measured for that baro, scaled up while the
 *        baros are known to be wrong: the whole transonic window, and after apogee only a
 *        baro that jumps away from the prediction (an ejection charge spike)
 */
void Dynamite::MeasureBaros(const SensorFrame& frame, const TareResult& tare)
{
	// Last frame's fused speed decides the transonic window
	const float speed = std::fabs(measurement.velocity);
	const bool transonic = hasFused && speed > DYNAMITE_TRANSONIC_LOW_MPS && speed < DYNAMITE_TRANSONIC_HIGH_MPS;
	const bool ejectionWindow = frame.timestampMs < baroDistrustUntilMs;
	const float predicted = measurement.altitude + measurement.velocity * frame.dt;
	status.baroDistrusted = transonic;

	for (uint8_t i = 0; i < 2; i++) {
		const uint8_t src = DYN_SRC_BARO0 + i;
		status.sourceValid[src] = tare.baroValid[i] && (frame.validMask & (SENSOR_FRAME_VALID_BARO0 << i));
		if (!status.sourceValid[src]) {
			continue;
		}

		status.sourceAltitude[src] = Barometric::PressureToAltitude(frame.baro[i].pressure) - tare.groundAltitude[i];

		// Distrusting a good baro after apogee left only the drifted IMU altitude (10-28 m errors
		// in the replay test flights), so only a baro that actually jumped is distrusted
		const bool spike = ejectionWindow && hasFused
						   && std::fabs(status.sourceAltitude[src] - predicted) > DYNAMITE_EJECTION_GATE_M;
		const bool distrust = transonic || spike;
		status.baroDistrusted = status.baroDistrusted || spike;
		status.sourceSigma[src] = (distrust ? DYNAMITE_BARO_DISTRUST_SCALE : 1.0f) * MaxF(tare.altitudeStd[i], DYNAMITE_MIN_SIGMA_M);
	}
}

/**
 * @brief Altitude above the pad from each IMU: Madgwick rotates the accel into the earth
 *        frame, gravity is removed, then it is integrated twice.
 *        On the pad the integration only runs while the accel is above its noise, otherwise
 *        the IMU altitude is held at 0. In flight it always integrates.
 *        Sigma grows with integration time (accel noise integrated twice).
 *        A clipped accel (at the +/-16 g limit) reads too low: that IMU integrates the other
 *        IMU's accel instead, or the last fused accel if both are clipped
 */
void Dynamite::MeasureImus(const SensorFrame& frame, const TareResult& tare, bool inFlight)
{
	const float dt = frame.dt;

	// Pass 1: attitude, vertical accel and clipping for every IMU
	for (uint8_t i = 0; i < 2; i++) {
		const uint8_t src = DYN_SRC_IMU0 + i;
		status.sourceValid[src] = tare.imuValid[i] && madgwick[i].IsInitialized()
								  && (frame.validMask & (SENSOR_FRAME_VALID_IMU0 << i));
		status.accelClipped[i] = false;
		status.gyroClipped[i] = false;
		if (!status.sourceValid[src]) {
			continue;
		}

		for (uint8_t axis = 0; axis < 3; axis++) {
			status.accelClipped[i] = status.accelClipped[i]
				|| (std::fabs(frame.imu[i].accel[axis]) >= DYNAMITE_CLIP_FRACTION * DYNAMITE_ACCEL_RANGE_G[i] * 9.80665f);
			status.gyroClipped[i] = status.gyroClipped[i] || (std::fabs(frame.imu[i].gyro[axis]) >= DYNAMITE_GYRO_CLIP_RADPS);
		}

		madgwick[i].Update(frame.imu[i].gyro, frame.imu[i].accel, dt);
		status.tiltDeg[i] = madgwick[i].TiltDeg();
		status.accelUsed[i] = madgwick[i].AccelUsed();
		status.imuAccel[i] = madgwick[i].VerticalAccel(frame.imu[i].accel);
	}

	// The IMU reading clearly less total force than the other one is clipped. Checked every frame,
	// not only after launch detection: clipping starts at ignition. On the pad both read ~1 g
	if (status.sourceValid[DYN_SRC_IMU0] && status.sourceValid[DYN_SRC_IMU1]) {
		float force[2];
		for (uint8_t i = 0; i < 2; i++) {
			const float* a = frame.imu[i].accel;
			force[i] = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
		}
		if (force[0] < force[1] - DYNAMITE_IMU_DISAGREE_MPS2) {
			status.accelClipped[0] = true;
		} else if (force[1] < force[0] - DYNAMITE_IMU_DISAGREE_MPS2) {
			status.accelClipped[1] = true;
		}
	}

	// Replace clipped vertical accels
	const bool ok0 = status.sourceValid[DYN_SRC_IMU0] && !status.accelClipped[0];
	const bool ok1 = status.sourceValid[DYN_SRC_IMU1] && !status.accelClipped[1];
	for (uint8_t i = 0; i < 2; i++) {
		if (status.sourceValid[DYN_SRC_IMU0 + i] && status.accelClipped[i]) {
			const bool otherOk = (i == 0) ? ok1 : ok0;
			status.imuAccel[i] = otherOk ? status.imuAccel[1 - i] : MaxF(status.imuAccel[i], measurement.acceleration);
		}
	}
	if (status.accelClipped[0] || status.accelClipped[1] || status.gyroClipped[0] || status.gyroClipped[1]) {
		status.clipFrames++;
	}

	// Pass 2: integrate
	for (uint8_t i = 0; i < 2; i++) {
		const uint8_t src = DYN_SRC_IMU0 + i;
		if (!status.sourceValid[src]) {
			continue;
		}

		status.imuMoving[i] = inFlight
							  || (std::fabs(status.imuAccel[i]) > DYNAMITE_ACC_MOTION_SIGMAS * accelStd[i]);

		if (status.imuMoving[i] && dt > 0.0f) {
			// Trapezoid rule on both integrals. Using only this frame's accel for the whole step
			// overestimated the velocity while thrust ramps up (~1.5 m/s, +3 m altitude by burnout)
			const float previousVelocity = status.imuVelocity[i];
			status.imuVelocity[i] += 0.5f * (imuPrevAccel[i] + status.imuAccel[i]) * dt;
			imuAltitude[i] += 0.5f * (previousVelocity + status.imuVelocity[i]) * dt;
			imuMovingTime[i] += dt;
		} else if (!status.imuMoving[i]) {
			// Still on the pad: zero velocity update
			status.imuVelocity[i] = 0.0f;
			imuAltitude[i] = 0.0f;
			imuMovingTime[i] = 0.0f;
		}
		imuPrevAccel[i] = status.imuMoving[i] ? status.imuAccel[i] : 0.0f;

		status.sourceAltitude[src] = imuAltitude[i];
		const float accStd = inFlight ? MaxF(accelStd[i], DYNAMITE_FLIGHT_ACC_STD_MPS2) : accelStd[i];
		status.sourceSigma[src] = 0.5f * accStd * imuMovingTime[i] * imuMovingTime[i] + DYNAMITE_MIN_SIGMA_M;
	}
}

/**
 * @brief dynamite(): predict the altitude from the last fused state, weight every source by
 *        1 / (|predicted - measured| + sigma + eps), normalise, and blend.
 *        A source that drifts away from the others loses trust the same frame
 */
void Dynamite::Fuse(float dt)
{
	// Net vertical accel: mean of the valid IMUs
	float accSum = 0.0f;
	uint8_t accCount = 0;
	for (uint8_t i = 0; i < 2; i++) {
		if (status.sourceValid[DYN_SRC_IMU0 + i]) {
			accSum += status.imuAccel[i];
			accCount++;
		}
	}
	const float acceleration = (accCount > 0) ? (accSum / accCount) : 0.0f;

	// Accel only drives the prediction while an IMU is moving (always in flight).
	// On the pad its noise would otherwise integrate into the velocity
	bool anyMoving = false;
	for (uint8_t i = 0; i < 2; i++) {
		anyMoving = anyMoving || (status.sourceValid[DYN_SRC_IMU0 + i] && status.imuMoving[i]);
	}
	const float feedForwardAcc = anyMoving ? acceleration : 0.0f;

	status.predicted = hasFused
					   ? (measurement.altitude + measurement.velocity * dt + 0.5f * feedForwardAcc * dt * dt)
					   : 0.0f;

	float weightSum = 0.0f;
	for (uint8_t i = 0; i < DYN_SRC_COUNT; i++) {
		status.weight[i] = 0.0f;
		if (status.sourceValid[i]) {
			status.weight[i] = 1.0f / (std::fabs(status.predicted - status.sourceAltitude[i])
									   + status.sourceSigma[i] + DYNAMITE_TRUST_EPS_M);
			weightSum += status.weight[i];
		}
	}

	float fused = 0.0f;
	float variance = 0.0f;
	for (uint8_t i = 0; i < DYN_SRC_COUNT; i++) {
		status.weight[i] /= weightSum;
		fused += status.weight[i] * status.sourceAltitude[i];
		variance += status.weight[i] * status.weight[i] * status.sourceSigma[i] * status.sourceSigma[i];
	}

	// Alpha-beta velocity: propagate with the accel, correct with a small share of the
	// altitude residual. A raw (fused - last) / dt turns 0.2 m of baro jitter into ~10 m/s at 50 Hz
	if (hasFused && dt > 0.0f) {
		const float residual = fused - status.predicted;
		const float scaled = DYNAMITE_VEL_BETA * dt / DYNAMITE_VEL_BETA_DT_S;
		const float beta = (scaled < 1.0f) ? scaled : 1.0f;
		measurement.velocity += feedForwardAcc * dt + (beta / dt) * residual;
	} else {
		measurement.velocity = 0.0f;
	}
	measurement.altitude = fused;
	measurement.acceleration = acceleration;
	measurement.sigma = std::sqrt(variance);
	hasFused = true;
}

void Dynamite::PrintLine() const
{
	SOAR_PRINT("DYN | B0=%ld B1=%ld I0=%ld I1=%ld | pred=%ld alt=%ld mm vel=%ld mm/s acc=%ld mm/s2\n",
			   ToMilli(status.sourceAltitude[DYN_SRC_BARO0]), ToMilli(status.sourceAltitude[DYN_SRC_BARO1]),
			   ToMilli(status.sourceAltitude[DYN_SRC_IMU0]), ToMilli(status.sourceAltitude[DYN_SRC_IMU1]),
			   ToMilli(status.predicted), ToMilli(measurement.altitude),
			   ToMilli(measurement.velocity), ToMilli(measurement.acceleration));
}

void Dynamite::PrintStatus() const
{
	SOAR_PRINT("Dynamite: valid=%d alt=%ld mm vel=%ld mm/s acc=%ld mm/s2 sigma=%ld mm predicted=%ld mm\n",
			   measurement.valid ? 1 : 0, ToMilli(measurement.altitude), ToMilli(measurement.velocity),
			   ToMilli(measurement.acceleration), ToMilli(measurement.sigma), ToMilli(status.predicted));

	for (uint8_t i = 0; i < DYN_SRC_COUNT; i++) {
		SOAR_PRINT("  %s valid=%d alt=%ld mm sigma=%ld mm trust (x1000)=%ld\n",
				   SOURCE_NAMES[i], status.sourceValid[i] ? 1 : 0,
				   ToMilli(status.sourceAltitude[i]), ToMilli(status.sourceSigma[i]),
				   ToMilli(status.weight[i]));
	}

	for (uint8_t i = 0; i < 2; i++) {
		SOAR_PRINT("  IMU%u acc=%ld mm/s2 vel=%ld mm/s moving=%d tilt=%ld mdeg accelUsed=%d clip(acc/gyro)=%d/%d\n",
				   static_cast<unsigned int>(i), ToMilli(status.imuAccel[i]), ToMilli(status.imuVelocity[i]),
				   status.imuMoving[i] ? 1 : 0, ToMilli(status.tiltDeg[i]), status.accelUsed[i] ? 1 : 0,
				   status.accelClipped[i] ? 1 : 0, status.gyroClipped[i] ? 1 : 0);
	}
	SOAR_PRINT("  clipped frames=%lu baroDistrusted=%d\n",
			   static_cast<unsigned long>(status.clipFrames), status.baroDistrusted ? 1 : 0);
}
