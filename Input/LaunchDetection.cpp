/*
 * LaunchDetection.cpp
 *
 *  Created on: Oct 17, 2026
 *      Author: Aly Masani
 */

/**
 ********************************************************************************
 * @file    LaunchDetection.cpp
 * @author  Aly Masani
 * @date    Oct 17, 2026
 * @brief   Internal launch detection, modelled on the legacy EverestTask dynamite() /
 *          recalculateGain(): every source gets a trust weight from how close it is to
 *          the predicted altitude, the weighted altitude is compared against the ground
 *          noise measured by the tare. Uses sensor data only, never the rocket state
 ********************************************************************************
 */

/************************************
 * INCLUDES
 ************************************/
#include "Inc/LaunchDetection.hpp"
#include "Inc/Barometric.hpp"
#include "Inc/PrintUtils.hpp"
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

const char* const SOURCE_NAMES[LAUNCH_SRC_COUNT] = {"BARO0", "BARO1", "IMU"};
}

/************************************
 * FUNCTION DEFINITIONS
 ************************************/
void LaunchDetection::Reset()
{
	imuAltitude = 0.0f;
	imuMovingTime = 0.0f;
	hasFused = false;
	status = LaunchDetectionStatus{};
}

bool LaunchDetection::Update(const SensorFrame& calibrated, const TareResult& tare)
{
	if (status.launched) {
		return true;
	}

	MeasureBaros(calibrated, tare);
	MeasureImu(calibrated, tare);

	bool anyValid = false;
	for (uint8_t i = 0; i < LAUNCH_SRC_COUNT; i++) {
		anyValid = anyValid || status.sourceValid[i];
	}
	if (!anyValid) {
		return false;
	}

	FuseAltitude(calibrated.dt);

	// Launch: fused altitude above the ground noise for a few frames in a row
	status.threshold = MaxF(LAUNCH_ALT_SIGMAS * status.sigma, LAUNCH_ALT_MIN_M);
	if (status.altitude > status.threshold) {
		status.confirmCount++;
	} else {
		status.confirmCount = 0;
	}

	if (status.confirmCount >= LAUNCH_CONFIRM_FRAMES) {
		Latch(calibrated.timestampMs);
	}

	return status.launched;
}

/**
 * @brief Altitude above the pad from each baro, each against its own tare reference.
 *        Sigma is the altitude noise the tare measured for that baro
 */
void LaunchDetection::MeasureBaros(const SensorFrame& frame, const TareResult& tare)
{
	for (uint8_t i = 0; i < 2; i++) {
		const uint8_t src = LAUNCH_SRC_BARO0 + i;
		status.sourceValid[src] = tare.baroValid[i] && (frame.validMask & (SENSOR_FRAME_VALID_BARO0 << i));
		if (!status.sourceValid[src]) {
			continue;
		}

		status.sourceAltitude[src] = Barometric::PressureToAltitude(frame.baro[i].pressure) - tare.groundAltitude[i];
		status.sourceSigma[src] = MaxF(tare.altitudeStd[i], LAUNCH_MIN_SIGMA_M);
	}
}

/**
 * @brief Altitude above the pad from the IMUs: net vertical accel (body +Z, pad gravity
 *        reading removed) integrated twice. Only integrates while the accel is above its
 *        noise, otherwise the rocket is sitting still and the IMU altitude is held at 0.
 *        Its sigma grows with integration time (accel noise integrated twice)
 */
void LaunchDetection::MeasureImu(const SensorFrame& frame, const TareResult& tare)
{
	float accelSum = 0.0f;
	float accelStd = 0.0f;
	uint8_t used = 0;

	for (uint8_t i = 0; i < 2; i++) {
		if (tare.imuValid[i] && (frame.validMask & (SENSOR_FRAME_VALID_IMU0 << i))) {
			accelSum += frame.imu[i].accel[2] - tare.gravityBody[i][2];
			// Largest axis noise, accelStd is in sensor axes
			for (uint8_t axis = 0; axis < 3; axis++) {
				accelStd = MaxF(accelStd, tare.accelStd[i][axis]);
			}
			used++;
		}
	}

	status.sourceValid[LAUNCH_SRC_IMU] = (used > 0);
	if (used == 0) {
		return;
	}

	status.imuAccel = accelSum / used;
	status.imuMoving = std::fabs(status.imuAccel) > LAUNCH_ACC_MOTION_SIGMAS * accelStd;

	const float dt = frame.dt;
	if (status.imuMoving && dt > 0.0f) {
		status.imuVelocity += status.imuAccel * dt;
		imuAltitude += status.imuVelocity * dt;
		imuMovingTime += dt;
	} else if (!status.imuMoving) {
		// Still on the pad: zero velocity update
		status.imuVelocity = 0.0f;
		imuAltitude = 0.0f;
		imuMovingTime = 0.0f;
	}

	status.sourceAltitude[LAUNCH_SRC_IMU] = imuAltitude;
	status.sourceSigma[LAUNCH_SRC_IMU] = 0.5f * accelStd * imuMovingTime * imuMovingTime + LAUNCH_MIN_SIGMA_M;
}

/**
 * @brief dynamite(): predict the altitude from the last fused state, weight every source by
 *        1 / (|predicted - measured| + sigma + eps), normalise, and fuse.
 *        A source that drifts away from the others loses trust the same frame
 */
void LaunchDetection::FuseAltitude(float dt)
{
	status.predicted = hasFused ? (status.altitude + status.velocity * dt) : 0.0f;

	float weightSum = 0.0f;
	for (uint8_t i = 0; i < LAUNCH_SRC_COUNT; i++) {
		status.weight[i] = 0.0f;
		if (status.sourceValid[i]) {
			status.weight[i] = 1.0f / (std::fabs(status.predicted - status.sourceAltitude[i])
									   + status.sourceSigma[i] + LAUNCH_TRUST_EPS_M);
			weightSum += status.weight[i];
		}
	}

	float fused = 0.0f;
	float variance = 0.0f;
	for (uint8_t i = 0; i < LAUNCH_SRC_COUNT; i++) {
		status.weight[i] /= weightSum;
		fused += status.weight[i] * status.sourceAltitude[i];
		variance += status.weight[i] * status.weight[i] * status.sourceSigma[i] * status.sourceSigma[i];
	}

	status.velocity = (hasFused && dt > 0.0f) ? (fused - status.altitude) / dt : 0.0f;
	status.altitude = fused;
	status.sigma = std::sqrt(variance);
	hasFused = true;
}

void LaunchDetection::Latch(uint32_t timestampMs)
{
	status.launched = true;
	status.launchTimeMs = timestampMs;

	SOAR_PRINT("LaunchDetection: LAUNCH detected at t=%lums\n", static_cast<unsigned long>(timestampMs));
	PrintStatus();
}

void LaunchDetection::PrintAltitudes() const
{
	SOAR_PRINT("ALT launch=%d | B0=%ld B1=%ld IMU=%ld | pred=%ld fused=%ld thr=%ld mm | conf=%u/%u\n",
			   status.launched ? 1 : 0,
			   ToMilli(status.sourceAltitude[LAUNCH_SRC_BARO0]),
			   ToMilli(status.sourceAltitude[LAUNCH_SRC_BARO1]),
			   ToMilli(status.sourceAltitude[LAUNCH_SRC_IMU]),
			   ToMilli(status.predicted), ToMilli(status.altitude), ToMilli(status.threshold),
			   static_cast<unsigned int>(status.confirmCount), static_cast<unsigned int>(LAUNCH_CONFIRM_FRAMES));
}

void LaunchDetection::PrintStatus() const
{
	SOAR_PRINT("  launch=%d t=%lums alt=%ld mm vel=%ld mm/s sigma=%ld mm thr=%ld mm confirm=%u/%u\n",
			   status.launched ? 1 : 0, static_cast<unsigned long>(status.launchTimeMs),
			   ToMilli(status.altitude), ToMilli(status.velocity), ToMilli(status.sigma),
			   ToMilli(status.threshold),
			   static_cast<unsigned int>(status.confirmCount), static_cast<unsigned int>(LAUNCH_CONFIRM_FRAMES));

	SOAR_PRINT("  predicted=%ld mm imuAcc=%ld mm/s2 imuVel=%ld mm/s moving=%d\n",
			   ToMilli(status.predicted), ToMilli(status.imuAccel), ToMilli(status.imuVelocity),
			   status.imuMoving ? 1 : 0);

	for (uint8_t i = 0; i < LAUNCH_SRC_COUNT; i++) {
		SOAR_PRINT("  %s valid=%d alt=%ld mm sigma=%ld mm trust (x1000)=%ld\n",
				   SOURCE_NAMES[i], status.sourceValid[i] ? 1 : 0,
				   ToMilli(status.sourceAltitude[i]), ToMilli(status.sourceSigma[i]),
				   ToMilli(status.weight[i]));
	}
}
