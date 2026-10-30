/*
 * ApogeeDetection.cpp
 *
 *  Created on: Oct 24, 2026
 *      Author: Aly Masani
 */

/**
 ********************************************************************************
 * @file    ApogeeDetection.cpp
 * @author  Aly Masani
 * @date    Oct 24, 2026
 * @brief   Apogee detection, ported from the legacy HALO::apogeeDetection() /
 *          calculateConfidence() / calculateVelocityConfidence() /
 *          calculateAccelerationConfidence() (Legacy/HALO.hpp), with a velocity backup
 ********************************************************************************
 */

/************************************
 * INCLUDES
 ************************************/
#include "Inc/ApogeeDetection.hpp"
#include "PrintUtils.hpp"
#include "SystemDefines.hpp"
#include <cmath>

/************************************
 * PRIVATE MACROS AND DEFINES
 ************************************/
namespace {
const char* SourceName(ApogeeSource s)
{
	switch (s)
	{
	case ApogeeSource::CONFIDENCE: return "confidence";
	case ApogeeSource::BACKUP:     return "velocity backup";
	case ApogeeSource::RESTORED:   return "restored after power loss";
	case ApogeeSource::MISSED:     return "missed during power loss";
	default:                       return "none";
	}
}
}

/************************************
 * FUNCTION DEFINITIONS
 ************************************/
void ApogeeDetection::Reset()
{
	head = 0;
	count = 0;
	maxAvgAltitude = -1e30f;
	altitudeAccumulator = 0.0f;
	status = ApogeeStatus{};
}

void ApogeeDetection::Restore(bool reached, float maxAltitude, uint32_t timestampMs)
{
	Reset();
	status.maxAltitude = maxAltitude;
	status.maxTimeMs = timestampMs;
	if (reached) {
		status.reached = true;
		status.source = ApogeeSource::RESTORED;
		status.detectTimeMs = timestampMs;
	}
}

bool ApogeeDetection::Update(const float x[3], const float sd[3], uint32_t timestampMs)
{
	if (status.reached) {
		return false;
	}

	if (x[0] > status.maxAltitude) {
		status.maxAltitude = x[0];
		status.maxTimeMs = timestampMs;
	}
	// Failsafe: never apogee while the measurement says the board has not moved (NoteMeasurement)
	const bool highEnough = (status.maxAltitude > APOGEE_MIN_ALT_M) && status.moved;

	// Legacy updateBuffer(): moving window of the state
	windowAlt[head] = x[0];
	windowVel[head] = x[1];
	windowAcc[head] = x[2];
	head = (head + 1) % APOGEE_WINDOW;
	if (count < APOGEE_WINDOW) {
		count++;
	}

	// Backup: velocity negative for a while
	status.negativeVelFrames = (x[1] < 0.0f && status.negativeVelFrames < 255) ? status.negativeVelFrames + 1 : 0;
	if (highEnough && status.negativeVelFrames >= APOGEE_BACKUP_FRAMES) {
		Latch(ApogeeSource::BACKUP, timestampMs);
		return true;
	}

	if (count < APOGEE_WINDOW) {
		return false;
	}

	float sumAlt = 0.0f, sumVel = 0.0f, sumAcc = 0.0f;
	for (uint8_t i = 0; i < APOGEE_WINDOW; i++) {
		sumAlt += windowAlt[i];
		sumVel += windowVel[i];
		sumAcc += windowAcc[i];
	}

	if (Confidence(sumAlt / APOGEE_WINDOW, sumVel / APOGEE_WINDOW, sumAcc / APOGEE_WINDOW, sd) && highEnough) {
		Latch(ApogeeSource::CONFIDENCE, timestampMs);
		return true;
	}
	return false;
}

/**
 * @brief Legacy apogeeDetection() body: scores the window averages, true when apogee
 */
bool ApogeeDetection::Confidence(float avgAlt, float avgVel, float avgAcc, const float sd[3])
{
	const float sdAlt = (sd[0] > APOGEE_CONF_MIN_SD[0]) ? sd[0] : APOGEE_CONF_MIN_SD[0];
	const float sdVel = (sd[1] > APOGEE_CONF_MIN_SD[1]) ? sd[1] : APOGEE_CONF_MIN_SD[1];
	const float sdAcc = (sd[2] > APOGEE_CONF_MIN_SD[2]) ? sd[2] : APOGEE_CONF_MIN_SD[2];

	// calculateConfidence(): how far below the highest average altitude we are
	if (avgAlt > maxAvgAltitude) {
		maxAvgAltitude = avgAlt;
	}
	const float difference = maxAvgAltitude - avgAlt;
	if (difference > 0.0f) {
		altitudeAccumulator += difference;
	}
	status.altConfidence = (altitudeAccumulator >= sdAlt)
						   ? 1.0f
						   : ((difference > altitudeAccumulator) ? difference : altitudeAccumulator) / sdAlt;

	// calculateVelocityConfidence(): only once the velocity is inside its sigma of zero
	status.velConfidence = 0.0f;
	if (avgVel < sdVel) {
		status.velConfidence = 1.0f - std::fabs(avgVel / (2.0f * sdVel));
	}

	// calculateAccelerationConfidence(): how close to free fall
	status.accConfidence = 1.0f - std::fabs(avgAcc - APOGEE_TARGET_ACC_MPS2) / (2.0f * sdAcc);

	// Same clamps as legacy (upper only)
	if (status.altConfidence > 1.0f) status.altConfidence = 1.0f;
	if (status.velConfidence > 1.0f) status.velConfidence = 1.0f;
	if (status.accConfidence > 1.0f) status.accConfidence = 1.0f;
	if (status.altConfidence == 0.0f && status.velConfidence == 0.0f) {
		status.accConfidence = 0.0f;
	}

	status.totalConfidence = APOGEE_ALT_WEIGHT * status.altConfidence
						   + APOGEE_VEL_WEIGHT * status.velConfidence
						   + APOGEE_ACC_WEIGHT * status.accConfidence;

	if (status.totalConfidence >= 1.0f) {
		status.hitOne = true;
	}

	// Rising: remember the peak
	if (status.totalConfidence > status.maxConfidence) {
		status.maxConfidence = status.totalConfidence;
		return false;
	}

	// Reached 1 and now falling
	if (status.hitOne && status.totalConfidence < status.maxConfidence) {
		return true;
	}

	// Never reached 1 but clearly peaked and dipped
	return (status.maxConfidence > APOGEE_DIP_MIN_PEAK)
		   && (status.totalConfidence <= status.maxConfidence * APOGEE_DIP_RATIO);
}

void ApogeeDetection::Latch(ApogeeSource source, uint32_t timestampMs)
{
	status.reached = true;
	status.source = source;
	status.detectTimeMs = timestampMs;
	SOAR_PRINT("ApogeeDetection: APOGEE (%s) %ld mm at t=%lums, confirmed t=%lums, scenarios -> post apogee\n",
			   SourceName(source), ToMilli(status.maxAltitude), static_cast<unsigned long>(status.maxTimeMs),
			   static_cast<unsigned long>(timestampMs));
}

void ApogeeDetection::PrintStatus() const
{
	SOAR_PRINT("Apogee: reached=%d by %s max=%ld mm at t=%lums detected t=%lums\n",
			   status.reached ? 1 : 0, SourceName(status.source), ToMilli(status.maxAltitude),
			   static_cast<unsigned long>(status.maxTimeMs), static_cast<unsigned long>(status.detectTimeMs));
	SOAR_PRINT("  confidence (x1000): alt=%ld vel=%ld acc=%ld total=%ld peak=%ld hitOne=%d negVelFrames=%u/%u moved=%d\n",
			   ToMilli(status.altConfidence), ToMilli(status.velConfidence), ToMilli(status.accConfidence),
			   ToMilli(status.totalConfidence), ToMilli(status.maxConfidence), status.hitOne ? 1 : 0,
			   static_cast<unsigned int>(status.negativeVelFrames), static_cast<unsigned int>(APOGEE_BACKUP_FRAMES), status.moved ? 1 : 0);
}
