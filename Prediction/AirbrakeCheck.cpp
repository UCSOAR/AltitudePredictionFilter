/*
 * AirbrakeCheck.cpp
 *
 *  Created on: Oct 30, 2026
 *      Author: Aly Masani
 */

/**
 ********************************************************************************
 * @file    AirbrakeCheck.cpp
 * @author  Aly Masani
 * @date    Oct 30, 2026
 * @brief   Airbrake level from the predicted apogee, logged only (nothing is actuated)
 ********************************************************************************
 */

/************************************
 * INCLUDES
 ************************************/
#include "Inc/AirbrakeCheck.hpp"
#include "PrintUtils.hpp"
#include "SystemDefines.hpp"
#include <cmath>

/************************************
 * FUNCTION DEFINITIONS
 ************************************/
bool AirbrakeCheck::Update(const ApogeePrediction& prediction, bool burnedOut, bool moved, bool afterApogee,
						   uint32_t timestampMs)
{
	if (status.flagged || afterApogee || !burnedOut || !moved) {
		return false;
	}
	if (!prediction.valid || prediction.afterApogee || (timestampMs - prediction.madeAtMs) > AIRBRAKE_MAX_PREDICTION_AGE_MS) {
		return false;
	}
	if (status.checks > 0 && (timestampMs - lastCheckMs) < AIRBRAKE_CHECK_PERIOD_MS) {
		return false;
	}
	lastCheckMs = timestampMs;
	if (status.checks < 0xFFFF) {
		status.checks++;
	}

	status.overshoot = prediction.apogeeAltitude - target;
	const float levels = std::floor(status.overshoot / AIRBRAKE_M_PER_LEVEL);
	status.level = (levels <= 0.0f) ? 0 : ((levels >= AIRBRAKE_MAX_LEVEL) ? AIRBRAKE_MAX_LEVEL : static_cast<uint8_t>(levels));

	// Low / high are -/+ PREDICT_SIGMAS sigma. Level 0: confidence that the apogee stays below the
	// level 1 overshoot; level >= 1: that it is above this level's overshoot
	float sigma = (prediction.apogeeAltitudeHigh - prediction.apogeeAltitudeLow) / (2.0f * PREDICT_SIGMAS);
	if (sigma < 1.0f) {
		sigma = 1.0f;
	}
	const float edge = target + AIRBRAKE_M_PER_LEVEL * ((status.level == 0) ? 1.0f : status.level);
	const float above = 0.5f * (1.0f + std::erf((prediction.apogeeAltitude - edge) / (sigma * 1.41421356f)));
	status.confidence = (status.level == 0) ? (1.0f - above) : above;

	if (status.level >= 1) {
		status.flagged = true;
		status.flaggedMs = timestampMs;
		SOAR_PRINT("AirbrakeCheck: FLAG level %u, confidence %ld/1000, overshoot %ld mm at t=%lums (log only)\n",
				   static_cast<unsigned int>(status.level), ToMilli(status.confidence), ToMilli(status.overshoot),
				   static_cast<unsigned long>(timestampMs));
	}
	return true;
}

int32_t AirbrakeCheck::EventValue() const
{
	float o = std::round(status.overshoot);
	o = (o < -32768.0f) ? -32768.0f : ((o > 32767.0f) ? 32767.0f : o);
	const uint16_t overshoot = static_cast<uint16_t>(static_cast<int16_t>(o));
	const uint8_t conf = static_cast<uint8_t>(std::round(status.confidence * 255.0f));
	return static_cast<int32_t>(static_cast<uint32_t>(status.level) | (static_cast<uint32_t>(conf) << 8)
								| (static_cast<uint32_t>(overshoot) << 16));
}

void AirbrakeCheck::PrintStatus() const
{
	SOAR_PRINT("Airbrakes (log only): flagged=%d level=%u confidence=%ld/1000 overshoot=%ld mm target=%ld mm "
			   "flagged t=%lums checks=%u\n",
			   status.flagged ? 1 : 0, static_cast<unsigned int>(status.level), ToMilli(status.confidence),
			   ToMilli(status.overshoot), ToMilli(target), static_cast<unsigned long>(status.flaggedMs),
			   static_cast<unsigned int>(status.checks));
}
