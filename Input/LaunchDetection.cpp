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
 * @brief   Internal launch detection: the dynamite fused altitude must leave the ground
 *          noise measured by the tare. Uses sensor data only, never the rocket state
 ********************************************************************************
 */

/************************************
 * INCLUDES
 ************************************/
#include "Inc/LaunchDetection.hpp"
#include "Inc/PrintUtils.hpp"
#include "SystemDefines.hpp"

/************************************
 * PRIVATE MACROS AND DEFINES
 ************************************/
namespace {
float MaxF(float a, float b)
{
	return (a > b) ? a : b;
}
}

/************************************
 * FUNCTION DEFINITIONS
 ************************************/
void LaunchDetection::Reset()
{
	status = LaunchDetectionStatus{};
}

bool LaunchDetection::Update(const EverestMeasurement& everest, uint32_t timestampMs)
{
	if (status.launched) {
		return true;
	}
	if (!everest.valid) {
		return false;
	}

	// Fused altitude above the ground noise for a few frames in a row
	status.altitude = everest.altitude;
	status.threshold = MaxF(LAUNCH_ALT_SIGMAS * everest.sigma, LAUNCH_ALT_MIN_M);
	if (status.altitude > status.threshold) {
		status.confirmCount++;
	} else {
		status.confirmCount = 0;
	}

	if (status.confirmCount >= LAUNCH_CONFIRM_FRAMES) {
		Latch(timestampMs);
	}

	return status.launched;
}

void LaunchDetection::Latch(uint32_t timestampMs)
{
	status.launched = true;
	status.launchTimeMs = timestampMs;

	SOAR_PRINT("LaunchDetection: LAUNCH detected at t=%lums\n", static_cast<unsigned long>(timestampMs));
	PrintStatus();
}

void LaunchDetection::PrintStatus() const
{
	SOAR_PRINT("  launch=%d t=%lums alt=%ld mm thr=%ld mm confirm=%u/%u\n",
			   status.launched ? 1 : 0, static_cast<unsigned long>(status.launchTimeMs),
			   ToMilli(status.altitude), ToMilli(status.threshold),
			   static_cast<unsigned int>(status.confirmCount), static_cast<unsigned int>(LAUNCH_CONFIRM_FRAMES));
}
