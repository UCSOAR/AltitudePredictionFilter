/*
 * FilterOutput.cpp
 *
 *  Created on: Oct 24, 2026
 *      Author: Aly Masani
 */

/**
 ********************************************************************************
 * @file    FilterOutput.cpp
 * @author  Aly Masani
 * @date    Oct 24, 2026
 * @brief   OUTPUT stage: builds the FilterData that is published and logged to flash
 ********************************************************************************
 */

/************************************
 * INCLUDES
 ************************************/
#include "Inc/FilterOutput.hpp"
#include "PrintUtils.hpp"
#include "SystemDefines.hpp"
#include <cmath>

/************************************
 * FUNCTION DEFINITIONS
 ************************************/
bool FilterOutput::Due(ApfState state, uint32_t timestampMs)
{
	const uint32_t period = (state == ApfState::LAUNCH) ? FILTER_OUTPUT_FLIGHT_PERIOD_MS : FILTER_OUTPUT_SLOW_PERIOD_MS;
	if (published == 0 || (timestampMs - lastMs) >= period) {
		lastMs = timestampMs;
		return true;
	}
	return false;
}

FilterData FilterOutput::Build(uint32_t timestampMs, ApfState state, float groundAltitude,
							   const Dynamite& dynamite, const LaunchDetection& launch, const Ukf& ukf,
							   const ApogeeDetection& apogee, const ApogeePredictor& predictor)
{
	FilterData d{};
	d.timestampMs = timestampMs;
	d.groundAltitude = groundAltitude;
	d.apfState = static_cast<uint8_t>(state);

	// Best estimate: the UKF once it runs (LAUNCH / LANDED hold), Dynamite before that
	if (ukf.IsInitialized()) {
		const float* x = ukf.State();
		d.altitude = x[0];
		d.velocity = x[1];
		d.acceleration = x[2];
		d.altitudeStd = std::sqrt(ukf.Variance(0));
		d.velocityStd = std::sqrt(ukf.Variance(1));
		d.flags |= FILTER_FLAG_UKF;
		if (ukf.GetStatus().inflation > 1.0f) {
			d.flags |= FILTER_FLAG_MODEL_MISFIT;
		}
		if (ukf.IsAfterApogee() && SCENARIOS[ukf.GetStatus().scenario.caseA].afterSynthetic) {
			d.flags |= FILTER_FLAG_SYNTHETIC_DESCENT;
		}
	} else {
		const EverestMeasurement& ev = dynamite.GetMeasurement();
		d.altitude = ev.altitude;
		d.velocity = ev.velocity;
		d.acceleration = ev.acceleration;
		d.altitudeStd = ev.sigma;
	}

	const DynamiteStatus& ds = dynamite.GetStatus();
	if (ds.accelClipped[0] || ds.accelClipped[1] || ds.gyroClipped[0] || ds.gyroClipped[1]) {
		d.flags |= FILTER_FLAG_IMU_CLIPPED;
	}
	if (ds.baroDistrusted) {
		d.flags |= FILTER_FLAG_BARO_DISTRUSTED;
	}
	// From the APF state, so it follows whatever decides launch (the board's FlightTask on DAQLite)
	(void)launch;
	if (state == ApfState::LAUNCH || state == ApfState::LANDED) {
		d.flags |= FILTER_FLAG_LAUNCHED;
	}
	if (apogee.IsReached()) {
		d.flags |= FILTER_FLAG_APOGEE;
	}

	const ApogeePrediction& p = predictor.Get();
	if (p.valid) {
		d.flags |= FILTER_FLAG_PREDICTION;
		// Time measured from when the prediction was made, so it stays right between predictions
		const float since = (timestampMs - p.madeAtMs) / 1000.0f;
		if (!p.afterApogee) {
			d.predictedApogee = p.apogeeAltitude;
			d.timeToEvent = p.timeToApogee - since;
		} else {
			d.predictedApogee = apogee.GetStatus().maxAltitude;
			d.timeToEvent = p.timeToLanding - since;
		}
	}
	return d;
}

void FilterOutput::PrintLine(const FilterData& d)
{
	SOAR_PRINT("OUT t=%lums state=%u alt=%ld mm vel=%ld mm/s acc=%ld mm/s2 sd=%ld mm | apogee=%ld mm in %ld ms | flags=0x%02X\n",
			   static_cast<unsigned long>(d.timestampMs), static_cast<unsigned int>(d.apfState),
			   ToMilli(d.altitude), ToMilli(d.velocity), ToMilli(d.acceleration), ToMilli(d.altitudeStd),
			   ToMilli(d.predictedApogee), ToMilli(d.timeToEvent), static_cast<unsigned int>(d.flags));
}

void FilterOutput::PrintLast() const
{
	if (!hasLast) {
		SOAR_PRINT("APF output: nothing published yet\n");
		return;
	}
	SOAR_PRINT("APF output: %lu published, ground %ld mm ASL\n",
			   static_cast<unsigned long>(published), ToMilli(last.groundAltitude));
	PrintLine(last);
	SOAR_PRINT("  flags: launched=%d apogee=%d prediction=%d ukf=%d imuClip=%d baroDistrust=%d misfit=%d syntheticDescent=%d\n",
			   (last.flags & FILTER_FLAG_LAUNCHED) ? 1 : 0, (last.flags & FILTER_FLAG_APOGEE) ? 1 : 0,
			   (last.flags & FILTER_FLAG_PREDICTION) ? 1 : 0, (last.flags & FILTER_FLAG_UKF) ? 1 : 0,
			   (last.flags & FILTER_FLAG_IMU_CLIPPED) ? 1 : 0, (last.flags & FILTER_FLAG_BARO_DISTRUSTED) ? 1 : 0,
			   (last.flags & FILTER_FLAG_MODEL_MISFIT) ? 1 : 0, (last.flags & FILTER_FLAG_SYNTHETIC_DESCENT) ? 1 : 0);
}
