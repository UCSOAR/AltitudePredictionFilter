/*
 * ApogeePredictor.cpp
 *
 *  Created on: Oct 24, 2026
 *      Author: Aly Masani
 */

/**
 ********************************************************************************
 * @file    ApogeePredictor.cpp
 * @author  Aly Masani
 * @date    Oct 24, 2026
 * @brief   Apogee / trajectory / landing prediction by projecting the UKF state forward
 *          with the scenario model (legacy HALO::predictNStates)
 ********************************************************************************
 */

/************************************
 * INCLUDES
 ************************************/
#include "Inc/ApogeePredictor.hpp"
#include "PrintUtils.hpp"
#include "SystemDefines.hpp"
#include <cmath>
#include <cstdio>

/************************************
 * FUNCTION DEFINITIONS
 ************************************/
bool ApogeePredictor::Update(const Ukf& ukf, uint32_t timestampMs)
{
	if (!ukf.IsInitialized()) {
		return false;
	}
	if (runs > 0 && (timestampMs - lastRunMs) < PREDICT_PERIOD_MS) {
		return false;
	}
	lastRunMs = timestampMs;
	runs++;

	const float* x = ukf.State();
	const ScenarioModel model = ukf.Scenarios();	// Copy: the projection never changes the filter

	prediction.madeAtMs = timestampMs;
	prediction.afterApogee = model.IsAfterApogee();

	if (!prediction.afterApogee) {
		const float sdVel = std::sqrt(ukf.Variance(1));
		float tLow = 0.0f, tHigh = 0.0f;
		prediction.apogeeAltitude = ProjectApogee(model, x[0], x[1], prediction.timeToApogee,
												  prediction.trajectory, &prediction.trajectoryCount);
		prediction.apogeeAltitudeLow = ProjectApogee(model, x[0], x[1] - PREDICT_SIGMAS * sdVel, tLow, nullptr, nullptr);
		prediction.apogeeAltitudeHigh = ProjectApogee(model, x[0], x[1] + PREDICT_SIGMAS * sdVel, tHigh, nullptr, nullptr);
		prediction.valid = prediction.timeToApogee < PREDICT_MAX_S;
	} else {
		ProjectLanding(model, x[0], x[1]);
		prediction.valid = prediction.timeToLanding < PREDICT_LANDING_MAX_S;
	}
	return true;
}

float ApogeePredictor::ProjectApogee(const ScenarioModel& model, float alt, float vel, float& timeToApogee,
									 float* trajectory, uint8_t* trajectoryCount) const
{
	float t = 0.0f;
	float nextPoint = 1.0f;
	if (trajectoryCount) {
		*trajectoryCount = 0;
	}

	while (vel > 0.0f && t < PREDICT_MAX_S) {
		// Midpoint step: acceleration from the scenarios at the start and at the half step
		const float a1 = model.Match(alt, vel).acceleration;
		const float midVel = vel + 0.5f * a1 * PREDICT_STEP_S;
		const float midAlt = alt + 0.5f * vel * PREDICT_STEP_S;
		const float a2 = model.Match(midAlt, midVel).acceleration;

		const float newVel = vel + a2 * PREDICT_STEP_S;
		if (newVel <= 0.0f) {
			// Apogee inside this step: stop exactly where the velocity crosses zero
			const float frac = vel / (vel - newVel);
			alt += 0.5f * vel * frac * PREDICT_STEP_S;
			t += frac * PREDICT_STEP_S;
			vel = 0.0f;
			break;
		}
		alt += 0.5f * (vel + newVel) * PREDICT_STEP_S;
		vel = newVel;
		t += PREDICT_STEP_S;

		if (trajectory && *trajectoryCount < PREDICT_TRAJECTORY_POINTS && t >= nextPoint - 1e-4f) {
			trajectory[(*trajectoryCount)++] = alt;
			nextPoint += 1.0f;
		}
	}

	timeToApogee = t;
	return alt;
}

void ApogeePredictor::ProjectLanding(const ScenarioModel& model, float alt, float vel)
{
	float t = 0.0f;
	float nextPoint = 1.0f;
	prediction.trajectoryCount = 0;

	while (alt > 0.0f && t < PREDICT_LANDING_MAX_S) {
		const float a = model.Match(alt, vel).acceleration;
		const float newVel = vel + a * PREDICT_LANDING_STEP_S;
		const float newAlt = alt + 0.5f * (vel + newVel) * PREDICT_LANDING_STEP_S;
		if (newAlt <= 0.0f) {
			const float frac = alt / (alt - newAlt);
			t += frac * PREDICT_LANDING_STEP_S;
			vel = vel + frac * (newVel - vel);
			alt = 0.0f;
			break;
		}
		alt = newAlt;
		vel = newVel;
		t += PREDICT_LANDING_STEP_S;

		while (prediction.trajectoryCount < PREDICT_TRAJECTORY_POINTS && t >= nextPoint - 1e-4f) {
			prediction.trajectory[prediction.trajectoryCount++] = alt;
			nextPoint += 1.0f;
		}
	}

	prediction.timeToLanding = t;
	prediction.landingVelocity = vel;
}

void ApogeePredictor::PrintStatus() const
{
	if (!prediction.valid) {
		SOAR_PRINT("Predictor: no prediction yet (runs in LAUNCH)\n");
		return;
	}
	if (!prediction.afterApogee) {
		SOAR_PRINT("Predictor: apogee %ld mm (low %ld / high %ld) in %ld ms, made t=%lums\n",
				   ToMilli(prediction.apogeeAltitude), ToMilli(prediction.apogeeAltitudeLow),
				   ToMilli(prediction.apogeeAltitudeHigh), ToMilli(prediction.timeToApogee),
				   static_cast<unsigned long>(prediction.madeAtMs));
	} else {
		SOAR_PRINT("Predictor: landing in %ld ms at %ld mm/s, made t=%lums\n",
				   ToMilli(prediction.timeToLanding), ToMilli(prediction.landingVelocity),
				   static_cast<unsigned long>(prediction.madeAtMs));
	}
	// One print for the whole trajectory (one SOAR_PRINT per point would flood the print queue)
	char line[220];
	int len = std::snprintf(line, sizeof(line), "  trajectory (m, every 1 s):");
	for (uint8_t i = 0; i < prediction.trajectoryCount && len > 0 && len < static_cast<int>(sizeof(line)) - 8; i++) {
		len += std::snprintf(line + len, sizeof(line) - len, " %ld", static_cast<long>(prediction.trajectory[i]));
	}
	SOAR_PRINT("%s\n", line);
}
