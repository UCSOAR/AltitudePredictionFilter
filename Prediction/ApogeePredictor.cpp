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
		prediction.apogeeAltitude = ProjectApogee(model, x[0], x[1], x[2], prediction.timeToApogee,
												  prediction.trajectory, &prediction.trajectoryCount);
		prediction.apogeeAltitudeLow = ProjectApogee(model, x[0], x[1] - PREDICT_SIGMAS * sdVel, x[2], tLow, nullptr, nullptr);
		prediction.apogeeAltitudeHigh = ProjectApogee(model, x[0], x[1] + PREDICT_SIGMAS * sdVel, x[2], tHigh, nullptr, nullptr);
		prediction.valid = prediction.timeToApogee < PREDICT_MAX_S;
	} else {
		ProjectLanding(model, x[0], x[1], x[2]);
		prediction.valid = prediction.timeToLanding < PREDICT_LANDING_MAX_S;
	}
	return true;
}

float ApogeePredictor::ProjectApogee(const ScenarioModel& model, float alt, float vel, float acc0, float& timeToApogee,
									 float* trajectory, uint8_t* trajectoryCount) const
{
#if APF_SCENARIO_MODEL == APF_SCENARIO_MODEL_LEGACY
	// Legacy HALO::predictNStates: step the state itself through the legacy model, with its own
	// gain memory starting fresh and no 50 / 50 reset (legacy commented it out there)
	LegacyGainMemory memory;
	memory.Reset();
	float state[3] = {alt, vel, acc0};
	float t = 0.0f;
	float nextPoint = 1.0f;
	if (trajectoryCount) {
		*trajectoryCount = 0;
	}
	while (state[1] > 0.0f && t < PREDICT_MAX_S) {
		const LegacyStep step = model.LegacyPredict(state, PREDICT_STEP_S, memory, false);
		if (step.state[1] <= 0.0f) {
			// Apogee inside this step: stop where the velocity crosses zero
			const float frac = state[1] / (state[1] - step.state[1]);
			state[0] += frac * (step.state[0] - state[0]);
			t += frac * PREDICT_STEP_S;
			break;
		}
		// The scenarios end before apogee (the beforeSimsF2 cases stop a few m/s short): their last
		// point is as far as they go, take it as the apogee
		if (step.state[0] == state[0] && step.state[1] == state[1]) {
			break;
		}
		for (uint8_t i = 0; i < 3; i++) {
			state[i] = step.state[i];
		}
		t += PREDICT_STEP_S;
		if (trajectory && *trajectoryCount < PREDICT_TRAJECTORY_POINTS && t >= nextPoint - 1e-4f) {
			trajectory[(*trajectoryCount)++] = state[0];
			nextPoint += 1.0f;
		}
	}
	timeToApogee = t;
	return state[0];
#else
	// The KD-tree matches on (alt, vel, acc): carry the projected acceleration step to step
	float acc = acc0;
	auto accel = [&](float a, float v) { return model.Match(a, v, acc).acceleration; };
	float t = 0.0f;
	float nextPoint = 1.0f;
	if (trajectoryCount) {
		*trajectoryCount = 0;
	}

	while (vel > 0.0f && t < PREDICT_MAX_S) {
		// Midpoint step: acceleration from the scenarios at the start and at the half step
		const float a1 = accel(alt, vel);
		const float midVel = vel + 0.5f * a1 * PREDICT_STEP_S;
		const float midAlt = alt + 0.5f * vel * PREDICT_STEP_S;
		const float a2 = accel(midAlt, midVel);

		acc = a2;
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
#endif
}

void ApogeePredictor::ProjectLanding(const ScenarioModel& model, float alt, float vel, float acc0)
{
	float t = 0.0f;
	float nextPoint = 1.0f;
	prediction.trajectoryCount = 0;
#if APF_SCENARIO_MODEL == APF_SCENARIO_MODEL_LEGACY
	// Legacy model, as ProjectApogee: the state steps through the after apogee scenarios
	LegacyGainMemory memory;
	memory.Reset();
	float state[3] = {alt, vel, acc0};
	while (state[0] > 0.0f && t < PREDICT_LANDING_MAX_S) {
		const LegacyStep step = model.LegacyPredict(state, PREDICT_LANDING_STEP_S, memory, false);
		if (step.state[0] <= 0.0f) {
			const float frac = state[0] / (state[0] - step.state[0]);
			t += frac * PREDICT_LANDING_STEP_S;
			state[1] += frac * (step.state[1] - state[1]);
			state[0] = 0.0f;
			break;
		}
		// The scenarios end before the ground (e.g. ascent only cases): nothing left to step through
		if (step.state[0] == state[0] && step.state[1] == state[1]) {
			t = PREDICT_LANDING_MAX_S;
			break;
		}
		for (uint8_t i = 0; i < 3; i++) {
			state[i] = step.state[i];
		}
		t += PREDICT_LANDING_STEP_S;
		while (prediction.trajectoryCount < PREDICT_TRAJECTORY_POINTS && t >= nextPoint - 1e-4f) {
			prediction.trajectory[prediction.trajectoryCount++] = state[0];
			nextPoint += 1.0f;
		}
	}
	prediction.timeToLanding = t;
	prediction.landingVelocity = state[1];
#else
	float a_prev = acc0;	// KD-tree matches on (alt, vel, acc): carry the acceleration

	while (alt > 0.0f && t < PREDICT_LANDING_MAX_S) {
		const float a = model.Match(alt, vel, a_prev).acceleration;
		a_prev = a;
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
#endif
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
