/*
 * ScenarioModel.cpp
 *
 *  Created on: Oct 24, 2026
 *      Author: Aly Masani
 */

/**
 ********************************************************************************
 * @file    ScenarioModel.cpp
 * @author  Aly Masani
 * @date    Oct 24, 2026
 * @brief   Expected acceleration from the simulated flights (Estimate/ScenarioData),
 *          used by the UKF predict step. Pre apogee uses each case's BEFORE table,
 *          post apogee its AFTER table (the scenario mode switch in the diagram)
 ********************************************************************************
 */

/************************************
 * INCLUDES
 ************************************/
#include "Inc/ScenarioModel.hpp"
#include <cmath>

/************************************
 * FUNCTION DEFINITIONS
 ************************************/
ScenarioMatch ScenarioModel::Match(float altitude, float velocity) const
{
	// Two closest cases by velocity distance at this altitude
	float bestDist[2] = {1e30f, 1e30f};
	float bestAcc[2] = {0.0f, 0.0f};
	float bestVel[2] = {0.0f, 0.0f};
	uint8_t bestIdx[2] = {0, 0};

	for (uint8_t i = 0; i < SCENARIO_COUNT; i++) {
		if (static_cast<int8_t>(i) == excludedCase) {
			continue;
		}
		float caseVel = 0.0f;
		float caseAcc = 0.0f;
		Lookup(SCENARIOS[i], altitude, caseVel, caseAcc);

		const float dist = std::fabs(velocity - caseVel);
		if (dist < bestDist[0]) {
			bestDist[1] = bestDist[0]; bestAcc[1] = bestAcc[0]; bestVel[1] = bestVel[0]; bestIdx[1] = bestIdx[0];
			bestDist[0] = dist; bestAcc[0] = caseAcc; bestVel[0] = caseVel; bestIdx[0] = i;
		} else if (dist < bestDist[1]) {
			bestDist[1] = dist; bestAcc[1] = caseAcc; bestVel[1] = caseVel; bestIdx[1] = i;
		}
	}

	ScenarioMatch match{};
	match.caseA = bestIdx[0];
	match.caseB = bestIdx[1];
	match.velocityError = velocity - bestVel[0];

	// Inverse distance blend
	const float wA = 1.0f / (bestDist[0] + SCENARIO_MATCH_EPS_MPS);
	const float wB = (SCENARIO_COUNT > 1) ? 1.0f / (bestDist[1] + SCENARIO_MATCH_EPS_MPS) : 0.0f;
	match.weightA = wA / (wA + wB);
	match.acceleration = match.weightA * bestAcc[0] + (1.0f - match.weightA) * bestAcc[1];

	return match;
}

void ScenarioModel::Lookup(const ScenarioCase& c, float altitude, float& velocity, float& acceleration) const
{
	const ScenarioPoint* rows = afterApogee ? c.after : c.before;
	const uint16_t count = afterApogee ? c.afterCount : c.beforeCount;

	if (count == 0) {
		velocity = 0.0f;
		acceleration = 0.0f;
		return;
	}

	// Before apogee altitude rises row by row, after apogee it falls.
	// Flip the comparison so one binary search handles both
	const bool rising = !afterApogee;
	auto below = [rising](float rowAlt, float alt) { return rising ? (rowAlt <= alt) : (rowAlt >= alt); };

	// Clamp to the ends of the table
	if (!below(rows[0].altitude, altitude)) {
		velocity = rows[0].velocity;
		acceleration = rows[0].acceleration;
		return;
	}
	if (below(rows[count - 1].altitude, altitude)) {
		velocity = rows[count - 1].velocity;
		acceleration = rows[count - 1].acceleration;
		return;
	}

	// rows[lo] is on the start side of altitude, rows[hi] past it
	uint16_t lo = 0;
	uint16_t hi = count - 1;
	while (hi - lo > 1) {
		const uint16_t mid = static_cast<uint16_t>((lo + hi) / 2);
		if (below(rows[mid].altitude, altitude)) {
			lo = mid;
		} else {
			hi = mid;
		}
	}

	const float span = rows[hi].altitude - rows[lo].altitude;
	const float w = (std::fabs(span) > 1e-6f) ? (altitude - rows[lo].altitude) / span : 0.0f;
	velocity = rows[lo].velocity + w * (rows[hi].velocity - rows[lo].velocity);
	acceleration = rows[lo].acceleration + w * (rows[hi].acceleration - rows[lo].acceleration);
}
