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
 *          used by the UKF predict step. Legacy HALO KD-tree scenario search: nearest
 *          (altitude, velocity, acceleration) point per case, two closest cases, legacy
 *          gains between them. Pre apogee uses the BEFORE trees, post apogee the AFTER trees
 ********************************************************************************
 */

/************************************
 * INCLUDES
 ************************************/
#include "Inc/ScenarioModel.hpp"
#include "Inc/KDTree.hpp"
#include <cmath>
#include <limits>

/************************************
 * PRIVATE FUNCTIONS
 ************************************/
namespace {
/**
 * @brief One KD-tree per case and phase (legacy Scenario::createTree), built once from
 *        the constant tables and shared by every ScenarioModel copy (the predictor copies it)
 */
KDTree& Tree(uint8_t caseIndex, bool after)
{
	static KDTree trees[SCENARIO_COUNT][2];
	static bool built = false;
	if (!built) {
		for (uint8_t c = 0; c < SCENARIO_COUNT; c++) {
			for (uint8_t phase = 0; phase < 2; phase++) {
				const ScenarioPoint* rows = phase ? SCENARIOS[c].after : SCENARIOS[c].before;
				const uint16_t count = phase ? SCENARIOS[c].afterCount : SCENARIOS[c].beforeCount;
				std::vector<std::array<float, 3>> points;
				points.reserve(count);
				for (uint16_t i = 0; i < count; i++) {
					points.push_back({rows[i].altitude, rows[i].velocity, rows[i].acceleration});
				}
				trees[c][phase] = KDTree(points);
			}
		}
		built = true;
	}
	return trees[caseIndex][after ? 1 : 0];
}

/**
 * @brief Legacy HALO::predictNextValues gain calculation, unchanged (including that `both`
 *        and `vector1Further` are not reset between axes). gainV1[i] is the share of
 *        vector1 on axis i, gainV2 = 1 - gainV1
 */
void LegacyGains(const float vector1[3], const float vector2[3], const float X_in[3], float gainV1[3], float gainV2[3])
{
	bool both = false;
	bool vector1Further = false;

	for (int i = 0; i < 3; i++) {
		float distance = std::abs(vector1[i] - vector2[i]);

		if (distance < 1) {
			// similar so defaulting
			gainV1[i] = 0.5;
			gainV2[i] = 0.5;
		} else {
			if (vector1[i] > X_in[i]) {
				// below vector1
				if (vector2[i] > X_in[i]) {
					// point below both lines
					both = true;
					if (vector1[i] > vector2[i]) {
						// vector1 on top
						vector1Further = true;
					} else {
						vector1Further = false;
					}
				}
			} else {
				if (vector2[i] > X_in[i]) {
					// point between lines
				} else {
					// point above both lines
					if (vector1[i] < vector2[i]) {
						vector1Further = true;
					} else {
						vector1Further = false;
					}
					both = true;
				}
			}

			if (both) {
				float longestDistance = 0;
				float distance1 = 0;

				if (vector1Further) {
					longestDistance = std::abs(vector1[i] - X_in[i]);
					distance1 = std::abs(vector2[i] - X_in[i]);
				} else {
					longestDistance = std::abs(vector2[i] - X_in[i]);
					distance1 = std::abs(vector1[i] - X_in[i]);
				}

				float relativeFactor = longestDistance / distance1;

				if (!vector1Further) {
					gainV1[i] = relativeFactor / (relativeFactor + 1);
				} else {
					gainV1[i] = 1 - (relativeFactor / (relativeFactor + 1));
				}

				gainV2[i] = 1 - gainV1[i];
			} else {
				// point between lines
				gainV1[i] = 1 - (std::abs(vector1[i] - X_in[i]) / distance);
				gainV2[i] = 1 - gainV1[i];
			}
		}
	}
}

/**
 * @brief Legacy HALO::findNearestScenarios: nearest point of every case (KD-tree), Euclidean
 *        distance, the two lowest kept (legacy loop: strict <, starting from int max)
 */
struct Nearest
{
	float distances[SCENARIO_COUNT];
	std::array<float, 3> point[SCENARIO_COUNT];
	size_t index[SCENARIO_COUNT];
	int lowest;
	int secondLowest;
};

void FindNearest(const std::array<float, 3>& measurement, bool afterApogee, int8_t excludedCase, Nearest& n)
{
	for (uint8_t s = 0; s < SCENARIO_COUNT; s++) {
		n.distances[s] = std::numeric_limits<float>::max();
		n.index[s] = 0;
		n.point[s] = {0.0f, 0.0f, 0.0f};
		if (static_cast<int8_t>(s) == excludedCase) {
			continue;
		}
		const pointIndex result = Tree(s, afterApogee).nearest_pointIndex(measurement);
		n.point[s] = result.first;
		n.index[s] = result.second;
		n.distances[s] = std::sqrt(std::pow(measurement[0] - n.point[s][0], 2) +
								   std::pow(measurement[1] - n.point[s][1], 2) +
								   std::pow(measurement[2] - n.point[s][2], 2));
	}

	float lowestDistance = std::numeric_limits<int>::max();
	float secondLowestDistance = std::numeric_limits<int>::max();
	n.lowest = 0;
	n.secondLowest = 0;
	for (int s = 0; s < SCENARIO_COUNT; s++) {
		if (n.distances[s] < lowestDistance) {
			secondLowestDistance = lowestDistance;
			n.secondLowest = n.lowest;
			lowestDistance = n.distances[s];
			n.lowest = s;
		} else if (n.distances[s] < secondLowestDistance) {
			secondLowestDistance = n.distances[s];
			n.secondLowest = s;
		}
	}
}

/**
 * @brief Where a state sits on one case's flight, at any rate: the KD-tree point (row k) refined
 *        to the closest point on the lines row k-1 -> k and k -> k+1 (same raw units as the match).
 *        Gives that point (alt, vel, acc) and its time along the case
 */
void Locate(const ScenarioPoint* rows, uint16_t count, size_t k, const float state[3], float point[3], float& time)
{
	point[0] = rows[k].altitude;
	point[1] = rows[k].velocity;
	point[2] = rows[k].acceleration;
	time = rows[k].time;
	float best = std::numeric_limits<float>::max();

	for (int side = -1; side <= 0; side++) {
		const int a = static_cast<int>(k) + side;	// line rows[a] -> rows[a + 1]
		if (a < 0 || a + 1 >= static_cast<int>(count)) {
			continue;
		}
		const float p[3] = {rows[a].altitude, rows[a].velocity, rows[a].acceleration};
		const float d[3] = {rows[a + 1].altitude - p[0], rows[a + 1].velocity - p[1], rows[a + 1].acceleration - p[2]};
		const float len2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
		float u = 0.0f;
		if (len2 > 0.0f) {
			u = ((state[0] - p[0]) * d[0] + (state[1] - p[1]) * d[1] + (state[2] - p[2]) * d[2]) / len2;
			u = (u < 0.0f) ? 0.0f : ((u > 1.0f) ? 1.0f : u);
		}
		float q[3], dist = 0.0f;
		for (uint8_t i = 0; i < 3; i++) {
			q[i] = p[i] + u * d[i];
			dist += (state[i] - q[i]) * (state[i] - q[i]);
		}
		if (dist < best) {
			best = dist;
			for (uint8_t i = 0; i < 3; i++) {
				point[i] = q[i];
			}
			time = rows[a].time + u * (rows[a + 1].time - rows[a].time);
		}
	}
}

/**
 * @brief A case's (alt, vel, acc) at a time, linear between rows, held at the ends
 *        (legacy evaluateVectorAtTime, interpolated instead of the next row)
 */
void AtTime(const ScenarioPoint* rows, uint16_t count, float time, float out[3])
{
	if (time <= rows[0].time || count == 1) {
		out[0] = rows[0].altitude;
		out[1] = rows[0].velocity;
		out[2] = rows[0].acceleration;
		return;
	}
	if (time >= rows[count - 1].time) {
		out[0] = rows[count - 1].altitude;
		out[1] = rows[count - 1].velocity;
		out[2] = rows[count - 1].acceleration;
		return;
	}
	uint16_t lo = 0, hi = count - 1;
	while (hi - lo > 1) {
		const uint16_t mid = static_cast<uint16_t>((lo + hi) / 2);
		if (rows[mid].time <= time) {
			lo = mid;
		} else {
			hi = mid;
		}
	}
	const float w = (time - rows[lo].time) / (rows[hi].time - rows[lo].time);
	out[0] = rows[lo].altitude + w * (rows[hi].altitude - rows[lo].altitude);
	out[1] = rows[lo].velocity + w * (rows[hi].velocity - rows[lo].velocity);
	out[2] = rows[lo].acceleration + w * (rows[hi].acceleration - rows[lo].acceleration);
}
}

/************************************
 * FUNCTION DEFINITIONS
 ************************************/
ScenarioMatch ScenarioModel::Match(float altitude, float velocity, float acceleration) const
{
	const std::array<float, 3> measurement = {altitude, velocity, acceleration};

	Nearest n;
	FindNearest(measurement, afterApogee, excludedCase, n);
	const int lowestDistanceIndex = n.lowest;
	const int secondLowestDistanceIndex = n.secondLowest;
	const float lowestDistance = n.distances[lowestDistanceIndex];
	const float secondLowestDistance = n.distances[secondLowestDistanceIndex];

	// Legacy gains between the two nearest points, used on the acceleration axis
	const float* v1 = n.point[lowestDistanceIndex].data();
	const float* v2 = n.point[secondLowestDistanceIndex].data();
	float gainV1[3] = {0, 0, 0}, gainV2[3] = {0, 0, 0};
	LegacyGains(v1, v2, measurement.data(), gainV1, gainV2);
	// Legacy divides by the distance to a point, 0 when the state sits exactly on it (legacy caught
	// the resulting NaN later, in the Kalman gain): fall back to its "similar" default, 0.5 / 0.5
	if (!std::isfinite(gainV1[2])) {
		gainV1[2] = 0.5f;
		gainV2[2] = 0.5f;
	}

	ScenarioMatch match{};
	match.caseA = static_cast<uint8_t>(lowestDistanceIndex);
	match.caseB = static_cast<uint8_t>(secondLowestDistanceIndex);
	match.distanceA = lowestDistance;
	match.distanceB = secondLowestDistance;
	match.pointA = static_cast<uint16_t>(n.index[lowestDistanceIndex]);
	match.pointB = static_cast<uint16_t>(n.index[secondLowestDistanceIndex]);
	match.weightA = gainV1[2];
	match.acceleration = gainV1[2] * v1[2] + gainV2[2] * v2[2];
	match.velocityError = velocity - v1[1];
	return match;
}

LegacyStep ScenarioModel::LegacyPredict(const float in[3], float dt, LegacyGainMemory& memory, bool resetOnPairChange) const
{
	const std::array<float, 3> measurement = {in[0], in[1], in[2]};
	Nearest n;
	FindNearest(measurement, afterApogee, excludedCase, n);

	LegacyStep step{};
	step.caseA = static_cast<uint8_t>(n.lowest);
	step.caseB = static_cast<uint8_t>(n.secondLowest);

	// Current point (vector1 / vector2) and the point dt later (vector1Future / vector2Future) of both cases
	float current[2][3], future[2][3];
	const uint8_t pair[2] = {step.caseA, step.caseB};
	for (uint8_t k = 0; k < 2; k++) {
		const ScenarioCase& c = SCENARIOS[pair[k]];
		const ScenarioPoint* rows = afterApogee ? c.after : c.before;
		const uint16_t count = afterApogee ? c.afterCount : c.beforeCount;
		float time = 0.0f;
		Locate(rows, count, n.index[pair[k]], in, current[k], time);
		AtTime(rows, count, time + dt, future[k]);
	}

	// Legacy predictNextValues
	float gainV1[3] = {0, 0, 0}, gainV2[3] = {0, 0, 0};
	LegacyGains(current[0], current[1], in, gainV1, gainV2);
	for (uint8_t i = 0; i < 3; i++) {
		// Legacy divides by the distance to a point, 0 when the state sits exactly on it (legacy
		// caught the NaN later, in the Kalman gain): its "similar" default, 0.5 / 0.5
		if (!std::isfinite(gainV1[i])) {
			gainV1[i] = 0.5f;
			gainV2[i] = 0.5f;
		}
	}

	if (memory.firstTime) {
		for (uint8_t i = 0; i < 3; i++) {
			memory.prevGain1[i] = gainV1[i];
			memory.prevGain2[i] = gainV2[i];
		}
		memory.firstTime = false;
	}

	// Pair changed since this state's last step: blend 50 / 50 this step
	if (resetOnPairChange && (memory.caseA != step.caseA || memory.caseB != step.caseB)) {
		for (uint8_t i = 0; i < 3; i++) {
			memory.prevGain1[i] = 0.5f;
			memory.prevGain2[i] = 0.5f;
		}
		step.gainsReset = true;
	}

	// Interpolate between the two scenarios with LAST step's gains, then keep this step's
	for (uint8_t i = 0; i < 3; i++) {
		step.state[i] = memory.prevGain1[i] * future[0][i] + memory.prevGain2[i] * future[1][i];
		memory.prevGain1[i] = gainV1[i];
		memory.prevGain2[i] = gainV2[i];
	}
	memory.caseA = step.caseA;
	memory.caseB = step.caseB;
	return step;
}

void ScenarioModel::Distances(float altitude, float velocity, float acceleration, float out[SCENARIO_COUNT]) const
{
	Nearest n;
	FindNearest({altitude, velocity, acceleration}, afterApogee, excludedCase, n);
	for (uint8_t s = 0; s < SCENARIO_COUNT; s++) {
		out[s] = n.distances[s];
	}
}
