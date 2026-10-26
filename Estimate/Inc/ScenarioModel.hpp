/*
 * ScenarioModel.hpp
 *
 *  Created on: Oct 24, 2026
 *      Author: Aly Masani
 */

#ifndef SCENARIOMODEL_SCENARIOMODEL_HPP_
#define SCENARIOMODEL_SCENARIOMODEL_HPP_

#include "ScenarioData.hpp"

/************************************
 * TYPEDEFS
 ************************************/
/**
 * @brief Result of matching a state against the simulated flights
 */
struct ScenarioMatch
{
	float acceleration;		// m/s^2, blended from the two closest cases
	uint8_t caseA;			// Closest case (index into SCENARIOS)
	uint8_t caseB;			// Second closest case
	float weightA;			// Share of caseA in the blend (0..1), legacy gain on the acceleration axis
	float distanceA;		// Euclidean (alt, vel, acc) distance to caseA's nearest point
	float distanceB;		// ... and to caseB's
	uint16_t pointA;		// Index of that nearest point in caseA's table
	uint16_t pointB;
	float velocityError;	// m/s, state velocity - caseA nearest point velocity
};

/**
 * @brief Legacy gain memory of one state being stepped (one per UKF sigma point, one for a
 *        prediction run). Legacy HALO listOfGainsSigmaPoints / firstTime / scenariosGainsList.
 *        Legacy kept the pair in a 6 slot list indexed by counter % 6 with 7 sigma points, so the
 *        slots drifted and points compared against another point's pair. Here every sigma point
 *        keeps its own pair (the fix)
 */
struct LegacyGainMemory
{
	float prevGain1[3];		// Gains of the last step, used for THIS step's blend (legacy one step lag)
	float prevGain2[3];
	uint8_t caseA;			// Pair used last step (legacy scenariosGainsList, starts {0, 0})
	uint8_t caseB;
	bool firstTime;			// First step: use this step's own gains

	void Reset()
	{
		for (uint8_t i = 0; i < 3; i++) {
			prevGain1[i] = 0.5f;
			prevGain2[i] = 0.5f;
		}
		caseA = 0;
		caseB = 0;
		firstTime = true;
	}
};

/**
 * @brief One legacy model step
 */
struct LegacyStep
{
	float state[3];			// Predicted altitude, velocity, acceleration
	uint8_t caseA;			// The two closest cases this step
	uint8_t caseB;
	bool gainsReset;		// Pair changed: blended 50 / 50 this step
};

/************************************
 * CLASS DEFINITIONS
 ************************************/
/**
 * @brief Scenario model for the UKF predict step, the legacy HALO KD-tree search
 *        (Legacy/HALO.cpp findNearestScenarios + Legacy/KDTree.cpp, unchanged):
 *          1. every case has a KD-tree of its (altitude, velocity, acceleration) points,
 *             one for before apogee, one for after
 *          2. per case: nearest point to the state (raw units, squared Euclidean)
 *          3. the two cases with the smallest Euclidean distance are kept
 *          4. legacy predictNextValues gains (position of the state between the two points,
 *             per axis) blend their accelerations
 *
 *        NOTE: the KD-trees use std::shared_ptr nodes on the heap (as legacy), one node per
 *        table row, built on the first Match(). Check the heap before flying this
 */
class ScenarioModel
{
	public:
		void SetAfterApogee(bool after) { afterApogee = after; }
		bool IsAfterApogee() const { return afterApogee; }

		// Testing only (leave-one-case-out): skip this case, -1 = use all cases
		void ExcludeCase(int8_t index) { excludedCase = index; }

		ScenarioMatch Match(float altitude, float velocity, float acceleration) const;

		/**
		 * @brief Legacy HALO dynamic model (findNearestScenarios + predictNextValues): the state dt
		 *        from now is the blend of the two closest scenarios dt after their matched point.
		 *        Rate independent: the matched point is refined onto the scenario's line between
		 *        rows and the future point is interpolated at +dt (legacy took the next table row at
		 *        or after +1/3 s, which only works at its 3 Hz)
		 * @param resetOnPairChange legacy resets to 50 / 50 in the filter, not in predictNStates
		 */
		LegacyStep LegacyPredict(const float in[3], float dt, LegacyGainMemory& memory, bool resetOnPairChange) const;

		// Distance from a state to every case's nearest point (what the cases are ranked by), for test plots
		void Distances(float altitude, float velocity, float acceleration, float out[SCENARIO_COUNT]) const;

	private:
		bool afterApogee = false;
		int8_t excludedCase = -1;
};

#endif /* SCENARIOMODEL_SCENARIOMODEL_HPP_ */
