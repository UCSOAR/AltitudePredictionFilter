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
 * MACROS AND DEFINES
 ************************************/
// Velocity distance below this counts as an exact match (stops a 1/0 weight)
constexpr float SCENARIO_MATCH_EPS_MPS = 0.5f;

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
	float weightA;			// Share of caseA in the blend (0..1)
	float velocityError;	// m/s, state velocity - caseA velocity at this altitude
};

/************************************
 * CLASS DEFINITIONS
 ************************************/
/**
 * @brief Scenario model for the UKF predict step (replaces the legacy HALO KD-tree).
 *        Every case is looked up at the state's altitude (altitude is monotonic within
 *        a phase: rising before apogee, falling after), the two cases whose velocity
 *        at that altitude is closest to the state are blended by inverse velocity distance
 */
class ScenarioModel
{
	public:
		void SetAfterApogee(bool after) { afterApogee = after; }
		bool IsAfterApogee() const { return afterApogee; }

		// Testing only (leave-one-case-out): skip this case, -1 = use all cases
		void ExcludeCase(int8_t index) { excludedCase = index; }

		ScenarioMatch Match(float altitude, float velocity) const;

	private:
		/**
		 * @brief Velocity and acceleration of one case at an altitude (linear interpolation,
		 *        clamped to the ends of the table)
		 */
		void Lookup(const ScenarioCase& c, float altitude, float& velocity, float& acceleration) const;

		bool afterApogee = false;
		int8_t excludedCase = -1;
};

#endif /* SCENARIOMODEL_SCENARIOMODEL_HPP_ */
