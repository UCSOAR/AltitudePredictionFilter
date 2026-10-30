/*
 * StateSelection.hpp
 *
 *  Created on: Oct 24, 2026
 *      Author: Aly Masani
 */

#ifndef STATESELECTION_STATESELECTION_HPP_
#define STATESELECTION_STATESELECTION_HPP_

#include <stdint.h>

/************************************
 * TYPEDEFS
 ************************************/
/**
 * @brief Internal state of the altitude prediction filter (state transitions diagram).
 *        Separate from the rocket state, which only gates the tare (pad) and landing
 */
enum class ApfState : uint8_t
{
	PRE_LAUNCH = 0,		// Taring on the pad
	DETECTION,			// Tared, watching for launch
	LAUNCH,				// Launch detected: Process -> Estimate -> Prediction -> Output every frame
	LANDED				// Rocket state TOUCHDOWN / RECOVERY: hold the last estimate
};

/************************************
 * CLASS DEFINITIONS
 ************************************/
/**
 * @brief Decides which stages run for every SensorFrame.
 *          PRE_LAUNCH -> DETECTION  tare done
 *          DETECTION  -> LAUNCH     launch detected
 *          LAUNCH     -> LANDED     rocket state TOUCHDOWN / RECOVERY
 *          any        -> PRE_LAUNCH re-tare (only accepted while the rocket state is a pad state)
 */
class StateSelection
{
	public:
		/**
		 * @param tared          PreProcess has a finished tare
		 * @param launched       LaunchDetection has latched
		 * @param rocketLanded   Rocket state is TOUCHDOWN or RECOVERY
		 * @param timestampMs    Frame time, used to stamp the transition
		 */
		void Update(bool tared, bool launched, bool rocketLanded, uint32_t timestampMs);

		// Back to PRE_LAUNCH, called together with the re-tare
		void Reset(uint32_t timestampMs);

		// Power loss: start in the saved state (LAUNCH or LANDED)
		void Restore(ApfState saved, uint32_t timestampMs) { Transition(saved, timestampMs); }

		ApfState GetState() const { return state; }
		uint32_t GetStateEntryMs() const { return stateEntryMs; }
		static const char* Name(ApfState s);

		void PrintStatus(uint32_t nowMs) const;

	private:
		void Transition(ApfState next, uint32_t timestampMs);

		ApfState state = ApfState::PRE_LAUNCH;
		uint32_t stateEntryMs = 0;
};

#endif /* STATESELECTION_STATESELECTION_HPP_ */
