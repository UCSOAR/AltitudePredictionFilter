/*
 * StateSelection.cpp
 *
 *  Created on: Oct 24, 2026
 *      Author: Aly Masani
 */

/**
 ********************************************************************************
 * @file    StateSelection.cpp
 * @author  Aly Masani
 * @date    Oct 24, 2026
 * @brief   Internal state machine of the altitude prediction filter
 ********************************************************************************
 */

/************************************
 * INCLUDES
 ************************************/
#include "Inc/StateSelection.hpp"
#include "SystemDefines.hpp"

/************************************
 * FUNCTION DEFINITIONS
 ************************************/
void StateSelection::Update(bool tared, bool launched, bool rocketLanded, uint32_t timestampMs)
{
	switch (state)
	{
	case ApfState::PRE_LAUNCH:
		if (tared) {
			Transition(ApfState::DETECTION, timestampMs);
		}
		break;

	case ApfState::DETECTION:
		if (!tared) {
			// Tare was thrown away (rocket left a pad state before the tare finished)
			Transition(ApfState::PRE_LAUNCH, timestampMs);
		} else if (launched) {
			Transition(ApfState::LAUNCH, timestampMs);
		}
		break;

	case ApfState::LAUNCH:
		if (rocketLanded) {
			Transition(ApfState::LANDED, timestampMs);
		}
		break;

	case ApfState::LANDED:
	default:
		break;
	}
}

void StateSelection::Reset(uint32_t timestampMs)
{
	if (state != ApfState::PRE_LAUNCH) {
		Transition(ApfState::PRE_LAUNCH, timestampMs);
	}
}

const char* StateSelection::Name(ApfState s)
{
	switch (s)
	{
	case ApfState::PRE_LAUNCH: return "PRE_LAUNCH";
	case ApfState::DETECTION:  return "DETECTION";
	case ApfState::LAUNCH:     return "LAUNCH";
	case ApfState::LANDED:     return "LANDED";
	default:                   return "UNKNOWN";
	}
}

void StateSelection::PrintStatus(uint32_t nowMs) const
{
	SOAR_PRINT("APF state: %s for %lums\n", Name(state),
			   static_cast<unsigned long>(nowMs - stateEntryMs));
}

void StateSelection::Transition(ApfState next, uint32_t timestampMs)
{
	SOAR_PRINT("APF state: %s -> %s at t=%lums\n", Name(state), Name(next),
			   static_cast<unsigned long>(timestampMs));
	state = next;
	stateEntryMs = timestampMs;
}
