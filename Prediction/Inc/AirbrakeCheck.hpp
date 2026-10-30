/*
 * AirbrakeCheck.hpp
 *
 *  Created on: Oct 30, 2026
 *      Author: Aly Masani
 */

#ifndef AIRBRAKECHECK_AIRBRAKECHECK_HPP_
#define AIRBRAKECHECK_AIRBRAKECHECK_HPP_

#include "ApogeePredictor.hpp"
#include <stdint.h>

/************************************
 * MACROS AND DEFINES
 ************************************/
// LOG ONLY: nothing here moves the airbrakes. The flag, level and confidence are logged so the
// check can be compared with the real apogee after the flight
constexpr float AIRBRAKE_TARGET_APOGEE_M = 3048.0f;	// m above the pad (10,000 ft), set to the flight's target
constexpr float AIRBRAKE_M_PER_LEVEL = 100.0f;		// Overshoot by 300 m -> level 3
constexpr uint8_t AIRBRAKE_MAX_LEVEL = 5;
constexpr uint32_t AIRBRAKE_CHECK_PERIOD_MS = 2000;	// Check every 2 s
constexpr uint32_t AIRBRAKE_MAX_PREDICTION_AGE_MS = 500;	// Only use a fresh prediction

/************************************
 * TYPEDEFS
 ************************************/
struct AirbrakeStatus
{
	bool flagged;			// Latched: a level >= 1 was flagged, no more checks
	uint8_t level;			// 0 (no airbrakes) .. AIRBRAKE_MAX_LEVEL, last check
	float confidence;		// 0..1, last check: chance the apogee really is above this level's overshoot
	float overshoot;		// m, predicted apogee - target, last check
	uint32_t flaggedMs;		// Frame time of the flag
	uint16_t checks;		// Checks run
};

/************************************
 * CLASS DEFINITIONS
 ************************************/
/**
 * @brief PREDICTION stage, airbrake level check (log only). After burnout and before apogee, every
 *        AIRBRAKE_CHECK_PERIOD_MS: overshoot = predicted apogee - target, level = overshoot / 100 m
 *        (0..5). Confidence = P(apogee > target + level * 100 m), from the predictor's low / high
 *        band (+-2 sigma). The first level >= 1 is latched and the check stops
 */
class AirbrakeCheck
{
	public:
		void Reset() { status = AirbrakeStatus{}; lastCheckMs = 0; }

		/**
		 * @param prediction  ApogeePredictor::Get()
		 * @param burnedOut   UKF acceleration below 0 (motor out)
		 * @param moved       ApogeeDetection::HasMoved() (the board really left the pad)
		 * @param afterApogee Apogee detected
		 * @return true when a check ran this frame (log it)
		 */
		bool Update(const ApogeePrediction& prediction, bool burnedOut, bool moved, bool afterApogee, uint32_t timestampMs);

		const AirbrakeStatus& GetStatus() const { return status; }

		// Target apogee, m above the pad (AIRBRAKE_TARGET_APOGEE_M unless changed, e.g. replay tests)
		void SetTarget(float targetM) { target = targetM; }
		float GetTarget() const { return target; }

		// Event value: level (bits 0-7), confidence /255 (bits 8-15), overshoot m int16 (bits 16-31)
		int32_t EventValue() const;

		void PrintStatus() const;

	private:
		AirbrakeStatus status{};
		uint32_t lastCheckMs = 0;
		float target = AIRBRAKE_TARGET_APOGEE_M;
};

#endif /* AIRBRAKECHECK_AIRBRAKECHECK_HPP_ */
