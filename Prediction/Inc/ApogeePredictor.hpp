/*
 * ApogeePredictor.hpp
 *
 *  Created on: Oct 24, 2026
 *      Author: Aly Masani
 */

#ifndef APOGEEPREDICTOR_APOGEEPREDICTOR_HPP_
#define APOGEEPREDICTOR_APOGEEPREDICTOR_HPP_

#include "../../Estimate/Inc/Ukf.hpp"

/************************************
 * MACROS AND DEFINES
 ************************************/
constexpr uint32_t PREDICT_PERIOD_MS = 200;		// Run every 200 ms (5 Hz), not every frame
constexpr float PREDICT_STEP_S = 0.1f;			// Projection step (same as the scenario tables)
constexpr float PREDICT_MAX_S = 60.0f;			// Give up projecting after this long (apogee)
constexpr float PREDICT_LANDING_MAX_S = 900.0f;	// Landing projection limit
constexpr float PREDICT_LANDING_STEP_S = 1.0f;	// Descent is slow, bigger steps
constexpr float PREDICT_SIGMAS = 2.0f;			// Low / high apogee from velocity -/+ this many sigma
constexpr uint8_t PREDICT_TRAJECTORY_POINTS = 30;	// Trajectory: altitude every second ahead

/************************************
 * TYPEDEFS
 ************************************/
struct ApogeePrediction
{
	bool valid;
	bool afterApogee;				// true: apogee passed, the landing is predicted instead
	uint32_t madeAtMs;				// Frame time the prediction was made

	// Before apogee
	float apogeeAltitude;			// m above the pad
	float apogeeAltitudeLow;		// velocity - 2 sigma
	float apogeeAltitudeHigh;		// velocity + 2 sigma
	float timeToApogee;				// s from madeAtMs

	// After apogee
	float timeToLanding;			// s from madeAtMs
	float landingVelocity;			// m/s at touchdown

	// Trajectory predictor: altitude 1 s, 2 s, ... ahead
	uint8_t trajectoryCount;
	float trajectory[PREDICT_TRAJECTORY_POINTS];
};

/************************************
 * CLASS DEFINITIONS
 ************************************/
/**
 * @brief PREDICTION stage: apogee + trajectory predictor (legacy HALO::predictNStates idea).
 *        Projects the UKF state forward with the same scenario model the UKF uses, on a copy,
 *        so the filter itself is never touched:
 *          before apogee: until the velocity reaches 0 -> predicted apogee height and time
 *          after apogee:  until the altitude reaches 0 -> predicted landing time and speed
 */
class ApogeePredictor
{
	public:
		void Reset() { prediction = ApogeePrediction{}; lastRunMs = 0; runs = 0; }

		/**
		 * @brief Runs at most every PREDICT_PERIOD_MS
		 * @return true when a new prediction was made this frame
		 */
		bool Update(const Ukf& ukf, uint32_t timestampMs);

		const ApogeePrediction& Get() const { return prediction; }
		void PrintStatus() const;

	private:
		/**
		 * @brief Steps alt / vel forward with the scenario acceleration until the velocity reaches 0
		 * @return apogee altitude, timeToApogee set, trajectory filled if not null
		 */
		float ProjectApogee(const ScenarioModel& model, float alt, float vel, float& timeToApogee,
							float* trajectory, uint8_t* trajectoryCount) const;
		void ProjectLanding(const ScenarioModel& model, float alt, float vel);

		ApogeePrediction prediction{};
		uint32_t lastRunMs = 0;
		uint32_t runs = 0;
};

#endif /* APOGEEPREDICTOR_APOGEEPREDICTOR_HPP_ */
