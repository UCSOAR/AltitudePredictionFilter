/*
 * FilterOutput.hpp
 *
 *  Created on: Oct 24, 2026
 *      Author: Aly Masani
 */

#ifndef FILTEROUTPUT_FILTEROUTPUT_HPP_
#define FILTEROUTPUT_FILTEROUTPUT_HPP_

#include "SensorDataTypes.hpp"
#include "StateSelection.hpp"
#include "LaunchDetection.hpp"
#include "../../Process/Inc/Dynamite.hpp"
#include "../../Estimate/Inc/Ukf.hpp"
#include "../../Prediction/Inc/ApogeeDetection.hpp"
#include "../../Prediction/Inc/ApogeePredictor.hpp"

/************************************
 * MACROS AND DEFINES
 ************************************/
// Publish / log rate: 10 Hz in flight, once a second otherwise (keeps the pad log small).
// Every frame (50 Hz) in flight overloaded LoggingTask on the bench and stalled polling
constexpr uint32_t FILTER_OUTPUT_FLIGHT_PERIOD_MS = 100;
constexpr uint32_t FILTER_OUTPUT_SLOW_PERIOD_MS = 1000;

// FilterData::flags
enum FilterFlags : uint8_t
{
	FILTER_FLAG_LAUNCHED          = (1U << 0),	// Launch detected
	FILTER_FLAG_APOGEE            = (1U << 1),	// Apogee detected
	FILTER_FLAG_PREDICTION        = (1U << 2),	// predictedApogee / timeToEvent are valid
	FILTER_FLAG_UKF               = (1U << 3),	// Estimate from the UKF (otherwise Dynamite)
	FILTER_FLAG_IMU_CLIPPED       = (1U << 4),	// An IMU clipped this frame
	FILTER_FLAG_BARO_DISTRUSTED   = (1U << 5),	// Baros distrusted this frame (transonic / ejection)
	FILTER_FLAG_MODEL_MISFIT      = (1U << 6),	// UKF inflated P this frame (flight unlike the sims)
	FILTER_FLAG_SYNTHETIC_DESCENT = (1U << 7)	// Post apogee scenarios are synthetic (placeholder data)
};

/************************************
 * CLASS DEFINITIONS
 ************************************/
/**
 * @brief OUTPUT stage: turns the filter state into FilterData for the DataBroker / flash log.
 *        The publish itself stays in SensorFrameTask (needs the RTOS), so this is PC testable
 */
class FilterOutput
{
	public:
		void Reset() { lastMs = 0; published = 0; hasLast = false; }

		// true when FilterData should be published this frame
		bool Due(ApfState state, uint32_t timestampMs);

		FilterData Build(uint32_t timestampMs, ApfState state, float groundAltitude,
						 const Dynamite& dynamite, const LaunchDetection& launch, const Ukf& ukf,
						 const ApogeeDetection& apogee, const ApogeePredictor& predictor);

		// Remember what was published (for APF_Output)
		void Published(const FilterData& data) { last = data; hasLast = true; published++; }

		void PrintLast() const;
		static void PrintLine(const FilterData& data);

	private:
		uint32_t lastMs = 0;
		uint32_t published = 0;
		bool hasLast = false;
		FilterData last{};
};

#endif /* FILTEROUTPUT_FILTEROUTPUT_HPP_ */
