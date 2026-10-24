/*
 * LaunchDetection.hpp
 *
 *  Created on: Oct 17, 2026
 *      Author: Aly Masani
 */

#ifndef LAUNCHDETECTION_LAUNCHDETECTION_HPP_
#define LAUNCHDETECTION_LAUNCHDETECTION_HPP_

#include "../../Process/Inc/Dynamite.hpp"
/************************************
 * MACROS AND DEFINES
 ************************************/
// All values are starting points, tune them against recorded pad / flight data

// Launch when the fused altitude is above max(LAUNCH_ALT_SIGMAS * fused sigma, LAUNCH_ALT_MIN_M)
// for LAUNCH_CONFIRM_FRAMES frames in a row. The floor stops slow pad pressure drift from triggering
constexpr float LAUNCH_ALT_SIGMAS = 5.0f;
constexpr float LAUNCH_ALT_MIN_M = 1.0f;
constexpr uint8_t LAUNCH_CONFIRM_FRAMES = 3;

/************************************
 * TYPEDEFS
 ************************************/
/**
 * @brief Latest detector values, kept for debug prints
 */
struct LaunchDetectionStatus
{
	float altitude;			// m above the pad, from the Everest measurement
	float threshold;		// m, launch threshold this frame
	uint8_t confirmCount;
	bool launched;
	uint32_t launchTimeMs;
};

/************************************
 * CLASS DEFINITIONS
 ************************************/
/**
 * @brief Internal launch detection of the altitude prediction filter, from sensor data only
 *        (never from the rocket state). Launch once the dynamite fused altitude stays
 *        above the ground noise. Once launched it stays launched until Reset()
 */
class LaunchDetection
{
	public:
		void Reset();

		/**
		 * @param everest     Fused measurement from Dynamite
		 * @param timestampMs Frame time, stamps the launch
		 * @return true once launch is detected (latched)
		 */
		bool Update(const EverestMeasurement& everest, uint32_t timestampMs);

		bool IsLaunched() const { return status.launched; }
		const LaunchDetectionStatus& GetStatus() const { return status; }
		void PrintStatus() const;

	private:
		void Latch(uint32_t timestampMs);

		LaunchDetectionStatus status{};
};

#endif /* LAUNCHDETECTION_LAUNCHDETECTION_HPP_ */
