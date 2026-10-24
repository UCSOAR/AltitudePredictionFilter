/*
 * LaunchDetection.hpp
 *
 *  Created on: Oct 17, 2026
 *      Author: Aly Masani
 */

#ifndef LAUNCHDETECTION_LAUNCHDETECTION_HPP_
#define LAUNCHDETECTION_LAUNCHDETECTION_HPP_

#include "TareStats.hpp"
/************************************
 * MACROS AND DEFINES
 ************************************/
// All values are starting points, tune them against recorded pad / flight data

// Altitude sources fused with trust weights (same idea as the legacy dynamite())
enum LaunchAltSource : uint8_t
{
	LAUNCH_SRC_BARO0 = 0,
	LAUNCH_SRC_BARO1,
	LAUNCH_SRC_IMU,
	LAUNCH_SRC_COUNT
};

// Trust weight = 1 / (|predicted - measured| + sigma + LAUNCH_TRUST_EPS_M).
// Larger eps flattens the weights (legacy used 100 m, which made them almost equal),
// smaller eps punishes a disagreeing sensor harder
constexpr float LAUNCH_TRUST_EPS_M = 1.0f;

// Smallest altitude noise used for any source, stops a near zero sigma from taking all the trust
constexpr float LAUNCH_MIN_SIGMA_M = 0.01f;

// IMU motion gate: net vertical accel above this many accel sigmas means the rocket is moving.
// Below it the IMU altitude is held at 0 (on the pad the IMU knows we have not moved)
constexpr float LAUNCH_ACC_MOTION_SIGMAS = 5.0f;

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
	// Per source, indexed by LaunchAltSource
	bool  sourceValid[LAUNCH_SRC_COUNT];
	float sourceAltitude[LAUNCH_SRC_COUNT];		// m above the pad
	float sourceSigma[LAUNCH_SRC_COUNT];		// m, expected noise of that source
	float weight[LAUNCH_SRC_COUNT];				// Normalised trust weights, sum to 1

	// Fused
	float predicted;		// m, fused altitude extrapolated from the last frame
	float altitude;			// m above the pad, trust weighted
	float velocity;			// m/s, up positive
	float sigma;			// m, noise of the fused altitude
	float threshold;		// m, launch threshold this frame

	// IMU channel
	float imuAccel;			// Net vertical accel (gravity removed), m/s^2
	float imuVelocity;		// m/s, integrated only while moving
	bool  imuMoving;

	uint8_t confirmCount;
	bool launched;
	uint32_t launchTimeMs;
};

/************************************
 * CLASS DEFINITIONS
 ************************************/
/**
 * @brief Internal launch detection of the altitude prediction filter, from sensor data only
 *        (never from the rocket state). Every frame:
 *          1. Altitude above the pad from each source (BARO0, BARO1, IMU double integration)
 *          2. Trust weight per source from how far it is from the predicted altitude (dynamite)
 *          3. Trust weighted fused altitude
 *          4. Launch once the fused altitude stays above the ground noise
 *        Once launched it stays launched until Reset()
 */
class LaunchDetection
{
	public:
		void Reset();

		/**
		 * @param calibrated Calibrated frame from PreProcess (IMUs in body axes)
		 * @param tare       Finished tare (ground references and noise)
		 * @return true once launch is detected (latched)
		 */
		bool Update(const SensorFrame& calibrated, const TareResult& tare);

		bool IsLaunched() const { return status.launched; }
		const LaunchDetectionStatus& GetStatus() const { return status; }
		void PrintStatus() const;
		void PrintAltitudes() const;	// One line: every source, predicted, fused, launch flag

	private:
		void MeasureBaros(const SensorFrame& frame, const TareResult& tare);
		void MeasureImu(const SensorFrame& frame, const TareResult& tare);
		void FuseAltitude(float dt);
		void Latch(uint32_t timestampMs);

		// IMU integration, only runs while the motion gate is open
		float imuAltitude = 0.0f;
		float imuMovingTime = 0.0f;		// s since the motion gate opened

		bool hasFused = false;			// false until the first fused altitude exists

		LaunchDetectionStatus status{};
};

#endif /* LAUNCHDETECTION_LAUNCHDETECTION_HPP_ */
