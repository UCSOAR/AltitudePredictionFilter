/*
 * ApogeeDetection.hpp
 *
 *  Created on: Oct 24, 2026
 *      Author: Aly Masani
 */

#ifndef APOGEEDETECTION_APOGEEDETECTION_HPP_
#define APOGEEDETECTION_APOGEEDETECTION_HPP_

#include <stdint.h>

/************************************
 * MACROS AND DEFINES
 ************************************/
// Legacy HALO::apogeeDetection() settings
constexpr uint8_t APOGEE_WINDOW = 10;				// Frames averaged ("window of 10 is best", legacy)
constexpr float APOGEE_ALT_WEIGHT = 0.5f;
constexpr float APOGEE_VEL_WEIGHT = 0.7f;
constexpr float APOGEE_ACC_WEIGHT = 0.4f;
constexpr float APOGEE_TARGET_ACC_MPS2 = -9.81f;	// Free fall at apogee
constexpr float APOGEE_DIP_MIN_PEAK = 0.7f;			// Trigger if the peak was above this...
constexpr float APOGEE_DIP_RATIO = 0.8f;			// ...and the confidence fell to this share of it

// Smallest sigmas used in the confidences (alt m, vel m/s, acc m/s^2). Legacy's P was large; ours
// is tight (vel sigma ~0.2 m/s), which made the "velocity near zero" window about one frame long
// and the confidence could miss its peak. Found with the replay test (Tools/replay)
constexpr float APOGEE_CONF_MIN_SD[3] = {1.0f, 2.0f, 1.0f};

// Ignore apogee below this altitude above the pad (launch transients, bench tests)
constexpr float APOGEE_MIN_ALT_M = 10.0f;

// Backup (not in legacy): velocity negative for this many frames in a row also means apogee,
// in case the confidence never peaks the way legacy expects (~0.5 s at 40 Hz)
constexpr uint8_t APOGEE_BACKUP_FRAMES = 20;

/************************************
 * TYPEDEFS
 ************************************/
enum class ApogeeSource : uint8_t
{
	NONE = 0,
	CONFIDENCE,		// Legacy confidence method
	BACKUP			// Velocity negative for APOGEE_BACKUP_FRAMES
};

struct ApogeeStatus
{
	bool reached;
	ApogeeSource source;
	float maxAltitude;		// m above the pad, highest estimate so far
	uint32_t maxTimeMs;		// Frame time of maxAltitude
	uint32_t detectTimeMs;	// Frame time apogee was confirmed

	// Legacy confidences (last frame)
	float altConfidence;
	float velConfidence;
	float accConfidence;
	float totalConfidence;
	float maxConfidence;
	bool hitOne;
	uint8_t negativeVelFrames;
};

/************************************
 * CLASS DEFINITIONS
 ************************************/
/**
 * @brief PREDICTION stage, apogee detection. Port of the legacy HALO::apogeeDetection():
 *        a 10 frame window of the UKF state is scored against the UKF uncertainty
 *          alt: how far it has dropped from its highest point, vs the altitude sigma
 *          vel: how close to zero (or below) the velocity is, vs the velocity sigma
 *          acc: how close the acceleration is to free fall (-9.81)
 *        total = 0.5 alt + 0.7 vel + 0.4 acc. Apogee when the total has reached 1 and starts
 *        falling, or peaked above 0.7 and fell to 80% of the peak.
 *        On apogee the scenario model switches to the post apogee tables (once)
 */
class ApogeeDetection
{
	public:
		void Reset();

		/**
		 * @param x   UKF state: altitude (m above pad), velocity (m/s), acceleration (m/s^2)
		 * @param sd  UKF standard deviations of x
		 * @param timestampMs Frame time
		 * @return true on the frame apogee is confirmed (only once)
		 */
		bool Update(const float x[3], const float sd[3], uint32_t timestampMs);

		bool IsReached() const { return status.reached; }
		const ApogeeStatus& GetStatus() const { return status; }
		void PrintStatus() const;

	private:
		bool Confidence(float avgAlt, float avgVel, float avgAcc, const float sd[3]);
		void Latch(ApogeeSource source, uint32_t timestampMs);

		// Moving window (legacy buffer)
		float windowAlt[APOGEE_WINDOW] = {};
		float windowVel[APOGEE_WINDOW] = {};
		float windowAcc[APOGEE_WINDOW] = {};
		uint8_t head = 0;
		uint8_t count = 0;

		// Legacy calculateConfidence() state
		float maxAvgAltitude = -1e30f;
		float altitudeAccumulator = 0.0f;

		ApogeeStatus status{};
};

#endif /* APOGEEDETECTION_APOGEEDETECTION_HPP_ */
