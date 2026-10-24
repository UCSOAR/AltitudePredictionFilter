/*
 * Dynamite.hpp
 *
 *  Created on: Oct 24, 2026
 *      Author: Aly Masani
 */

#ifndef DYNAMITE_DYNAMITE_HPP_
#define DYNAMITE_DYNAMITE_HPP_

#include "TareStats.hpp"
#include "Madgwick.hpp"

/************************************
 * MACROS AND DEFINES
 ************************************/
// All values are starting points, tune them against recorded pad / flight data

// Altitude sources fused with trust weights
enum DynamiteSource : uint8_t
{
	DYN_SRC_BARO0 = 0,
	DYN_SRC_BARO1,
	DYN_SRC_IMU0,
	DYN_SRC_IMU1,
	DYN_SRC_COUNT
};

// Trust weight = 1 / (|predicted - measured| + sigma + DYNAMITE_TRUST_EPS_M).
// Larger eps flattens the weights (legacy used 100 m, which made them almost equal),
// smaller eps punishes a disagreeing sensor harder
constexpr float DYNAMITE_TRUST_EPS_M = 1.0f;

// Smallest altitude noise used for any source, stops a near zero sigma from taking all the trust
constexpr float DYNAMITE_MIN_SIGMA_M = 0.01f;

// IMU motion gate (pad only): net vertical accel above this many accel sigmas means moving.
// Below it the IMU altitude is held at 0, on the pad the IMU knows we have not moved
constexpr float DYNAMITE_ACC_MOTION_SIGMAS = 5.0f;

// IMU clipping, ranges as configured in the drivers (keep in sync):
//   IMU0 LSM6DSO   +/-16 g, IMU1 LSM6DSO32 +/-32 g, both gyros +/-2000 dps (data in 0.1 dps)
// Readings this close to a limit are clipped
constexpr float DYNAMITE_ACCEL_RANGE_G[2] = {16.0f, 32.0f};
constexpr float DYNAMITE_GYRO_RANGE_DPS = 2000.0f;
constexpr float DYNAMITE_CLIP_FRACTION = 0.97f;
constexpr float DYNAMITE_GYRO_CLIP_RADPS = DYNAMITE_CLIP_FRACTION * DYNAMITE_GYRO_RANGE_DPS * 3.14159265f / 180.0f;

// In flight, an IMU reading this much less total force than the other one is treated as clipped
// (clipping only ever makes a reading smaller, e.g. an IMU left on a lower range)
constexpr float DYNAMITE_IMU_DISAGREE_MPS2 = 5.0f;

// In flight the accel noise is motor / airframe vibration, far above the pad tare noise.
// Used as the floor of the IMU accel std once in flight
constexpr float DYNAMITE_FLIGHT_ACC_STD_MPS2 = 0.5f;

// Baro distrust: sigma multiplier while the baros are known to be wrong
//  - transonic: shock waves over the static ports (Mach ~0.7 - 1.2 at the pad altitude)
//  - ejection:  charge pressure spike, a few seconds after apogee (DistrustBarosUntil)
constexpr float DYNAMITE_TRANSONIC_LOW_MPS = 240.0f;
constexpr float DYNAMITE_TRANSONIC_HIGH_MPS = 410.0f;
constexpr float DYNAMITE_BARO_DISTRUST_SCALE = 20.0f;
constexpr uint32_t DYNAMITE_EJECTION_DISTRUST_MS = 3000;	// After apogee detection
constexpr float DYNAMITE_EJECTION_GATE_M = 10.0f;			// In that window, a baro this far from the prediction is a spike

// Velocity smoothing (alpha-beta tracker on the fused altitude, IMU accel as feed forward).
// Each frame the velocity is corrected by DYNAMITE_VEL_BETA * (fused - predicted) / dt.
// Smaller = smoother but slower to correct, larger = follows the baro noise. Tune at the flight poll rate
constexpr float DYNAMITE_VEL_BETA = 0.01f;

/************************************
 * TYPEDEFS
 ************************************/
/**
 * @brief Everest measurement handed to Launch Detection and the UKF.
 *        Altitude is above the pad, up positive
 */
struct EverestMeasurement
{
	bool  valid;
	float altitude;			// m, trust weighted
	float velocity;			// m/s, smoothed (alpha-beta tracker on the fused altitude)
	float acceleration;		// m/s^2, net vertical (earth frame, gravity removed), IMU mean
	float sigma;			// m, noise of the fused altitude
};

/**
 * @brief Latest values of every source, kept for debug prints
 */
struct DynamiteStatus
{
	bool  sourceValid[DYN_SRC_COUNT];
	float sourceAltitude[DYN_SRC_COUNT];	// m above the pad
	float sourceSigma[DYN_SRC_COUNT];		// m
	float weight[DYN_SRC_COUNT];			// Normalised trust weights, sum to 1
	float predicted;						// m

	// Per IMU
	float imuAccel[2];		// Net vertical accel, m/s^2
	float imuVelocity[2];	// m/s
	bool  imuMoving[2];
	float tiltDeg[2];		// Madgwick tilt from the pad attitude
	bool  accelUsed[2];		// Madgwick accel correction active (false under thrust / coast)
	bool  accelClipped[2];	// Accel at the range limit, the other IMU / last fused accel is used instead
	bool  gyroClipped[2];	// Gyro at the data limit, the attitude can no longer be trusted
	uint32_t clipFrames;	// Frames with any clipped IMU since Init

	bool  baroDistrusted;	// Baro sigma scaled up this frame (transonic or ejection)
};

/************************************
 * CLASS DEFINITIONS
 ************************************/
/**
 * @brief PROCESS stage: Madgwick attitude per IMU + baro altitude, blended with dynamite
 *        trust weights (legacy EverestTask dynamite() / recalculateGain()). Every frame:
 *          1. Madgwick per IMU -> earth frame vertical accel -> IMU altitude (double integration)
 *          2. Baro altitude above the pad, each against its own tare reference
 *          3. Predict the altitude from the last fused state
 *          4. Trust weight per source = 1 / (|predicted - measured| + sigma + eps), normalise
 *          5. Fused Everest measurement (alt, vel, acc, sigma)
 */
class Dynamite
{
	public:
		// Start from the finished tare (Madgwick initial attitude from the tare gravity)
		void Init(const TareResult& tare);
		void Reset();

		/**
		 * @param calibrated Calibrated frame (gyro bias removed, body axes)
		 * @param tare       Finished tare
		 * @param inFlight   true once launch is detected: the pad zero velocity hold stops
		 *                   and the IMUs integrate freely
		 */
		const EverestMeasurement& Update(const SensorFrame& calibrated, const TareResult& tare, bool inFlight);

		// Distrust the baros until this frame time (ejection charges after apogee)
		void DistrustBarosUntil(uint32_t timestampMs) { baroDistrustUntilMs = timestampMs; }

		bool IsInitialized() const { return initialized; }
		const EverestMeasurement& GetMeasurement() const { return measurement; }
		const DynamiteStatus& GetStatus() const { return status; }

		void PrintStatus() const;
		void PrintLine() const;		// One line: every source, predicted, fused

	private:
		void MeasureBaros(const SensorFrame& frame, const TareResult& tare);
		void MeasureImus(const SensorFrame& frame, const TareResult& tare, bool inFlight);
		void Fuse(float dt);

		Madgwick madgwick[2];

		// IMU integration
		float imuAltitude[2] = {0.0f, 0.0f};
		float imuPrevAccel[2] = {0.0f, 0.0f};	// Last frame's vertical accel, for the trapezoid rule
		float imuMovingTime[2] = {0.0f, 0.0f};	// s since the integration started
		float accelStd[2] = {0.0f, 0.0f};		// Largest axis accel noise from the tare

		uint32_t baroDistrustUntilMs = 0;

		bool initialized = false;
		bool hasFused = false;

		EverestMeasurement measurement{};
		DynamiteStatus status{};
};

#endif /* DYNAMITE_DYNAMITE_HPP_ */
