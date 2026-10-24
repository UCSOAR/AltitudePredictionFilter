/*
 * TareStats.hpp
 *
 *  Created on: Oct 17, 2026
 *      Author: Aly Masani
 */

#ifndef TARESTATS_TARESTATS_HPP_
#define TARESTATS_TARESTATS_HPP_

#include "SensorFrame.hpp"
/************************************
 * MACROS AND DEFINES
 ************************************/
// Tare window: both the time and the sample count must be reached
constexpr uint32_t TARE_MIN_DURATION_MS = 5000;
constexpr uint32_t TARE_MIN_SAMPLES = 20;

// Max standard deviation allowed during the tare window, above this the rocket
// is considered moving and the tare restarts. Tune against pad data
constexpr float TARE_MAX_ACCEL_STD_MPS2 = 0.3f;
constexpr float TARE_MAX_GYRO_STD_RADPS = 0.05f;

/************************************
 * TYPEDEFS
 ************************************/
/**
 * @brief Results of the pad tare, used to calibrate frames, initialize the filter
 *        and as the reference for launch detection.
 *        A sensor is only valid if it delivered data for most of the tare window
 */
struct TareResult
{
	// Barometers
	bool  baroValid[2];
	float groundPressure[2];	// Mean pressure per baro (Pa), each baro keeps its own reference
	float pressureStd[2];		// Pa
	float groundAltitude[2];	// Standard atmosphere altitude of groundPressure (m)
	float altitudeStd[2];		// Pressure noise converted to altitude noise (m)

	// IMUs
	bool  imuValid[2];
	float gravity[2][3];		// Mean accel at rest, sensor axes (m/s^2), includes gravity
	float gravityBody[2][3];	// gravity rotated into rocket body axes (m/s^2)
	float gyroBias[2][3];		// Mean gyro at rest, sensor axes (rad/s)
	float accelStd[2][3];		// m/s^2, sensor axes
	float gyroStd[2][3];		// rad/s, sensor axes
	float padTiltDeg[2];		// Angle between body +Z and measured "up"
	bool  alignmentOk;			// Gravity on the expected body axis for every valid IMU

	uint32_t sampleCount;
	uint32_t durationMs;
};

/************************************
 * CLASS DEFINITIONS
 ************************************/
/**
 * @brief Running mean / variance (Welford's algorithm), no sample buffer needed
 */
class RunningStat
{
	public:
		void Reset() { n = 0; mean = 0.0f; m2 = 0.0f; }
		void Add(float x);
		uint32_t Count() const { return n; }
		float Mean() const { return mean; }
		float Variance() const { return (n > 1) ? (m2 / (n - 1)) : 0.0f; }
		float Std() const;

	private:
		uint32_t n = 0;
		float mean = 0.0f;
		float m2 = 0.0f;
};

/**
 * @brief Accumulates the statistics of every sensor channel over the tare window
 */
class TareStats
{
	public:
		void Reset();
		void Add(const SensorFrame& in);

		uint32_t SampleCount() const { return sampleCount; }
		uint32_t ElapsedMs() const { return windowEndMs - windowStartMs; }
		bool WindowComplete() const;

		// True if no valid IMU axis is noisier than the tare limits
		bool IsStill() const;

		// Fills the statistics fields of TareResult (baro, gravity, bias, noise, validity)
		void FillResult(TareResult& out) const;

	private:
		bool ImuValid(uint8_t imu) const;
		bool BaroValid(uint8_t baro) const;

		RunningStat pressureStat[2];
		RunningStat accelStat[2][3];
		RunningStat gyroStat[2][3];
		uint32_t sampleCount = 0;
		uint32_t windowStartMs = 0;
		uint32_t windowEndMs = 0;
};

#endif /* TARESTATS_TARESTATS_HPP_ */
