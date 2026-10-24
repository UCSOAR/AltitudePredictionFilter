/*
 * TareStats.cpp
 *
 *  Created on: Oct 17, 2026
 *      Author: Aly Masani
 */

/**
 ********************************************************************************
 * @file    TareStats.cpp
 * @author  Aly Masani
 * @date    Oct 17, 2026
 * @brief   Tare window statistics: mean and noise of every IMU axis and baro
 ********************************************************************************
 */

/************************************
 * INCLUDES
 ************************************/
#include "Inc/TareStats.hpp"
#include <cmath>

/************************************
 * FUNCTION DEFINITIONS
 ************************************/
void RunningStat::Add(float x)
{
	n++;
	const float delta = x - mean;
	mean += delta / n;
	m2 += delta * (x - mean);
}

float RunningStat::Std() const
{
	return std::sqrt(Variance());
}

void TareStats::Reset()
{
	for (uint8_t i = 0; i < 2; i++) {
		pressureStat[i].Reset();
		for (uint8_t axis = 0; axis < 3; axis++) {
			accelStat[i][axis].Reset();
			gyroStat[i][axis].Reset();
		}
	}

	sampleCount = 0;
	windowStartMs = 0;
	windowEndMs = 0;
}

void TareStats::Add(const SensorFrame& in)
{
	if (sampleCount == 0) {
		windowStartMs = in.timestampMs;
	}
	windowEndMs = in.timestampMs;

	for (uint8_t i = 0; i < 2; i++) {
		if (in.validMask & (SENSOR_FRAME_VALID_IMU0 << i)) {
			for (uint8_t axis = 0; axis < 3; axis++) {
				accelStat[i][axis].Add(in.imu[i].accel[axis]);
				gyroStat[i][axis].Add(in.imu[i].gyro[axis]);
			}
		}

		if (in.validMask & (SENSOR_FRAME_VALID_BARO0 << i)) {
			pressureStat[i].Add(in.baro[i].pressure);
		}
	}

	sampleCount++;
}

bool TareStats::WindowComplete() const
{
	return (sampleCount >= TARE_MIN_SAMPLES) && (ElapsedMs() >= TARE_MIN_DURATION_MS);
}

bool TareStats::IsStill() const
{
	for (uint8_t i = 0; i < 2; i++) {
		if (!ImuValid(i)) {
			continue;
		}
		for (uint8_t axis = 0; axis < 3; axis++) {
			if (accelStat[i][axis].Std() > TARE_MAX_ACCEL_STD_MPS2 ||
				gyroStat[i][axis].Std() > TARE_MAX_GYRO_STD_RADPS) {
				return false;
			}
		}
	}
	return true;
}

/**
 * @brief A sensor is valid if it delivered data for at least half of the window
 */
bool TareStats::ImuValid(uint8_t imu) const
{
	return accelStat[imu][0].Count() * 2 >= sampleCount && sampleCount > 0;
}

bool TareStats::BaroValid(uint8_t baro) const
{
	return pressureStat[baro].Count() * 2 >= sampleCount && sampleCount > 0;
}

void TareStats::FillResult(TareResult& out) const
{
	for (uint8_t i = 0; i < 2; i++) {
		out.baroValid[i] = BaroValid(i);
		out.groundPressure[i] = pressureStat[i].Mean();
		out.pressureStd[i] = pressureStat[i].Std();

		out.imuValid[i] = ImuValid(i);
		for (uint8_t axis = 0; axis < 3; axis++) {
			out.gravity[i][axis]  = accelStat[i][axis].Mean();
			out.accelStd[i][axis] = accelStat[i][axis].Std();
			out.gyroBias[i][axis] = gyroStat[i][axis].Mean();
			out.gyroStd[i][axis]  = gyroStat[i][axis].Std();
		}
	}

	out.sampleCount = sampleCount;
	out.durationMs = ElapsedMs();
}
