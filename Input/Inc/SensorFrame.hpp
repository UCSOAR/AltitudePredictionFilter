/*
 * SensorFrame.hpp
 *
 *  Created on: Oct 17, 2026
 *      Author: Aly Masani
 */

#ifndef SENSORFRAME_SENSORFRAME_HPP_
#define SENSORFRAME_SENSORFRAME_HPP_

#include <stdint.h>

/************************************
 * TYPEDEFS
 ************************************/
/**
 * @brief Validity bits for SensorFrame::validMask, one bit per sensor
 */
enum SensorFrameValidBits : uint8_t
{
	SENSOR_FRAME_VALID_IMU0  = (1U << 0),
	SENSOR_FRAME_VALID_IMU1  = (1U << 1),
	SENSOR_FRAME_VALID_BARO0 = (1U << 2),
	SENSOR_FRAME_VALID_BARO1 = (1U << 3),
	SENSOR_FRAME_VALID_MAG   = (1U << 4),
	SENSOR_FRAME_VALID_ALL   = 0b11111
};

/**
 * @param accel Acceleration in m/s^2
 * @param gyro  Angular rate in rad/s
 * @param temp  Temperature in degC
 */
struct SensorFrameIMU
{
	float accel[3];
	float gyro[3];
	float temp;
};

/**
 * @param pressure Pressure in Pa
 * @param temp     Temperature in degC
 */
struct SensorFrameBaro
{
	float pressure;
	float temp;
};

/**
 * @brief One complete poll of all sensors, converted to SI units.
 *        Index 0/1 matches the sensor id (IMU0 = LSM6DSO, IMU1 = LSM6DSO32,
 *        BARO0 = MS5607, BARO1 = MS5611).
 *        Raw frames are in sensor axes; calibrated frames (from PreProcess) have the
 *        IMUs in rocket body axes (+Z toward the nose)
 * @param sequence    Increments every frame, gaps mean dropped frames
 * @param timestampMs Time the poll arrived (ms since boot)
 * @param dt          Seconds since the previous frame, 0 on the first frame
 * @param mag         Magnetic field in uT (sensor axes)
 * @param validMask   SensorFrameValidBits, a cleared bit means do not use that sensor
 */
struct SensorFrame
{
	uint32_t sequence;
	uint32_t timestampMs;
	float dt;
	SensorFrameIMU imu[2];
	SensorFrameBaro baro[2];
	float mag[3];
	uint8_t validMask;
};

#endif /* SENSORFRAME_SENSORFRAME_HPP_ */
