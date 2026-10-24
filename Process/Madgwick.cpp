/*
 * Madgwick.cpp
 *
 *  Created on: Oct 24, 2026
 *      Author: Aly Masani
 */

/**
 ********************************************************************************
 * @file    Madgwick.cpp
 * @author  Aly Masani
 * @date    Oct 24, 2026
 * @brief   6 axis Madgwick AHRS (S. Madgwick 2010, IMU gradient descent update).
 *          Replaces the legacy Infusion port: no mag, no gravity removed before the
 *          update (legacy subtracted the full pad accel, gravity included)
 ********************************************************************************
 */

/************************************
 * INCLUDES
 ************************************/
#include "Inc/Madgwick.hpp"
#include <cmath>

/************************************
 * PRIVATE MACROS AND DEFINES
 ************************************/
namespace {
constexpr float RAD_TO_DEG = 180.0f / 3.14159265f;

float Norm3(const float v[3])
{
	return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}
}

/************************************
 * FUNCTION DEFINITIONS
 ************************************/
/**
 * @brief Shortest rotation that takes the measured "up" (gravity reading) onto earth +Z
 */
void Madgwick::Init(const float gravity[3])
{
	const float n = Norm3(gravity);
	if (n < 1e-3f) {
		initialized = false;
		return;
	}

	gravityMagnitude = n;
	for (uint8_t i = 0; i < 3; i++) {
		padUp[i] = gravity[i] / n;
	}

	// q = [1 + up.z, up x z] normalised (up x z = (up.y, -up.x, 0))
	const float dot = padUp[2];
	if (dot < -0.9999f) {
		// Upside down: 180 deg about X
		q[0] = 0.0f; q[1] = 1.0f; q[2] = 0.0f; q[3] = 0.0f;
	} else {
		q[0] = 1.0f + dot;
		q[1] = padUp[1];
		q[2] = -padUp[0];
		q[3] = 0.0f;
		const float qn = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
		for (float& c : q) {
			c /= qn;
		}
	}

	initialized = true;
}

void Madgwick::Update(const float gyro[3], const float accel[3], float dt)
{
	if (!initialized || dt <= 0.0f) {
		return;
	}

	float q0 = q[0], q1 = q[1], q2 = q[2], q3 = q[3];
	const float gx = gyro[0], gy = gyro[1], gz = gyro[2];

	// Rate of change of the quaternion from the gyro
	float qDot0 = 0.5f * (-q1 * gx - q2 * gy - q3 * gz);
	float qDot1 = 0.5f * ( q0 * gx + q2 * gz - q3 * gy);
	float qDot2 = 0.5f * ( q0 * gy - q1 * gz + q3 * gx);
	float qDot3 = 0.5f * ( q0 * gz + q1 * gy - q2 * gx);

	// Accel correction, only while the accel is mostly gravity
	const float aNorm = Norm3(accel);
	accelUsed = std::fabs(aNorm - gravityMagnitude) < MADGWICK_ACCEL_REJECT_MPS2;
	if (accelUsed) {
		const float ax = accel[0] / aNorm, ay = accel[1] / aNorm, az = accel[2] / aNorm;

		const float _2q0 = 2.0f * q0, _2q1 = 2.0f * q1, _2q2 = 2.0f * q2, _2q3 = 2.0f * q3;
		const float _4q0 = 4.0f * q0, _4q1 = 4.0f * q1, _4q2 = 4.0f * q2;
		const float _8q1 = 8.0f * q1, _8q2 = 8.0f * q2;
		const float q0q0 = q0 * q0, q1q1 = q1 * q1, q2q2 = q2 * q2, q3q3 = q3 * q3;

		// Gradient of the error between estimated and measured gravity direction
		float s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay;
		float s1 = _4q1 * q3q3 - _2q3 * ax + 4.0f * q0q0 * q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az;
		float s2 = 4.0f * q0q0 * q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az;
		float s3 = 4.0f * q1q1 * q3 - _2q1 * ax + 4.0f * q2q2 * q3 - _2q2 * ay;

		const float sn = std::sqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3);
		if (sn > 1e-9f) {
			qDot0 -= MADGWICK_BETA * s0 / sn;
			qDot1 -= MADGWICK_BETA * s1 / sn;
			qDot2 -= MADGWICK_BETA * s2 / sn;
			qDot3 -= MADGWICK_BETA * s3 / sn;
		}
	}

	q0 += qDot0 * dt;
	q1 += qDot1 * dt;
	q2 += qDot2 * dt;
	q3 += qDot3 * dt;

	const float qn = std::sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
	q[0] = q0 / qn; q[1] = q1 / qn; q[2] = q2 / qn; q[3] = q3 / qn;
}

/**
 * @brief v_earth = R(q) * v_imu
 */
void Madgwick::Rotate(const float in[3], float out[3]) const
{
	const float w = q[0], x = q[1], y = q[2], z = q[3];

	out[0] = (1.0f - 2.0f * (y * y + z * z)) * in[0] + 2.0f * (x * y - w * z) * in[1] + 2.0f * (x * z + w * y) * in[2];
	out[1] = 2.0f * (x * y + w * z) * in[0] + (1.0f - 2.0f * (x * x + z * z)) * in[1] + 2.0f * (y * z - w * x) * in[2];
	out[2] = 2.0f * (x * z - w * y) * in[0] + 2.0f * (y * z + w * x) * in[1] + (1.0f - 2.0f * (x * x + y * y)) * in[2];
}

float Madgwick::VerticalAccel(const float accel[3]) const
{
	float earth[3];
	Rotate(accel, earth);
	return earth[2] - gravityMagnitude;
}

/**
 * @brief Angle between earth +Z and where the pad "up" axis points now
 */
float Madgwick::TiltDeg() const
{
	float upNow[3];
	Rotate(padUp, upNow);
	float c = upNow[2];
	c = (c > 1.0f) ? 1.0f : ((c < -1.0f) ? -1.0f : c);
	return std::acos(c) * RAD_TO_DEG;
}
