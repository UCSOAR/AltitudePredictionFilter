/*
 * AxisAlignment.cpp
 *
 *  Created on: Oct 17, 2026
 *      Author: Aly Masani
 */

/**
 ********************************************************************************
 * @file    AxisAlignment.cpp
 * @author  Aly Masani
 * @date    Oct 17, 2026
 * @brief   Known IMU mounting transforms and the pad gravity check
 ********************************************************************************
 */

/************************************
 * INCLUDES
 ************************************/
#include "Inc/AxisAlignment.hpp"
#include <cmath>

/************************************
 * PRIVATE MACROS AND DEFINES
 ************************************/
namespace {
constexpr float RAD_TO_DEG = 180.0f / 3.14159265f;
}

/************************************
 * FUNCTION DEFINITIONS
 ************************************/
void AxisAlignment::Apply(const AxisMap& map, const float sensor[3], float body[3])
{
	for (uint8_t row = 0; row < 3; row++) {
		body[row] = map.m[row][0] * sensor[0]
				  + map.m[row][1] * sensor[1]
				  + map.m[row][2] * sensor[2];
	}
}

void AxisAlignment::ApplyToImu(uint8_t imu, const float sensor[3], float body[3])
{
	Apply(IMU_MOUNTING[imu], sensor, body);
}

void AxisAlignment::CheckGravity(TareResult& tare)
{
	tare.alignmentOk = true;

	for (uint8_t i = 0; i < 2; i++) {
		ApplyToImu(i, tare.gravity[i], tare.gravityBody[i]);

		if (!tare.imuValid[i]) {
			tare.padTiltDeg[i] = 0.0f;
			continue;
		}

		const float* g = tare.gravityBody[i];
		const float magnitude = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
		if (magnitude <= 0.0f) {
			tare.padTiltDeg[i] = 180.0f;
			tare.alignmentOk = false;
			continue;
		}

		// Clamp so float rounding can never push acos outside [-1, 1]
		float cosTilt = g[2] / magnitude;
		cosTilt = (cosTilt > 1.0f) ? 1.0f : ((cosTilt < -1.0f) ? -1.0f : cosTilt);

		tare.padTiltDeg[i] = std::acos(cosTilt) * RAD_TO_DEG;
		if (tare.padTiltDeg[i] > PAD_MAX_TILT_DEG) {
			tare.alignmentOk = false;
		}
	}
}
