/*
 * AxisAlignment.hpp
 *
 *  Created on: Oct 17, 2026
 *      Author: Aly Masani
 */

#ifndef AXISALIGNMENT_AXISALIGNMENT_HPP_
#define AXISALIGNMENT_AXISALIGNMENT_HPP_

#include "TareStats.hpp"

/**
 * @brief Rotates IMU readings from sensor axes into rocket body axes.
 *        Body frame: +Z toward the nose, X/Y across the rocket.
 *        On the pad an accelerometer reads +1 g along "up", so a correctly mounted
 *        IMU sees its tare gravity on body +Z
 */
namespace AxisAlignment
{
/**
 * @brief Signed axis permutation: body[row] = sum(m[row][col] * sensor[col]).
 *        Each row has exactly one entry of +1 or -1
 */
struct AxisMap
{
	int8_t m[3][3];
};

inline constexpr AxisMap IDENTITY = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};

// Mounting of each IMU in the rocket, index = IMU id.
#ifdef APF_PC_REPLAY
// PC replay / Taber test: the generated flights are already in body axes
inline constexpr AxisMap IMU_MOUNTING[2] = {IDENTITY, IDENTITY};
#else
// DAQLite, measured on the bench (APF_RawStream, board flat then held nose up):
//   nose (body +Z) = IMU0 +X = IMU1 -Y
//   board normal (sensor +Z, the same for both IMUs) = body +X, body Y completes the right hand set
// Nose up both IMUs see gravity on body +Z (10 deg tilt, hand held), flat on body +X,
// and the two IMUs agree within ~20 mg in both positions
inline constexpr AxisMap IMU0_DAQLITE = {{{0, 0, 1}, {0, -1, 0}, {1, 0, 0}}};	// LSM6DSO
inline constexpr AxisMap IMU1_DAQLITE = {{{0, 0, 1}, {-1, 0, 0}, {0, -1, 0}}};	// LSM6DSO32 "A"
inline constexpr AxisMap IMU_MOUNTING[2] = {IMU0_DAQLITE, IMU1_DAQLITE};
#endif

// Max angle between body +Z and measured "up" on the pad before the alignment is flagged
constexpr float PAD_MAX_TILT_DEG = 30.0f;

// body = map * sensor
void Apply(const AxisMap& map, const float sensor[3], float body[3]);

// Rotate a reading of IMU imu (0/1) into body axes
void ApplyToImu(uint8_t imu, const float sensor[3], float body[3]);

/**
 * @brief Fills tare.gravityBody, tare.padTiltDeg and tare.alignmentOk from tare.gravity.
 *        alignmentOk is false if any valid IMU sees "up" more than PAD_MAX_TILT_DEG
 *        away from body +Z, which usually means IMU_MOUNTING is wrong
 */
void CheckGravity(TareResult& tare);
}

#endif /* AXISALIGNMENT_AXISALIGNMENT_HPP_ */
