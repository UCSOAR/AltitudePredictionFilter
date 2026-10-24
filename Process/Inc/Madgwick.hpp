/*
 * Madgwick.hpp
 *
 *  Created on: Oct 24, 2026
 *      Author: Aly Masani
 */

#ifndef MADGWICK_MADGWICK_HPP_
#define MADGWICK_MADGWICK_HPP_

#include <stdint.h>

/************************************
 * MACROS AND DEFINES
 ************************************/
// Gradient descent gain (Madgwick 2010). Larger trusts the accel more, smaller trusts the gyro more
constexpr float MADGWICK_BETA = 0.1f;

// Accel correction only runs while |accel| is within this of the pad gravity reading.
// Under thrust (|a| >> g) or in free fall / coast (|a| ~ 0) the accel is not gravity,
// so the attitude runs on the gyro alone (same idea as the legacy Fusion acceleration rejection)
constexpr float MADGWICK_ACCEL_REJECT_MPS2 = 1.0f;

/************************************
 * CLASS DEFINITIONS
 ************************************/
/**
 * @brief 6 axis Madgwick AHRS (gyro + accel, no mag) for one IMU.
 *        Earth frame: +Z up. The quaternion rotates IMU frame vectors into the earth frame,
 *        so the result does not depend on how the IMU is mounted
 */
class Madgwick
{
	public:
		/**
		 * @brief Start the attitude from the tare gravity reading, no convergence time needed
		 * @param gravity    Mean accel at rest in this IMU's frame (m/s^2), from the tare
		 */
		void Init(const float gravity[3]);

		/**
		 * @param gyro  rad/s, bias already removed
		 * @param accel m/s^2, includes gravity
		 * @param dt    s since the last update
		 */
		void Update(const float gyro[3], const float accel[3], float dt);

		/**
		 * @brief Net vertical acceleration: accel rotated into the earth frame, minus the
		 *        pad gravity magnitude (cancels this IMU's scale error at 1 g). Up positive
		 */
		float VerticalAccel(const float accel[3]) const;

		void Rotate(const float in[3], float out[3]) const;	// IMU frame -> earth frame
		float TiltDeg() const;								// Angle between IMU "up" at the pad and now
		bool IsInitialized() const { return initialized; }
		bool AccelUsed() const { return accelUsed; }		// false while accel correction is rejected
		const float* Quaternion() const { return q; }

	private:
		float q[4] = {1.0f, 0.0f, 0.0f, 0.0f};	// w, x, y, z
		float gravityMagnitude = 9.80665f;
		float padUp[3] = {0.0f, 0.0f, 1.0f};	// IMU frame "up" at the pad
		bool initialized = false;
		bool accelUsed = false;
};

#endif /* MADGWICK_MADGWICK_HPP_ */
