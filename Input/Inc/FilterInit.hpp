/*
 * FilterInit.hpp
 *
 *  Created on: Oct 17, 2026
 *      Author: Aly Masani
 */

#ifndef FILTERINIT_FILTERINIT_HPP_
#define FILTERINIT_FILTERINIT_HPP_

#include "TareStats.hpp"

/************************************
 * TYPEDEFS
 ************************************/
/**
 * @brief Fixed size 3x3 matrix, avoids pulling Eigen into the flight code
 */
struct Mat3
{
	float m[3][3];
};

/**
 * @brief Initial conditions for the UKF. State order: altitude, velocity, acceleration.
 *        Altitude is relative to the pad (x0 altitude = 0), up is positive
 * @param groundAltitude Absolute pad altitude (m, standard atmosphere), for reference only
 * @param x0 Initial state [m, m/s, m/s^2]
 * @param P0 Initial state covariance, built from the measured tare noise
 * @param Q  Process noise spectral density (per second, the UKF scales it by dt)
 * @param R  Measurement noise covariance floors (alt, vel, acc)
 */
struct FilterInitState
{
	bool valid;
	float groundAltitude;
	float x0[3];
	Mat3 P0;
	Mat3 Q;
	Mat3 R;
};

/************************************
 * FUNCTION DECLARATIONS
 ************************************/
namespace FilterInit
{
// Velocity on the pad is known to be 0, but a zero variance would make the filter
// ignore every velocity correction, so keep a small floor on every P0 diagonal
constexpr float P0_VEL_VAR = 0.01f;		// (0.1 m/s)^2
constexpr float P0_MIN_ALT_VAR = 0.01f;	// (0.1 m)^2
constexpr float P0_MIN_ACC_VAR = 0.01f;	// (0.1 m/s^2)^2

/**
 * @brief Builds x0, P0, Q and R from a finished tare
 */
FilterInitState Build(const TareResult& tare);
}

#endif /* FILTERINIT_FILTERINIT_HPP_ */
