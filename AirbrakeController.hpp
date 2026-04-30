/*
 * Airbreakcontroller.hpp
 *
 *  Created on: Apr 30, 2026
 *      Author: Andrey
 */

#ifndef ALTITUDEPREDICTIONFILTER_AIRBRAKECONTROLLER_HPP_
#define ALTITUDEPREDICTIONFILTER_AIRBRAKECONTROLLER_HPP_

#pragma once

#include <vector>
#include <array>
#include "HALO.hpp"

// ============================================================
// Compile-time mission constants
// ============================================================
#define TARGET_APOGEE_FT        10000.0f   // feet — set per mission
#define AIRBRAKE_OVERSHOOT_BAND 400.0f     // ft above target before brakes engage
#define AIRBRAKE_LEVELS         10

// ============================================================
// Per-level deployment window descriptor.
//
// Each entry describes the LAST viable {altitude, velocity,
// acceleration} state at which deploying that brake level can
// still slow the rocket enough to hit TARGET_APOGEE.
//
// These vectors are extracted from the airbrake sim datasets:
// for level L, find the earliest time-step at which the sim
// trajectory converges to TARGET_APOGEE and record that state.
//
// Format matches Scenario data: {altitude, velocity, accel}
// ============================================================
struct AirbrakeDeploymentWindow {
    float altitude;     // ft  — maximum altitude still inside window
    float velocity;     // ft/s
    float acceleration; // ft/s^2
};

// ============================================================
// AirbrakeController
//
// Lives inside HALO (or alongside it). After burnout:
//   1. projectApogee()   — uses existing gains + airbrake scenarios
//   2. selectLevel()     — finds lowest level whose window fits
//   3. ratchet enforced  — level can only increase
//   4. actuate()         — sends level over CAN (stub comment)
// ============================================================
class AirbrakeController {
public:
    // Call once at init with the 10 sets of airbrake scenarios
    // (same Scenario structure as nominal flight, one set per level)
    void init(std::vector<std::vector<Scenario>>& brakeScenariosPerLevel,
              const std::array<AirbrakeDeploymentWindow, AIRBRAKE_LEVELS>& windows);

    // Call every filter cycle after burnout is confirmed.
    // X0         — current filtered state {alt, vel, accel}
    // prevGain1  — gain from sigma-point 0 nearest-scenario interpolation
    // prevGain2  — same, second scenario
    // scenarios  — the TWO nearest airbrake scenarios already found by HALO
    //              (same pair used in predictNextValues)
    // Returns the commanded brake level (0 = retracted, 10 = full deploy).
    int update(const VectorXf& X0,
               const std::vector<float>& prevGain1,
               const std::vector<float>& prevGain2,
               const std::vector<std::vector<float>>& nearestVectors);

    int currentLevel() const { return currentLevel_; }

private:
    // Projects apogee altitude using the two nearest scenario future-vectors
    // and the interpolation gains — identical math to predictNextValues but
    // walks forward until velocity crosses zero.
    float projectApogee(const std::vector<float>& prevGain1,
                        const std::vector<float>& prevGain2,
                        const std::vector<std::vector<float>>& nearestVectors) const;

    // Finds the lowest brake level whose deployment window we are still inside.
    // Returns 0 if no level is needed (not overshooting) or all windows passed.
    int selectLevel(const VectorXf& X0) const;

    // 10 sets of Scenario objects — index 0 = level 1, index 9 = level 10.
    // Each set has the same structure as the nominal HALO scenarios but
    // represents flight with airbrakes deployed at that level.
    std::vector<std::vector<Scenario>>* brakeScenariosPerLevel_ = nullptr;

    // Pre-computed deployment windows, one per level (ascending aggressiveness)
    std::array<AirbrakeDeploymentWindow, AIRBRAKE_LEVELS> windows_{};

    // Ratchet state — level only ever increases
    int currentLevel_ = 0;

    bool initialized_ = false;
};



#endif /* ALTITUDEPREDICTIONFILTER_AIRBRAKECONTROLLER_HPP_ */
