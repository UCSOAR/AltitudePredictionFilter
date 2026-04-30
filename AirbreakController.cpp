#include <AirbrakeController.hpp>
#include <cmath>
#include <limits>

// ============================================================
// init
// ============================================================
void AirbrakeController::init(
    std::vector<std::vector<Scenario>>& brakeScenariosPerLevel,
    const std::array<AirbrakeDeploymentWindow, AIRBRAKE_LEVELS>& windows)
{
    brakeScenariosPerLevel_ = &brakeScenariosPerLevel;
    windows_     = windows;
    currentLevel_ = 0;
    initialized_  = true;
}

// ============================================================
// update  — called every filter cycle post-burnout
// ============================================================
int AirbrakeController::update(
    const VectorXf&                         X0,
    const std::vector<float>&               prevGain1,
    const std::vector<float>&               prevGain2,
    const std::vector<std::vector<float>>&  nearestVectors)
{
    if (!initialized_) return 0;

    // --- 1. Project apogee using existing gains + scenario future-vectors ---
    float projectedApogee = projectApogee(prevGain1, prevGain2, nearestVectors);

    SOAR_PRINT("AirbrakeController - projected apogee: %.1f ft, target: %.1f ft\n",
               projectedApogee, static_cast<float>(TARGET_APOGEE_FT));

    // --- 2. Are we overshooting beyond the deadband? ---
    if (projectedApogee <= TARGET_APOGEE_FT + AIRBRAKE_OVERSHOOT_BAND) {
        // No action needed — within acceptable band
        return currentLevel_;
    }

    // --- 3. Select the minimum viable level for our current state ---
    int needed = selectLevel(X0);

    // --- 4. Ratchet: only go up, never retract ---
    if (needed > currentLevel_) {
        currentLevel_ = needed;

        SOAR_PRINT("AirbrakeController - escalating to level %d\n", currentLevel_);

        // [CANBUS] Publish airbrake level command:
        //   e.g. CanBus::Publish(AIRBRAKE_LEVEL_CMD, currentLevel_);
        // The receiving actuator task maps 1-10 → servo angle and drives brakes.
    }

    return currentLevel_;
}

// ============================================================
// projectApogee
//
// Reuses prevGain1/prevGain2 (already computed by HALO's
// predictNextValues for sigma-point 0) to interpolate between
// the two nearest scenario future-vectors, then extrapolates
// forward using kinematic integration until velocity ≤ 0.
//
// This is deliberately lightweight — no KDTree lookup, no new
// scenario search. We are just extending what HALO already did.
// ============================================================
float AirbrakeController::projectApogee(
    const std::vector<float>&               prevGain1,
    const std::vector<float>&               prevGain2,
    const std::vector<std::vector<float>>&  nearestVectors) const
{
    // nearestVectors layout (same as findNearestScenarios output):
    //   [0] currentVector1  {alt, vel, accel, time, ...}
    //   [1] futureVector1
    //   [2] currentVector2
    //   [3] futureVector2

    if (nearestVectors.size() < 4) return 0.0f;

    const auto& v1 = nearestVectors[1];  // future step of scenario 1
    const auto& v2 = nearestVectors[3];  // future step of scenario 2

    // Interpolated next-step state using the gains HALO already computed
    float alt   = prevGain1[0] * v1[0] + prevGain2[0] * v2[0];
    float vel   = prevGain1[1] * v1[1] + prevGain2[1] * v2[1];
    float accel = prevGain1[2] * v1[2] + prevGain2[2] * v2[2];

    // Kinematic forward integration until velocity flips sign (apogee)
    // Uses a fixed timestep matching HALO's deltaTime (1/3 s default).
    // This is a straight-line extrapolation — cheap, no scenario lookup.
    const float dt       = 1.0f / 3.0f;
    const int   maxSteps = 500;           // ~167 s ceiling, well past any apogee

    for (int i = 0; i < maxSteps; i++) {
        if (vel <= 0.0f) break;

        float newVel = vel + accel * dt;
        float newAlt = alt + (vel + newVel) * 0.5f * dt;

        // Acceleration decays toward gravity as drag decreases with speed
        // (simple linear drag model — matches what the sims use post-burnout)
        accel = accel + (-9.81f - accel) * 0.05f;

        alt = newAlt;
        vel = newVel;
    }

    return alt;
}

// ============================================================
// selectLevel
//
// Walks brake levels 1→10 (ascending aggressiveness).
// For each level, checks whether the current filtered state
// {altitude, velocity} is still inside that level's deployment
// window — i.e. we haven't yet passed the last viable moment
// to deploy level L and still hit TARGET_APOGEE.
//
// The first level whose window we are still inside is the
// minimum viable level. Return it so the ratchet can compare
// against the previously commanded level.
//
// "Inside the window" means:
//   current_altitude  <= window.altitude   (haven't passed that altitude yet)
//   current_velocity  >= window.velocity   (still fast enough that level L works)
//
// If we've passed every window (all levels' opportunities gone),
// return AIRBRAKE_LEVELS (maximum — best we can do).
// ============================================================
int AirbrakeController::selectLevel(const VectorXf& X0) const
{
    // X0 layout in HALO: {altitude, velocity, acceleration}
    float curAlt = X0(0);
    float curVel = X0(1);

    for (int lvl = 0; lvl < AIRBRAKE_LEVELS; lvl++) {
        const auto& w = windows_[lvl];

        bool altStillInWindow = (curAlt <= w.altitude);
        bool velStillInWindow = (curVel >= w.velocity);

        if (altStillInWindow && velStillInWindow) {
            // Level lvl+1 (1-indexed) is the minimum viable level
            return lvl + 1;
        }
    }

    // Passed all windows — command maximum brakes
    return AIRBRAKE_LEVELS;
}
