#ifndef TESTING_BUILD
#include "SystemDefines.hpp"
#endif

#ifdef TESTING_BUILD
#define SOAR_PRINT(...) printf(__VA_ARGS__)
#include <cstdio>
#endif

#include <cmath>
#include <limits>
#include "AirbrakeController.hpp"

// ============================================================
// init
// ============================================================
void AirbrakeController::init(
    const std::vector<float>& apogees)
{
    // copy apogees used for projection
    this->apogees = apogees;

    // point at the compiled-in LUT from AirbrakeSims.hpp
    this->airbrake_sims = airbrake_lut;
}

// ============================================================
// calculate_level - returns airbrake level computed from measurements
// ============================================================
int AirbrakeController::calculate_level(
    const std::vector<float>& gain_1, const std::vector<float>& gain_2,
    int scenario_index_1, int scenario_index_2,
    uint32_t curr_alt
    )
{
    // guard
    if (this->airbrake_sims == nullptr || this->apogees.empty()){
        return this->currentLevel_;
    }

    // Use the altitude component of the gain vectors (index 0) as scalar weights.
    float w1 = (gain_1.size() > 0) ? gain_1[0] : 1.0f;
    float w2 = (gain_2.size() > 0) ? gain_2[0] : 1.0f;

    // Step 1: Check projected apogee > target
    float alt_projected = 0.0f;
    if (scenario_index_1 >= 0 && scenario_index_1 < (int)this->apogees.size())
        alt_projected += this->apogees[scenario_index_1] * w1;
    if (scenario_index_2 >= 0 && scenario_index_2 < (int)this->apogees.size())
        alt_projected += this->apogees[scenario_index_2] * w2;

    int curr_lvl_local = 0;
    float curr_alt_airbrakes = TARGET_APOGEE_M + 1.0f;
    float prev_alt_airbrakes = curr_alt_airbrakes;

    // find next 1k checkpoint (1500 -> 2k) and convert to LUT index (0..ALT_LEVELS-1)
    uint32_t checkpoint_level = (curr_alt + 999u) / 1000u; // 1..N
    if (checkpoint_level == 0) checkpoint_level = 1;
    if (checkpoint_level > ALT_LEVELS) checkpoint_level = ALT_LEVELS;
    size_t alt_idx = (size_t)(checkpoint_level - 1);

    if (alt_projected > TARGET_APOGEE_M){
        // increment brake level until predicted apogee reaches <= target or we hit max level
        while (curr_alt_airbrakes > TARGET_APOGEE_M && curr_lvl_local < BRAKE_LEVELS){
            prev_alt_airbrakes = curr_alt_airbrakes;
            uint16_t sim1 = this->airbrake_sims[scenario_index_1][alt_idx][curr_lvl_local];
            uint16_t sim2 = this->airbrake_sims[scenario_index_2][alt_idx][curr_lvl_local];
            curr_alt_airbrakes = sim1 * w1 + sim2 * w2;
            curr_lvl_local += 1;
        }

        // curr_lvl_local was incremented after computing curr_alt_airbrakes; candidate level is curr_lvl_local-1
        int candidate_lvl = curr_lvl_local - 1;
        if (candidate_lvl < 0) candidate_lvl = 0;
        if (candidate_lvl >= BRAKE_LEVELS) candidate_lvl = BRAKE_LEVELS - 1;

        // choose the level (current candidate vs previous) that gives smaller absolute error to target
        float diff_curr = std::fabs(TARGET_APOGEE_M - curr_alt_airbrakes);
        float diff_prev = std::fabs(prev_alt_airbrakes - TARGET_APOGEE_M);
        if (diff_prev < diff_curr) candidate_lvl = std::max(0, candidate_lvl - 1);

        if (candidate_lvl != this->currentLevel_){
            this->currentLevel_ = candidate_lvl;
            // TODO: alert level change
        }
    }

    return this->currentLevel_;
}

