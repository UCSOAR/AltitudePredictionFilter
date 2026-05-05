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
    // used for projection
    this->apogees = apogees;

    // used for airbrake level [SIM_I][ALT_deployed][BRAKE_LEVEL]
    this->airbrake_sims = airbrake_lut;

}

// ============================================================
// calculate_level - returns airbrake level computed from measurements
// ============================================================
int AirbrakeController::calculate_level(
    uint16_t gain_1, uint16_t gain_2,
    int scenario_index_1, int scenario_index_2,
    uint32_t curr_alt
    )
{
    // Step 1:Check Projected Apogee > Target
    uint32_t alt_projected = this->apogees[scenario_index_1] * gain_1 + this->apogees[scenario_index_2] * gain_2;
    int curr_lvl = 0;
    uint32_t curr_alt_airbrakes = TARGET_APOGEE_M + 1;
    uint32_t prev_alt_airbrakes = 0.0f;

    // find next 1k checkpoint (1500 -> 2k)
    uint32_t checkpoint_alt = curr_alt % 1000 + 1;

    if (alt_projected > TARGET_APOGEE_M){
        while(curr_alt_airbrakes > TARGET_APOGEE_M){
            prev_alt_airbrakes = curr_alt_airbrakes;
            curr_alt_airbrakes = this->airbrake_lut[scenario_index_1][checkpoint_alt][curr_lvl] * gain_1 + this->airbrake_lut[scenario_index_2][checkpoint_alt][curr_lvl] + gain_2;
        }
        // altitude_airbrakes can be undershot (5.3k at lvl 2 = z - 1, 4.7k at lvl = 3 = z <---)
        if (TARGET_APOGEE_M  - curr_alt_airbrakes > prev_alt_airbrakes - TARGET_APOGEE_M){
            // use z - 1 since it has smaller diff
            return curr_lvl - 1;
        }else{
            return curr_lvl;
        }
    }
}

