#pragma once

// ---------------------------------------------------------------
// INCLUDE ORDER REQUIREMENT:
//   Do NOT include this file directly. It is included at the
//   bottom of HALO.hpp, after Scenario and VectorXf are defined.
//   Any TU that needs AirbrakeController should just include HALO.hpp.
// ---------------------------------------------------------------

#include <array>
#include <vector>

#include "CommonDefinesFilter.hpp"
#include <vector>
#include "AirbrakeSims.hpp"

#include "Eigen\Dense"
#include "Eigen\Cholesky"
using namespace Eigen;


#ifndef TARGET_APOGEE_M
#define TARGET_APOGEE_M        5000.0f
#endif
#ifndef AIRBRAKE_OVERSHOOT_BAND
#define AIRBRAKE_OVERSHOOT_BAND 200.0f
#endif
#define AIRBRAKE_LEVELS 10

// ============================================================
// AirbrakeController
// ============================================================
class AirbrakeController {
public:

    void init(const std::vector<float>& apogees);

    int calculate_level(const std::vector<float>& gain_1, const std::vector<float>& gain_2,
    int scenario_index_1, int scenario_index_2,
    uint32_t curr_alt);

private:
    int  currentLevel_ = 0;
    std::vector<float> apogees;
    // pointer to the global LUT declared in AirbrakeSims.hpp
    uint16_t (*airbrake_sims)[ALT_LEVELS][BRAKE_LEVELS] = nullptr;
};
