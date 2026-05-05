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

    void init();

    int calculate_level();

private:
    int  currentLevel_ = 0;
    std::vector<float>* apogees;
    uint16_t* airbrake_lut[SIM_COUNT][ALT_LEVELS][BRAKE_LEVELS];
};
