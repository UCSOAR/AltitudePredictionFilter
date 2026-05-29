/*
 * AirbrakeSims.hpp
 *
 *  Created on: May 3, 2026
 *      Author: Andrey
 */

#ifndef ALTITUDEPREDICTIONFILTER_AIRBRAKESIMS_HPP_
#define ALTITUDEPREDICTIONFILTER_AIRBRAKESIMS_HPP_

#pragma once
#include <cstdint>
#include <vector>

#define SIM_COUNT 3
#define ALT_LEVELS 10    // 1k–10k ft
#define BRAKE_LEVELS 10  // 0–9 airbrake levels

// [sim][deployed_alt][airbrake_lvl] = max_alt (ft)

// Declare the LUT as extern here; define it in AirbrakeSims.cpp to avoid
// multiple-definition linker errors when this header is included by many TUs.
extern uint16_t airbrake_lut[SIM_COUNT][ALT_LEVELS][BRAKE_LEVELS];


#endif /* ALTITUDEPREDICTIONFILTER_AIRBRAKESIMS_HPP_ */
