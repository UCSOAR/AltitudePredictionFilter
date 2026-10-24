/*
 * Barometric.cpp
 *
 *  Created on: Oct 17, 2026
 *      Author: Aly Masani
 */

/**
 ********************************************************************************
 * @file    Barometric.cpp
 * @author  Aly Masani
 * @date    Oct 17, 2026
 * @brief   h = 44330 * (1 - (p / p0)^(1 / 5.2558))
 ********************************************************************************
 */

/************************************
 * INCLUDES
 ************************************/
#include "Inc/Barometric.hpp"
#include <cmath>

/************************************
 * PRIVATE MACROS AND DEFINES
 ************************************/
namespace {
constexpr float FORMULA_SCALE_M = 44330.0f;
constexpr float FORMULA_EXPONENT = 1.0f / 5.2558f;
}

/************************************
 * FUNCTION DEFINITIONS
 ************************************/
float Barometric::PressureToAltitude(float pressurePa)
{
	if (pressurePa < MIN_VALID_PRESSURE_PA) {
		return 0.0f;
	}

	return FORMULA_SCALE_M * (1.0f - std::pow(pressurePa / SEA_LEVEL_PRESSURE_PA, FORMULA_EXPONENT));
}

float Barometric::AltitudePerPascal(float pressurePa)
{
	if (pressurePa < MIN_VALID_PRESSURE_PA) {
		return 0.0f;
	}

	return FORMULA_SCALE_M * FORMULA_EXPONENT
			* std::pow(pressurePa / SEA_LEVEL_PRESSURE_PA, FORMULA_EXPONENT - 1.0f)
			/ SEA_LEVEL_PRESSURE_PA;
}
