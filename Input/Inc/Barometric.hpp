/*
 * Barometric.hpp
 *
 *  Created on: Oct 17, 2026
 *      Author: Aly Masani
 */

#ifndef BAROMETRIC_BAROMETRIC_HPP_
#define BAROMETRIC_BAROMETRIC_HPP_

/**
 * @brief Pressure to altitude conversion (standard atmosphere / barometric formula),
 *        ported from convertToAltitude() in the legacy EverestTask
 */
namespace Barometric
{
constexpr float SEA_LEVEL_PRESSURE_PA = 101325.0f;
constexpr float MIN_VALID_PRESSURE_PA = 10000.0f;	// Below this the reading is treated as invalid

/**
 * @brief Altitude above sea level for a pressure (standard atmosphere)
 * @return altitude in m, 0 for invalid pressures
 */
float PressureToAltitude(float pressurePa);

/**
 * @brief Magnitude of d(altitude)/d(pressure) at a pressure, used to convert
 *        pressure noise (Pa) into altitude noise (m). About 0.09 m/Pa near the ground
 * @return m/Pa, 0 for invalid pressures
 */
float AltitudePerPascal(float pressurePa);
}

#endif /* BAROMETRIC_BAROMETRIC_HPP_ */
