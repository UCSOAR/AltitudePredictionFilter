/*
 * PrintUtils.hpp
 *
 *  Created on: Oct 17, 2026
 *      Author: Aly Masani
 */

#ifndef PRINTUTILS_PRINTUTILS_HPP_
#define PRINTUTILS_PRINTUTILS_HPP_

/**
 * @brief SOAR_PRINT has no %f support (newlib-nano), so floats are printed as
 *        integer milli-units with %ld
 */
inline long ToMilli(float x)
{
	return static_cast<long>(x * 1000.0f + ((x >= 0.0f) ? 0.5f : -0.5f));
}

#endif /* PRINTUTILS_PRINTUTILS_HPP_ */
