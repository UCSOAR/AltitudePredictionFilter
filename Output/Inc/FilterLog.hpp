/*
 * FilterLog.hpp
 *
 *  Created on: Oct 24, 2026
 *      Author: Aly Masani
 */

#ifndef FILTERLOG_FILTERLOG_HPP_
#define FILTERLOG_FILTERLOG_HPP_

#include "SensorDataTypes.hpp"
#include <cstring>

/************************************
 * FUNCTION DECLARATIONS
 ************************************/
/**
 * @brief Packs FilterData into the 15 data bytes of one 20 byte flash log slot
 *        (slot = 1 type byte + 4 timestamp bytes + these 15). Header only so the
 *        LoggingTask, the flash dump and the PC replay all use the exact same layout.
 *
 *   offset  size  field                 unit
 *    0      4     altitude              int32  cm
 *    4      2     velocity              int16  0.1 m/s   (+/-3276 m/s)
 *    6      2     acceleration          int16  0.01 m/s^2 (+/-327 m/s^2)
 *    8      2     altitudeStd           uint16 cm
 *    10     2     predictedApogee       uint16 m
 *    12     1     timeToEvent           int8   s (clamped to +/-127)
 *    13     1     apfState
 *    14     1     flags
 */
namespace FilterLog
{
constexpr uint8_t PAYLOAD_BYTES = 15;

inline int32_t Round(float v) { return static_cast<int32_t>(v + ((v >= 0.0f) ? 0.5f : -0.5f)); }
inline int32_t Clamp(int32_t v, int32_t lo, int32_t hi) { return (v < lo) ? lo : ((v > hi) ? hi : v); }

inline void Pack(const FilterData& d, uint8_t* out)
{
	const int32_t alt = Round(d.altitude * 100.0f);
	const int16_t vel = static_cast<int16_t>(Clamp(Round(d.velocity * 10.0f), -32768, 32767));
	const int16_t acc = static_cast<int16_t>(Clamp(Round(d.acceleration * 100.0f), -32768, 32767));
	const uint16_t sd = static_cast<uint16_t>(Clamp(Round(d.altitudeStd * 100.0f), 0, 65535));
	const uint16_t apo = static_cast<uint16_t>(Clamp(Round(d.predictedApogee), 0, 65535));
	const int8_t tte = static_cast<int8_t>(Clamp(Round(d.timeToEvent), -127, 127));

	std::memcpy(out + 0, &alt, 4);
	std::memcpy(out + 4, &vel, 2);
	std::memcpy(out + 6, &acc, 2);
	std::memcpy(out + 8, &sd, 2);
	std::memcpy(out + 10, &apo, 2);
	std::memcpy(out + 12, &tte, 1);
	out[13] = d.apfState;
	out[14] = d.flags;
}

inline FilterData Unpack(const uint8_t* in, uint32_t timestampMs)
{
	int32_t alt;
	int16_t vel, acc;
	uint16_t sd, apo;
	int8_t tte;
	std::memcpy(&alt, in + 0, 4);
	std::memcpy(&vel, in + 4, 2);
	std::memcpy(&acc, in + 6, 2);
	std::memcpy(&sd, in + 8, 2);
	std::memcpy(&apo, in + 10, 2);
	std::memcpy(&tte, in + 12, 1);

	FilterData d{};
	d.timestampMs = timestampMs;
	d.altitude = alt / 100.0f;
	d.velocity = vel / 10.0f;
	d.acceleration = acc / 100.0f;
	d.altitudeStd = sd / 100.0f;
	d.predictedApogee = static_cast<float>(apo);
	d.timeToEvent = static_cast<float>(tte);
	d.apfState = in[13];
	d.flags = in[14];
	return d;
}
}

#endif /* FILTERLOG_FILTERLOG_HPP_ */
