/*
 * FilterLog.hpp
 *
 *  Created on: Oct 24, 2026
 *      Author: Aly Masani
 */

#ifndef FILTERLOG_FILTERLOG_HPP_
#define FILTERLOG_FILTERLOG_HPP_

#include "SensorDataTypes.hpp"
#include <cstdint>
#include <cstring>

/************************************
 * FUNCTION DECLARATIONS
 ************************************/
/**
 * @brief Flash log layouts of the altitude prediction filter. Header only so the LoggingTask,
 *        the flash dump (LoggingService) and the PC tools all use the exact same layout.
 *        A slot is 20 bytes: 1 type byte + 4 timestamp bytes + 15 payload bytes (these).
 *
 *   FILTER (every published FilterData: every frame in flight, 10 Hz otherwise)
 *   offset  size  field                 unit
 *    0      3     altitude              int24  cm (+/-83 km)
 *    3      2     velocity              int16  0.1 m/s   (+/-3276 m/s)
 *    5      2     acceleration          int16  0.01 m/s^2 (+/-327 m/s^2)
 *    7      2     altitudeStd           uint16 cm
 *    9      2     predictedApogee       uint16 m
 *    11     2     timeToEvent           int16  0.1 s (+/-3276 s, a whole descent fits)
 *    13     1     apfState
 *    14     1     flags                 FILTER_FLAG_*
 *
 *   FILTER2 (diagnostics, 1 Hz)
 *    0      2     velocityStd           uint16 cm/s
 *    2      1     caseA                 int8   nearest scenario (-1 none)
 *    3      1     caseB                 int8
 *    4      1     weightA               uint8  /255
 *    5      2     innovationAlt         int16  cm
 *    7      1     gainAlt               uint8  /255 (0 = model, 1 = sensors)
 *    8      1     gainVel               uint8  /255
 *    9      1     gainAcc               uint8  /255
 *    10     1     baroWeight0           uint8  /255 (Dynamite trust)
 *    11     1     baroWeight1           uint8  /255
 *    12     2     dynamiteSigma         uint16 cm
 *    14     1     reserved
 *
 *   EVENT (written through to flash at once)
 *    0      1     event                 FilterLog::Event
 *    1      4     value                 int32  (altitude cm for flight events)
 *    5      1     apfState
 */
namespace FilterLog
{
constexpr uint8_t PAYLOAD_BYTES = 15;

// Flight events, FilterData::event. Logged as an EVENT record and flushed to flash at once
enum Event : uint8_t
{
	EVENT_NONE = 0,
	EVENT_LAUNCH,				// APF DETECTION -> LAUNCH (board FlightTask launch)
	EVENT_APOGEE,				// Apogee detected, value = apogee altitude
	EVENT_APOGEE_MISSED,		// Power came back after apogee: marked reached, value = saved max altitude
	EVENT_LANDED,				// APF LAUNCH -> LANDED
	EVENT_BARO_DISTRUST_ON,		// Baros distrusted (transonic / ejection)
	EVENT_BARO_DISTRUST_OFF,
	EVENT_IMU_CLIP,				// An IMU clipped (at most once a second)
	EVENT_STATE_RESTORED,		// Booted into a saved FLIGHT / LANDED state, value = saved phase
	EVENT_FALSE_LAUNCH,			// FlightTask went back to GROUND without the launch ever climbing: back to the pad, re-tare
	EVENT_AIRBRAKE_CHECK,		// Airbrake check (log only), every 2 s from burnout, value = AirbrakeCheck::EventValue
	EVENT_AIRBRAKE_FLAG			// First check with level >= 1 (latched, record time = time flagged), same value
};

inline int32_t Round(float v) { return static_cast<int32_t>(v + ((v >= 0.0f) ? 0.5f : -0.5f)); }
inline int32_t Clamp(int32_t v, int32_t lo, int32_t hi) { return (v < lo) ? lo : ((v > hi) ? hi : v); }
inline uint8_t Unit8(float v) { return static_cast<uint8_t>(Clamp(Round(v * 255.0f), 0, 255)); }

inline void Pack(const FilterData& d, uint8_t* out)
{
	const int32_t alt = Clamp(Round(d.altitude * 100.0f), -8388608, 8388607);
	const int16_t vel = static_cast<int16_t>(Clamp(Round(d.velocity * 10.0f), -32768, 32767));
	const int16_t acc = static_cast<int16_t>(Clamp(Round(d.acceleration * 100.0f), -32768, 32767));
	const uint16_t sd = static_cast<uint16_t>(Clamp(Round(d.altitudeStd * 100.0f), 0, 65535));
	const uint16_t apo = static_cast<uint16_t>(Clamp(Round(d.predictedApogee), 0, 65535));
	const int16_t tte = static_cast<int16_t>(Clamp(Round(d.timeToEvent * 10.0f), -32768, 32767));
	out[0] = static_cast<uint8_t>(alt & 0xFF);			// int24, little endian
	out[1] = static_cast<uint8_t>((alt >> 8) & 0xFF);
	out[2] = static_cast<uint8_t>((alt >> 16) & 0xFF);
	std::memcpy(out + 3, &vel, 2);
	std::memcpy(out + 5, &acc, 2);
	std::memcpy(out + 7, &sd, 2);
	std::memcpy(out + 9, &apo, 2);
	std::memcpy(out + 11, &tte, 2);
	out[13] = d.apfState;
	out[14] = d.flags;
}

inline FilterData Unpack(const uint8_t* in, uint32_t timestampMs)
{
	int32_t alt = static_cast<int32_t>(in[0]) | (static_cast<int32_t>(in[1]) << 8) | (static_cast<int32_t>(in[2]) << 16);
	if (alt & 0x800000) {
		alt -= 0x1000000;		// sign extend the int24
	}
	int16_t vel, acc, tte;
	uint16_t sd, apo;
	std::memcpy(&vel, in + 3, 2);
	std::memcpy(&acc, in + 5, 2);
	std::memcpy(&sd, in + 7, 2);
	std::memcpy(&apo, in + 9, 2);
	std::memcpy(&tte, in + 11, 2);
	FilterData d{};
	d.timestampMs = timestampMs;
	d.altitude = alt / 100.0f;
	d.velocity = vel / 10.0f;
	d.acceleration = acc / 100.0f;
	d.altitudeStd = sd / 100.0f;
	d.predictedApogee = static_cast<float>(apo);
	d.timeToEvent = tte / 10.0f;
	d.apfState = in[13];
	d.flags = in[14];
	return d;
}

inline void PackDiagnostics(const FilterData& d, uint8_t* out)
{
	const uint16_t vsd = static_cast<uint16_t>(Clamp(Round(d.velocityStd * 100.0f), 0, 65535));
	const int16_t inn = static_cast<int16_t>(Clamp(Round(d.innovationAlt * 100.0f), -32768, 32767));
	const uint16_t sig = static_cast<uint16_t>(Clamp(Round(d.dynamiteSigma * 100.0f), 0, 65535));
	std::memcpy(out + 0, &vsd, 2);
	out[2] = static_cast<uint8_t>(d.caseA);
	out[3] = static_cast<uint8_t>(d.caseB);
	out[4] = Unit8(d.weightA);
	std::memcpy(out + 5, &inn, 2);
	out[7] = Unit8(d.gainAlt);
	out[8] = Unit8(d.gainVel);
	out[9] = Unit8(d.gainAcc);
	out[10] = Unit8(d.baroWeight[0]);
	out[11] = Unit8(d.baroWeight[1]);
	std::memcpy(out + 12, &sig, 2);
	out[14] = 0;
}

inline void PackEvent(const FilterData& d, uint8_t* out)
{
	std::memset(out, 0, PAYLOAD_BYTES);
	out[0] = d.event;
	std::memcpy(out + 1, &d.eventValue, 4);
	out[5] = d.apfState;
}

/**
 * @brief What survives a power loss (LoggingService keeps it in its own flash sector).
 *        Saved on every phase change, at apogee and every few seconds in flight
 */
enum Phase : uint8_t
{
	PHASE_PAD = 0,			// Before launch: normal boot, tare, wait for launch
	PHASE_FLIGHT = 1,		// Board FlightTask in FLIGHT
	PHASE_LANDED = 2		// Back on the ground after a flight: stays here until flight_reset
};

constexpr uint16_t TARE_BYTES = 192;	// TareResult copy (static_assert in SensorFrameTask)

struct PersistedState
{
	uint32_t bootCount;
	uint8_t phase;				// Phase
	uint8_t apfState;
	uint8_t apogeeReached;
	uint8_t tareValid;			// tare[] holds a finished TareResult
	float maxAltitude;			// m above the pad
	float lastAltitude;			// m, last estimate when saved
	float lastVelocity;			// m/s
	float predictedApogee;		// m, apogee predictor before apogee (0 = none): a missed apogee uses it
	uint8_t tare[TARE_BYTES];
};
}

#endif /* FILTERLOG_FILTERLOG_HPP_ */
