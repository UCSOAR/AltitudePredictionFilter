#pragma once
#include <cstdint>

enum class FILTER_STATE : uint8_t
{
    PRE_START = 0,
    STARTED = 1,
    TAREING = 2,
    TARED = 3,
    LAUNCHED = 4,
    PRE_APOGEE_AIRBRAKES = 5,
    POST_APOGEE_AIRBRAKES = 6,
    POST_APOGEE_NOAIRBRAKES = 7,
    RESTARTING = 10,
    FAILED = 255
};