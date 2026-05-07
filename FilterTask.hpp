/*
 * filterTask.hpp
 *
 *  Created on: Apr 25, 2026
 *      Author: Andrey
 */

#ifndef ALTITUDEPREDICTIONFILTER_FILTERTASK_HPP_
#define ALTITUDEPREDICTIONFILTER_FILTERTASK_HPP_

/************************************
 * INCLUDES
 ************************************/
#include "Task.hpp"
#include "SystemDefines.hpp"

/************************************
 * MACROS AND DEFINES
 ************************************/
#define CMD_SET_FILTER_STATE      0x10
#define CMD_REQUEST_FILTER_STATE  0x11
#define CMD_FILTER_STATE          0x12

/************************************
 * TYPEDEFS
 ************************************/
enum class FilterEvent : uint8_t {
    NONE = 0,
    TARE_INITIATED,
	LAUNCH_COMMAND,
    LAUNCH_DETECTED,
    BURNOUT_DETECTED,
    ENABLE_AIRBRAKE,
    APOGEE_DETECTED,
    RESET,

    // CAN requests
    CAN_SET_STATE,
    CAN_REQUEST_STATE
};

enum class FilterState : uint8_t {
    TARE = 0,
    LAUNCH,
    BOOST,
    BURNOUT_DETECTED,
    AIRBRAKE_CONTROL,
    APOGEE_DETECTED
};

/************************************
 * CLASS DEFINITIONS
 ************************************/
class FilterTask : public Task
{
public:
    static FilterTask& Inst() {
        static FilterTask inst;
        return inst;
    }

    void InitTask();
    // Refresh rate API (milliseconds)
    void SetRefreshMs(uint32_t ms);
    uint32_t GetRefreshMs();

protected:
    static void RunTask(void* pvParams) { FilterTask::Inst().Run(pvParams); } // Static Task Interface, passes control to the instance Run();
    void Run(void * pvParams); // Main run code
    void HandleCommand(Command& cm);

private:
    // Private Functions
    FilterTask();        // Private constructor
    FilterTask(const FilterTask&);                        // Prevent copy-construction
    FilterTask& operator=(const FilterTask&);            // Prevent assignment
    // Refresh interval in milliseconds used for the periodic filter step
    uint32_t refreshMs_;

    // STATE MACHINE
    FilterState state_;
    FilterState requestedState_;   // from CAN
    bool stateOverride_;           // allow external override

    void ProcessState();
    void HandleEvent(FilterEvent evt);
    void TransitionTo(FilterState newState);

    // CAN helpers
    void SendStateCAN(FilterState state);
    void HandleCANCommand(Command& cm);

    void RequestFilterState(FilterState newState)
    {
        this->TransitionTo(newState);
    }
};

/************************************
 * FUNCTION DECLARATIONS
 ************************************/

#endif /* ALTITUDEPREDICTIONFILTER_FILTERTASK_HPP_ */
