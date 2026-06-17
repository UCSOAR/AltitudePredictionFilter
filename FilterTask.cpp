/**
 ********************************************************************************
 * @file    ${file_name}
 * @author  ${user}
 * @date    ${date}
 * @brief   This is a template source file to create a new task in our firmware
 *
 * Setup Steps
 * 1. Define the Task Queue Depth in SystemDefines.hpp
 * 2. Define the Task Stack Depth in SystemDefines.hpp
 * 3. Define the Task Priority in SystemDefines.hpp
 * 4. Replace all placeholders marked with a $ sign
 ********************************************************************************
 */

/************************************
 * INCLUDES
 ************************************/
#include <FilterTask.hpp>
#include "SystemDefines.hpp"
#include "DataBroker.hpp"
#include "SensorDataTypes.hpp"
#include "Command.hpp"
#include "everest.hpp"
#include "cmsis_os.h"
#include "FreeRTOS.h"
#include <cstring>

/************************************
 * PRIVATE MACROS AND DEFINES
 ************************************/

/************************************
 * VARIABLES
 ************************************/
// (no file-scope variables required)

/************************************
 * FUNCTION DECLARATIONS
 ************************************/

/************************************
 * FUNCTION DEFINITIONS
 ************************************/

/**
 * @brief Constructor for filterTask
 */
FilterTask::FilterTask() : Task(TASK_FILTER_QUEUE_DEPTH_OBJS), refreshMs_(500)
{
}

void FilterTask::SetRefreshMs(uint32_t ms)
{
    taskENTER_CRITICAL();
    refreshMs_ = ms;
    taskEXIT_CRITICAL();
}

uint32_t FilterTask::GetRefreshMs()
{
    taskENTER_CRITICAL();
    uint32_t v = refreshMs_;
    taskEXIT_CRITICAL();
    return v;
}

/**
 * @brief Initialize the filterTask
 *        Do not modify this function aside from adding the task name
 */
void FilterTask::InitTask()
{
    // Make sure the task is not already initialized
    SOAR_ASSERT(rtTaskHandle == nullptr, "Cannot initialize FilterTask task twice");

    BaseType_t rtValue =
        xTaskCreate((TaskFunction_t)FilterTask::RunTask,
            (const char*)"FilterTask",
            (uint16_t)TASK_FILTER_STACK_DEPTH_WORDS,
            (void*)this,
            (UBaseType_t)TASK_FILTER_PRIORITY,
            (TaskHandle_t*)&rtTaskHandle);

    SOAR_ASSERT(rtValue == pdPASS, "FilterTask::InitTask() - xTaskCreate() failed");
    // Subscribe to sensor data so this task receives DataBroker messages
    DataBroker::Subscribe<IMUData>(this);
    DataBroker::Subscribe<BaroData>(this);
    DataBroker::Subscribe<MagData>(this);
    DataBroker::Subscribe<GPSData>(this);
    DataBroker::Subscribe<TimeStampData>(this);
}

/**
 * @brief Instance Run loop for the Task, runs on scheduler start as long as the task is initialized.
 * @param pvParams RTOS Passed void parameters, contains a pointer to the object instance, should not be used
 */
void FilterTask::Run(void * pvParams)
{
    // Get singleton instance of EverestTask to process incoming data
    EverestTask &everest = EverestTask::getEverest();

    // track last run time so the filter block executes only at the configured refresh rate
    TickType_t lastRunTick = xTaskGetTickCount();
    uint32_t lastRunMs = static_cast<uint32_t>(lastRunTick * portTICK_PERIOD_MS);

    // store latest timestamp published by LoggingTask (if any)
    static uint32_t latestLoggedTimestampMs = 0;

    while (1) {

        /* Calculate time until next scheduled filter run */
        uint32_t refresh = GetRefreshMs();
        if (refresh == 0) refresh = 1;

        TickType_t now = xTaskGetTickCount();
        uint32_t timestamp = static_cast<uint32_t>(now * portTICK_PERIOD_MS);
        uint32_t elapsedMs = timestamp - lastRunMs;
        uint32_t waitMs = (elapsedMs >= refresh) ? 0 : (refresh - elapsedMs);
//       SOAR_PRINT("filterTask - Elapsed: %u ms, waiting for %u ms until next filter step\n", elapsedMs, waitMs);

        // Wait for either a command or the next scheduled run
        Command cm;
        bool res = qEvtQueue->Receive(cm, waitMs);

        // If a command arrived, handle it (but do NOT run the filter block here unless scheduled)
        if (res) {
            if (cm.GetCommand() == DATA_BROKER_COMMAND) {
                DataBrokerMessageTypes mt = DataBroker::getMessageType(cm);
                if (mt == DataBrokerMessageTypes::TIME_DATA) {
                    TimeStampData ts = DataBroker::ExtractData<TimeStampData>(cm);
                    latestLoggedTimestampMs = ts.timestamp_ms;
                } else {
                    everest.Extract(cm);

                    // If Everest not initialized yet, drive its tare/initialization
                    if (everest.everestInitialized == 0) {
                        // choose authoritative time if available
                        uint32_t nowMsLocal = latestLoggedTimestampMs != 0 ? latestLoggedTimestampMs : static_cast<uint32_t>(xTaskGetTickCount() * portTICK_PERIOD_MS);
                        everest.updateDeltaTime(static_cast<float>(nowMsLocal));
                        everest.initEverest();
                    }
                }
                cm.Reset();
            } else {
                HandleCommand(cm);
            }
        }

        // After waiting (either due to message or timeout), check if it's time to run the filter step
        now = xTaskGetTickCount();
        timestamp = static_cast<uint32_t>(now * portTICK_PERIOD_MS);
        elapsedMs = timestamp - lastRunMs;
        if (elapsedMs >= refresh) {
            // Prefer the authoritative timestamp from LoggingTask if available
            float currentTime = static_cast<float>(latestLoggedTimestampMs != 0 ? latestLoggedTimestampMs : timestamp);
            std::vector<float> halo = everest.QueueEverest(currentTime);

            // If HALO returned a valid state, publish FilterData
            if (!halo.empty()) {
                FilterData fd{};
                uint32_t ts = static_cast<uint32_t>(currentTime);
                static_assert(sizeof(float) == sizeof(uint32_t), "float must be 32-bit");
                memcpy(&fd.altPredicted, &halo[0], sizeof(float));
                if (halo.size() > 1) memcpy(&fd.veloPredicted, &halo[1], sizeof(float));
                if (halo.size() > 2) memcpy(&fd.accelPredicted, &halo[2], sizeof(float));
                fd.timePredicted = ts;

                SOAR_PRINT("filterTask - Publishing FilterData: Altitude: %f, Velocity: %f, Acceleration: %f, Time: %u\n",
                    halo[0], halo.size() > 1 ? halo[1] : 0.0f, halo.size() > 2 ? halo[2] : 0.0f, ts);

                DataBroker::Publish<FilterData>(&fd);
            }

            // reset lastRunTick to now so next interval measures from here
            // Only update when the filter actually produced output (halo non-empty)
            if (!halo.empty()) {
                lastRunTick = now;
                lastRunMs = timestamp;
            }
        }
    }
}

/**
 * @brief Handles a command
 * @param cm Command reference to handle
 */
void FilterTask::HandleCommand(Command& cm)
{
    // listen to CAN commands for state
    switch (cm.GetCommand()) {

    default:
        SOAR_PRINT("filterTask - Received Unsupported Command {%d}\n", cm.GetCommand());
        break;
    }

    //No matter what we happens, we must reset allocated data
    cm.Reset();
}
