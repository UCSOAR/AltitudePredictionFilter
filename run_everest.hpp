#ifndef RUN_EVEREST_HPP
#define RUN_EVEREST_HPP

// Minimal declarations so other translation units can start the test task
// without including the .cpp implementation.

#include "cmsis_os.h"

#define TESTING_BUILD

// RunEverestTask: provides a proper Task wrapper for the test publisher used by
// the filter during development. Call `StartRunEverestTask()` from system
// initialization to create the RTOS task

/**
 * @brief Starts the Everest data injection task.
 *
 * Creates a one-shot FreeRTOS task that publishes test sensor data
 * (from taberLaunch / baroData tables) to the DataBroker at
 * AVAILABLE_MEAS_REFRESH_MS intervals, using live FreeRTOS tick
 * timestamps. The task self-deletes when all data has been published.
 *
 * Call once from system init — before starting the RTOS scheduler.
 * FilterTask must already be initialized so it can receive the published data.
 */
void StartRunEverestInjection();

/**
 * @brief PC/desktop test entry point (not compiled into firmware).
 *
 * Runs the full tare → stationary → flight sequence synchronously,
 * feeding taberLaunch data directly into EverestTask without RTOS.
 * Only available when TESTING_BUILD or _WIN32 is defined.
 */
#if defined(TESTING_BUILD) || defined(_WIN32)
int main_test_run_everest();
#endif

#endif // RUN_EVEREST_HPP
