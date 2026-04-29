#ifndef RUN_EVEREST_HPP
#define RUN_EVEREST_HPP

// Minimal declarations so other translation units can start the test task
// without including the .cpp implementation.

#include "cmsis_os.h"

// RunEverestTask: provides a proper Task wrapper for the test publisher used by
// the filter during development. Call `StartRunEverestTask()` from system
// initialization to create the RTOS task.

class RunEverestTask;

class RunEverestTask {
public:
	static RunEverestTask& Inst();
	void InitTask();

private:
	RunEverestTask();
	void Run(void* pvParams);
	static void RunTask(void* pvParams);
};

// Convenience wrapper called from `main_system.cpp` to start the task.
void StartRunEverestTask();

// Optional: lightweight entry used for host unit testing.
int main_test_run_everest();

#endif // RUN_EVEREST_HPP
