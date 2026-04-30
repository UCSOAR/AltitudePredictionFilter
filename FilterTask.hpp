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

/************************************
 * TYPEDEFS
 ************************************/

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
};

/************************************
 * FUNCTION DECLARATIONS
 ************************************/

#endif /* ALTITUDEPREDICTIONFILTER_FILTERTASK_HPP_ */
