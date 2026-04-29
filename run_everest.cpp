
#define LOGMETRICS
#define MAX_LINE_LENGTH 1024

#include "everest.hpp"
#include "input_data.hpp"
#include "gpsData.hpp"
#include "DataBroker.hpp"
#include "SensorDataTypes.hpp"
#if defined(TESTING_BUILD) || defined(_WIN32)
#include <direct.h>
#endif
#include "cmsis_os.h"
#include "FreeRTOS.h"
#include "task.h"
#include "run_everest.hpp"

// Internal state and helpers (internal linkage)
namespace {
bool cycleSensors = 0;
bool hasGPS = 1;
bool hasIMU1 = 1;
bool hasIMU2 = 1;
bool hasMag = 1;
bool hasBaro = 0;
float deltaTime = 0.333f;
float stationaryTime = 0.0f;

int sensorThisCycle = 0;

float everest_time = 0.0f;
float timestamp = 0.0f;

float findClosestTime(float timestamp) {
  float altitude = 1000;
  for (int i = 0; i < gpsData1.size(); i++) {
    if (gpsData1[i][0] > timestamp) {
      altitude = gpsData1[i][1];
      break;
    }
  }
  return altitude;
}

int getSensorData(EverestTask* everest, int i, int stationary) {
  IMUData_Everest imuE1{};
  IMUData_Everest imuE2{};
  BarosData baro1{};
  BarosData baro2{};
  float gps = 0;

  float accelX = taberLaunch[i][1];
  float accelY = taberLaunch[i][2];
  float accelZ = taberLaunch[i][3];

  float gyroX = taberLaunch[i][4];
  float gyroY = taberLaunch[i][5];
  float gyroZ = taberLaunch[i][6];

  float magX = taberLaunch[i][7];
  float magY = taberLaunch[i][8];
  float magZ = taberLaunch[i][9];

  float pressure = baroData[i][1];

  if (stationary) {
    pressure = baroData[0][1];
    if (hasIMU1) {
      imuE1 = {everest_time, 0, 0, 0, 0, 0, 1, 0, 0, 0};
    } else {
      imuE1 = {everest_time, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    }

    if (hasIMU2) {
      imuE2 = {everest_time, 0, 0, 0, 0, 0, 1, 0, 0, 0};
    } else {
      imuE2 = {everest_time, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    }

    if (hasBaro) {
      baro1 = {everest_time, pressure, 0, 0};
      baro2 = {everest_time, pressure, 0, 0};
    } else {
      baro1 = {everest_time, 0, 0, 0};
      baro2 = {everest_time, 0, 0, 0};
    }

    gps = hasGPS ? findClosestTime(everest_time) : 0;

  } else {
    timestamp = taberLaunch[i][0];

    if (hasIMU1) {
      if (hasMag) {
        imuE1 = {timestamp, gyroX, gyroY, gyroZ, accelX, accelY, accelZ, magX, magY, magZ};
      } else {
        imuE1 = {timestamp, gyroX, gyroY, gyroZ, accelX, accelY, accelZ, 0, 0, 0};
      }
    } else {
      imuE1 = {timestamp, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    }

    if (hasIMU2) {
      if (hasMag) {
        imuE2 = {timestamp, gyroX, gyroY, gyroZ, accelX, accelY, accelZ, magX, magY, magZ};
      } else {
        imuE2 = {timestamp, gyroX, gyroY, gyroZ, accelX, accelY, accelZ, 0, 0, 0};
      }
    } else {
      imuE2 = {timestamp, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    }

    if (hasBaro) {
      baro1 = {timestamp, pressure, 0, 0};
      baro2 = {timestamp, pressure, 0, 0};
    } else {
      baro1 = {timestamp, 0, 0, 0};
      baro2 = {timestamp, 0, 0, 0};
    }

    gps = hasGPS ? findClosestTime(timestamp) : 0;
  }

  // publish measurements to DataBroker (simulate sensors)
  if (cycleSensors) {
   switch (sensorThisCycle) {
     case 0:
      if (!hasMag) {
        imuE1.magX = 0;
        imuE1.magY = 0;
        imuE1.magZ = 0;
      }
      if (!hasIMU1) {
        imuE1 = {
            timestamp, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        };
      }
      {
        // Publish IMU1 to DataBroker for subscribers (e.g., filterTask)
        IMUData imuPub{};
        imuPub.gyro.x = static_cast<int16_t>(imuE1.gyroX);
        imuPub.gyro.y = static_cast<int16_t>(imuE1.gyroY);
        imuPub.gyro.z = static_cast<int16_t>(imuE1.gyroZ);
        imuPub.accel.x = static_cast<int16_t>(imuE1.accelX);
        imuPub.accel.y = static_cast<int16_t>(imuE1.accelY);
        imuPub.accel.z = static_cast<int16_t>(imuE1.accelZ);
        imuPub.id = 0;
        SOAR_PRINT("run_everest - Publishing IMU1 at t=%f\n", imuE1.time);
        DataBroker::Publish<IMUData>(&imuPub);
      }
       sensorThisCycle++;
       break;
     case 1:
       if (hasIMU2) {
         imuE2 = {
             timestamp + (1.0f / 15.0f),  // (1/3) / 5
             gyroX,
             gyroY,
             gyroZ,
             accelX,
             accelY,
             accelZ,
             magX,
             magY,
             magZ,
         };
       }
       if (!hasMag) {
         imuE2.magX = 0;
         imuE2.magY = 0;
         imuE2.magZ = 0;
       }
       if (!hasIMU2) {
         imuE2 = {
             timestamp + (1.0f / 15.0f),  // (1/3) / 5
             0,
             0,
             0,
             0,
             0,
             0,
             0,
             0,
             0,
         };
       }
      {
        IMUData imuPub2{};
        imuPub2.gyro.x = static_cast<int16_t>(imuE2.gyroX);
        imuPub2.gyro.y = static_cast<int16_t>(imuE2.gyroY);
        imuPub2.gyro.z = static_cast<int16_t>(imuE2.gyroZ);
        imuPub2.accel.x = static_cast<int16_t>(imuE2.accelX);
        imuPub2.accel.y = static_cast<int16_t>(imuE2.accelY);
        imuPub2.accel.z = static_cast<int16_t>(imuE2.accelZ);
        imuPub2.id = 1;
        SOAR_PRINT("run_everest - Publishing IMU2 at t=%f\n", imuE2.time);
        DataBroker::Publish<IMUData>(&imuPub2);
      }
       sensorThisCycle++;
       break;
     case 2:
       if (hasBaro) {
         baro1 = {timestamp + (2.0f / 15.0f), pressure, 0, 0};
       } else {
         baro1 = {0, pressure, 0, 0};
       }
      {
        BaroData baroPub1{};
        baroPub1.pressure = static_cast<uint32_t>(baro1.pressure);
        baroPub1.temp = 0;
        baroPub1.id = 0;
        SOAR_PRINT("run_everest - Publishing BARO1 at t=%f pressure=%u\n", baro1.time, baro1.pressure);
        DataBroker::Publish<BaroData>(&baroPub1);
      }
       sensorThisCycle++;
       break;
     case 3:
       if (hasBaro) {
         baro2 = {timestamp + (3.0f / 15.0f), pressure, 0, 0};
       } else {
         baro2 = {timestamp + (3.0f / 15.0f), 0, 0, 0};
       }
      {
        BaroData baroPub2{};
        baroPub2.pressure = static_cast<uint32_t>(baro2.pressure);
        baroPub2.temp = 0;
        baroPub2.id = 1;
        SOAR_PRINT("run_everest - Publishing BARO2 at t=%f pressure=%u\n", baro2.time, baro2.pressure);
        DataBroker::Publish<BaroData>(&baroPub2);
      }
       sensorThisCycle++;
       break;
     case 4:
       gps = hasGPS ? findClosestTime(float(timestamp + (4.0f / 15.0f))) : 0;
      {
        GPSData gpsPub{};
        gpsPub.antennaAltitude_.altitude_ = static_cast<int32_t>(gps);
        SOAR_PRINT("run_everest - Publishing GPS altitude=%f\n", gps);
        DataBroker::Publish<GPSData>(&gpsPub);
      }
       sensorThisCycle = 0;
       break;
   }
 } else {
   if (hasIMU1) {
    IMUData imuPub{};
    imuPub.gyro.x = static_cast<int16_t>(imuE1.gyroX);
    imuPub.gyro.y = static_cast<int16_t>(imuE1.gyroY);
    imuPub.gyro.z = static_cast<int16_t>(imuE1.gyroZ);
    imuPub.accel.x = static_cast<int16_t>(imuE1.accelX);
    imuPub.accel.y = static_cast<int16_t>(imuE1.accelY);
    imuPub.accel.z = static_cast<int16_t>(imuE1.accelZ);
    imuPub.id = 0;
    SOAR_PRINT("run_everest - Publishing IMU1 at t=%f\n", imuE1.time);
    DataBroker::Publish<IMUData>(&imuPub);
   }
   if (hasIMU2) {
    IMUData imuPub2{};
    imuPub2.gyro.x = static_cast<int16_t>(imuE2.gyroX);
    imuPub2.gyro.y = static_cast<int16_t>(imuE2.gyroY);
    imuPub2.gyro.z = static_cast<int16_t>(imuE2.gyroZ);
    imuPub2.accel.x = static_cast<int16_t>(imuE2.accelX);
    imuPub2.accel.y = static_cast<int16_t>(imuE2.accelY);
    imuPub2.accel.z = static_cast<int16_t>(imuE2.accelZ);
    imuPub2.id = 1;
    SOAR_PRINT("run_everest - Publishing IMU2 at t=%f\n", imuE2.time);
    DataBroker::Publish<IMUData>(&imuPub2);
   }
   if (hasBaro) {
    BaroData baroPub1{};
    baroPub1.pressure = static_cast<uint32_t>(baro1.pressure);
    baroPub1.temp = 0;
    baroPub1.id = 0;
    BaroData baroPub2{};
    baroPub2.pressure = static_cast<uint32_t>(baro2.pressure);
    baroPub2.temp = 0;
    baroPub2.id = 1;
    SOAR_PRINT("run_everest - Publishing BAROs at t=%f\n", baro1.time);
    DataBroker::Publish<BaroData>(&baroPub1);
    DataBroker::Publish<BaroData>(&baroPub2);
   }
   if (hasGPS) {
     GPSData gpsPub{};
     gpsPub.antennaAltitude_.altitude_ = static_cast<int32_t>(gps);
    SOAR_PRINT("run_everest - Publishing GPS altitude=%f\n", gps);
     DataBroker::Publish<GPSData>(&gpsPub);
   }
 }
}

/**
* Serves to just initialize structs
*/
int main_test_run_everest() {
 EverestTask everest = EverestTask();
 // open files. moved here and out of madgwick setup.
#if defined(TESTING_BUILD)
 everest.openFiles();
#endif
 // read first line and preset the deltaTime to timestamp
 char line[MAX_LINE_LENGTH];
 std::clock_t start;
 float totalTime = 0;

 /**** TARE PHASE ****/
 while (!everest.everestInitialized) {
   getSensorData(&everest, 0, 1);

   everest_time += 0.333f;
   if (everest.everestInitialized == 1) {
     break;
   }
   if (everest.everestInitialized == 0) {
     everest.updateDeltaTime(everest_time);
     everest.initEverest();
   }
 }

 /**** STATIONARY PHASE ****/
 while (stationaryTime > 0 && everest_time <= stationaryTime) {
   // Tokenize the line using strtok
   // Parse accelerometer readings (X, Y, Z)
   everest_time += deltaTime;

   // get measurements (this will be done from subscribing to the sensors)
   getSensorData(&everest, 0, 1);
 }

 for (int i = 0; i < taberLaunch.size(); i++) {
   // start timer for iteration
   start = std::clock();

   if (cycleSensors) {  // If not using a time window, we poll and call everest
                        // as fast as possible. here we will poll a sensor then
                        // run 5 times every 1/3rd sec.
     for (int j = 0; j < 5; j++) {
       getSensorData(&everest, i, 0);
       std::vector<float> haloData =
           everest.QueueEverest(timestamp + ((float)j / 15.0f));
     }
   } else {  // if using a time window, we will collect 5 measurements in 1/3rd
             // sec then run the filter.
     getSensorData(&everest, i, 0);
     std::vector<float> haloData = everest.QueueEverest(timestamp);
   }

   clock_t endTime = std::clock();

   totalTime += endTime - start;

   if (i == taberLaunch.size() - 13) {
     std::cout << "Overall time:\t\t\t\t\t\t\t\t\t"
               << totalTime / (double)CLOCKS_PER_SEC << std::endl;
     break;
   }
 }

#if defined(LOGON) || defined(LOGMETRICS)
 // exit to close files. No idea if this is a good idea on a board.
 exit(0);
#endif

 return 0;
}

// RunEverestTask implementation (Task wrapper)
RunEverestTask& RunEverestTask::Inst()
{
    static RunEverestTask inst;
    return inst;
}

RunEverestTask::RunEverestTask()
    : Task(TASK_FILTER_QUEUE_DEPTH_OBJS)
{
}

void RunEverestTask::InitTask()
{
    SOAR_ASSERT(rtTaskHandle == nullptr, "Cannot initialize RunEverestTask twice");

    BaseType_t rtValue = xTaskCreate((TaskFunction_t)RunEverestTask::RunTask,
        (const char*)"runEverest",
        (uint16_t)TASK_FILTER_STACK_DEPTH_WORDS,
        (void*)this,
        (UBaseType_t)TASK_FILTER_PRIORITY,
        (TaskHandle_t*)&rtTaskHandle);

    SOAR_ASSERT(rtValue == pdPASS, "RunEverestTask::InitTask() - xTaskCreate() failed");
}

void RunEverestTask::RunTask(void* pvParams)
{
    RunEverestTask* self = static_cast<RunEverestTask*>(pvParams);
    if (self) self->Run(pvParams);
    // Should never return
    vTaskDelete(NULL);
}

void RunEverestTask::Run(void* pvParams)
{
    (void)pvParams;
    // simple periodic publisher based on deltaTime
    while (1) {
        for (int i = 0; i < taberLaunch.size(); i++) {
            getSensorData(nullptr, i, 0);
            uint32_t delayMs = static_cast<uint32_t>(deltaTime * 1000.0f);
            if (delayMs == 0) delayMs = 1;
            vTaskDelay(pdMS_TO_TICKS(delayMs));
        }
    }
}

// Start the run_everest RTOS task. Call this from system init to start publishing test data.
void StartRunEverestTask()
{
    RunEverestTask::Inst().InitTask();
}
