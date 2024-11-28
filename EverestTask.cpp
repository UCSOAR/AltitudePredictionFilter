// Altitude estimation using multiple sensors
#include "everestTaskHPP.hpp"
#include <stdio.h>
#include <ctime>
#include <string.h>

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Data.cpp"
#include "input_data.cpp"

#define LOGON
#define TIMERON

// #define printf(...) ;

// task specific defines
// #include "main.h"
// #include "Data.h"
// #include "DebugTask.hpp"
// #include "Task.hpp"
// #include "DMBProtocolTask.hpp"
// #include "TelemetryMessage.hpp"
// #include "FlashTask.hpp"
// #include <string.h>

/**
 * To run:  g++ Infusion.cpp EverestTask.cpp -o Everest
 *          ./Everest
 */

using namespace std;

// SETTINGS (mostly for debugging, keep default for run)
enum debug_level {
  RAW = 0,        // raw data
  Secondary = 1,  // all operations before dynamite
  Dynamite = 2,   // everything during dynamite
  Third = 3,      // after dynamite
  ALL = 4,        // all
  NONE = 5,       // none
  HAL0 = 6,       // HALO
};
bool isTared = false;
debug_level debug = ALL;
bool firstSampleAfterCalibration = true;
bool useSTD = false;

// INTERNAL VARIABLES
// double timeInSeconds = 2;
double theTime = CALIBRATION_TIME * RATE_BARO;
double sum = 0;
double pressureSum = 0;
static float previousTimestamp = 0;
std::vector<double> sumZeroOffsetAccel;
std::vector<double> sumZeroOffsetAccel2;
std::vector<double> sumZeroOffsetGyro;
std::vector<double> sumZeroOffsetGyro2;

// Instantiate Everest
madAhrs* ahrs;
Infusion* infusion;

EverestTask everest = EverestTask::getEverest();
kinematics* Kinematics = everest.getKinematics();  // tare to ground

madAhrsFlags flags;
madAhrsInternalStates internalStates;
EverestData everestData;

//----------------------------------Task
// Integration----------------------------------//
// EverestTask::EverestTask() : Task(TASK_EVEREST_QUEUE_DEPTH_OBJS)
// {
//     // Initialize the task
//     MadgwickSetup();
//     everestData = (EverestData*)soar_malloc(sizeof(EverestData));
// }

// /**
//  * @brief Creates a task for the FreeRTOS Scheduler
//  */
// void EverestTask::InitTask()
// {
//     // Make sure the task is not already initialized
//     SOAR_ASSERT(rtTaskHandle == nullptr, "Cannot initialize Everest task
//     twice");

//     // Start the task
//     BaseType_t rtValue =
//         xTaskCreate((TaskFunction_t)EverestTask::RunTask,
//             (const char*)"EverestTask",
//             (uint16_t)TASK_EVEREST_TASK_STACK_DEPTH_WORDS,
//             (void*)this,
//             (UBaseType_t)TASK_EVEREST_TASK_PRIORITY,
//             (TaskHandle_t*)&rtTaskHandle);

//     //Ensure creation succeded
//     SOAR_ASSERT(rtValue == pdPASS, "EverestTask::InitTask() - xTaskCreate()
//     failed");
// }

// /**
//  * @brief KalmanFilterTask run loop
//  * @param pvParams Currently unused task context
//  */
// void EverestTask::Run(void* pvParams)
// {

//     //Task run loop
//     while (1) {

//         Command cm;

//         //Wait forever for a command
//         qEvtQueue->ReceiveWait(cm);

//         //Process the command
//         HandleCommand(cm);

//         cm.Reset();
//     }

// }

// /**
//  * @brief Handles a command
//  * @param cm Command reference to handle
//  */
// void EverestTask::HandleCommand(Command& cm)
// {
//     //TODO: Since this task will stall for a few milliseconds, we may need a
//     way to eat the whole queue (combine similar eg. REQUEST commands and eat
//     to WDG command etc)
//     //TODO: Maybe a HandleEvtQueue instead that takes in the whole queue and
//     eats the whole thing in order of non-blocking to blocking

//     //Switch for the GLOBAL_COMMAND
//     switch (cm.GetCommand()) {
//     case TASK_SPECIFIC_COMMAND: {
//     	if(cm.GetTaskCommand()==COPY_DATA){
//     		*everestData = *(EverestData*) cm.GetDataPointer();
//     	}else{
//             *everestData = *(EverestData*) cm.GetDataPointer();
// 			HandleRequestCommand(cm.GetTaskCommand());
//     	}
//         break;
//     }
//     default:
//         SOAR_PRINT("EverestTask - Received Unsupported Command {%d}\n",
//         cm.GetCommand()); break;
//     }

//     //No matter what happens, we must reset allocated data
//     cm.Reset();
// }

// /**
//  * @brief Handles a Request Command
//  * @param taskCommand The command to handle
//  */
// void EverestTask::HandleRequestCommand(uint16_t taskCommand)
// {
//     //Switch for task specific command within DATA_COMMAND
//     switch (taskCommand) {
//     case UPDATE:
//         TaskWrapper(everestData, MadAxesAlignmentPXPYNZ,
//         MadAxesAlignmentPXPYNZ); break;
//     case TEST:
//     	break;
//     case RETARE:
//         setIsTared(false);
//     default:
//         SOAR_PRINT("EverestTask - Received Unsupported REQUEST_COMMAND
//         {%d}\n", taskCommand); break;
//     }
// }

/**
 * @brief Calls finalWrapper with data and alignment
 * IMPORTANT: PASS 0s FOR NOT UPDATED MEASUREMENTS (BAROS)
 */
double EverestTask::TaskWrapper(EverestData everestData,
                                MadAxesAlignment alignment,
                                MadAxesAlignment alignment2) {
  return this->finalWrapper(
      everestData.accelX1, everestData.accelY1, everestData.accelZ1,
      everestData.gyroX1, everestData.gyroY1, everestData.gyroZ1,
      everestData.magX1, everestData.magY1, everestData.magZ1,
      everestData.accelX2, everestData.accelY2, everestData.accelZ2,
      everestData.gyroX2, everestData.gyroY2, everestData.gyroZ2,
      everestData.magX2, everestData.magY2, everestData.magZ2,
      everestData.pressure1, everestData.pressure2, everestData.timeIMU1,
      everestData.timeIMU2, everestData.timeBaro1, everestData.timeBaro2,
      alignment, alignment2);
}

//----------------------------------EVEREST-------------------------------------------//
/**
 * @brief Only done once. Sets pointers for Madgwick
 *     Internal
 */
void EverestTask::MadgwickSetup() {
  // Attaches Madgwick to Everest
  infusion = everest.ExternalInitialize();
  ahrs = infusion->getMadAhrs();

  // calculateSTDCoefficients();

  // Define calibration (replace with actual calibration data if available)
  const madMatrix gyroscopeMisalignment = {1.0f, 0.0f, 0.0f, 0.0f, 1.0f,
                                           0.0f, 0.0f, 0.0f, 1.0f};
  const madVector gyroscopeSensitivity = {1.0f, 1.0f, 1.0f};
  const madVector gyroscopeOffset = {0.0f, 0.0f, 0.0f};
  const madMatrix accelerometerMisalignment = {1.0f, 0.0f, 0.0f, 0.0f, 1.0f,
                                               0.0f, 0.0f, 0.0f, 1.0f};
  const madVector accelerometerSensitivity = {1.0f, 1.0f, 1.0f};
  const madVector accelerometerOffset = {0.0f, 0.0f, 0.0f};
  const madMatrix softIronMatrix = {1.0f, 0.0f, 0.0f, 0.0f, 1.0f,
                                    0.0f, 0.0f, 0.0f, 1.0f};

  internalStates = infusion->madAhrsGetInternalStates(ahrs);
  flags = infusion->madAhrsGetFlags(ahrs);

  const madVector hardIronOffset = {0.0f, 0.0f, 0.0f};

  // Initialise algorithms
  madOffset offset = infusion->getOffset();

  infusion->madOffsetInitialise(&offset, SAMPLE_RATE);
  infusion->madAhrsInitialise(ahrs);

  // Set AHRS algorithm settings
  madAhrsSettings settings = {
      EarthConventionEnu,
      0.5f,
      2000.0f, /* replace this with actual gyroscope range in degrees/s */
      10.0f,
      10.0f,
      5 * SAMPLE_RATE, /* 5 seconds */
  };

  infusion->madAhrsSetSettings(ahrs, &settings);
}

/**
 * @brief Wrapper for Madgwick, does offset calc and passes
 *       data to Madgwick
 *
 *      Internal
 * @param data IMUData struct
 *
 */
void EverestTask::MadgwickWrapper(IMUData data) {
  // Infusion infusion = infusion;
  const float timestamp = data.time;
  madVector gyroscope = {
      data.gyroX, data.gyroY,
      data.gyroZ};  // replace this with actual gyroscope data in degrees/s
  madVector accelerometer = {
      data.accelX, data.accelY,
      data.accelZ};  // replace this with actual accelerometer data in g
  madVector mag = {data.magX, data.magY, data.magZ};

  // Update gyroscope offset correction algorithm
  madOffset offset = infusion->getOffset();
  gyroscope = infusion->madOffsetUpdate(&offset, gyroscope);

  // printf("Roll %0.3f, Pitch %0.3f, Yaw %0.3f, X %0.3f, Y %0.3f, Z %0.3f\n",
  //        euler.angle.roll, euler.angle.pitch, euler.angle.yaw,
  //        earth.axis.x, earth.axis.y, earth.axis.z);

  // Calculate delta time (in seconds) to account for gyroscope sample clock
  // error static float previousTimestamp;
  float deltaTime = (float)(timestamp - previousTimestamp);
  previousTimestamp = timestamp;

  this->state.deltaTimeIMU = deltaTime;

  if (debug == Secondary || debug == ALL) {
    // SOAR_PRINT("Averaged: (%.6f, %.6f, %.6f) deg/s, Accel: (%.6f, %.6f, %.6f)
    // g Time: %f\n",
    //     data.gyroX, data.gyroY, data.gyroZ, data.accelX, data.accelY,
    //     data.accelZ, deltaTime);
    // printf("Mag: (%.6f, %.6f, %.6f) uT\n", mag.axis.x, mag.axis.y,
    // mag.axis.z);
  }

  // Update gyroscope AHRS algorithm
  // infusion->madAhrsUpdateNoMagnetometer(ahrs, gyroscope, accelerometer,
  // deltaTime);
  infusion->madAhrsUpdate(ahrs, gyroscope, accelerometer, mag, deltaTime);

  // madAhrsInternalStates internal;
  // madAhrsFlags flags;

  madEuler euler = infusion->getEuler(ahrs);
  madVector earth = infusion->madAhrsGetEarthAcceleration(ahrs);

  internalStates = infusion->madAhrsGetInternalStates(infusion->getMadAhrs());
  flags = infusion->madAhrsGetFlags(infusion->getMadAhrs());

  // write to file
  //    fprintf(file, "%f,", timestamp);
  //
  //    fprintf(file, "%f,%f,%f,", euler.angle.roll, euler.angle.pitch,
  //    euler.angle.yaw);
  //
  //    fprintf(file, "%f,%d,%.0f,%.0f,%d,%.0f,%d,%d,%d,%d,%f",
  //    internalStates.accelerationError, internalStates.accelerometerIgnored,
  //    internalStates.accelerationRecoveryTrigger,
  //    internalStates.magneticError, internalStates.magnetometerIgnored,
  //    internalStates.magneticRecoveryTrigger, flags.initialising,
  //    flags.angularRateRecovery, flags.accelerationRecovery,
  //    flags.magneticRecovery, earth.axis.z);

  everest.state.earthAcceleration = earth.axis.z;

  // fprintf(file, "\n");

  //    if(debug == Secondary || debug == ALL){
  //        SOAR_PRINT("%f,%d,%.0f,%.0f,%d,%.0f,%d,%d,%d,%d\n",
  //        internalStates.accelerationError,
  //        internalStates.accelerometerIgnored,
  //        internalStates.accelerationRecoveryTrigger,
  //        internalStates.magneticError, internalStates.magnetometerIgnored,
  //        internalStates.magneticRecoveryTrigger, flags.initialising,
  //        flags.angularRateRecovery, flags.accelerationRecovery,
  //        flags.magneticRecovery);
  //    }
}

/**
 * @brief Averages IMUs and feeds them to Madgwick wrapper
 *      Should be called every time IMU data is updated
 *
 *    Internal
 */
void EverestTask::IMU_Update(const IMUData& imu1, const IMUData& imu2) {
  int numberOfSamples = 2;
  // Update IMU1
  this->internalIMU_1.time = imu1.time;
  this->internalIMU_1.gyroX = imu1.gyroX;
  this->internalIMU_1.gyroZ = imu1.gyroZ;
  this->internalIMU_1.gyroY = imu1.gyroY;

  this->internalIMU_1.accelX = imu1.accelX;
  this->internalIMU_1.accelY = imu1.accelY;
  this->internalIMU_1.accelZ = imu1.accelZ;

  this->internalIMU_1.magX = imu1.magX;
  this->internalIMU_1.magY = imu1.magY;
  this->internalIMU_1.magZ = imu1.magZ;

  // Update IMU2
  this->internalIMU_2.time = imu2.time;

  this->internalIMU_2.gyroX = imu2.gyroX;
  this->internalIMU_2.gyroY = imu2.gyroY;
  this->internalIMU_2.gyroZ = imu2.gyroZ;

  this->internalIMU_2.accelX = imu2.accelX;
  this->internalIMU_2.accelY = imu2.accelY;
  this->internalIMU_2.accelZ = imu2.accelZ;

  this->internalIMU_2.magX = imu2.magX;
  this->internalIMU_2.magY = imu2.magY;
  this->internalIMU_2.magZ = imu2.magZ;

  if (isinf(internalIMU_1.accelX)) {
    numberOfSamples -= 1;
    // printf("%d",numberOfSamples);
    this->internalIMU_1.gyroX = 0;
    this->internalIMU_1.gyroY = 0;
    this->internalIMU_1.gyroZ = 0;

    this->internalIMU_1.accelX = 0;
    this->internalIMU_1.accelY = 0;
    this->internalIMU_1.accelZ = 0;

    this->internalIMU_1.magX = 0;
    this->internalIMU_1.magY = 0;
    this->internalIMU_1.magZ = 0;
  } else {
    // Apply calibration
    // printf("uncalibrated: %f, %f, %f ", this->internalIMU_1.accelX,
    // this->internalIMU_1.accelY, this->internalIMU_1.accelZ);

    this->internalIMU_1.accelX =
        this->internalIMU_1.accelX - this->zeroOffsetAccel[0];
    this->internalIMU_1.accelY =
        this->internalIMU_1.accelY - this->zeroOffsetAccel[1];
    this->internalIMU_1.accelZ =
        this->internalIMU_1.accelZ - this->zeroOffsetAccel[2];

    // printf("-> offset (%f, %f, %f) = calibrated accel (%f, %f, %f)\n",
    // this->zeroOffsetAccel[0], this->zeroOffsetAccel[1],
    // this->zeroOffsetAccel[2], this->internalIMU_1.accelX,
    // this->internalIMU_1.accelY, this->internalIMU_1.accelZ);

    // printf("uncalibrated: %f, %f, %f ", this->internalIMU_1.gyroX,
    // this->internalIMU_1.gyroY, this->internalIMU_1.gyroZ);

    this->internalIMU_1.gyroX =
        this->internalIMU_1.gyroX - this->zeroOffsetGyro[0];
    this->internalIMU_1.gyroY =
        this->internalIMU_1.gyroY - this->zeroOffsetGyro[1];
    this->internalIMU_1.gyroZ =
        this->internalIMU_1.gyroZ - this->zeroOffsetGyro[2];

    // printf("-> offset (%f, %f, %f) = calibrated gyro (%f, %f, %f)\n",
    // this->zeroOffsetGyro[0], this->zeroOffsetGyro[1],
    // this->zeroOffsetGyro[2], this->internalIMU_1.gyroX,
    // this->internalIMU_1.gyroY, this->internalIMU_1.gyroZ);
  }

  if (isinf(internalIMU_2.accelX)) {
    numberOfSamples -= 1;
    // printf("%d",numberOfSamples);
    this->internalIMU_2.gyroX = 0;
    this->internalIMU_2.gyroY = 0;
    this->internalIMU_2.gyroZ = 0;

    this->internalIMU_2.accelX = 0;
    this->internalIMU_2.accelY = 0;
    this->internalIMU_2.accelZ = 0;

    this->internalIMU_2.magX = 0;
    this->internalIMU_2.magY = 0;
    this->internalIMU_2.magZ = 0;
  } else {
    // Apply calibration
    // printf("uncalibrated: %f, %f, %f ", this->internalIMU_2.accelX,
    // this->internalIMU_2.accelY, this->internalIMU_2.accelZ);

    this->internalIMU_2.accelX =
        this->internalIMU_2.accelX - this->zeroOffsetAccel2[0];
    this->internalIMU_2.accelY =
        this->internalIMU_2.accelY - this->zeroOffsetAccel2[1];
    this->internalIMU_2.accelZ =
        this->internalIMU_2.accelZ - this->zeroOffsetAccel2[2];

    // printf("-> offset (%f, %f, %f) = calibrated accel2 (%f, %f, %f)\n",
    // this->zeroOffsetAccel2[0], this->zeroOffsetAccel2[1],
    // this->zeroOffsetAccel2[2], this->internalIMU_2.accelX,
    // this->internalIMU_2.accelY, this->internalIMU_2.accelZ);

    // printf("uncalibrated: %f, %f, %f ", this->internalIMU_2.gyroX,
    // this->internalIMU_2.gyroY, this->internalIMU_2.gyroZ);

    this->internalIMU_2.gyroX =
        this->internalIMU_2.gyroX - this->zeroOffsetGyro2[0];
    this->internalIMU_2.gyroY =
        this->internalIMU_2.gyroY - this->zeroOffsetGyro2[1];
    this->internalIMU_2.gyroZ =
        this->internalIMU_2.gyroZ - this->zeroOffsetGyro2[2];

    // printf("-> offset (%f, %f, %f) = calibrated gyro2 (%f, %f, %f)\n",
    // this->zeroOffsetGyro[0], this->zeroOffsetGyro[1],
    // this->zeroOffsetGyro[2], this->internalIMU_2.gyroX,
    // this->internalIMU_2.gyroY, this->internalIMU_2.gyroZ);
  }

// Calculate average of IMU parameters
#define averageIMU this->state.avgIMU

  averageIMU.gyroX =
      (this->internalIMU_1.gyroX + this->internalIMU_2.gyroX) / numberOfSamples;
  averageIMU.gyroY =
      (this->internalIMU_1.gyroY + this->internalIMU_2.gyroY) / numberOfSamples;
  averageIMU.gyroZ =
      (this->internalIMU_1.gyroZ + this->internalIMU_2.gyroZ) / numberOfSamples;

  averageIMU.accelX =
      (this->internalIMU_1.accelX + this->internalIMU_2.accelX) /
      numberOfSamples;
  averageIMU.accelY =
      (this->internalIMU_1.accelY + this->internalIMU_2.accelY) /
      numberOfSamples;
  averageIMU.accelZ =
      (this->internalIMU_1.accelZ + this->internalIMU_2.accelZ) /
      numberOfSamples;

  averageIMU.magX =
      (this->internalIMU_1.magX + this->internalIMU_2.magX) / numberOfSamples;
  averageIMU.magY =
      (this->internalIMU_1.magY + this->internalIMU_2.magY) / numberOfSamples;
  averageIMU.magZ =
      (this->internalIMU_1.magX + this->internalIMU_2.magZ) / numberOfSamples;

  averageIMU.time =
      (this->internalIMU_1.time + this->internalIMU_2.time) / numberOfSamples;

#undef averageIMU

  if (numberOfSamples == 0) {
    this->state.avgIMU = {imu1.time, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  }
  // feed to Madgwick
  // MadgwickWrapper(state.avgIMU);
  this->MadgwickWrapper(state.avgIMU);
}

/**
 * @brief updates baro and delta time
 *     Should be called every time baro data is updated
 *      Internal
 */
void EverestTask::Baro_Update(const BarosData& Baro1, const BarosData& Baro2) {
  // Update Baros
  this->baro1.time = Baro1.time;
  this->baro1.pressure = Baro1.pressure;
  this->baro1.deltaTime = Baro1.time - this->baro1.previousTime;
  this->baro1.previousTime = Baro1.time;

  this->baro2.time = Baro2.time;
  this->baro2.pressure = Baro2.pressure;
  this->baro2.deltaTime = Baro2.time - this->baro2.previousTime;
  this->baro2.previousTime = Baro2.time;

  if (debug == RAW || debug == ALL) {
    // SOAR_PRINT("Baro1: %.f Pa, Baro2: %.f Pa, Baro3: %.f Pa, RealBaro: %.f
    // m\n",
    //     baro1.pressure, baro2.pressure, baro3.pressure, realBaro.altitude);
  }
}

/**
 * @brief Calls IMU and Baro update functions and calculates altitude
 *      calls Dynamite and updates altitude list
 *
 * @return calculated altitude
 *
 *    External (only function that should be called after instantiation of
 * Everest to pass sensor data to Everest for altitude calculation)
 */
double EverestTask::ExternalUpdate(IMUData imu1, IMUData imu2, BarosData baro1,
                                   BarosData baro2) {
  if (!isTared) {
    // printf("Taring in progress\n");
    everest.tare(imu1, imu2, baro1, baro2);
    return 0;
  }

  everest.IMU_Update(imu1, imu2);

  if (debug == Third || debug == ALL) {
    // SOAR_PRINT("After IMU Update IMU Altitude: %f\n",
    // everest.state.avgIMU.altitude);
  }

  everest.Baro_Update(baro1, baro2);

  double finalAlt = everest.dynamite();

  //    if(debug == Dynamite || debug == ALL){
  //        printf("After Dynamite: %f\n", finalAlt);
  //    }

  // Update altitude list
  this->AltitudeList.secondLastAltitude = this->AltitudeList.lastAltitude;
  this->AltitudeList.lastAltitude = finalAlt;

  //    fprintf(file, ",%f\n", finalAlt); // write to file

  return finalAlt;
}

/**
 * @brief (Currently does not work, use final wrapper) Wraps External Update
 * with alignment, returns External Update with aligned data
 */
double EverestTask::AlignedExternalUpdate(IMUData imu1, IMUData imu2,
                                          BarosData baro1, BarosData baro2,
                                          MadAxesAlignment alignment) {
  // align
  madVector alignedIMU1 =
      infusion->AxesSwitch({imu1.accelX, imu1.accelY, imu1.accelZ}, alignment);
  madVector alignedIMUGyro1 =
      infusion->AxesSwitch({imu1.gyroX, imu1.gyroY, imu1.gyroZ}, alignment);

  madVector alignedIMU2 =
      infusion->AxesSwitch({imu2.accelX, imu2.accelY, imu2.accelZ}, alignment);
  madVector alignedIMUGyro2 =
      infusion->AxesSwitch({imu2.gyroX, imu2.gyroY, imu2.gyroZ}, alignment);

  if (debug == Secondary || debug == ALL) {
    // SOAR_PRINT("Unaligned IMU1: (%.6f, %.6f, %.6f) g, (%.6f, %.6f, %.6f)
    // deg/s\n",
    //     imu1.accelX, imu1.accelY, imu1.accelZ, imu1.gyroX, imu1.gyroY,
    //     imu1.gyroZ);

    // SOAR_PRINT("Unaligned IMU2: (%.6f, %.6f, %.6f) g, (%.6f, %.6f, %.6f)
    // deg/s\n",
    //     imu2.accelX, imu2.accelY, imu2.accelZ, imu2.gyroX, imu2.gyroY,
    //     imu2.gyroZ);

    // SOAR_PRINT("Alignment: %d\n", alignment);
  }

  // put aligned data into IMUData struct
  imu1.accelX = alignedIMU1.axis.x;
  imu1.accelY = alignedIMU1.axis.y;
  imu1.accelZ = alignedIMU1.axis.z;

  imu1.gyroX = alignedIMUGyro1.axis.x;
  imu1.gyroY = alignedIMUGyro1.axis.y;
  imu1.gyroZ = alignedIMUGyro1.axis.z;

  // IMU 2
  imu2.accelX = alignedIMU2.axis.x;
  imu2.accelY = alignedIMU2.axis.y;
  imu2.accelZ = alignedIMU2.axis.z;

  imu2.gyroX = alignedIMUGyro2.axis.x;
  imu2.gyroY = alignedIMUGyro2.axis.y;
  imu2.gyroZ = alignedIMUGyro2.axis.z;

  if (debug == Secondary || debug == ALL) {
    // SOAR_PRINT("Aligned IMU1: (%.6f, %.6f, %.6f) g, (%.6f, %.6f, %.6f)
    // deg/s\n",
    //     imu1.accelX, imu1.accelY, imu1.accelZ, imu1.gyroX, imu1.gyroY,
    //     imu1.gyroZ);

    // SOAR_PRINT("Aligned IMU2: (%.6f, %.6f, %.6f) g, (%.6f, %.6f, %.6f)
    // deg/s\n",
    //     imu2.accelX, imu2.accelY, imu2.accelZ, imu2.gyroX, imu2.gyroY,
    //     imu2.gyroZ);
  }

  return ExternalUpdate(imu1, imu2, baro1, baro2);
}

/**
 * Calculates altitude using IMU sensor data and kinematic equations.
 *
 * @param avgIMU with the average sensor data from the IMU
 *
 * @category Internal | ASYNCHRONOUS
 *
 * @return calculated altitude
 */
double EverestTask::deriveForAltitudeIMU(IMUData avgIMU) {
  // double accelerationZ = avgIMU.accelX * -9.81;
  double accelerationZ = everest.state.earthAcceleration * -9.81;
  double initialVelocity = this->getKinematics()->initialVelo;
  double initialAltitude = this->Kinematics.initialAlt;
  double deltaTime = this->state.deltaTimeIMU;

  // Derive altitude from IMU
  double finalVelocity = initialVelocity + accelerationZ * deltaTime;
  double altitude =
      initialAltitude + (initialVelocity + finalVelocity) * deltaTime / 2.0;

  if (debug == Secondary || debug == ALL) {
    // SOAR_PRINT("\nKinematics\n");
    // SOAR_PRINT("IMU Initial Altitude: %f\n", initialAltitude);
    // SOAR_PRINT("IMU Velocity: %f\n", initialVelocity);
    // SOAR_PRINT("IMU Acceleration: %f\n", accelerationZ);
    // SOAR_PRINT("IMU Delta Time: %f\n", deltaTime);
    // SOAR_PRINT("Derived Altitude: %f\n", altitude);
  }

  // return altitude
  return altitude;
}

/**
 * Calculates the altitude based on the given pressure using the barometric
 * formula
 *
 * @param pressure pressure from baros in Pa
 *
 * @return altitude in meters
 *
 * @category Internal
 */
double convertToAltitude(double pressure) {
  double seaLevelPressure = 1013.25;  // sea level pressure in hPa
  pressure = pressure / 100.0;        // convert to hPa
  double altitude = 44330.0 * (1.0 - pow(pressure / seaLevelPressure,
                                         1 / 5.2558));  // barometric formula

  // If pressure is less than 100, altitude is 0
  if (pressure < 100) {
    altitude = 0;
  }

  if (debug == Dynamite || debug == ALL) {
    // SOAR_PRINT("\nConversion \n");
    // SOAR_PRINT("Pressure: %.f hPa, Altitude: %.f m\n", pressure, altitude);
  }

  return altitude;
}

/**
 * @brief Multi-system trust algorithm. Assumes measurements are updated
 * @returns normalised altitude
 *
 * @category Internal | Asynchronous
 */
double EverestTask::dynamite() {
  double IMUAltitude = deriveForAltitudeIMU(everest.state.avgIMU);
  this->state.avgIMU.altitude = IMUAltitude;

  double BaroAltitude1 = convertToAltitude(this->baro1.pressure);
  this->baro1.altitude = BaroAltitude1;

  double BaroAltitude2 = convertToAltitude(everest.baro2.pressure);
  this->baro2.altitude = BaroAltitude2;

  if (debug == Dynamite || debug == ALL) {
    // SOAR_PRINT("\nDynamite\n");
    // SOAR_PRINT("Baro1 Altitude: %f\n", BaroAltitude1);
    // SOAR_PRINT("Baro2 Altitude: %f\n", BaroAltitude2);
    // SOAR_PRINT("Baro3 Altitude: %f\n", BaroAltitude3);
    // SOAR_PRINT("Real Baro Altitude: %f\n", RealBaroAltitude);
    // SOAR_PRINT("IMU Altitude: %f\n", IMUAltitude);
  }

  // // distributing measurement
  // double distributed_IMU_Altitude = (IMUAltitude *
  // everest.state.gain_IMU)/pow(everest.state.std_IMU, 2); double
  // distributed_Baro_Altitude1 = (BaroAltitude1 *
  // everest.state.gain_Baro1)/pow(everest.state.std_Baro1, 2); double
  // distributed_Baro_Altitude2 = (BaroAltitude2 *
  // everest.state.gain_Baro2)/pow(everest.state.std_Baro2, 2); double
  // distributed_Baro_Altitude3 = (BaroAltitude3 *
  // everest.state.gain_Baro3)/pow(everest.state.std_Baro3,2); double
  // distributed_RealBaro_Altitude = (RealBaroAltitude *
  // everest.state.gain_Real_Baro)/pow(everest.state.std_Real_Baro,2);

  // if pressure is zero, set gain to zero
  if (everest.baro1.pressure == 0) {
    everest.state.gain_Baro1 = 0;
  } else if (everest.state.gain_Baro1 == 0) {
    // if not zero, set gain to previous gain
    everest.state.gain_Baro1 = everest.state.prev_gain_Baro1;
  }

  // if not zero, set gain to zero
  if (everest.baro2.pressure == 0) {
    everest.state.gain_Baro2 = 0;
  } else if (everest.state.gain_Baro2 == 0) {
    // if not zero, set gain to previous gain
    everest.state.gain_Baro2 = everest.state.prev_gain_Baro2;
  }

  // distribute measurements based on gain
  double distributed_IMU_Altitude = IMUAltitude * everest.state.gain_IMU;
  double distributed_Baro_Altitude1 =
      (BaroAltitude1 * everest.state.gain_Baro1);
  double distributed_Baro_Altitude2 =
      (BaroAltitude2 * everest.state.gain_Baro2);

  if (debug == Dynamite || debug == ALL) {
    // SOAR_PRINT("\nDistributed\n");
    // SOAR_PRINT("Distributed IMU Altitude: %f\n", distributed_IMU_Altitude);
    // SOAR_PRINT("Gain IMU: %f\n", everest.state.gain_IMU);
    // SOAR_PRINT("STD IMU: %f\n\n", everest.state.std_IMU);

    // SOAR_PRINT("Distributed Baro1 Altitude: %f\n",
    // distributed_Baro_Altitude1); SOAR_PRINT("Gain Baro1: %f\n",
    // everest.state.gain_Baro1); SOAR_PRINT("STD Baro1: %f\n\n",
    // everest.state.std_Baro1);

    // SOAR_PRINT("Distributed Baro2 Altitude: %f\n",
    // distributed_Baro_Altitude2); SOAR_PRINT("Gain Baro2: %f\n",
    // everest.state.gain_Baro2); SOAR_PRINT("STD Baro2: %f\n\n",
    // everest.state.std_Baro2);
  }

  // summation of distributed measurements
  double distributed_Sum = distributed_IMU_Altitude +
                           distributed_Baro_Altitude1 +
                           distributed_Baro_Altitude2;

  if (debug == Dynamite || debug == ALL) {
    // SOAR_PRINT("Distributed Sum: %f\n\n", distributed_Sum);
  }

  //    if(debug == Dynamite || debug == ALL){
  //        // printf("Sum STD: %f\n\n", sumSTD1);
  //    }

  // summation of gains
  double sumGain = everest.state.gain_IMU + everest.state.gain_Baro1 +
                   everest.state.gain_Baro2;

  if (debug == Dynamite || debug == ALL) {
    // SOAR_PRINT("Sum Gain: %f\n\n", sumGain);
  }

  // normalised altitude
  double normalised_Altitude = (distributed_Sum) / sumGain;

  if (debug == Dynamite || debug == ALL) {
    // SOAR_PRINT("Normalised Altitude: %f\n\n", normalised_Altitude);
  }

  // Update Kinematics
  Kinematics.finalAltitude = normalised_Altitude;

  if (debug == Dynamite || debug == ALL) {
    // SOAR_PRINT("Final Altitude: %f\n\n", Kinematics.finalAltitude);
  }

  // update velocity
  Kinematics.initialVelo = (Kinematics.finalAltitude - Kinematics.initialAlt) /
                           (this->state.deltaTimeIMU);

  if (debug == Dynamite || debug == ALL) {
    // SOAR_PRINT("Initial Velocity: %f\n", Kinematics.initialVelo);
  }

  // update altitude
  Kinematics.initialAlt = Kinematics.finalAltitude;

  recalculateGain(normalised_Altitude);

  // Save the gains that are not zero as previous gains
  // so once we have recovery phase these old gains are used
  if (everest.state.gain_IMU != 0) {
    everest.state.prev_gain_IMU = everest.state.gain_IMU;
  }
  if (everest.state.gain_Baro1 != 0) {
    everest.state.prev_gain_Baro1 = everest.state.gain_Baro1;
  }
  if (everest.state.gain_Baro2 != 0) {
    everest.state.prev_gain_Baro2 = everest.state.gain_Baro2;
  }

  if (debug == Dynamite || debug == ALL) {
    // SOAR_PRINT("Previous Gains\n");
    // SOAR_PRINT("Prev Gain IMU: %f\n", everest.state.prev_gain_IMU);
    // SOAR_PRINT("Prev Gain Baro1: %f\n", everest.state.prev_gain_Baro1);
    // SOAR_PRINT("Prev Gain Baro2: %f\n", everest.state.prev_gain_Baro2);
    // SOAR_PRINT("Prev Gain Baro3: %f\n", everest.state.prev_gain_Baro3);
    // SOAR_PRINT("Prev Gain Real Baro: %f\n\n",
    // everest.state.prev_gain_Real_Baro);
  }

  return normalised_Altitude;
}

// @brief calculation - new gain = 1 / abs(estimate - measurement)
void EverestTask::recalculateGain(double estimate) {
  double gainedEstimate = deriveChangeInVelocityToGetAltitude(
      estimate);  // pre integrated for altitude

  double gain_IMU =
      1 / fabsf(gainedEstimate -
                this->state.avgIMU.altitude);  // change to previous trusts
  double gain_Baro1 = 1 / fabsf(gainedEstimate - this->baro1.altitude);
  double gain_Baro2 = 1 / fabsf(gainedEstimate - this->baro2.altitude);

  //    if(debug == Third || debug == ALL){
  //        printf("\nRecalculate Gain - Before normalization\n");
  //        printf("Gain IMU: %f\n", gain_IMU);
  //        printf("Gain Baro1: %f\n", gain_Baro1);
  //        printf("Gain Baro2: %f\n", gain_Baro2);
  //        printf("Gain Baro3: %f\n", gain_Baro3);
  //        printf("Gain Real Baro: %f\n", gain_Real_Baro);
  //        printf("Gained Estimate: %f\n", gainedEstimate);

  //        printf("Altitude: %f\n", estimate);
  //        printf("Baro1: %f\n", this->baro1.altitude);
  //        printf("Baro2: %f\n", this->baro2.altitude);
  //        printf("Baro3: %f\n", this->baro3.altitude);
  //        printf("Real Baro: %f\n\n", this->realBaro.altitude);
  //    }

  // normalise
  this->state.gain_IMU = gain_IMU / (gain_IMU + gain_Baro1 + gain_Baro2);
  this->state.gain_Baro1 = gain_Baro1 / (gain_IMU + gain_Baro1 + gain_Baro2);
  this->state.gain_Baro2 = gain_Baro2 / (gain_IMU + gain_Baro1 + gain_Baro2);

  if (debug == Dynamite || debug == ALL) {
    // SOAR_PRINT("\nRecalculate Gain\n");
    // SOAR_PRINT("New Gain IMU: %f\n", this->state.gain_IMU);
    // SOAR_PRINT("New Gain Baro1: %f\n", this->state.gain_Baro1);
    // SOAR_PRINT("New Gain Baro2: %f\n", this->state.gain_Baro2);
    // SOAR_PRINT("New Gain Baro3: %f\n", this->state.gain_Baro3);
    // SOAR_PRINT("New Gain Real Baro: %f\n\n", this->state.gain_Real_Baro);
  }
}

/**
 * @brief Converts the STDs to coefficients
 */
void EverestTask::calculateSTDCoefficients() {
  // calculate standard deviation coefficients
  double std_IMU = this->state.gain_IMU;
  double std_Baro1 = this->state.gain_Baro1;
  double std_Baro2 = this->state.gain_Baro2;

  double sumSTD1 = pow(everest.state.gain_IMU, 2) +
                   pow(everest.state.gain_Baro1, 2) +
                   pow(everest.state.gain_Baro2, 2);

  // normalise
  this->state.std_IMU = pow(std_IMU, 2) / sumSTD1;
  this->state.std_Baro1 = pow(std_Baro1, 2) / sumSTD1;
  this->state.std_Baro2 = pow(std_Baro2, 2) / sumSTD1;

  //    if(debug == Dynamite || debug == ALL){
  //        printf("\nStandard Deviation Coefficients\n");
  //        printf("STD IMU: %f\n", this->state.std_IMU);
  //        printf("STD Baro1: %f\n", this->state.std_Baro1);
  //        printf("STD Baro2: %f\n", this->state.std_Baro2);
  //        printf("STD Baro3: %f\n", this->state.std_Baro3);
  //        printf("STD Real Baro: %f\n\n", this->state.std_Real_Baro);
  //    }
}

/**
 * @brief Calculates the derivative of the altitude
 *
 * @param estimate the estimated altitude
 *
 * @return velocity
 *
 *   Internal
 */
double EverestTask::deriveChangeInVelocityToGetAltitude(double estimate) {
  double deltaTimeAverage = (this->baro1.deltaTime + this->baro2.deltaTime +
                             this->state.deltaTimeIMU) /
                            3.0;

  double velocityZ = (this->AltitudeList.secondLastAltitude -
                      4 * this->AltitudeList.lastAltitude + 3 * estimate) /
                     (2.0 * deltaTimeAverage);

  double newAltitude =
      this->AltitudeList.lastAltitude + velocityZ * deltaTimeAverage;

  if (debug == Dynamite || debug == ALL) {
    // SOAR_PRINT("\nDerivative for new gain\n");
    // SOAR_PRINT("Velocity: %f\n", velocityZ);
    // SOAR_PRINT("New Altitude: %f\n", newAltitude);
    // SOAR_PRINT("Delta Time Average: %f\n\n", deltaTimeAverage);
  }

  return newAltitude;
}

/**
 * @brief Getter for finalAltitude
 *
 * @return kinematics struct
 *
 */
double getFinalAltitude() { return Kinematics->finalAltitude; }

/**
 * @brief Tares the altitude to the ground and calibrates zero ground offset for
 * IMU
 *
 * @category Internal | Asynchronous
 *
 * Call function 10*RefreshRate times to get the initial altitude
 *
 * Once finished will print the tared altitude and set it as the initial
 * altitude
 */
void EverestTask::tare(IMUData& imu1, IMUData& imu2, BarosData baro1,
                       BarosData baro2) {
  // for 10 seconds collect baro
  // decrement the time

  // average pressures
  // have to do since function is weird
  double average = 0;
  int numberOfSamples = 0;

  if (baro1.pressure != 0) {
    average = average + convertToAltitude(baro1.pressure);
    numberOfSamples++;

    if (debug == Secondary || debug == ALL) {
      // SOAR_PRINT("average: %f number: %d \n", average, numberOfSamples);
    }
  }

  if (baro2.pressure != 0) {
    average = average + convertToAltitude(baro2.pressure);
    numberOfSamples++;

    if (debug == Secondary || debug == ALL) {
      // SOAR_PRINT("average: %f number: %d \n", average, numberOfSamples);
    }
  }

  if (numberOfSamples != 0) {
    sum += average / numberOfSamples;
  }

  // if(imu1.accelX > 99999){
  //     printf("imu1.accelX: %f", imu1.accelX);
  // }else{
  this->zeroOffsetAccel = {this->zeroOffsetAccel[0] + imu1.accelX,
                           this->zeroOffsetAccel[1] + imu1.accelY,
                           this->zeroOffsetAccel[2] + imu1.accelZ};
  this->zeroOffsetGyro = {this->zeroOffsetGyro[0] + imu1.gyroX,
                          this->zeroOffsetGyro[1] + imu1.gyroY,
                          this->zeroOffsetGyro[2] + imu1.gyroZ};

  // printf("zeroOffsetAccel[0]: %f, zeroOffsetAccel[1]: %f, zeroOffsetAccel[2]:
  // %f \n", this->zeroOffsetAccel[0], this->zeroOffsetAccel[1],
  // this->zeroOffsetAccel[2]);

  // printf("zeroOffsetGyro[0]: %f, zeroOffsetGyro[1]: %f, zeroOffsetGyro[2]: %f
  // \n", this->zeroOffsetGyro[0], this->zeroOffsetGyro[1],
  // this->zeroOffsetGyro[2]);

  if (debug == Secondary || debug == ALL) {
    // SOAR_PRINT("average: %f number: %d \n", average, numberOfSamples);
  }

  // }

  // if(!isinf(imu2.accelX)){
  this->zeroOffsetAccel2 = {this->zeroOffsetAccel2[0] + imu2.accelX,
                            this->zeroOffsetAccel2[1] + imu2.accelY,
                            this->zeroOffsetAccel2[2] + imu2.accelZ};
  this->zeroOffsetGyro2 = {this->zeroOffsetGyro2[0] + imu2.gyroX,
                           this->zeroOffsetGyro2[1] + imu2.gyroY,
                           this->zeroOffsetGyro2[2] + imu2.gyroZ};

  // printf("zeroOffsetAccel2[0]: %f, zeroOffsetAccel2[1]: %f,
  // zeroOffsetAccel2[2]: %f \n", this->zeroOffsetAccel2[0],
  // this->zeroOffsetAccel2[1], this->zeroOffsetAccel2[2]);

  // printf("zeroOffsetGyro2[0]: %f, zeroOffsetGyro2[1]: %f, zeroOffsetGyro2[2]:
  // %f \n", this->zeroOffsetGyro2[0], this->zeroOffsetGyro2[1],
  // this->zeroOffsetGyro2[2]);

  if (debug == Secondary || debug == ALL) {
    // SOAR_PRINT("average: %f number: %d \n", average, numberOfSamples);
  }
  // }

  if (debug == Secondary || debug == ALL) {
    // SOAR_PRINT("Tare Sum: %f\n", sum);
    // SOAR_PRINT("Number of samples %f\n", numberOfSamples);
  }

  if (theTime == 0) {
    sum = sum / (CALIBRATION_TIME * RATE_BARO);
    this->Kinematics.initialAlt = sum;

    this->zeroOffsetAccel = {
        this->zeroOffsetAccel[0] / (CALIBRATION_TIME * RATE_BARO),
        this->zeroOffsetAccel[1] / (CALIBRATION_TIME * RATE_BARO),
        this->zeroOffsetAccel[2] / (CALIBRATION_TIME * RATE_BARO)};

    this->zeroOffsetGyro = {
        this->zeroOffsetGyro[0] / (CALIBRATION_TIME * RATE_BARO),
        this->zeroOffsetGyro[1] / (CALIBRATION_TIME * RATE_BARO),
        this->zeroOffsetGyro[2] / (CALIBRATION_TIME * RATE_BARO)};

    this->zeroOffsetAccel2 = {
        this->zeroOffsetAccel2[0] / (CALIBRATION_TIME * RATE_BARO),
        this->zeroOffsetAccel2[1] / (CALIBRATION_TIME * RATE_BARO),
        this->zeroOffsetAccel2[2] / (CALIBRATION_TIME * RATE_BARO)};

    this->zeroOffsetGyro2 = {
        this->zeroOffsetGyro2[0] / (CALIBRATION_TIME * RATE_BARO),
        this->zeroOffsetGyro2[1] / (CALIBRATION_TIME * RATE_BARO),
        this->zeroOffsetGyro2[2] / (CALIBRATION_TIME * RATE_BARO)};

    // SOAR_PRINT("Tare Initial Altitude: %f\n", this->Kinematics.initialAlt);
    isTared = true;

    // printf("\nCalibration offsets:\n  accel1(%f,%f,%f),\n accel2(%f,%f,%f),\n
    // gyro(%f,%f,%f),\n  gyro2(%f,%f,%f)\n\n", this->zeroOffsetAccel[0],
    // this->zeroOffsetAccel[1], this->zeroOffsetAccel[2],
    // this->zeroOffsetGyro[0],  this->zeroOffsetGyro[1],
    // this->zeroOffsetGyro[2],
    // this->zeroOffsetAccel2[0],this->zeroOffsetAccel2[1],
    // this->zeroOffsetAccel2[2], this->zeroOffsetGyro2[0],
    // this->zeroOffsetGyro2[1], this->zeroOffsetGyro2[2]
    // );

    // ExternalUpdate(imu1, imu2, baro1, baro2, baro3, realBaro);
  }

  // call to update time and offsets for these structs
  IMU_Update(imu1, imu2);

  theTime -= 1;
}

/**
 * @brief accel is in m/s -> gs, gyro is passed in dps, pressure is in Pa, real
 * is for altitude ONLY in m Aligns before sending to update
 */
double EverestTask::finalWrapper(
    float accelX1, float accelY1, float accelZ1, float gyroX1, float gyroY1,
    float gyroZ1, float magX1, float magY1, float magZ1, float accelX2,
    float accelY2, float accelZ2, float gyroX2, float gyroY2, float gyroZ2,
    float magX2, float magY2, float magZ2, float pressure1, float pressure2,
    float timeIMU1, float timeIMU2, float timeBaro1, float timeBaro2,
    MadAxesAlignment alignment, MadAxesAlignment alignment2) {
  IMUData sensorData = {
      timeIMU1,
      gyroX1,
      gyroY1,
      gyroZ1,
      (float)(accelX1 / 9.81),
      (float)(accelY1 / 9.81),
      (float)(accelZ1 / 9.81),
      magX1,
      magY1,
      magZ1,
  };

  IMUData sensorData2 = {
      timeIMU2,
      gyroX2,
      gyroY2,
      gyroZ2,
      (float)(accelX2 / 9.81),
      (float)(accelY2 / 9.81),
      (float)(accelZ2 / 9.81),
      magX2,
      magY2,
      magZ2,
  };

  BarosData baro1 = {timeBaro1, pressure1, 0, 0};

  BarosData baro2 = {timeBaro2, pressure2, 0, 0};

  // align
  madVector imu1Gyro = {sensorData.gyroX, sensorData.gyroY, sensorData.gyroZ};
  madVector imu1Accel = {sensorData.accelX, sensorData.accelY,
                         sensorData.accelZ};
  madVector imu1Mag = {sensorData.magX, sensorData.magY, sensorData.magZ};

  madVector imu1GyroAligned = infusion->AxesSwitch(imu1Gyro, alignment);
  madVector imu1AccelAligned = infusion->AxesSwitch(imu1Accel, alignment);
  madVector imu1MagAligned = infusion->AxesSwitch(imu1Mag, alignment);

  madVector imu2Gyro = {sensorData2.gyroX, sensorData2.gyroY,
                        sensorData2.gyroZ};
  madVector imu2Accel = {sensorData2.accelX, sensorData2.accelY,
                         sensorData2.accelZ};
  madVector imu2Mag = {sensorData2.magX, sensorData2.magY, sensorData2.magZ};

  madVector imu2GyroAligned = infusion->AxesSwitch(imu2Gyro, alignment2);
  madVector imu2AccelAligned = infusion->AxesSwitch(imu2Accel, alignment2);
  madVector imu2MagAligned = infusion->AxesSwitch(imu2Mag, alignment2);

  //    if(debug == Secondary || debug == ALL){
  //        printf("Aligned: Gyro: (%.6f, %.6f, %.6f) deg/s, Accel: (%.6f, %.6f,
  //        %.6f) g\n",
  //            imu1GyroAligned.axis.x, imu1GyroAligned.axis.y,
  //            imu1GyroAligned.axis.z, imu1AccelAligned.axis.x,
  //            imu1AccelAligned.axis.y, imu1AccelAligned.axis.z);
  //    }

  // feed vectors into sensorData structs
  sensorData.gyroX = imu1GyroAligned.axis.x;
  sensorData.gyroY = imu1GyroAligned.axis.y;
  sensorData.gyroZ = imu1GyroAligned.axis.z;

  sensorData.accelX = imu1AccelAligned.axis.x;
  sensorData.accelY = imu1AccelAligned.axis.y;
  sensorData.accelZ = imu1AccelAligned.axis.z;

  sensorData.magX = imu1MagAligned.axis.x;
  sensorData.magY = imu1MagAligned.axis.y;
  sensorData.magZ = imu1MagAligned.axis.z;

  // second IMU
  sensorData2.gyroX = imu2GyroAligned.axis.x;
  sensorData2.gyroY = imu2GyroAligned.axis.y;
  sensorData2.gyroZ = imu2GyroAligned.axis.z;

  sensorData2.accelX = imu2AccelAligned.axis.x;
  sensorData2.accelY = imu2AccelAligned.axis.y;
  sensorData2.accelZ = imu2AccelAligned.axis.z;

  sensorData2.magX = imu2MagAligned.axis.x;
  sensorData2.magY = imu2MagAligned.axis.y;
  sensorData2.magZ = imu2MagAligned.axis.z;

  double eAltitude =
      everest.ExternalUpdate(sensorData, sensorData2, baro1, baro2);

  // SOAR_PRINT("Altitude: %f\n", eAltitude);

  return eAltitude;
}

/**
 * @brief Resets isTared flag to re-initialize the tare
 */
void setIsTare(bool isTare) { isTared = isTare; }

// --------------------------------------------------- END OF EVEREST
// ---------------------------------------------------//
void monteCarloSims() {
  // Monte Carlo simulations
  // 1. Generate random scenarios
  // 2. Run the simulations
  // 3. Find the nearest scenarios
  // 4. Run the simulations for the nearest scenarios
  // 5. Find the nearest scenarios
  // 6. Repeat until convergence

  // 1. Generate random scenarios
  // 2. Run the simulations
  // 3. Find the nearest scenarios
  // 4. Run the simulations for the nearest scenarios
  // 5. Find the nearest scenarios
  // 6. Repeat until convergence

  // everest.MadgwickSetup();

  // HALO halo = HALO();

  // MatrixXf Q(3,3);
  // Q <<    0, 0, 0,
  //         0, 0, 0,
  //         0, 0, 0;

  // // everest covariance matrix
  // MatrixXf R0(3,3);
  // R0 <<   0, 0, 0,
  //         0, 0, 0,
  //         0, 0, 0;

  // MatrixXf P0(3,3);
  // P0 <<   0, 0, 0,
  //         0, 0, 0,
  //         0, 0, 0;
}

#define MAX_LINE_LENGTH 1024

/**
 * Serves to just initialize structs
 */
int main() {
  // Setup Madgwick
  // Attach Madgwick to Everest
  everest.MadgwickSetup();

  HALO halo = HALO();

  MatrixXf Q(3, 3);
  Q << 100, 0, 0, 0, 40, 0, 0, 0, 8;

  // everest covariance matrix
  MatrixXf R0(3, 3);
  R0 << 200, 0.5, 0.5, 0.5, 100, 1, 0.5, 1, 10;

  MatrixXf P0(3, 3);
  P0 << 50, 0, 0, 0, 0, 0, 0, 0, 0;

// clear log files
#ifdef LOGON
  if (remove("gains.txt") == 0) {
    std::cout << "File deleted successfully: " << "gains.txt" << std::endl;
  } else {
    std::perror("Error deleting file");
  }

  if (remove("predictedValues.txt") == 0) {
    std::cout << "File deleted successfully: " << "predictedValues.txt"
              << std::endl;
  } else {
    std::perror("Error deleting file");
  }

  if (remove("log.txt") == 0) {
    std::cout << "File deleted successfully: " << "log.txt" << std::endl;
  } else {
    std::perror("Error deleting file");
  }

  if (remove("nearestScenarios.txt") == 0) {
    std::cout << "File deleted successfully: " << "nearestScenarios.txt"
              << std::endl;
  } else {
    std::perror("Error deleting file");
  }

  if (remove("nearestScenariosFormatted.txt") == 0) {
    std::cout << "File deleted successfully: "
              << "nearestScenariosFormatted.txt" << std::endl;
  } else {
    std::perror("Error deleting file");
  }

  if (remove("sigmaPoints.txt") == 0) {
    std::cout << "File deleted successfully: " << "sigmaPoints.txt"
              << std::endl;
  } else {
    std::perror("Error deleting file");
  }

  if (remove("sigmaPoints1.txt") == 0) {
    std::cout << "File deleted successfully: " << "sigmaPoints1.txt"
              << std::endl;
  } else {
    std::perror("Error deleting file");
  }

  if (remove("sigmaPoints2.txt") == 0) {
    std::cout << "File deleted successfully: " << "sigmaPoints2.txt"
              << std::endl;
  } else {
    std::perror("Error deleting file");
  }

  if (remove("sigmaPoints3.txt") == 0) {
    std::cout << "File deleted successfully: " << "sigmaPoints3.txt"
              << std::endl;
  } else {
    std::perror("Error deleting file");
  }

  if (remove("sigmaPoints4.txt") == 0) {
    std::cout << "File deleted successfully: " << "sigmaPoints4.txt"
              << std::endl;
  } else {
    std::perror("Error deleting file");
  }

  if (remove("sigmaPoints5.txt") == 0) {
    std::cout << "File deleted successfully: " << "sigmaPoints5.txt"
              << std::endl;
  } else {
    std::perror("Error deleting file");
  }

  if (remove("sigmaPoints6.txt") == 0) {
    std::cout << "File deleted successfully: " << "sigmaPoints6.txt"
              << std::endl;
  } else {
    std::perror("Error deleting file");
  }

  FILE* predictedValues = fopen("predictedValues.txt", "a+");

  if (!predictedValues) {
    fprintf(stderr, "Error opening HALO.txt...exiting\n");
    exit(1);
  }

  fprintf(
      predictedValues,
      "s1_alt,s1_velo,s1_acc,s2_alt,s2_velo,s2_acc,s3_alt,s3_velo,s3_acc,s4_"
      "alt,s4_velo,s4_acc,s5_alt,s5_velo,s5_acc,s6_alt,s6_velo,s6_acc\n");

  fclose(predictedValues);

  FILE* gains = fopen("gains.txt", "a+");
  if (!gains) {
    fprintf(stderr, "Error opening HALO.txt...exiting\n");
    exit(1);
  }

  fprintf(gains, "gain_IMU,gain_Baro1,gain_Baro2\n");

  fclose(gains);

  FILE* sigmaPoints = fopen("sigmaPoints.txt", "a+");
  if (!sigmaPoints) {
    fprintf(stderr, "Error opening HALO.txt...exiting\n");
    exit(1);
  }

  fprintf(sigmaPoints, "alt,velo,acc\n");

  fclose(sigmaPoints);

  FILE* sigmaPoints1 = fopen("sigmaPoints1.txt", "a+");
  if (!sigmaPoints1) {
    fprintf(stderr, "Error opening HALO.txt...exiting\n");
    exit(1);
  }

  fprintf(sigmaPoints1, "alt,velo,acc\n");

  fclose(sigmaPoints1);

  FILE* sigmaPoints2 = fopen("sigmaPoints2.txt", "a+");

  if (!sigmaPoints2) {
    fprintf(stderr, "Error opening HALO.txt...exiting\n");
    exit(1);
  }

  fprintf(sigmaPoints2, "alt,velo,acc\n");

  fclose(sigmaPoints2);

  FILE* sigmaPoints3 = fopen("sigmaPoints3.txt", "a+");

  if (!sigmaPoints3) {
    fprintf(stderr, "Error opening HALO.txt...exiting\n");
    exit(1);
  }

  fprintf(sigmaPoints3, "alt,velo,acc\n");

  fclose(sigmaPoints3);

  FILE* sigmaPoints4 = fopen("sigmaPoints4.txt", "a+");

  if (!sigmaPoints4) {
    fprintf(stderr, "Error opening HALO.txt...exiting\n");
    exit(1);
  }

  fprintf(sigmaPoints4, "alt,velo,acc\n");

  fclose(sigmaPoints4);

  FILE* sigmaPoints5 = fopen("sigmaPoints5.txt", "a+");

  if (!sigmaPoints5) {
    fprintf(stderr, "Error opening HALO.txt...exiting\n");
    exit(1);
  }

  fprintf(sigmaPoints5, "alt,velo,acc\n");

  fclose(sigmaPoints5);

  FILE* sigmaPoints6 = fopen("sigmaPoints6.txt", "a+");

  if (!sigmaPoints6) {
    fprintf(stderr, "Error opening HALO.txt...exiting\n");
    exit(1);
  }

  fprintf(sigmaPoints6, "alt,velo,acc\n");

  fclose(sigmaPoints6);

  FILE* nearestScenarios = fopen("nearestScenarios.txt", "a+");

  if (!nearestScenarios) {
    fprintf(stderr, "Error opening HALO.txt...exiting\n");
    exit(1);
  }

  fprintf(nearestScenarios,
          "lowestDistance,secondLowestDistance,firstScenario,SecondScenario\n");

  fclose(nearestScenarios);

  // test purposes
  FILE* file = fopen(
      "HALO.txt",
      "w+");  // Open the file for appending or create it if it doesn't exist
  halo.file = file;
  if (!file) {
    fprintf(stderr, "Error opening HALO.txt...exiting\n");
    exit(1);
  }

  fprintf(file,
          "Time,Everest_Alt,Everest_Velo,Everest_Accel,Halo_ALt, Halo_Velo, "
          "Halo_accel\n");

#endif

/////// home
#ifdef HOME

  FILE* file1 = fopen(
      "C:/Users/andin/OneDrive/Documents/AllRepos/UnscentedKalmanFilter/"
      "EverestLibrary_HALO/EverestL/EverestLibrary/Imu_Baro.csv",
      "r");

  FILE* simsFile = fopen(
      "C:/Users/andin/OneDrive/Documents/AllRepos/UnscentedKalmanFilter/"
      "EverestLibrary_HALO/EverestL/EverestLibrary/beforeSimsF2_short.csv",
      "r");

  FILE* simsAfterFile = fopen(
      "C:/Users/andin/OneDrive/Documents/AllRepos/UnscentedKalmanFilter/"
      "EverestLibrary_HALO/EverestL/EverestLibrary/beforeSimsF2_short.csv",
      "r");

  FILE* taberFile = fopen(
      "C:/Users/andin/OneDrive/Documents/AllRepos/UnscentedKalmanFilter/"
      "EverestLibrary_HALO/EverestL/EverestLibrary/taber_launch_formattedF.csv",
      "r");

#endif

////// away
#ifndef HOME

  FILE* file1 = fopen(
      "C:/Users/Andrey/Documents/UKFRepo/UnscentedKalmanFilter/"
      "EverestLibrary_HALO/EverestL/EverestLibrary/Imu_Baro.csv",
      "r");

  FILE* simsFile = fopen(
      "C:/Users/Andrey/Documents/UKFRepo/UnscentedKalmanFilter/"
      "EverestLibrary_HALO/EverestL/EverestLibrary/beforeSimsF2_Short.csv",
      "r");

  FILE* simsAfterFile = fopen(
      "C:/Users/Andrey/Documents/UKFRepo/UnscentedKalmanFilter/"
      "EverestLibrary_HALO/EverestL/EverestLibrary/beforeSimsF2_Short.csv",
      "r");

  FILE* taberFile = fopen(
      "C:/Users/Andrey/Documents/UKFRepo/UnscentedKalmanFilter/"
      "EverestLibrary_HALO/EverestL/EverestLibrary/taber_launch_formattedF.csv",
      "r");

#endif

  //////

  if (!file1) {
    perror("Error opening Imu_Baro.csv");
    return 1;
  }

  if (!simsFile) {
    perror("Error opening simsFile.csv");
    return 1;
  }

  if (!simsAfterFile) {
    perror("Error opening simsAfterFile.csv");
    return 1;
  }

  if (!taberFile) {
    perror("Error opening taberFile.csv");
    return 1;
  }

  // read first line and preset the deltaTime to timestamp
  char line[MAX_LINE_LENGTH];
  std::clock_t start;
  double duration;

  int howMany = 1;
  int i = 0;
  int numberOfScenarios = 6;

  // std::vector<std::vector<std::vector<float>>> scenarioListofVectorsBefore =
  // {{{}}, {{}}, {{}}, {{}}, {{}}, {{}}}; for(int i = 0; i < numberOfScenarios;
  // i++){
  //     scenarioListofVectorsBefore. = howMany;
  // }
  float alt;
  float velo;
  float acc;
  float time1;
  std::vector<float> temp = {0, 0, 0, 0};

  // // skip header
  // fgets(line, sizeof(line), simsFile);

  // while (fgets(line, sizeof(line), simsFile)) {
  //     char *token = strtok(line, ",");

  //     for(int i = 0; i < numberOfScenarios; i++){
  //         alt = atof(token); // Convert the time value to float
  //         token = strtok(NULL, ",");
  //         velo = atof(token);
  //         token = strtok(NULL, ",");
  //         acc = atof(token);
  //         token = strtok(NULL, ",");
  //         time1 = atof(token);
  //         token = strtok(NULL, ",");

  //         temp = {alt, velo, acc, time1};
  //         // printf("temp[%d] = (%f,%f,%f,%f)\n", i, temp[0], temp[1],
  //         temp[2], temp[3]); if(temp[0] != 0 && temp[1] != 0 && temp[2] != 0
  //         && temp[3] != 0){
  //             scenarioListofVectorsBefore[i].push_back(temp);
  //         }

  //     }

  //     // printf("\n");

  // }

  // std::vector<std::vector<std::vector<float>>> scenarioListofVectorsAfter =
  // {{{}}, {{}}, {{}}, {{}}, {{}}, {{}}};

  // // skip header
  // fgets(line, sizeof(line), simsAfterFile);

  // while (fgets(line, sizeof(line), simsAfterFile)) {
  //     char *token = strtok(line, ",");

  //     for(int i = 0; i < numberOfScenarios; i++){
  //         alt = atof(token); // Convert the time value to float
  //         token = strtok(NULL, ",");
  //         velo = atof(token);
  //         token = strtok(NULL, ",");
  //         acc = atof(token);
  //         token = strtok(NULL, ",");
  //         time1 = atof(token);
  //         token = strtok(NULL, ",");

  //         std::vector<float> temp = {alt, velo, acc, time1};
  //         // printf("after temp[%d] = (%f,%f,%f,%f)\n", i, temp[0], temp[1],
  //         temp[2], temp[3]); if(temp[0] != 0 && temp[1] != 0 && temp[2] != 0
  //         && temp[3] != 0){
  //             scenarioListofVectorsAfter[i].push_back(temp);
  //         }

  //     }

  //     // printf("\n");
  // }

  // scenarioListofVectorsBefore[0].erase(scenarioListofVectorsBefore[0].begin());
  // scenarioListofVectorsBefore[1].erase(scenarioListofVectorsBefore[1].begin());
  // scenarioListofVectorsBefore[2].erase(scenarioListofVectorsBefore[2].begin());
  // scenarioListofVectorsBefore[3].erase(scenarioListofVectorsBefore[3].begin());
  // scenarioListofVectorsBefore[4].erase(scenarioListofVectorsBefore[4].begin());
  // scenarioListofVectorsBefore[5].erase(scenarioListofVectorsBefore[5].begin());

  // scenarioListofVectorsAfter[0].erase(scenarioListofVectorsAfter[0].begin());
  // scenarioListofVectorsAfter[1].erase(scenarioListofVectorsAfter[1].begin());
  // scenarioListofVectorsAfter[2].erase(scenarioListofVectorsAfter[2].begin());
  // scenarioListofVectorsAfter[3].erase(scenarioListofVectorsAfter[3].begin());
  // scenarioListofVectorsAfter[4].erase(scenarioListofVectorsAfter[4].begin());
  // scenarioListofVectorsAfter[5].erase(scenarioListofVectorsAfter[5].begin());

  std::vector<std::vector<float>> megaList;

  // Helper lambda to add scenario index to each point
  auto addScenarioIndex = [](const std::vector<std::vector<float>>& sim,
                             int index) {
    std::vector<std::vector<float>> result;
    for (const auto& point : sim) {
      std::vector<float> newPoint = {point.at(0), point.at(1), point.at(2)};
      newPoint.push_back(static_cast<float>(index));
      result.push_back(newPoint);
      printf("newPoint: %f, %f, %f, %f\n", newPoint[0], newPoint[1],
             newPoint[2], newPoint[3]);
    }

    return result;
  };

  // Combine all simulation data into one list with scenario index
  auto sim1WithIndex = addScenarioIndex(sim1, 1);
  auto sim2WithIndex = addScenarioIndex(sim2, 2);
  auto sim3WithIndex = addScenarioIndex(sim3, 3);
  auto sim4WithIndex = addScenarioIndex(sim4, 4);
  auto sim5WithIndex = addScenarioIndex(sim5, 5);
  auto sim6WithIndex = addScenarioIndex(sim6, 6);

  megaList.insert(megaList.end(), sim1WithIndex.begin(), sim1WithIndex.end());
  megaList.insert(megaList.end(), sim2WithIndex.begin(), sim2WithIndex.end());
  megaList.insert(megaList.end(), sim3WithIndex.begin(), sim3WithIndex.end());
  megaList.insert(megaList.end(), sim4WithIndex.begin(), sim4WithIndex.end());
  megaList.insert(megaList.end(), sim5WithIndex.begin(), sim5WithIndex.end());
  megaList.insert(megaList.end(), sim6WithIndex.begin(), sim6WithIndex.end());

  // Create the mega tree
  KDTree megaTree(megaList);

  // create scenarios with before and after lists
  Scenario scenario1 = Scenario{sim1, sim1, 1};
  // scenario1.createTree();
  scenario1.passMegaTree(megaTree);
  Scenario scenario2 = Scenario{sim2, sim2, 2};
  // scenario2.createTree();
  scenario2.passMegaTree(megaTree);
  Scenario scenario3 = Scenario{sim3, sim3, 3};
  // scenario3.createTree();
  scenario3.passMegaTree(megaTree);
  Scenario scenario4 = Scenario{sim4, sim4, 4};
  // scenario4.createTree();
  scenario4.passMegaTree(megaTree);
  Scenario scenario5 = Scenario{sim5, sim5, 5};
  // scenario5.createTree();
  scenario5.passMegaTree(megaTree);
  Scenario scenario6 = Scenario{sim6, sim6, 6};
  // scenario6.createTree();
  scenario6.passMegaTree(megaTree);

  std::vector<Scenario> scenarios = {scenario1, scenario2, scenario3,
                                     scenario4, scenario5, scenario6};

  float totalTime = 0;

  // add to a list of scenarios
  // for(int i = 0; i < numberOfScenarios; i++){
  //     Scenario scenario = Scenario{scenarioListofVectorsBefore[i],
  //     scenarioListofVectorsAfter[i], "Scenario " + std::to_string(i)};
  //     scenarios.push_back(scenario);
  // }

  halo.setScenarios(scenarios);
  halo.deltaTime = 1 / 3;

  // printf("size of vector: %d\n", taberLaunch.size());

  for (int i = 0; i < 555; i++) {
    // Tokenize the line using strtok
    // Parse accelerometer readings (X, Y, Z)
    float time = taberLaunch.at(i).at(0);
    // printf("time: %f\n", time);
    float accelX = taberLaunch[i][1];
    float accelY = taberLaunch[i][2];
    float accelZ = taberLaunch[i][3];

    // Parse gyroscope readings (X, Y, Z)
    float gyroX = taberLaunch[i][4];
    float gyroY = taberLaunch[i][5];
    float gyroZ = taberLaunch[i][6];

    // Parse magnetometer readings (X, Y, Z)
    float magX = taberLaunch[i][7];
    float magY = taberLaunch[i][8];
    float magZ = taberLaunch[i][9];

    // printf("magX: %f, magY: %f, magZ: %f\n", magX, magY, magZ);
    // printf("gyroX: %f, gyroY: %f, gyroZ: %f\n", gyroX, gyroY, gyroZ);
    // printf("accelX: %f, accelY: %f, accelZ: %f\n", accelX, accelY, accelZ);

    // Parse pressure readings
    float pressure = baroData[i][1];

    // printf("pressure: %f\n", pressure);

    IMUData sensorData = {
        time, gyroX, gyroY, gyroZ, accelX, accelY, accelZ, magX, magY, magZ,
    };

    IMUData sensorData2 = {
        time, gyroX, gyroY, gyroZ, accelX, accelY, accelZ, magX, magY, magZ,
    };

    start = std::clock();

    BarosData baro1 = {time, pressure, 0, 0};

    BarosData baro2 = {time, pressure, 0, 0};

    // if(howMany <= 13){

    // printf("\n#%d
    // Sample--------------------------------------------------------------------------\n\n",
    // howMany);

    // // Example: Print all sensor readings
    // if(debug == RAW || debug == ALL){
    //     printf("Raw Time: %.6f s, Gyro: (%.6f, %.6f, %.6f) deg/s, Accel:
    //     (%.6f, %.6f, %.6f) g Pressure: (%.f, %.f, %.f, %.f)\n",
    //         time, sensorData.gyroX, sensorData.gyroY, sensorData.gyroZ,
    //         sensorData.accelX, sensorData.accelY, sensorData.accelZ,
    //         baro1.pressure, baro2.pressure, baro3.pressure,
    //         realBaro.pressure);
    // }

    // madVector imu1Gyro = {sensorData.gyroX, sensorData.gyroY,
    // sensorData.gyroZ}; madVector imu1Accel = {sensorData.accelX,
    // sensorData.accelY, sensorData.accelZ};

    // madVector imu1GyroAligned = infusion->AxesSwitch(imu1Gyro,
    // MadAxesAlignmentPXPYNZ); madVector imu1AccelAligned =
    // infusion->AxesSwitch(imu1Accel, MadAxesAlignmentPXPYNZ);

    // madVector imu2Gyro = {sensorData2.gyroX, sensorData2.gyroY,
    // sensorData2.gyroZ}; madVector imu2Accel = {sensorData2.accelX,
    // sensorData2.accelY, sensorData2.accelZ};

    // madVector imu2GyroAligned = infusion->AxesSwitch(imu2Gyro,
    // MadAxesAlignmentPXPYNZ); madVector imu2AccelAligned =
    // infusion->AxesSwitch(imu2Accel, MadAxesAlignmentPXPYNZ);

    // if(debug == Secondary || debug == ALL){
    //     printf("Aligned: Gyro: (%.6f, %.6f, %.6f) deg/s, Accel: (%.6f, %.6f,
    //     %.6f) g\n",
    //         imu1GyroAligned.axis.x, imu1GyroAligned.axis.y,
    //         imu1GyroAligned.axis.z, imu1AccelAligned.axis.x,
    //         imu1AccelAligned.axis.y, imu1AccelAligned.axis.z);
    // }

    // // feed vectors into sensorData structs
    // sensorData.gyroX = imu1GyroAligned.axis.x;
    // sensorData.gyroY = imu1GyroAligned.axis.y;
    // sensorData.gyroZ = imu1GyroAligned.axis.z;

    // sensorData.accelX = imu1AccelAligned.axis.x;
    // sensorData.accelY = imu1AccelAligned.axis.y;
    // sensorData.accelZ = imu1AccelAligned.axis.z;

    // // second IMU
    // sensorData2.gyroX = imu2GyroAligned.axis.x;
    // sensorData2.gyroY = imu2GyroAligned.axis.y;
    // sensorData2.gyroZ = imu2GyroAligned.axis.z;

    // sensorData2.accelX = imu2AccelAligned.axis.x;
    // sensorData2.accelY = imu2AccelAligned.axis.y;
    // sensorData2.accelZ = imu2AccelAligned.axis.z;

    // everest.IMU_Update(sensorData, sensorData2);

    // double eAltitude = everest.AlignedExternalUpdate(sensorData, sensorData2,
    // baro1, baro2, baro3, realBaro, MadAxesAlignmentPXPYNZ);
    // double eAltitude = everest.ExternalUpdate(sensorData, sensorData2, baro1,
    // baro2, baro3, realBaro);
    //    double eAltitude = finalWrapper(accelX, accelY, accelZ, gyroX, gyroY,
    //    gyroZ, pressure, pressure, pressure, pressure, time, time, time, time,
    //    time, time, MadAxesAlignmentPXPYNZ, MadAxesAlignmentPXPYNZ);
    EverestData everestData = {
        sensorData.time,    sensorData2.time,   baro1.time,
        baro2.time,

        baro1.pressure,     baro2.pressure,

        sensorData.accelX,  sensorData.accelY,  sensorData.accelZ,
        sensorData.gyroX,   sensorData.gyroY,   sensorData.gyroZ,
        sensorData.magX,    sensorData.magY,    sensorData.magZ,

        sensorData2.accelX, sensorData2.accelY, sensorData2.accelZ,
        sensorData2.gyroX,  sensorData2.gyroY,  sensorData2.gyroZ,
        sensorData2.magX,   sensorData2.magY,   sensorData2.magZ,
    };

    double eAltitude = everest.TaskWrapper(everestData, MadAxesAlignmentPXPYNZ,
                                           MadAxesAlignmentPXPYNZ);
    double eVelocity = everest.getKinematics()->initialVelo;
    double eAccelerationZ = (everest.state.earthAcceleration - 1) * -9.81;

    if (howMany == 7) {
      VectorXf X0(3);
      X0 << everest.Kinematics.initialAlt, 0, 0;

      // Initialize with tare / GPS values
      halo.init(X0, P0, Q, R0);
    }

    // start timer for iteration
    start = std::clock();

    if (howMany > 7) {
      halo.setTime(time);

      halo.setStateVector(eAccelerationZ, eVelocity, eAltitude);

      std::vector<double> unitedStates = {halo.X0[0], halo.X0[1], halo.X0[2]};

      // printf("\nFinal Measurements - time, eAltitude, HAltitude, HVelo,
      // Haccel:\n %f,%f,%f,%f,%f\n", time, eAltitude, unitedStates[0],
      // unitedStates[1], unitedStates[2]);

#ifdef LOGON
      fprintf(file, "%f,%f,%f,%f,%f,%f,%f\n", time, eAltitude, eVelocity,
              eAccelerationZ, unitedStates[0], unitedStates[1],
              unitedStates[2]);
#endif
    }

    // update Kinematics

    // }

    howMany++;

    clock_t endTime = std::clock();

    totalTime += endTime - start;

    // printf("Time for one more (seconds): %f\n", duration/CLOCKS_PER_SEC);
  }
  // }

  // printf("Overall for %d samples: %f", howMany, totalTime/CLOCKS_PER_SEC);
  howMany = howMany - 7;

  std::cout << "Overall for " << howMany << " samples:\t\t\t\t\t\t\t\t\t"
            << totalTime / (double)CLOCKS_PER_SEC << std::endl;

#ifdef TIMERON
  std::cout << "Update time:\t\t\t\t\t\t\t\t\t\t" << halo.updateTime.count()
            << std::endl;
  std::cout << "Predict time:\t\t\t\t\t\t\t\t\t\t" << halo.predictTime.count()
            << std::endl;

  std::cout << "\tTriangulationTime:\t\t\t\t\t\t\t"
            << halo.triangulationTime.count() << std::endl;
  std::cout << "\tdModeltime:\t\t\t\t\t\t\t\t" << halo.dynamicModelTime.count()
            << std::endl;

  std::cout << "\t\tgetScenarioTime:\t\t\t\t\t" << halo.getScenarioTime.count()
            << std::endl;
  std::cout << "\t\tpPredictionTime:\t\t\t\t\t" << halo.PpredictionTime.count()
            << std::endl;
  std::cout << "\t\tprojErrorTime:\t\t\t\t\t\t" << halo.projErrorTime.count()
            << std::endl;
  std::cout << "\t\tpreMeanTime:\t\t\t\t\t\t" << halo.preMeanTime.count()
            << std::endl;
  std::cout << "\t\tsPointTime:\t\t\t\t\t\t" << halo.sPointTime.count()
            << std::endl;
  std::cout << "\t\tpredictLoopTime:\t\t\t\t\t" << halo.predictLoopTime.count()
            << std::endl;
  std::cout << "\t\tendPredictLoopTime:\t\t\t\t\t"
            << halo.endPredictLoopTime.count() << std::endl;
  std::cout << "\t\tnearestScenariosTime:\t\t\t\t\t"
            << halo.nearestScenariosTime.count() << std::endl;

  std::cout << "\t\t\t->loopScenariosTime:\t\t\t"
            << halo.loopScenariosTime.count() << std::endl;
  std::cout << "\t\t\t\t->getListsTime:\t\t" << halo.getListsTime.count()
            << std::endl;
  std::cout << "\t\t\t\t->othersTime:\t\t" << halo.othersTime.count()
            << std::endl;
  std::cout << "\t\t\t\t->KDTreeTime:\t\t" << halo.KDTreeTime.count()
            << std::endl;
  std::cout << "\t\t\t\t->twoDistancesTime:\t" << halo.twoDistancesTime.count()
            << std::endl;
  std::cout << "\t\t\t\t->emplaceBackTime:\t" << halo.emplaceBackTime.count()
            << std::endl;

  std::cout << "\t\t\t->vectorsTime:\t\t\t\t" << halo.vectorsTime.count()
            << std::endl;
  std::cout << "\t\t\t->push_backTime:\t\t\t" << halo.push_backTime.count()
            << std::endl;
#endif

  fclose(file1);

#ifdef LOGON
  fclose(file);
#endif

  return 0;
}