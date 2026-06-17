// Altitude estimation using multiple sensors
#include <stdio.h>
#include <ctime>
#include <string.h>

#include <cstdio>
#include <cerrno>

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include "everest.hpp"
#include "DataBroker.hpp"
#include "cmsis_os.h"

// #define LOGON
// #define LOGMETRICS

// #define TIMERON

// #define TESTING_BUILD

#ifdef TESTING_BUILD
#define SOAR_PRINT(...) printf(__VA_ARGS__)
#endif

#define printf(...) ;

#ifdef LOGMETRICS
static FILE* everestGains = NULL;
#endif

#ifdef TESTING_BUILD
#include <direct.h>
#include "input_data.hpp"
// #include "pre_flight_data.cpp"
#include "gpsData.hpp"
#endif

FILE* haloFile;
FILE* everestFile;
bool tared = false;

// singleton accessor
EverestTask& EverestTask::getEverest() {
  static EverestTask inst;
  return inst;
}

// Extract data from a DataBroker command and update measurements
void EverestTask::Extract(const Command& cm) {
  // Periodically refresh which measurements are "available" by checking
  // how long since the last sample from each publisher. This keeps
  // `availableMeasurements` up to date even if no new DataBroker messages
  // arrive for a short while.
  uint32_t nowMs = static_cast<uint32_t>(xTaskGetTickCount() * portTICK_PERIOD_MS);
  if ((nowMs - lastAvailableRefreshMs) >= AVAILABLE_MEAS_REFRESH_MS) {
    for (int i = 0; i < 5; ++i) {
      uint32_t last = lastSampleTimeMs[i];
      // mark available if we have a recent sample (within refresh window)
      availableMeasurements[i] = ((nowMs > last) && ((nowMs - last) <= AVAILABLE_MEAS_REFRESH_MS)) ? 1 : 0;
    }
    lastAvailableRefreshMs = nowMs;
  }

  auto msgType = DataBroker::getMessageType(cm);
  switch(msgType) {
    case DataBrokerMessageTypes::IMU_DATA: {
      IMUData imu = DataBroker::ExtractData<IMUData>(cm);
//      SOAR_PRINT("Everest::Extract - IMU_DATA received id=%d gyro=(%d,%d,%d) accel=(%d,%d,%d)\n",
//                 imu.id, imu.gyro.x, imu.gyro.y, imu.gyro.z, imu.accel.x, imu.accel.y, imu.accel.z);
      IMUData_Everest iev{};
      float ts = static_cast<float>(nowMs);
      iev.time = ts;
      iev.gyroX = static_cast<float>(imu.gyro.x);
      iev.gyroY = static_cast<float>(imu.gyro.y);
      iev.gyroZ = static_cast<float>(imu.gyro.z);
      iev.accelX = static_cast<float>(imu.accel.x);
      iev.accelY = static_cast<float>(imu.accel.y);
      iev.accelZ = static_cast<float>(imu.accel.z);
      iev.magX = 0.0f;
      iev.magY = 0.0f;
      iev.magZ = 0.0f;
      if (imu.id == 0) { IMU1_Measurements(iev); availableMeasurements[0] = 1; lastSampleTimeMs[0] = nowMs; }
      else { IMU2_Measurements(iev); availableMeasurements[1] = 1; lastSampleTimeMs[1] = nowMs; }
      break;
    }
    case DataBrokerMessageTypes::BARO_DATA: {
      BaroData b = DataBroker::ExtractData<BaroData>(cm);
//      SOAR_PRINT("Everest::Extract - BARO_DATA received id=%d pressure=%u temp=%d\n", b.id, b.pressure, b.temp);
      BarosData bd{};
      bd.time = static_cast<float>(nowMs);
      bd.pressure = b.pressure;
      bd.altitude = 0;
      if (b.id == 0) { Baro1_Measurements(bd); availableMeasurements[2] = 1; lastSampleTimeMs[2] = nowMs; }
      else { Baro2_Measurements(bd); availableMeasurements[3] = 1; lastSampleTimeMs[3] = nowMs; }
      break;
    }
    case DataBrokerMessageTypes::GPS_DATA: {
      GPSData g = DataBroker::ExtractData<GPSData>(cm);
//      SOAR_PRINT("Everest::Extract - GPS_DATA received altitude=%d\n", g.antennaAltitude_.altitude_);
      GPS_Measurements(static_cast<float>(g.antennaAltitude_.altitude_));
      availableMeasurements[4] = 1;
      lastSampleTimeMs[4] = nowMs;
      break;
    }
    case DataBrokerMessageTypes::MAG_DATA: {
      MagData m = DataBroker::ExtractData<MagData>(cm);
//      SOAR_PRINT("Everest::Extract - MAG_DATA received (%d,%d,%d)\n", m.magX, m.magY, m.magZ);
      // Assign to imu1 mag fields by default; tasks can set as needed
      this->everestData.magX1 = static_cast<float>(m.magX);
      this->everestData.magY1 = static_cast<float>(m.magY);
      this->everestData.magZ1 = static_cast<float>(m.magZ);
      break;
    }
    default:
      break;
  }
}

void EverestTask::initialize1(systemState& state) {
  this->state.gain_IMU = 4 / 10.0;  // change to actual initial trusts
  this->state.gain_Baro1 = 3 / 10.0;
  state.gain_Baro2 = 3 / 10.0;

  Kinematics.initialVelo = 0;
  Kinematics.initialAlt = 0;
  Kinematics.finalAltitude = 0;
}

Infusion* EverestTask::ExternalInitialize() {
  initialize1(state);
  return &madgwick;
}

/**
 * To run:  g++ Infusion.cpp EverestTask.cpp -o Everest
 *          ./Everest
 */

using namespace std;

// SETTINGS (mostly for debugging, keep default for run)
enum debug_level {
  RAW = 0,         // raw data
  Secondary = 1,   // all operations before dynamite
  Dynamite = 2,    // everything during dynamite
  Third = 3,       // after dynamite
  ALL = 4,         // all
  NONE = 5,        // none
  HAL0 = 6,        // HALO
  Calibration = 7  // Calibration
};

debug_level debug = ALL;

/**
 * @brief Calls finalWrapper with data and alignment
 * IMPORTANT: PASS 0s FOR NOT UPDATED MEASUREMENTS (BAROS)
 */
float EverestTask::TaskWrapper(EverestData everestData,
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

/**
 * @brief Only done once. Sets pointers for Madgwick
 *     Internal
 */
void EverestTask::MadgwickSetup() {
  // Attaches Madgwick to Everest
  infusion = this->ExternalInitialize();
  ahrs = infusion->getMadAhrs();

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

  infusion->madOffsetInitialise(&offset, REFRESH_RATE);
  infusion->madAhrsInitialise(ahrs);

  // Set AHRS algorithm settings
  madAhrsSettings settings = {
      EarthConventionEnu,
      0.5f,
      2000.0f, /* replace this with actual gyroscope range in degrees/s */
      10.0f,
      10.0f,
      5 * REFRESH_RATE, /* 5 seconds, we ought to use adaptive deltatime here.
                         */
  };

  infusion->madAhrsSetSettings(ahrs, &settings);

  madgwickInitialized = 1;
}

/**
 * @brief Wrapper for Madgwick, does offset calc and passes
 *       data to Madgwick
 *
 *      Internal
 * @param data IMUData_Everest struct
 *
 */
void EverestTask::MadgwickWrapper(IMUData_Everest data) {
  const float timestamp = data.time;
  madVector gyroscope = {data.gyroX, data.gyroY,
                         data.gyroZ};  // data in degrees/s
  madVector accelerometer = {data.accelX, data.accelY,
                             data.accelZ};  // data in g
  madVector mag = {data.magX, data.magY, data.magZ};

  // Update gyroscope offset correction algorithm
  madOffset offset = infusion->getOffset();
  gyroscope = infusion->madOffsetUpdate(&offset, gyroscope);

  if (debug == Secondary || debug == ALL) {
    SOAR_PRINT(
        "Averaged: (%.6f, %.6f, %.6f) deg/s, Accel: (%.6f, %.6f, %.6f)g Time: "
        "%f\n",
        data.gyroX, data.gyroY, data.gyroZ, data.accelX, data.accelY,
        data.accelZ, deltaTime);
    SOAR_PRINT("Mag: (%.6f, %.6f, %.6f) uT\n", mag.axis.x, mag.axis.y,
               mag.axis.z);
  }

  // Update gyroscope AHRS algorithm
  if (mag.axis.x == 0.0f && mag.axis.y == 0.0f && mag.axis.z == 0.0f) {
    infusion->madAhrsUpdateNoMagnetometer(ahrs, gyroscope, accelerometer,
                                          deltaTime);
  } else {
    infusion->madAhrsUpdate(ahrs, gyroscope, accelerometer, mag, deltaTime);
  }

  madEuler euler = infusion->getEuler(ahrs);
  madVector earth = infusion->madAhrsGetEarthAcceleration(ahrs);

  if (debug == Calibration || debug == ALL) {
    SOAR_PRINT("EARTH XYZ (adjusted for gravity): %f %f %f\n", earth.axis.x,
               earth.axis.y, earth.axis.z);
  }

  internalStates = infusion->madAhrsGetInternalStates(infusion->getMadAhrs());
  flags = infusion->madAhrsGetFlags(infusion->getMadAhrs());

  this->state.earthAcceleration = earth.axis.z;

  if (debug == Secondary || debug == ALL) {
    SOAR_PRINT("%f,%d,%.0f,%.0f,%d,%.0f,%d,%d,%d,%d\n",
               internalStates.accelerationError,
               internalStates.accelerometerIgnored,
               internalStates.accelerationRecoveryTrigger,
               internalStates.magneticError, internalStates.magnetometerIgnored,
               internalStates.magneticRecoveryTrigger, flags.initialising,
               flags.angularRateRecovery, flags.accelerationRecovery,
               flags.magneticRecovery);
  }
}

/**
 * @brief Averages IMUs and feeds them to Madgwick wrapper
 *      Should be called every time IMU data is updated
 *
 *    Internal
 */
void EverestTask::IMU_Update(const IMUData_Everest& imu1,
                             const IMUData_Everest& imu2) {
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

  if (isinf(internalIMU_1.accelX) && isinf(internalIMU_1.accelY) &&
      isinf(internalIMU_1.accelZ)) {
    numberOfSamples -= 1;
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

    if (debug == Calibration || debug == ALL) {
      SOAR_PRINT("uncalibrated: %f, %f, %f ", this->internalIMU_1.accelX,
                 this->internalIMU_1.accelY, this->internalIMU_1.accelZ);
    }

    if (debug == Calibration || debug == ALL) {
      SOAR_PRINT("uncalibrated: %f, %f, %f ", this->internalIMU_1.gyroX,
                 this->internalIMU_1.gyroY, this->internalIMU_1.gyroZ);
    }

    this->internalIMU_1.gyroX =
        this->internalIMU_1.gyroX - this->zeroOffsetGyro[0];
    this->internalIMU_1.gyroY =
        this->internalIMU_1.gyroY - this->zeroOffsetGyro[1];
    this->internalIMU_1.gyroZ =
        this->internalIMU_1.gyroZ - this->zeroOffsetGyro[2];

    if (debug == Calibration || debug == ALL) {
      SOAR_PRINT("-> offset (%f, %f, %f) = calibrated gyro (%f, %f, %f)\n",
                 this->zeroOffsetGyro[0], this->zeroOffsetGyro[1],
                 this->zeroOffsetGyro[2], this->internalIMU_1.gyroX,
                 this->internalIMU_1.gyroY, this->internalIMU_1.gyroZ);
    }
  }

  if (isinf(internalIMU_2.accelX) && isinf(internalIMU_2.accelY) &&
      isinf(internalIMU_2.accelZ)) {
    numberOfSamples -= 1;
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

    if (debug == Calibration || debug == ALL) {
      SOAR_PRINT("uncalibrated: %f, %f, %f ", this->internalIMU_2.accelX,
                 this->internalIMU_2.accelY, this->internalIMU_2.accelZ);
    }

    if (debug == Calibration || debug == ALL) {
      SOAR_PRINT("uncalibrated: %f, %f, %f ", this->internalIMU_2.gyroX,
                 this->internalIMU_2.gyroY, this->internalIMU_2.gyroZ);
    }

    this->internalIMU_2.gyroX =
        this->internalIMU_2.gyroX - this->zeroOffsetGyro2[0];
    this->internalIMU_2.gyroY =
        this->internalIMU_2.gyroY - this->zeroOffsetGyro2[1];
    this->internalIMU_2.gyroZ =
        this->internalIMU_2.gyroZ - this->zeroOffsetGyro2[2];

    if (debug == Calibration || debug == ALL) {
      SOAR_PRINT("-> offset (%f, %f, %f) = calibrated gyro2 (%f, %f, %f)\n",
                 this->zeroOffsetGyro2[0], this->zeroOffsetGyro2[1],
                 this->zeroOffsetGyro2[2], this->internalIMU_2.gyroX,
                 this->internalIMU_2.gyroY, this->internalIMU_2.gyroZ);
    }
  }

  if (numberOfSamples == 0) {
    this->state.avgIMU = {imu1.time, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    return;
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
      (this->internalIMU_1.magZ + this->internalIMU_2.magZ) / numberOfSamples;

  averageIMU.time =
      (this->internalIMU_1.time + this->internalIMU_2.time) / numberOfSamples;

#undef averageIMU

  // feed to Madgwick
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
    SOAR_PRINT("Baro1: %f Pa, Baro2: %f Pa\n", baro1.pressure, baro2.pressure);
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
float EverestTask::ExternalUpdate(IMUData_Everest imu1, IMUData_Everest imu2,
                                  BarosData baro1, BarosData baro2) {
  this->IMU_Update(imu1, imu2);

  if (debug == Third || debug == ALL) {
    SOAR_PRINT("After IMU Update IMU Altitude: %f\n",
               this->state.avgIMU.altitude);
  }

  this->Baro_Update(baro1, baro2);

  float finalAlt = this->dynamite();

  if (debug == Dynamite || debug == ALL) {
    SOAR_PRINT("After Dynamite: %f\n", finalAlt);
  }

  // Update altitude list
  this->AltitudeList.secondLastAltitude = this->AltitudeList.lastAltitude;
  this->AltitudeList.lastAltitude = finalAlt;

  return finalAlt;
}

/**
 * @brief (Currently does not work, use final wrapper) Wraps External Update
 * with alignment, returns External Update with aligned data
 */
float EverestTask::AlignedExternalUpdate(IMUData_Everest imu1,
                                         IMUData_Everest imu2, BarosData baro1,
                                         BarosData baro2,
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
    SOAR_PRINT("Unaligned IMU1:(%.6f, %.6f, %.6f)g,(%.6f, %.6f, %.6f)deg/s\n",
               imu1.accelX, imu1.accelY, imu1.accelZ, imu1.gyroX, imu1.gyroY,
               imu1.gyroZ);

    SOAR_PRINT("Unaligned IMU2:(%.6f, %.6f, %.6f)g,(%.6f, %.6f, %.6f)deg/s\n",
               imu2.accelX, imu2.accelY, imu2.accelZ, imu2.gyroX, imu2.gyroY,
               imu2.gyroZ);

    SOAR_PRINT("Alignment: %d\n", alignment);
  }

  // put aligned data into IMUData_Everest struct
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
    SOAR_PRINT("Aligned IMU1:(%.6f, %.6f, %.6f)g,(%.6f, %.6f, %.6f)deg/s\n",
               imu1.accelX, imu1.accelY, imu1.accelZ, imu1.gyroX, imu1.gyroY,
               imu1.gyroZ);

    SOAR_PRINT("Aligned IMU2:(%.6f, %.6f, %.6f)g,(%.6f, %.6f, %.6f)deg/s\n",
               imu2.accelX, imu2.accelY, imu2.accelZ, imu2.gyroX, imu2.gyroY,
               imu2.gyroZ);
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
float EverestTask::deriveForAltitudeIMU(IMUData_Everest avgIMU) {
  float accelerationZ = this->state.earthAcceleration * -9.81;
  float initialVelocity = this->Kinematics.initialVelo;
  float initialAltitude = this->Kinematics.initialAlt;

  // Derive altitude from IMU
  float finalVelocity = initialVelocity + accelerationZ * deltaTime;

  float altitude =
      initialAltitude + (initialVelocity + finalVelocity) * deltaTime / 2.0;

  if (debug == Secondary || debug == ALL) {
    SOAR_PRINT("\nKinematics\n");
    SOAR_PRINT("IMU Initial Altitude: %f\n", initialAltitude);
    SOAR_PRINT("IMU Velocity: %f\n", initialVelocity);
    SOAR_PRINT("IMU Acceleration: %fm/s^2\n", accelerationZ);
    SOAR_PRINT("IMU Delta Time: %f\n", deltaTime);
    SOAR_PRINT("Derived Altitude: %f\n", altitude);
  }

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
float convertToAltitude(float pressure) {
  float seaLevelPressure = 1013.25;  // sea level pressure in hPa
  pressure = pressure / 100.0;       // convert to hPa
  float altitude = 44330.0 * (1.0 - pow(pressure / seaLevelPressure,
                                        1 / 5.2558));  // barometric formula

  // If pressure is less than 100, altitude is 0
  if (pressure < 100) {
    altitude = 0;
  }

  if (debug == Dynamite || debug == ALL) {
    SOAR_PRINT("\nConversion \n");
    SOAR_PRINT("Pressure: %.f hPa, Altitude: %.f m\n", pressure, altitude);
  }

  return altitude;
}

/**
 * @brief Multi-system trust algorithm. Assumes measurements are updated
 * @returns normalised altitude
 *
 * @category Internal | Asynchronous
 */
float EverestTask::dynamite() {
  float IMUAltitude = deriveForAltitudeIMU(this->state.avgIMU);
  this->state.avgIMU.altitude = IMUAltitude;

  if (std::isnan(IMUAltitude) || std::isinf(IMUAltitude)) {
    if (debug == Dynamite || debug == ALL) {
      SOAR_PRINT("WARNING: IMU Altitude is NaN/Inf! Dropping IMU weight.\n");
    }
    this->state.gain_IMU = 0;
    IMUAltitude = 0;
  }

  float BaroAltitude1 = convertToAltitude(this->baro1.pressure);
  this->baro1.altitude = BaroAltitude1;

  float BaroAltitude2 = convertToAltitude(this->baro2.pressure);
  this->baro2.altitude = BaroAltitude2;

  float GPSAltitude = this->everestData.altitudeGPS;

  if (debug == Dynamite || debug == ALL) {
    SOAR_PRINT("\nDynamite\n");
    SOAR_PRINT("Baro1 Altitude: %f\n", BaroAltitude1);
    SOAR_PRINT("Baro2 Altitude: %f\n", BaroAltitude2);
    // SOAR_PRINT("Baro3 Altitude: %f\n", BaroAltitude3);
    // SOAR_PRINT("Real Baro Altitude: %f\n", RealBaroAltitude);
    SOAR_PRINT("IMU Altitude: %f\n", IMUAltitude);
    SOAR_PRINT("GPSAltitude: %f\n", GPSAltitude);
  }

  // if pressure is zero, set gain to zero
  if (this->baro1.pressure == 0) {
    this->state.gain_Baro1 = 0;
  } else if (this->state.gain_Baro1 == 0) {
    // if not zero, set gain to previous gain
    this->state.gain_Baro1 = this->state.prev_gain_Baro1;
  }

  // if not zero, set gain to zero
  if (this->baro2.pressure == 0) {
    this->state.gain_Baro2 = 0;
  } else if (this->state.gain_Baro2 == 0) {
    // if not zero, set gain to previous gain
    this->state.gain_Baro2 = this->state.prev_gain_Baro2;
  }

  // distribute measurements based on gain
  float distributed_IMU_Altitude = IMUAltitude * this->state.gain_IMU;
  float distributed_Baro_Altitude1 = (BaroAltitude1 * this->state.gain_Baro1);
  float distributed_Baro_Altitude2 = (BaroAltitude2 * this->state.gain_Baro2);

  // summation of distributed measurements
  float distributed_Sum = distributed_IMU_Altitude +
                          distributed_Baro_Altitude1 +
                          distributed_Baro_Altitude2;

  if (debug == Dynamite || debug == ALL) {
    SOAR_PRINT("Distributed Sum: %f\n\n", distributed_Sum);
  }

  // summation of gains
  float sumGain =
      this->state.gain_IMU + this->state.gain_Baro1 + this->state.gain_Baro2;

  if (debug == Dynamite || debug == ALL) {
    SOAR_PRINT("Sum Gain: %f\n\n", sumGain);
  }

  // normalised altitude
  float normalised_Altitude = (distributed_Sum) / sumGain + 0.0001;

  if (debug == Dynamite || debug == ALL) {
    SOAR_PRINT("Normalised Altitude: %f\n\n", normalised_Altitude);
  }

  // overrides normalised altitude with GPS altitude if available, Don't!!
  /*if (this->availableMeasurements[4] == 1) {
    std::cout << "GPS Altitude: " << this->everestData.altitudeGPS
              << std::endl;
    std::cout << "Normalised Altitude before: " << normalised_Altitude
              << std::endl;
    normalised_Altitude = this->everestData.altitudeGPS;
    std::cout << "Normalised Altitude after: " << normalised_Altitude
              << std::endl;
  }*/

  // Update Kinematics
  Kinematics.finalAltitude = normalised_Altitude;

  if (debug == Dynamite || debug == ALL) {
    SOAR_PRINT("Final Altitude: %f\n\n", Kinematics.finalAltitude);
  }

  // update velocity
  if (deltaTime > 0.01) {
    Kinematics.initialVelo =
        (Kinematics.finalAltitude - Kinematics.initialAlt) /
        ((deltaTime) + 0.0001);
  }

  if (debug == Dynamite || debug == ALL) {
    SOAR_PRINT("Initial Velocity: %f\n", Kinematics.initialVelo);
  }

  // update altitude
  Kinematics.initialAlt = Kinematics.finalAltitude;

  recalculateGain(normalised_Altitude);

  /*
  if (this->availableMeasurements[4] == 1) {
    updateGainsWithGPS();
  }*/

  // Save the gains that are not zero as previous gains
  // so once we have recovery phase these old gains are used
  if (this->state.gain_IMU != 0) {
    this->state.prev_gain_IMU = this->state.gain_IMU;
  }
  if (this->state.gain_Baro1 != 0) {
    this->state.prev_gain_Baro1 = this->state.gain_Baro1;
  }
  if (this->state.gain_Baro2 != 0) {
    this->state.prev_gain_Baro2 = this->state.gain_Baro2;
  }
  if (debug == Dynamite || debug == ALL) {
    SOAR_PRINT("Previous Gains\n");
    SOAR_PRINT("Prev Gain IMU: %f\n", this->state.prev_gain_IMU);
    SOAR_PRINT("Prev Gain Baro1: %f\n", this->state.prev_gain_Baro1);
    SOAR_PRINT("Prev Gain Baro2: %f\n", this->state.prev_gain_Baro2);
    // SOAR_PRINT("Prev Gain Baro3: %f\n", this->state.prev_gain_Baro3);
    // SOAR_PRINT("Prev Gain Real Baro: %f\n\n",
    // this->state.prev_gain_Real_Baro);
  }

  return normalised_Altitude;
}

/**
 * @brief Updates the gains with GPS data
 *
 */
void EverestTask::updateGainsWithGPS() {
  float gpsAltitude = this->everestData.altitudeGPS;
  float imuGPSDiff = fabsf(this->state.avgIMU.altitude - gpsAltitude);
  float baro1GPSDiff = fabsf(this->baro1.altitude - gpsAltitude);
  float baro2GPSDiff = fabsf(this->baro2.altitude - gpsAltitude);

  float totalDiff = imuGPSDiff + baro1GPSDiff + baro2GPSDiff;

  float gain_IMU = this->state.gain_IMU * (1 - imuGPSDiff / totalDiff);
  float gain_Baro1 = this->state.gain_Baro1 * (1 - baro1GPSDiff / totalDiff);
  float gain_Baro2 = this->state.gain_Baro2 * (1 - baro2GPSDiff / totalDiff);

  if (debug == Dynamite || debug == ALL) {
    SOAR_PRINT("\nUpdate Gains with GPS\n");
    SOAR_PRINT("IMU GPS Diff: %f\n", imuGPSDiff);
    SOAR_PRINT("Baro1 GPS Diff: %f\n", baro1GPSDiff);
    SOAR_PRINT("Baro2 GPS Diff: %f\n", baro2GPSDiff);
    SOAR_PRINT("Total Diff: %f\n", totalDiff);
    SOAR_PRINT("New Gain IMU: %f\n", gain_IMU);
    SOAR_PRINT("New Gain Baro1: %f\n", gain_Baro1);
    SOAR_PRINT("New Gain Baro2: %f\n", gain_Baro2);
  }

  this->state.gain_IMU = gain_IMU;
  this->state.gain_Baro1 = gain_Baro1;
  this->state.gain_Baro2 = gain_Baro2;
}

/**
 * @brief calculation: new gain = 1 / abs(estimate - measurement)
 *
 */
void EverestTask::recalculateGain(float estimate) {
  float gainedEstimate = deriveChangeInVelocityToGetAltitude(
      estimate);  // pre-integrated for altitude

  // epsilon prevents gains from approaching 0 or infinity.
  float epsilon = 100;

  float gain_IMU = 1 / (fabsf(gainedEstimate - this->state.avgIMU.altitude) +
                        epsilon);  // change to previous trusts
  float gain_Baro1 =
      1 / (fabsf(gainedEstimate - this->baro1.altitude) + epsilon);
  float gain_Baro2 =
      1 / (fabsf(gainedEstimate - this->baro2.altitude) + epsilon);

  if (debug == Third || debug == ALL) {
    SOAR_PRINT("\nRecalculate Gain - Before normalization\n");
    SOAR_PRINT("Gain IMU: %f\n", gain_IMU);
    SOAR_PRINT("Gain Baro1: %f\n", gain_Baro1);
    SOAR_PRINT("Gain Baro2: %f\n", gain_Baro2);
    SOAR_PRINT("Gained Estimate: %f\n", gainedEstimate);

    SOAR_PRINT("Altitude (Estimate): %f\n", estimate);
    SOAR_PRINT("IMU Altitude: %f\n", this->state.avgIMU.altitude);
    SOAR_PRINT("Baro1: %f\n", this->baro1.altitude);
    SOAR_PRINT("Baro2: %f\n", this->baro2.altitude);
    SOAR_PRINT("GPS Altitude: %f\n", this->everestData.altitudeGPS);
  }

  // normalise
  this->state.gain_IMU = gain_IMU / (gain_IMU + gain_Baro1 + gain_Baro2);
  this->state.gain_Baro1 = gain_Baro1 / (gain_IMU + gain_Baro1 + gain_Baro2);
  this->state.gain_Baro2 = gain_Baro2 / (gain_IMU + gain_Baro1 + gain_Baro2);

  if (debug == Dynamite || debug == ALL) {
    SOAR_PRINT("\nRecalculate Gain\n");
    SOAR_PRINT("New Gain IMU: %f\n", this->state.gain_IMU);
    SOAR_PRINT("New Gain Baro1: %f\n", this->state.gain_Baro1);
    SOAR_PRINT("New Gain Baro2: %f\n", this->state.gain_Baro2);
    // SOAR_PRINT("New Gain Baro3: %f\n", this->state.gain_Baro3);
    // SOAR_PRINT("New Gain Real Baro: %f\n\n", this->state.gain_Real_Baro);
  }
}

/**
 * @brief Converts the STDs to coefficients
 */
void EverestTask::calculateSTDCoefficients() {
  // calculate standard deviation coefficients
  float std_IMU = this->state.gain_IMU;
  float std_Baro1 = this->state.gain_Baro1;
  float std_Baro2 = this->state.gain_Baro2;

  float sumSTD1 = pow(this->state.gain_IMU, 2) +
                  pow(this->state.gain_Baro1, 2) +
                  pow(this->state.gain_Baro2, 2);

  // normalise
  this->state.std_IMU = pow(std_IMU, 2) / sumSTD1;
  this->state.std_Baro1 = pow(std_Baro1, 2) / sumSTD1;
  this->state.std_Baro2 = pow(std_Baro2, 2) / sumSTD1;

  if (debug == Dynamite || debug == ALL) {
    SOAR_PRINT("\nStandard Deviation Coefficients\n");
    SOAR_PRINT("STD IMU: %f\n", this->state.std_IMU);
    SOAR_PRINT("STD Baro1: %f\n", this->state.std_Baro1);
    SOAR_PRINT("STD Baro2: %f\n", this->state.std_Baro2);
  }
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
float EverestTask::deriveChangeInVelocityToGetAltitude(float estimate) {
  double deltaTimeAverage =
      (this->baro1.deltaTime + this->baro2.deltaTime + deltaTime) / 3.0;

  float velocityZ = (this->AltitudeList.secondLastAltitude -
                     4 * this->AltitudeList.lastAltitude + 3 * estimate) /
                    (2.0 * deltaTimeAverage);

  float newAltitude =
      this->AltitudeList.lastAltitude + velocityZ * deltaTimeAverage;

  if (debug == Dynamite || debug == ALL) {
    SOAR_PRINT("\nDerivative for new gain\n");
    SOAR_PRINT("Velocity: %f\n", velocityZ);
    SOAR_PRINT("New Altitude: %f\n", newAltitude);
    SOAR_PRINT("Delta Time Average: %f\n\n", deltaTimeAverage);
  }

  return newAltitude;
}

/**
 * @brief Getter for finalAltitude
 *
 * @return kinematics struct
 *
 */
float EverestTask::getFinalAltitude() { return Kinematics.finalAltitude; }

/**
 * @brief Average IMUs to feed into alignment function
 */
int EverestTask::averageIMU(IMUData_Everest& imu1, IMUData_Everest& imu2) {
  // average IMU data
  this->counterEverest += 1;

  // average IMU1
  if (counterEverest == 7) {
    avgIMU1Align.gyroX = avgIMU1Align.gyroX / 7;
    avgIMU1Align.gyroY = avgIMU1Align.gyroY / 7;
    avgIMU1Align.gyroZ = avgIMU1Align.gyroZ / 7;

    avgIMU1Align.accelX = avgIMU1Align.accelX / 7;
    avgIMU1Align.accelY = avgIMU1Align.accelY / 7;
    avgIMU1Align.accelZ = avgIMU1Align.accelZ / 7;

    avgIMU1Align.magX = avgIMU1Align.magX / 7;
    avgIMU1Align.magY = avgIMU1Align.magY / 7;
    avgIMU1Align.magZ = avgIMU1Align.magZ / 7;

    // average IMU2
    avgIMU2Align.gyroX = avgIMU2Align.gyroX / 7;
    avgIMU2Align.gyroY = avgIMU2Align.gyroY / 7;
    avgIMU2Align.gyroZ = avgIMU2Align.gyroZ / 7;

    avgIMU2Align.accelX = avgIMU2Align.accelX / 7;
    avgIMU2Align.accelY = avgIMU2Align.accelY / 7;
    avgIMU2Align.accelZ = avgIMU2Align.accelZ / 7;

    avgIMU2Align.magX = avgIMU2Align.magX / 7;
    avgIMU2Align.magY = avgIMU2Align.magY / 7;
    avgIMU2Align.magZ = avgIMU2Align.magZ / 7;

    counterEverest = 0;
    return 1;
  }

  avgIMU1Align.gyroX += imu1.gyroX;
  avgIMU1Align.gyroY += imu1.gyroY;
  avgIMU1Align.gyroZ += imu1.gyroZ;

  avgIMU1Align.accelX += imu1.accelX;
  avgIMU1Align.accelY += imu1.accelY;
  avgIMU1Align.accelZ += imu1.accelZ;

  avgIMU1Align.magX += imu1.magX;
  avgIMU1Align.magY += imu1.magY;
  avgIMU1Align.magZ += imu1.magZ;

  // average IMU2
  avgIMU2Align.gyroX += imu2.gyroX;
  avgIMU2Align.gyroY += imu2.gyroY;
  avgIMU2Align.gyroZ += imu2.gyroZ;

  avgIMU2Align.accelX += imu2.accelX;
  avgIMU2Align.accelY += imu2.accelY;
  avgIMU2Align.accelZ += imu2.accelZ;

  avgIMU2Align.magX += imu2.magX;
  avgIMU2Align.magY += imu2.magY;
  avgIMU2Align.magZ += imu2.magZ;

  return 0;
}

int EverestTask::findAlignment(IMUData_Everest& imu1, IMUData_Everest& imu2) {
  // check if averageIMU is ready
  if (averageIMU(imu1, imu2) == 0) {
    return 0;
  }

  MadAxesAlignment alignment1;
  if (avgIMU1Align.accelX < -0.9) {
    // -x -> z (PZ,PY, PX)
    alignment1 = MadAxesAlignmentPZNYPX;
  } else if (avgIMU1Align.accelX > 0.9) {
    // x -> z (PZ,PY, NX)
    // alignment1 = MadAxesAlignmentPZNYNX;
  } else if (avgIMU1Align.accelY < -0.9) {
    // -y -> z (PZ, PX, NY)
    // alignment1 = MadAxesAlignmentPZPXNY;
  } else if (avgIMU1Align.accelY > 0.9) {
    // y -> z (PZ, NX, PY)
    // alignment1 = MadAxesAlignmentPZNXPY;
  }

  MadAxesAlignment alignment2;

  if (avgIMU2Align.accelX < -0.9) {
    // -x -> z (PZ,PY, PX)
    alignment2 = MadAxesAlignmentPZNYPX;
  } else if (avgIMU2Align.accelX > 0.9) {
    // x -> z (PZ,PY, NX)
    // alignment2 = MadAxesAlignmentPZNYNX;
  } else if (avgIMU2Align.accelY < -0.9) {
    // -y -> z (PZ, PX, NY)
    // alignment2 = MadAxesAlignmentPZPXNY;
  } else if (avgIMU2Align.accelY > 0.9) {
    // y -> z (PZ, NX, PY)
    // alignment2 = MadAxesAlignmentPZNXPY;
  }

  this->alignment1 = alignment1;
  this->alignment2 = alignment2;

  return 1;
}

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
 *
 */
void EverestTask::tare(const IMUData_Everest& imu1, const IMUData_Everest& imu2,
                       const BarosData& baro1, const BarosData& baro2) {

  // average pressures
  if (baro1.pressure != 0) {
    sum += convertToAltitude(baro1.pressure);
    numberOfSamples++;

    if (debug == Secondary || debug == ALL) {
      SOAR_PRINT("sum: %f number: %d \n", sum, numberOfSamples);
    }
  }

  if (baro2.pressure != 0) {
    sum += convertToAltitude(baro2.pressure);
    numberOfSamples++;

    if (debug == Secondary || debug == ALL) {
      SOAR_PRINT("sum: %f number: %d \n", sum, numberOfSamples);
    }
  }

  if (!isinf(imu1.accelX) || !isinf(imu1.accelY) || !isinf(imu1.accelZ)) {
    this->zeroOffsetGyro = {this->zeroOffsetGyro[0] + imu1.gyroX,
                            this->zeroOffsetGyro[1] + imu1.gyroY,
                            this->zeroOffsetGyro[2] + imu1.gyroZ};

    imu1SampleCount++;
    if (debug == Calibration || debug == ALL) {
      SOAR_PRINT(
          "zeroOffsetGyro[0]:%f,zeroOffsetGyro[1]:%f,zeroOffsetGyro[2]:%f\n",
          this->zeroOffsetGyro[0], this->zeroOffsetGyro[1],
          this->zeroOffsetGyro[2]);
    }
  }

  if (!isinf(imu2.accelX) || !isinf(imu2.accelY) || !isinf(imu2.accelZ)) {
    this->zeroOffsetGyro2 = {this->zeroOffsetGyro2[0] + imu2.gyroX,
                             this->zeroOffsetGyro2[1] + imu2.gyroY,
                             this->zeroOffsetGyro2[2] + imu2.gyroZ};

    imu2SampleCount++;
    if (debug == Calibration || debug == ALL) {
      SOAR_PRINT(
          "zeroOffsetGyro2[0]:%f,zeroOffsetGyro2[1]:%f,zeroOffsetGyro2[2]:%f\n",
          this->zeroOffsetGyro2[0], this->zeroOffsetGyro2[1],
          this->zeroOffsetGyro2[2]);
    }
  }

  if (debug == Calibration || debug == ALL) {
    SOAR_PRINT("Tare Sum: %f\n", sum);
    SOAR_PRINT("Number of samples %d\n", numberOfSamples);
  }

  if (calibrationTimeRemaining == 0) {
    if (numberOfSamples != 0) {
      this->Kinematics.initialAlt = sum / numberOfSamples;
    } else {
      SOAR_PRINT("No samples collected, tare failed.\n");
    }

    if (imu1SampleCount != 0) {
      this->zeroOffsetGyro = {this->zeroOffsetGyro[0] / imu1SampleCount,
                              this->zeroOffsetGyro[1] / imu1SampleCount,
                              this->zeroOffsetGyro[2] / imu1SampleCount};
    }

    if (imu2SampleCount != 0) {
      this->zeroOffsetGyro2 = {this->zeroOffsetGyro2[0] / imu2SampleCount,
                               this->zeroOffsetGyro2[1] / imu2SampleCount,
                               this->zeroOffsetGyro2[2] / imu2SampleCount};
    }

    tared = true;

    if (debug == Calibration || debug == ALL) {
      SOAR_PRINT("Tare Initial Altitude: %f\n", this->Kinematics.initialAlt);
      SOAR_PRINT(
          "\nCalibration offsets:"
          "gyro(%f,%f,%f),\n  gyro2(%f,%f,%f)\n\n",

          this->zeroOffsetGyro[0], this->zeroOffsetGyro[1],
          this->zeroOffsetGyro[2], this->zeroOffsetGyro2[0],
          this->zeroOffsetGyro2[1], this->zeroOffsetGyro2[2]);
    }
  }

  // call to update time and offsets for these structs
  // IMU_Update(imu1, imu2);

  // keeps track of remaining time for tare
  calibrationTimeRemaining -= 1;
}

/**
 * @brief accel is in m/s -> gs, gyro is passed in dps, pressure is in Pa, real
 * is for altitude ONLY in m Aligns before sending to update
 */
float EverestTask::finalWrapper(
    float accelX1, float accelY1, float accelZ1, float gyroX1, float gyroY1,
    float gyroZ1, float magX1, float magY1, float magZ1, float accelX2,
    float accelY2, float accelZ2, float gyroX2, float gyroY2, float gyroZ2,
    float magX2, float magY2, float magZ2, float pressure1, float pressure2,
    float timeIMU1, float timeIMU2, float timeBaro1, float timeBaro2,
    MadAxesAlignment alignment, MadAxesAlignment alignment2) {
  // converts from m/s to gs
  IMUData_Everest sensorData = {timeIMU1,
                                gyroX1,
                                gyroY1,
                                gyroZ1,
                                (float)(accelX1),
                                (float)(accelY1),
                                (float)(accelZ1),
                                magX1,
                                magY1,
                                magZ1,
                                0};

  IMUData_Everest sensorData2 = {timeIMU2,
                                 gyroX2,
                                 gyroY2,
                                 gyroZ2,
                                 (float)(accelX2),
                                 (float)(accelY2),
                                 (float)(accelZ2),
                                 magX2,
                                 magY2,
                                 magZ2,
                                 0};

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

  if (debug == Secondary || debug == ALL) {
    SOAR_PRINT(
        "Aligned: Gyro: (%.6f, %.6f, %.6f) deg/s, Accel: (%.6f, %.6f,"
        "%.6f) g\n",
        imu1GyroAligned.axis.x, imu1GyroAligned.axis.y, imu1GyroAligned.axis.z,
        imu1AccelAligned.axis.x, imu1AccelAligned.axis.y,
        imu1AccelAligned.axis.z);
  }

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

  float eAltitude = this->ExternalUpdate(sensorData, sensorData2, baro1, baro2);

  return eAltitude;
}

/**
 * @brief initialized Halo, and passes Everest filtered values to HALO
 */
std::vector<float> EverestTask::EverestToHalo(EverestData everestData) {
  float eAltitude =
      this->TaskWrapper(everestData, this->alignment1, this->alignment2);
  float eVelocity = this->Kinematics.initialVelo;
  float eAccelerationZ = (this->state.earthAcceleration - 1) * -9.81;

  std::vector<float> haloData = {0, 0, 0};

  if (tared) {
    // TODO: Shouldn't everestTime be the source of truth?  At some point we
    // should decide. Update HALO
    haloData = halo.Halo_Input(&halo, haloInitialized, eAccelerationZ,
                               eVelocity, eAltitude, everestData.altitudeGPS,
                               timeEverest, deltaTime);
    SOAR_PRINT("%f,%f,%f\n", haloData[0], haloData[1], haloData[2]);
  }

  // Update HALO
  return haloData;
}

// tiny epsilon to prevent NaN.

// update IMU1
void EverestTask::IMU1_Measurements(IMUData_Everest imu1) {
  this->everestData.accelX1 = imu1.accelX;
  this->everestData.accelY1 = imu1.accelY;
  this->everestData.accelZ1 = imu1.accelZ;
  this->everestData.gyroX1 = imu1.gyroX;
  this->everestData.gyroY1 = imu1.gyroY;
  this->everestData.gyroZ1 = imu1.gyroZ;
  this->everestData.magX1 = imu1.magX;
  this->everestData.magY1 = imu1.magY;
  this->everestData.magZ1 = imu1.magZ;
  this->everestData.timeIMU1 = imu1.time + 1e-1f;
  this->availableMeasurements[0] = 1;
}

// update IMU2
void EverestTask::IMU2_Measurements(IMUData_Everest imu2) {
  this->everestData.accelX2 = imu2.accelX;
  this->everestData.accelY2 = imu2.accelY;
  this->everestData.accelZ2 = imu2.accelZ;
  this->everestData.gyroX2 = imu2.gyroX;
  this->everestData.gyroY2 = imu2.gyroY;
  this->everestData.gyroZ2 = imu2.gyroZ;
  this->everestData.magX2 = imu2.magX;
  this->everestData.magY2 = imu2.magY;
  this->everestData.magZ2 = imu2.magZ;
  this->everestData.timeIMU2 = imu2.time + 1e-1f;
  this->availableMeasurements[1] = 1;
}

// update Baro1
void EverestTask::Baro1_Measurements(BarosData baro1) {
  this->everestData.pressure1 = baro1.pressure;
  this->everestData.timeBaro1 = baro1.time + 1e-1f;
  this->availableMeasurements[2] = 1;
}

// update Baro2
void EverestTask::Baro2_Measurements(BarosData baro2) {
  this->everestData.pressure2 = baro2.pressure;
  this->everestData.timeBaro2 = baro2.time + 1e-1f;
  this->availableMeasurements[3] = 1;
}

// update GPS
void EverestTask::GPS_Measurements(float altitude) {
  this->everestData.altitudeGPS = altitude;
  this->availableMeasurements[4] = 1;
}

// Custom rounding function
float roundToDecimalPlaces(double value, int decimalPlaces) {
  double scale = std::pow(10.0, decimalPlaces);
  return std::round(value * scale) / scale;
}

std::vector<float> EverestTask::QueueEverest(float currentTime) {
  if (everestInitialized == 0) {
    // size 0 float indicates not ready. I doubt a union return type would be a
    // good solution here.
    return std::vector<float>();
  }

  updateDeltaTime(currentTime);

  // this code doesn't make sense since deltaTime is defined using timethis->
  // deltaTime will always be equal.
  // TODO: 0 here will be a threshold value for updates. If the filter is
  // updating too fast (doubtful) then we can limit it here.

  if (deltaTime > 0) {
    halo.gpsAvailable = this->availableMeasurements[4];
    if (this->availableMeasurements[0] == 1 &&
        this->availableMeasurements[1] == 1 &&
        this->availableMeasurements[2] == 1 &&
        this->availableMeasurements[3] == 1) {
      std::vector<float> haloData = this->EverestToHalo(this->everestData);
      // do not reset measurements, but leave stale data.

      return haloData;
    } else {
      // available[0] = IMU1, available[1] = IMU2, available[2] = Baro1,
      // available[3] = Baro2
      if (this->availableMeasurements[0] == 0) {
        // set to infinity
        this->everestData.accelX1 = std::numeric_limits<float>::infinity();
        this->everestData.accelY1 = std::numeric_limits<float>::infinity();
        this->everestData.accelZ1 = std::numeric_limits<float>::infinity();
      }

      if (this->availableMeasurements[1] == 0) {
        // set to infinity
        this->everestData.accelX2 = std::numeric_limits<float>::infinity();
        this->everestData.accelY2 = std::numeric_limits<float>::infinity();
        this->everestData.accelZ2 = std::numeric_limits<float>::infinity();
      }

      if (this->availableMeasurements[2] == 0) {
        this->everestData.pressure1 = 0;
      }

      if (this->availableMeasurements[3] == 0) {
        this->everestData.pressure2 = 0;
      }

      std::vector<float> haloData = this->EverestToHalo(this->everestData);

      // do not reset measurements, but leave stale data.

      return haloData;
    }
  }
  // If deltaTime <= 0 fall through: return empty vector to avoid undefined behavior
  return std::vector<float>();
}

void EverestTask::updateDeltaTime(float currentTime) {
  oldTime = timeEverest;
  deltaTime = currentTime - oldTime;
  timeEverest = currentTime;
}

void EverestTask::initEverest() {
  if (madgwickInitialized == 0) this->MadgwickSetup();

  if (this->isAligned == 0) {
    IMUData_Everest imu1 = {this->everestData.timeIMU1,
                            this->everestData.gyroX1,
                            this->everestData.gyroY1,
                            this->everestData.gyroZ1,
                            this->everestData.accelX1,
                            this->everestData.accelY1,
                            this->everestData.accelZ1,
                            this->everestData.magX1,
                            this->everestData.magY1,
                            this->everestData.magZ1,
                            0.0f};

    IMUData_Everest imu2 = {this->everestData.timeIMU2,
                            this->everestData.gyroX2,
                            this->everestData.gyroY2,
                            this->everestData.gyroZ2,
                            this->everestData.accelX2,
                            this->everestData.accelY2,
                            this->everestData.accelZ2,
                            this->everestData.magX2,
                            this->everestData.magY2,
                            this->everestData.magZ2,
                            0.0f};

    this->isAligned = this->findAlignment(imu1, imu2);
  }

  if (!tared) {
    IMUData_Everest imu1 = {this->everestData.timeIMU1,
                            this->everestData.gyroX1,
                            this->everestData.gyroY1,
                            this->everestData.gyroZ1,
                            this->everestData.accelX1,
                            this->everestData.accelY1,
                            this->everestData.accelZ1,
                            this->everestData.magX1,
                            this->everestData.magY1,
                            this->everestData.magZ1,
                            0.0f};

    IMUData_Everest imu2 = {this->everestData.timeIMU2,
                            this->everestData.gyroX2,
                            this->everestData.gyroY2,
                            this->everestData.gyroZ2,
                            this->everestData.accelX2,
                            this->everestData.accelY2,
                            this->everestData.accelZ2,
                            this->everestData.magX2,
                            this->everestData.magY2,
                            this->everestData.magZ2,
                            0.0f};

    BarosData baro1 = {this->everestData.timeBaro1, this->everestData.pressure1,
                       0, 0};
    BarosData baro2 = {this->everestData.timeBaro2, this->everestData.pressure2,
                       0, 0};

    this->tare(imu1, imu2, baro1, baro2);
  }

  if (haloInitialized == false) {
    // Done tareing, initialize once
    if (tared) {
      // Initialize HALO
      halo = HALO();
      // Set initial altitude to 1000 if it is zero (no baros in tareing)
      if (this->Kinematics.initialAlt == 0) {
        this->Kinematics.initialAlt = 1000;
      }
      halo.initializeHALO(this->Kinematics.initialAlt, &halo);
      haloInitialized = true;

      this->everestInitialized = 1;
    }
  }

  return;
}
