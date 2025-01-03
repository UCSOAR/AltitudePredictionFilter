// Altitude estimation using multiple sensors
#include "everestTaskHPP.hpp"
#include <stdio.h>
#include <ctime>
#include <string.h>
#include <direct.h>
#include <cstdio>
#include <cerrno>

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

#include "input_data.cpp"

#define LOGON
#define TIMERON
#define printf(...) ;

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
bool isTared = false;
debug_level debug = ALL;
bool firstSampleAfterCalibration = true;
bool useSTD = false;

// INTERNAL VARIABLES
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
  const float timestamp = data.time;
  madVector gyroscope = {data.gyroX, data.gyroY,
                         data.gyroZ};  // data in degrees/s
  madVector accelerometer = {data.accelX, data.accelY,
                             data.accelZ};  // data in g
  madVector mag = {data.magX, data.magY, data.magZ};

  // Update gyroscope offset correction algorithm
  madOffset offset = infusion->getOffset();
  gyroscope = infusion->madOffsetUpdate(&offset, gyroscope);

  // Calculate delta time (in seconds) to account for gyroscope sample clock
  // error
  float deltaTime = (float)(timestamp - previousTimestamp);
  previousTimestamp = timestamp;
  this->state.deltaTimeIMU = deltaTime;

  if (debug == Secondary || debug == ALL) {
    printf(
        "Averaged: (%.6f, %.6f, %.6f) deg/s, Accel: (%.6f, %.6f, %.6f)g Time: "
        "%f\n",
        data.gyroX, data.gyroY, data.gyroZ, data.accelX, data.accelY,
        data.accelZ, deltaTime);
    printf("Mag: (%.6f, %.6f, %.6f) uT\n", mag.axis.x, mag.axis.y, mag.axis.z);
  }

  // Update gyroscope AHRS algorithm
  infusion->madAhrsUpdate(ahrs, gyroscope, accelerometer, mag, deltaTime);

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
  // fprintf(file, "\n");

  everest.state.earthAcceleration = earth.axis.z;

  if (debug == Secondary || debug == ALL) {
    printf("%f,%d,%.0f,%.0f,%d,%.0f,%d,%d,%d,%d\n",
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
      printf("uncalibrated: %f, %f, %f ", this->internalIMU_1.accelX,
             this->internalIMU_1.accelY, this->internalIMU_1.accelZ);
    }

    this->internalIMU_1.accelX =
        this->internalIMU_1.accelX - this->zeroOffsetAccel[0];
    this->internalIMU_1.accelY =
        this->internalIMU_1.accelY - this->zeroOffsetAccel[1];
    this->internalIMU_1.accelZ =
        this->internalIMU_1.accelZ - this->zeroOffsetAccel[2];

    if (debug == Calibration || debug == ALL) {
      printf("-> offset (%f, %f, %f) = calibrated accel (%f, %f, %f)\n",
             this->zeroOffsetAccel[0], this->zeroOffsetAccel[1],
             this->zeroOffsetAccel[2], this->internalIMU_1.accelX,
             this->internalIMU_1.accelY, this->internalIMU_1.accelZ);
    }

    if (debug == Calibration || debug == ALL) {
      printf("uncalibrated: %f, %f, %f ", this->internalIMU_1.gyroX,
             this->internalIMU_1.gyroY, this->internalIMU_1.gyroZ);
    }

    this->internalIMU_1.gyroX =
        this->internalIMU_1.gyroX - this->zeroOffsetGyro[0];
    this->internalIMU_1.gyroY =
        this->internalIMU_1.gyroY - this->zeroOffsetGyro[1];
    this->internalIMU_1.gyroZ =
        this->internalIMU_1.gyroZ - this->zeroOffsetGyro[2];

    if (debug == Calibration || debug == ALL) {
      printf("-> offset (%f, %f, %f) = calibrated gyro (%f, %f, %f)\n",
             this->zeroOffsetGyro[0], this->zeroOffsetGyro[1],
             this->zeroOffsetGyro[2], this->internalIMU_1.gyroX,
             this->internalIMU_1.gyroY, this->internalIMU_1.gyroZ);
    }
  }

  if (isinf(internalIMU_2.accelX)) {
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
      printf("uncalibrated: %f, %f, %f ", this->internalIMU_2.accelX,
             this->internalIMU_2.accelY, this->internalIMU_2.accelZ);
    }

    this->internalIMU_2.accelX =
        this->internalIMU_2.accelX - this->zeroOffsetAccel2[0];
    this->internalIMU_2.accelY =
        this->internalIMU_2.accelY - this->zeroOffsetAccel2[1];
    this->internalIMU_2.accelZ =
        this->internalIMU_2.accelZ - this->zeroOffsetAccel2[2];

    if (debug == Calibration || debug == ALL) {
      printf("-> offset (%f, %f, %f) = calibrated accel2 (%f, %f, %f)\n",
             this->zeroOffsetAccel2[0], this->zeroOffsetAccel2[1],
             this->zeroOffsetAccel2[2], this->internalIMU_2.accelX,
             this->internalIMU_2.accelY, this->internalIMU_2.accelZ);
    }

    if (debug == Calibration || debug == ALL) {
      printf("uncalibrated: %f, %f, %f ", this->internalIMU_2.gyroX,
             this->internalIMU_2.gyroY, this->internalIMU_2.gyroZ);
    }

    this->internalIMU_2.gyroX =
        this->internalIMU_2.gyroX - this->zeroOffsetGyro2[0];
    this->internalIMU_2.gyroY =
        this->internalIMU_2.gyroY - this->zeroOffsetGyro2[1];
    this->internalIMU_2.gyroZ =
        this->internalIMU_2.gyroZ - this->zeroOffsetGyro2[2];

    if (debug == Calibration || debug == ALL) {
      printf("-> offset (%f, %f, %f) = calibrated gyro2 (%f, %f, %f)\n",
             this->zeroOffsetGyro2[0], this->zeroOffsetGyro2[1],
             this->zeroOffsetGyro2[2], this->internalIMU_2.gyroX,
             this->internalIMU_2.gyroY, this->internalIMU_2.gyroZ);
    }
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
    printf("Baro1: %f Pa, Baro2: %f Pa\n", baro1.pressure, baro2.pressure);
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
    printf("After IMU Update IMU Altitude: %f\n",
           everest.state.avgIMU.altitude);
  }

  everest.Baro_Update(baro1, baro2);

  double finalAlt = everest.dynamite();

  if (debug == Dynamite || debug == ALL) {
    printf("After Dynamite: %f\n", finalAlt);
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
    // printf("Unaligned IMU1: (%.6f, %.6f, %.6f) g, (%.6f, %.6f, %.6f)
    // deg/s\n",
    //     imu1.accelX, imu1.accelY, imu1.accelZ, imu1.gyroX, imu1.gyroY,
    //     imu1.gyroZ);

    // printf("Unaligned IMU2: (%.6f, %.6f, %.6f) g, (%.6f, %.6f, %.6f)
    // deg/s\n",
    //     imu2.accelX, imu2.accelY, imu2.accelZ, imu2.gyroX, imu2.gyroY,
    //     imu2.gyroZ);

    // printf("Alignment: %d\n", alignment);
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

  // summation of distributed measurements
  double distributed_Sum = distributed_IMU_Altitude +
                           distributed_Baro_Altitude1 +
                           distributed_Baro_Altitude2;

  if (debug == Dynamite || debug == ALL) {
    // SOAR_PRINT("Distributed Sum: %f\n\n", distributed_Sum);
  }

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

/**
 * @brief calculation: new gain = 1 / abs(estimate - measurement)
 *
 */
void EverestTask::recalculateGain(double estimate) {
  double gainedEstimate = deriveChangeInVelocityToGetAltitude(
      estimate);  // pre-integrated for altitude

  double gain_IMU =
      1 / fabsf(gainedEstimate -
                this->state.avgIMU.altitude);  // change to previous trusts
  double gain_Baro1 = 1 / fabsf(gainedEstimate - this->baro1.altitude);
  double gain_Baro2 = 1 / fabsf(gainedEstimate - this->baro2.altitude);

  if (debug == Third || debug == ALL) {
    printf("\nRecalculate Gain - Before normalization\n");
    printf("Gain IMU: %f\n", gain_IMU);
    printf("Gain Baro1: %f\n", gain_Baro1);
    printf("Gain Baro2: %f\n", gain_Baro2);
    printf("Gained Estimate: %f\n", gainedEstimate);

    printf("Altitude: %f\n", estimate);
    printf("Baro1: %f\n", this->baro1.altitude);
    printf("Baro2: %f\n", this->baro2.altitude);
  }

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

  if (debug == Dynamite || debug == ALL) {
    printf("\nStandard Deviation Coefficients\n");
    printf("STD IMU: %f\n", this->state.std_IMU);
    printf("STD Baro1: %f\n", this->state.std_Baro1);
    printf("STD Baro2: %f\n", this->state.std_Baro2);
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
  // average pressures
  double average = 0;  // have to do since function is weird
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

  this->zeroOffsetAccel = {this->zeroOffsetAccel[0] + imu1.accelX,
                           this->zeroOffsetAccel[1] + imu1.accelY,
                           this->zeroOffsetAccel[2] + imu1.accelZ};
  this->zeroOffsetGyro = {this->zeroOffsetGyro[0] + imu1.gyroX,
                          this->zeroOffsetGyro[1] + imu1.gyroY,
                          this->zeroOffsetGyro[2] + imu1.gyroZ};

  if (debug == Calibration | debug == ALL) {
    // printf("zeroOffsetAccel[0]: %f, zeroOffsetAccel[1]: %f,
    // zeroOffsetAccel[2]: %f \n", this->zeroOffsetAccel[0],
    // this->zeroOffsetAccel[1], this->zeroOffsetAccel[2]);

    // printf("zeroOffsetGyro[0]: %f, zeroOffsetGyro[1]: %f, zeroOffsetGyro[2]:
    // %f \n", this->zeroOffsetGyro[0], this->zeroOffsetGyro[1],
    // this->zeroOffsetGyro[2]);
  }

  if (debug == Secondary || debug == ALL) {
    // SOAR_PRINT("average: %f number: %d \n", average, numberOfSamples);
  }

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

    if (debug == Calibration || debug == ALL) {
      printf(
          "\nCalibration offsets:\n  accel1(%f,%f,%f),\n accel2(%f,%f,%f),\n"
          "gyro(%f,%f,%f),\n  gyro2(%f,%f,%f)\n\n",
          this->zeroOffsetAccel[0], this->zeroOffsetAccel[1],
          this->zeroOffsetAccel[2], this->zeroOffsetGyro[0],
          this->zeroOffsetGyro[1], this->zeroOffsetGyro[2],
          this->zeroOffsetAccel2[0], this->zeroOffsetAccel2[1],
          this->zeroOffsetAccel2[2], this->zeroOffsetGyro2[0],
          this->zeroOffsetGyro2[1], this->zeroOffsetGyro2[2]);
    }
  }

  // call to update time and offsets for these structs
  IMU_Update(imu1, imu2);

  // keeps track of remaining time for tare
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

  if (debug == Secondary || debug == ALL) {
    printf(
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
#define MAX_LINE_LENGTH 1024

/**
 * Serves to just initialize structs
 */
int main() {
  // Setup Madgwick and attach Madgwick to Everest
  everest.MadgwickSetup();

  HALO halo = HALO();

  MatrixXf Q(3, 3);
  Q << 100, 0, 0, 0, 40, 0, 0, 0, 8;

  // Covariance matrix
  MatrixXf R0(3, 3);
  R0 << 200, 0.5, 0.5, 0.5, 100, 1, 0.5, 1, 10;

  MatrixXf P0(3, 3);
  P0 << 50, 0, 0, 0, 0, 0, 0, 0, 0;

// create directory for results
#ifdef LOGON
  // Define the directory path
  std::string directoryPath = "testSuite/results";

  // Create the directory if it doesn't exist
  if (_mkdir("testSuite") == -1) {
    if (errno != EEXIST) {
      std::cerr << "Error creating directory: testSuite" << std::endl;
      return 1;
    }
  }

  if (_mkdir(directoryPath.c_str()) == -1) {
    if (errno != EEXIST) {
      std::cerr << "Error creating directory: " << directoryPath << std::endl;
      return 1;
    }
  }

  // File names
  std::vector<std::string> fileNames = {"gains.txt",
                                        "predictedValues.txt",
                                        "log.txt",
                                        "nearestScenarios.txt",
                                        "nearestScenariosFormatted.txt",
                                        "sigmaPoints.txt",
                                        "sigmaPoints1.txt",
                                        "sigmaPoints2.txt",
                                        "sigmaPoints3.txt",
                                        "sigmaPoints4.txt",
                                        "sigmaPoints5.txt",
                                        "sigmaPoints6.txt"};

  // Deleting files
  for (size_t i = 0; i < fileNames.size(); ++i) {
    std::string filePath = directoryPath + "/" + fileNames[i];
    if (std::remove(filePath.c_str()) == 0) {
      std::cout << "File deleted successfully: " << fileNames[i] << std::endl;
    } else {
      std::perror("Error deleting file");
    }
  }

  // Creating and writing to files
  std::vector<std::pair<std::string, std::string>> filesToCreate = {
      std::make_pair("predictedValues.txt",
                     "s1_alt,s1_velo,s1_acc,s2_alt,s2_velo,s2_acc,s3_alt,s3_"
                     "velo,s3_acc,s4_alt,s4_velo,s4_acc,s5_alt,s5_velo,s5_acc,"
                     "s6_alt,s6_velo,s6_acc\n"),
      std::make_pair("gains.txt", "gain_IMU,gain_Baro1,gain_Baro2\n"),
      std::make_pair("sigmaPoints.txt", "alt,velo,acc\n"),
      std::make_pair("sigmaPoints1.txt", "alt,velo,acc\n"),
      std::make_pair("sigmaPoints2.txt", "alt,velo,acc\n"),
      std::make_pair("sigmaPoints3.txt", "alt,velo,acc\n"),
      std::make_pair("sigmaPoints4.txt", "alt,velo,acc\n"),
      std::make_pair("sigmaPoints5.txt", "alt,velo,acc\n"),
      std::make_pair("sigmaPoints6.txt", "alt,velo,acc\n"),
      std::make_pair("nearestScenarios.txt",
                     "lowestDistance,secondLowestDistance,firstScenario,"
                     "SecondScenario\n")};

  for (size_t i = 0; i < filesToCreate.size(); ++i) {
    std::string filePath = directoryPath + "/" + filesToCreate[i].first;
    FILE* file = fopen(filePath.c_str(), "a+");
    if (file) {
      fputs(filesToCreate[i].second.c_str(), file);
      fclose(file);
      std::cout << "File written: " << filesToCreate[i].first << std::endl;
    } else {
      fprintf(stderr, "Error opening %s...exiting\n",
              filesToCreate[i].first.c_str());
      exit(1);
    }
  }

  FILE* file = fopen((directoryPath + "/HALO.txt").c_str(),
                     "w+");  // Open the file for writing
  if (!file) {
    fprintf(stderr, "Error opening HALO.txt...exiting\n");
    exit(1);
  }
  fprintf(file,
          "Time,Everest_Alt,Everest_Velo,Everest_Accel,Halo_Alt,Halo_Velo,Halo_"
          "Accel\n");

#endif

  // read first line and preset the deltaTime to timestamp
  char line[MAX_LINE_LENGTH];
  std::clock_t start;
  double duration;

  int howMany = 1;
  int i = 0;

  // create scenarios
  halo.createScenarios(&halo);

  float alt;
  float velo;
  float acc;
  float time1;
  std::vector<float> temp = {0, 0, 0, 0};

  halo.deltaTime = 1 / 3;

  float totalTime = 0;
  std::cout << "taberLaunch size: " << taberLaunch.size() << std::endl;

  for (int i = 0; i < taberLaunch.size(); i++) {
    std::cout << "Loop iteration: " << i << std::endl;
    // Tokenize the line using strtok
    // Parse accelerometer readings (X, Y, Z)
    float time = taberLaunch[i][0];
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

    // Parse pressure readings
    float pressure = baroData[i][1];

    IMUData sensorData = {
        time, gyroX, gyroY, gyroZ, accelX, accelY, accelZ, magX, magY, magZ,
    };

    IMUData sensorData2 = {
        time, gyroX, gyroY, gyroZ, accelX, accelY, accelZ, magX, magY, magZ,
    };

    BarosData baro1 = {time, pressure, 0, 0};

    BarosData baro2 = {time, pressure, 0, 0};

    printf(
        "\n#%d "
        "Sample----------------------------------------------------------------"
        "----------\n\n",
        howMany);

    // Print all sensor readings
    if (debug == RAW || debug == ALL) {
      printf(
          "Raw Time: %.6f s, Gyro: (%.6f, %.6f, %.6f) deg/s, Accel: (%.6f, "
          "%.6f, %.6f) g Pressure: (%.f, %.f, %.f, %.f)\n",
          time, sensorData.gyroX, sensorData.gyroY, sensorData.gyroZ,
          sensorData.accelX, sensorData.accelY, sensorData.accelZ,
          baro1.pressure, baro2.pressure);
    }

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

    if (i == 7) {
      VectorXf X0(3);
      X0 << everest.Kinematics.initialAlt, 0, 0;

      // Initialize with tare / GPS values
      halo.init(X0, P0, Q, R0);
    }

    // start timer for iteration
    start = std::clock();

    if (i > 7) {
      // enter measurements from Everest after tare
      halo.setTime(time);
      halo.setStateVector(eAccelerationZ, eVelocity, eAltitude);

      std::vector<double> unitedStates = {halo.X0[0], halo.X0[1], halo.X0[2]};

#ifdef LOGON
      fprintf(file, "%f,%f,%f,%f,%f,%f,%f\n", time, eAltitude, eVelocity,
              eAccelerationZ, unitedStates[0], unitedStates[1],
              unitedStates[2]);
#endif
    }

    clock_t endTime = std::clock();

    totalTime += endTime - start;

    if (i == taberLaunch.size() - 13) {
      std::cout << "Overall for " << howMany << " samples:\t\t\t\t\t\t\t\t\t"
                << totalTime / (double)CLOCKS_PER_SEC << std::endl;

#ifdef TIMERON
      std::cout << "Update time:\t\t\t\t\t\t\t\t\t\t" << halo.updateTime.count()
                << std::endl;
      std::cout << "Predict time:\t\t\t\t\t\t\t\t\t\t"
                << halo.predictTime.count() << std::endl;

      std::cout << "\tTriangulationTime:\t\t\t\t\t\t\t"
                << halo.triangulationTime.count() << std::endl;
      std::cout << "\tdModeltime:\t\t\t\t\t\t\t\t"
                << halo.dynamicModelTime.count() << std::endl;

      std::cout << "\t\tgetScenarioTime:\t\t\t\t\t"
                << halo.getScenarioTime.count() << std::endl;
      std::cout << "\t\tpPredictionTime:\t\t\t\t\t"
                << halo.PpredictionTime.count() << std::endl;
      std::cout << "\t\tprojErrorTime:\t\t\t\t\t\t"
                << halo.projErrorTime.count() << std::endl;
      std::cout << "\t\tpreMeanTime:\t\t\t\t\t\t" << halo.preMeanTime.count()
                << std::endl;
      std::cout << "\t\tsPointTime:\t\t\t\t\t\t" << halo.sPointTime.count()
                << std::endl;
      std::cout << "\t\tpredictLoopTime:\t\t\t\t\t"
                << halo.predictLoopTime.count() << std::endl;
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
      std::cout << "\t\t\t\t->twoDistancesTime:\t"
                << halo.twoDistancesTime.count() << std::endl;
      std::cout << "\t\t\t\t->emplaceBackTime:\t"
                << halo.emplaceBackTime.count() << std::endl;

      std::cout << "\t\t\t->vectorsTime:\t\t\t\t" << halo.vectorsTime.count()
                << std::endl;
      std::cout << "\t\t\t->push_backTime:\t\t\t" << halo.push_backTime.count()
                << std::endl;

      std::cout << "\t\t\t->->treeCreationTime:\t\t\t\t"
                << (halo.treeCreationTime).count() << std::endl;

#endif
    }
  }

  std::cout << "End of program" << std::endl;

#ifdef LOGON
  fclose(file);
#endif

  return 0;
}
