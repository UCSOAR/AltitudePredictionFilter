/**
 * @file everest.hpp
 * @brief
 */
#ifndef EVEREST_TASK_HPP
#define EVEREST_TASK_HPP

#include "infusion.hpp"
#include "KDTree.hpp"
#include "FilterState.hpp"
#include "HALO.hpp"
#include "Command.hpp"
#include <stdio.h>
#include <ctime>
#include <string.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

// #ifndef HOME
// #include "C:\Users\harry\Desktop\soar\extra\HALO.hpp"
// #endif

// Definitions
// CHANGE
#define REFRESH_RATE (3)  // replace this with actual sample rate
#define RATE_BARO (3)
#define CALIBRATION_TIME (20)
// How often (ms) to refresh `availableMeasurements` from latest sample times
#ifndef AVAILABLE_MEAS_REFRESH_MS
#define AVAILABLE_MEAS_REFRESH_MS 500
#endif

/* Macros/Enums
   ------------------------------------------------------------*/
enum EVEREST_TASK_COMMANDS { EVEREST_NONE = 0, UPDATE, TEST, RETARE };

/*Defines------------------------------------------------------------------*/
typedef struct {
  float initialVelo;
  float initialAlt;
  float finalAltitude;
} kinematics;

/**
 * @brief Keeps the data from the IMU sensor #1 time, gyroXYZ, accelXYZ, magXYZ
 */
typedef struct {
  float time;
  float gyroX, gyroY, gyroZ;
  float accelX, accelY, accelZ;
  float magX, magY, magZ;
  float altitude;
} IMUData_Everest;

/**
 * @brief Keeps whole system's states, including apogee detection results and
 * confidence values for each system
 */
typedef struct {
  float gain_IMU;
  float gain_Baro1;
  float gain_Baro2;

  // prev gains
  float prev_gain_IMU;
  float prev_gain_Baro1;
  float prev_gain_Baro2;

  float std_IMU;
  float std_Baro1;
  float std_Baro2;
  float std_GPS;

  IMUData_Everest avgIMU;
  float earthAcceleration;
} systemState;

typedef struct {
  float secondLastAltitude;
  float lastAltitude;
} altitudeList;

typedef struct {
  float time;
  float pressure;
  float altitude;

  float deltaTime;
  float previousTime;
} BarosData;

/**
 * @brief Just to carry all data between tasks
 */
typedef struct {
  float timeIMU1;
  float timeIMU2;
  float timeBaro1;
  float timeBaro2;

  float pressure1;
  float pressure2;

  float accelX1;
  float accelY1;
  float accelZ1;
  float gyroX1;
  float gyroY1;
  float gyroZ1;
  float magX1;
  float magY1;
  float magZ1;

  float accelX2;
  float accelY2;
  float accelZ2;
  float gyroX2;
  float gyroY2;
  float gyroZ2;
  float magX2;
  float magY2;
  float magZ2;

  float altitudeGPS;

} EverestData;

class EverestTask {
 public:
  void setFilterState(FILTER_STATE filterState);
  FILTER_STATE getFilterState();

  void IMU_Update(const IMUData_Everest& imu1, const IMUData_Everest& imu2);

  int averageIMU(IMUData_Everest& imu1, IMUData_Everest& imu2);

  int counterEverest = 0;
  IMUData_Everest avgIMU1Align = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  IMUData_Everest avgIMU2Align = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int isAligned = 0;

  int madgwickInitialized = 0;
  int everestInitialized = 0;

  bool firstSampleAfterCalibration = true;
  bool useSTD = false;

// don't use this on a board!
#ifdef TESTING_BUILD
  int openFiles();
#endif

  // INTERNAL VARIABLES
  double calibrationTimeRemaining = CALIBRATION_TIME * RATE_BARO;
  float sum = 0;
  float pressureSum = 0;
  float previousTimestamp = 0;
  bool haloInitialized = false;
  std::vector<float> sumZeroOffsetAccel;
  std::vector<float> sumZeroOffsetAccel2;
  std::vector<float> sumZeroOffsetGyro;
  std::vector<float> sumZeroOffsetGyro2;

  // Instantiate Everest
  madAhrs* ahrs;
  Infusion* infusion;

  float oldTime = 0;
  HALO halo;

  int imu1SampleCount = 0;
  int imu2SampleCount = 0;
  int numberOfSamples = 0;

  madAhrsFlags flags;
  madAhrsInternalStates internalStates;

  Infusion* ExternalInitialize();

  void initEverest();

  static EverestTask& getEverest();

  // Initialize system state
  systemState state = {};

  void Baro_Update(const BarosData& baro1, const BarosData& baro2);

  float dynamite();

  kinematics Kinematics;

  Infusion madgwick;

  altitudeList AltitudeList;

  void recalculateGain(float estimate);

  float deriveChangeInVelocityToGetAltitude(float estimate);

  void MadgwickWrapper(IMUData_Everest data);

  float ExternalUpdate(IMUData_Everest imu1, IMUData_Everest imu2,
                       BarosData baro1, BarosData baro2);

  float deriveForAltitudeIMU(IMUData_Everest avgIMU);

  float AlignedExternalUpdate(IMUData_Everest imu1, IMUData_Everest imu2,
                              BarosData baro1, BarosData baro2,
                              MadAxesAlignment alignment);

  void tare(const IMUData_Everest& imu1, const IMUData_Everest& imu2,
            const BarosData& baro1, const BarosData& baro2);

  void MadgwickSetup();

  void initialize1(systemState& state);

  void calculateSTDCoefficients();

  float TaskWrapper(EverestData everestData, MadAxesAlignment alignment,
                    MadAxesAlignment alignment2);

  float finalWrapper(float accelX1, float accelY1, float accelZ1, float gyroX1,
                     float gyroY1, float gyroZ1, float accelX2, float magX1,
                     float magY1, float magZ1, float accelY2, float accelZ2,
                     float gyroX2, float gyroY2, float gyroZ2, float magX2,
                     float magY2, float magZ2, float pressure1, float pressure2,
                     float timeIMU1, float timeIMU2, float timeBaro1,
                     float timeBaro2, MadAxesAlignment alignment,
                     MadAxesAlignment alignment2);

  bool getIsTared();

  void setIsTare(bool isTare);

  std::vector<float> EverestToHalo(EverestData everestData);

  // Extract data from a DataBroker command and update internal measurement buffers
  void Extract(const Command& cm);

  std::vector<float> QueueEverest(float currentTime);

  // calculate deltaTime and adjust everestTime and oldTime
  void updateDeltaTime(float currentTime);

  std::vector<int> availableMeasurements = {0, 0, 0, 0, 0};

  // last sample times (ms since scheduler start) for IMU1, IMU2, BARO1, BARO2, GPS
  uint32_t lastSampleTimeMs[5] = {0, 0, 0, 0, 0};
  // last time we refreshed `availableMeasurements` (ms)
  uint32_t lastAvailableRefreshMs = 0;

  EverestData everestData{};

  float timeEverest = 0;

  // the change between oldTime and timeEverest. Start value is dummy.
  float deltaTime = 0.333;

  void IMU1_Measurements(IMUData_Everest imu1);
  void IMU2_Measurements(IMUData_Everest imu2);
  void Baro1_Measurements(BarosData baro1);
  void Baro2_Measurements(BarosData baro2);

  int findAlignment(IMUData_Everest& imu1, IMUData_Everest& imu2);

  MadAxesAlignment alignment1, alignment2;

  void GPS_Measurements(float altitude);

  BarosData baro1, baro2;

  void updateGainsWithGPS();

  float getFinalAltitude();

 protected:
  FILTER_STATE filterState = FILTER_STATE::PRE_START;
  IMUData_Everest internalIMU_1, internalIMU_2;

  // our accelerometers are already normalized.
  // std::vector<float> zeroOffsetAccel = {0, 0, 0};
  // std::vector<float> zeroOffsetAccel2 = {0, 0, 0};
  std::vector<float> zeroOffsetGyro = {0, 0, 0};
  std::vector<float> zeroOffsetGyro2 = {0, 0, 0};

 private:
};

#endif
