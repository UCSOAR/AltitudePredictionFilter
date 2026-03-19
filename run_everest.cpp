#define LOGMETRICS
#define MAX_LINE_LENGTH 1024

#include "EverestTaskHPP.hpp"
#include "input_data.hpp"
#include "gpsData.hpp"
#include <direct.h>

// test config //
bool cycleSensors = 0;
bool hasGPS = 1;
bool hasIMU1 = 1;
bool hasIMU2 = 1;
bool hasBaro = 1;
float deltaTime = 0.333f;  // note this is still from taberLaunch. maybe we can
                           // find a way to fix it.
float stationaryTime = 0.0f;
///////////////////////////

int sensorThisCycle = 0;

float everest_time = 0.0f;
float timestamp = 0.0f;

float findClosestTime(float timestamp) {
  // cycle through the times until you find one bigger and return one or after
  // before it
  float altitude = 1000;
  int once = 0;
  for (int i = 0; i < gpsData1.size(); i++) {
    if (gpsData1[i][0] > timestamp) {
      altitude = gpsData1[i][1];
      once = 1;
      break;
    }
  }
  return altitude;
}

int getSensorData(EverestTask* everest, int i, int stationary) {
  IMUData_Everest sensorData;
  IMUData_Everest sensorData2;
  BarosData baro1;
  BarosData baro2;
  float gps;

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

  if (stationary) {
    pressure = baroData[0][1];

    sensorData = {
        everest_time, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    };

    sensorData2 = {
        everest_time, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    };

    baro1 = {everest_time, pressure, 0, 0};

    baro2 = {everest_time, pressure, 0, 0};

    gps = findClosestTime((float)everest_time);

  } else {
    // Tokenize the line using strtok
    // Parse accelerometer readings (X, Y, Z)
    timestamp = taberLaunch[i][0];

    sensorData = {
        timestamp, gyroX,  gyroY, gyroZ, accelX,
        accelY,    accelZ, magX,  magY,  magZ,
    };

    sensorData2 = {
        timestamp, gyroX,  gyroY, gyroZ, accelX,
        accelY,    accelZ, magX,  magY,  magZ,
    };

    baro1 = {timestamp, pressure, 0, 0};

    baro2 = {timestamp, pressure, 0, 0};

    gps = findClosestTime(float(timestamp));
  }

  // get measurements (this will be done from subscribing to the sensors)
  // simulates 5 sensor readings being made per 1/3rd sec.
  if (cycleSensors) {
    switch (sensorThisCycle) {
      case 0:
        everest->IMU1_Measurements(sensorData);
        sensorThisCycle++;
        break;
      case 1:
        sensorData2 = {
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
        everest->IMU2_Measurements(sensorData2);
        sensorThisCycle++;
        break;
      case 2:
        baro1 = {timestamp + (2.0f / 15.0f), pressure, 0, 0};
        everest->Baro1_Measurements(baro1);
        sensorThisCycle++;
        break;
      case 3:
        baro2 = {timestamp + (3.0f / 15.0f), pressure, 0, 0};
        everest->Baro2_Measurements(baro2);
        sensorThisCycle++;
        break;
      case 4:
        gps = findClosestTime(float(timestamp + (4.0f / 15.0f)));
        everest->GPS_Measurements(gps);
        sensorThisCycle = 0;
        break;
    }
  } else {
    everest->IMU1_Measurements(sensorData);
    everest->IMU2_Measurements(sensorData2);
    everest->Baro1_Measurements(baro1);
    everest->Baro2_Measurements(baro2);
    everest->GPS_Measurements(gps);
  }
}

/**
 * Serves to just initialize structs
 */
int main() {
  EverestTask everest = EverestTask();
  // open files. moved here and out of madgwick setup.
#if defined(LOGON) || defined(LOGMETRICS)
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
