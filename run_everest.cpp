#define LOGMETRICS
#define MAX_LINE_LENGTH 1024

#include "everest.hpp"
#include "input_data.hpp"
#include "gpsData.hpp"
#include "DataBroker.hpp"
#include "SensorDataTypes.hpp"
#include "run_everest.hpp"
#if defined(TESTING_BUILD) || defined(_WIN32)
//#include <direct.h>
#endif
#include "cmsis_os.h"
#include "FreeRTOS.h"
#include "task.h"


// Internal state and helpers (internal linkage)
namespace {

bool cycleSensors  = false;
bool hasGPS        = true;
bool hasIMU1       = true;
bool hasIMU2       = true;
bool hasMag        = true;
bool hasBaro       = false;

float stationaryTime = 0.0f;
float everest_time   = 0.0f;   // used only in PC test path
float timestamp      = 0.0f;   // used only in PC test path

int sensorThisCycle = 0;

// ---------------------------------------------------------------------------
// Helper: returns altitude from GPS table closest to the given timestamp
// ---------------------------------------------------------------------------
float findClosestTime(float ts)
{
    float altitude = 1000.0f;
    for (int i = 0; i < (int)gpsData1.size(); i++) {
        if (gpsData1[i][0] > ts) {
            altitude = gpsData1[i][1];
            break;
        }
    }
    return altitude;
}

// ---------------------------------------------------------------------------
// Helper: returns current time in seconds from the FreeRTOS tick counter.
// Falls back to everest_time when running outside an RTOS context (PC tests).
// ---------------------------------------------------------------------------
static inline float nowSeconds()
{
#if defined(TESTING_BUILD) || defined(_WIN32)
    return everest_time;
#else
    return static_cast<float>(xTaskGetTickCount()) * portTICK_PERIOD_MS / 1000.0f;
#endif
}

// ---------------------------------------------------------------------------
// Publish one set of sensor readings to the DataBroker.
//   nowMs   – FreeRTOS-derived timestamp (ms) attached to every reading
//   i       – index into taberLaunch / baroData arrays
//   stationary – use zeroed / initial values instead of flight data
// ---------------------------------------------------------------------------
void publishSensorData(uint32_t nowMs, int i, int stationary)
{
    const float ts = static_cast<float>(nowMs) / 1000.0f;

    // ---- raw flight data from tables ----
    float accelX  = taberLaunch[i][1];
    float accelY  = taberLaunch[i][2];
    float accelZ  = taberLaunch[i][3];
    float gyroX   = taberLaunch[i][4];
    float gyroY   = taberLaunch[i][5];
    float gyroZ   = taberLaunch[i][6];
    float magX    = taberLaunch[i][7];
    float magY    = taberLaunch[i][8];
    float magZ    = taberLaunch[i][9];
    float pressure = baroData[i][1];

    if (stationary) {
        accelX = accelY = gyroX = gyroY = gyroZ = 0.0f;
        accelZ    = 1.0f;   // gravity
        magX = magY = magZ = 0.0f;
        pressure  = baroData[0][1];
    }

    // ---- Build publish structs ----
    IMUData imu1{};
    imu1.id      = 0;
    imu1.gyro.x  = hasIMU1 ? static_cast<int16_t>(gyroX)  : 0;
    imu1.gyro.y  = hasIMU1 ? static_cast<int16_t>(gyroY)  : 0;
    imu1.gyro.z  = hasIMU1 ? static_cast<int16_t>(gyroZ)  : 0;
    imu1.accel.x = hasIMU1 ? static_cast<int16_t>(accelX) : 0;
    imu1.accel.y = hasIMU1 ? static_cast<int16_t>(accelY) : 0;
    imu1.accel.z = hasIMU1 ? static_cast<int16_t>(accelZ) : 0;

    IMUData imu2{};
    imu2.id      = 1;
    imu2.gyro.x  = hasIMU2 ? static_cast<int16_t>(gyroX)  : 0;
    imu2.gyro.y  = hasIMU2 ? static_cast<int16_t>(gyroY)  : 0;
    imu2.gyro.z  = hasIMU2 ? static_cast<int16_t>(gyroZ)  : 0;
    imu2.accel.x = hasIMU2 ? static_cast<int16_t>(accelX) : 0;
    imu2.accel.y = hasIMU2 ? static_cast<int16_t>(accelY) : 0;
    imu2.accel.z = hasIMU2 ? static_cast<int16_t>(accelZ) : 0;

    BaroData baro1{};
    baro1.id       = 0;
    baro1.pressure = hasBaro ? static_cast<uint32_t>(pressure) : 0;
    baro1.temp     = 0;

    BaroData baro2{};
    baro2.id       = 1;
    baro2.pressure = hasBaro ? static_cast<uint32_t>(pressure) : 0;
    baro2.temp     = 0;

    float gpsAlt = hasGPS ? findClosestTime(ts) : 0.0f;
    GPSData gps{};
    gps.antennaAltitude_.altitude_ = static_cast<int32_t>(gpsAlt);

    // ---- Publish ----
    if (cycleSensors) {
        // One sensor group published per call; caller invokes repeatedly
        switch (sensorThisCycle) {
            case 0:
                if (hasIMU1) {
                    SOAR_PRINT("run_everest [t=%ums] Publishing IMU1 accel=(%d,%d,%d) gyro=(%d,%d,%d)\n",
                               nowMs,
                               imu1.accel.x, imu1.accel.y, imu1.accel.z,
                               imu1.gyro.x,  imu1.gyro.y,  imu1.gyro.z);
                    DataBroker::Publish<IMUData>(&imu1);
                }
                sensorThisCycle++;
                break;
            case 1:
                if (hasIMU2) {
                    SOAR_PRINT("run_everest [t=%ums] Publishing IMU2 accel=(%d,%d,%d) gyro=(%d,%d,%d)\n",
                               nowMs,
                               imu2.accel.x, imu2.accel.y, imu2.accel.z,
                               imu2.gyro.x,  imu2.gyro.y,  imu2.gyro.z);
                    DataBroker::Publish<IMUData>(&imu2);
                }
                sensorThisCycle++;
                break;
            case 2:
                if (hasBaro) {
                    SOAR_PRINT("run_everest [t=%ums] Publishing BARO1 pressure=%u\n", nowMs, baro1.pressure);
                    DataBroker::Publish<BaroData>(&baro1);
                }
                sensorThisCycle++;
                break;
            case 3:
                if (hasBaro) {
                    SOAR_PRINT("run_everest [t=%ums] Publishing BARO2 pressure=%u\n", nowMs, baro2.pressure);
                    DataBroker::Publish<BaroData>(&baro2);
                }
                sensorThisCycle++;
                break;
            case 4:
                if (hasGPS) {
                    SOAR_PRINT("run_everest [t=%ums] Publishing GPS altitude=%d\n", nowMs, gps.antennaAltitude_.altitude_);
                    DataBroker::Publish<GPSData>(&gps);
                }
                sensorThisCycle = 0;
                break;
            default:
                sensorThisCycle = 0;
                break;
        }
    } else {
        // Publish all sensors at once
        if (hasIMU1) {
            SOAR_PRINT("run_everest [t=%ums] Publishing IMU1 accel=(%d,%d,%d) gyro=(%d,%d,%d)\n",
                       nowMs,
                       imu1.accel.x, imu1.accel.y, imu1.accel.z,
                       imu1.gyro.x,  imu1.gyro.y,  imu1.gyro.z);
            DataBroker::Publish<IMUData>(&imu1);
        }
        if (hasIMU2) {
            SOAR_PRINT("run_everest [t=%ums] Publishing IMU2 accel=(%d,%d,%d) gyro=(%d,%d,%d)\n",
                       nowMs,
                       imu2.accel.x, imu2.accel.y, imu2.accel.z,
                       imu2.gyro.x,  imu2.gyro.y,  imu2.gyro.z);
            DataBroker::Publish<IMUData>(&imu2);
        }
        if (hasBaro) {
            SOAR_PRINT("run_everest [t=%ums] Publishing BAROs pressure=%u\n", nowMs, baro1.pressure);
            DataBroker::Publish<BaroData>(&baro1);
            DataBroker::Publish<BaroData>(&baro2);
        }
        if (hasGPS) {
            SOAR_PRINT("run_everest [t=%ums] Publishing GPS altitude=%d\n", nowMs, gps.antennaAltitude_.altitude_);
            DataBroker::Publish<GPSData>(&gps);
        }
    }
}

// ---------------------------------------------------------------------------
// PC / desktop test harness — not compiled into firmware
// ---------------------------------------------------------------------------
#if defined(TESTING_BUILD) || defined(_WIN32)

int main_test_run_everest()
{
    EverestTask everest = EverestTask();
//    everest.openFiles();

    float deltaTime = AVAILABLE_MEAS_REFRESH_MS / 1000.0f;

    /**** TARE PHASE ****/
    while (!everest.everestInitialized) {
        uint32_t nowMs = static_cast<uint32_t>(everest_time * 1000.0f);
        publishSensorData(nowMs, 0, 1);
        everest_time += deltaTime;
        if (everest.everestInitialized == 1) break;
        everest.updateDeltaTime(everest_time);
        everest.initEverest();
    }

    /**** STATIONARY PHASE ****/
    while (stationaryTime > 0.0f && everest_time <= stationaryTime) {
        everest_time += deltaTime;
        uint32_t nowMs = static_cast<uint32_t>(everest_time * 1000.0f);
        publishSensorData(nowMs, 0, 1);
    }

    /**** FLIGHT DATA ****/
    std::clock_t start;
    float totalTime = 0.0f;

    for (int i = 0; i < (int)taberLaunch.size(); i++) {
        start = std::clock();

        uint32_t nowMs = static_cast<uint32_t>(taberLaunch[i][0] * 1000.0f);

        if (cycleSensors) {
            for (int j = 0; j < 5; j++) {
                publishSensorData(nowMs + static_cast<uint32_t>(j * (deltaTime * 1000.0f / 5.0f)), i, 0);
                everest.QueueEverest(static_cast<float>(nowMs) / 1000.0f + (float)j / 15.0f);
            }
        } else {
            publishSensorData(nowMs, i, 0);
            everest.QueueEverest(static_cast<float>(nowMs) / 1000.0f);
        }

        totalTime += static_cast<float>(std::clock() - start);

        if (i == (int)taberLaunch.size() - 13) {
            std::cout << "Overall time: " << totalTime / (double)CLOCKS_PER_SEC << std::endl;
            break;
        }
    }

#if defined(LOGON) || defined(LOGMETRICS)
    exit(0);
#endif
    return 0;
}

#endif // TESTING_BUILD || _WIN32

// ---------------------------------------------------------------------------
// RTOS injection task — publishes test data at the DataBroker refresh rate,
// yielding between each publish so the scheduler can run other tasks.
// ---------------------------------------------------------------------------
static void RunEverestInjectionTask(void* /*pvParams*/)
{
    const uint32_t delayMs = static_cast<uint32_t>(AVAILABLE_MEAS_REFRESH_MS);

    for (int i = 0; i < (int)taberLaunch.size(); i++) {
        // Timestamp from the live FreeRTOS clock, not from the data table
        uint32_t nowMs = static_cast<uint32_t>(xTaskGetTickCount()) * portTICK_PERIOD_MS;

        publishSensorData(nowMs, i, 0);

        // Yield to scheduler; resumes after delayMs
        osDelay(delayMs);
    }

    // Data exhausted — delete self, FilterTask continues running normally
    SOAR_PRINT("run_everest - injection complete, deleting task\n");
    vTaskDelete(nullptr);
}

} // namespace

// ---------------------------------------------------------------------------
// Public API — call once from system init to start injecting test data
// ---------------------------------------------------------------------------
void StartRunEverestInjection()
{
    BaseType_t ret = xTaskCreate(
        RunEverestInjectionTask,
        "everestInject",
        TASK_FILTER_STACK_DEPTH_WORDS,
        nullptr,
        TASK_FILTER_PRIORITY,
        nullptr
    );
    SOAR_ASSERT(ret == pdPASS, "StartRunEverestInjection - xTaskCreate failed");
}
