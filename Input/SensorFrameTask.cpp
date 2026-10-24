/*
 * SensorFrameTask.cpp
 *
 *  Created on: Oct 14, 2026
 *      Author: Aly Masani
 */


/**
 ********************************************************************************
 * @file    SensorFrameTask.cpp
 * @author  Aly Masani
 * @date    Oct 14, 2026
 * @brief   Collects one poll of IMU x2, Baro x2 and Mag DataBroker messages,
 *          converts them to SI units, timestamps them and passes the SensorFrame
 *          to the PreProcess stage
 ********************************************************************************
 */

/************************************
 * INCLUDES
 ************************************/
#include "Inc/SensorFrameTask.hpp"
#include "SystemDefines.hpp"
#include "Command.hpp"

#include "DataBroker.hpp"
#include "Task.hpp"
#include "PollingTask.hpp"
#include "Inc/PrintUtils.hpp"
#include <cstring>

/************************************
 * PRIVATE MACROS AND DEFINES
 ************************************/
namespace {
/**
 * @brief The tare may only run while the rocket is on the pad, never in flight
 */
bool IsPadState(RocketState state)
{
	switch (state)
	{
	case RS_PRELAUNCH:
	case RS_FILL:
	case RS_ARM:
	case RS_TEST:
		return true;

	default:
		return false;
	}
}

// Unit conversions from the driver units to SI
constexpr float MG_TO_MPS2         = 9.80665f / 1000.0f;                  // milli-g -> m/s^2
constexpr float MDPS_TO_RADPS      = (3.14159265f / 180.0f) / 1000.0f;    // milli-deg/s -> rad/s
constexpr float CENTI_DEGC_TO_DEGC = 1.0f / 100.0f;                       // 0.01 degC -> degC
constexpr float MGAUSS_TO_UT       = 0.1f;                                // milligauss -> microtesla

// Minimum time between raw debug prints so the UART can keep up at launch poll rates
constexpr uint32_t RAW_PRINT_MIN_PERIOD_MS = 250;

/**
 * @brief UART debug commands. Names must stay <= 16 characters (DebugTask buffer)
 *        and are case sensitive
 */
struct ApfCommand
{
	const char* name;
	SENSORFRAME_TASK_COMMANDS command;
};

constexpr ApfCommand APF_COMMANDS[] = {
	{"APF_Raw",        SF_CMD_PRINT_RAW},
	{"APF_Calibrated", SF_CMD_PRINT_CALIBRATED},
	{"APF_Tare",       SF_CMD_PRINT_TARE},
	{"APF_Alignment",  SF_CMD_PRINT_ALIGNMENT},
	{"APF_FilterInit", SF_CMD_PRINT_FILTER_INIT},
	{"APF_Launch",     SF_CMD_PRINT_LAUNCH},
	{"APF_Retare",     SF_CMD_RETARE},
	{"APF_RawStream",  SF_CMD_TOGGLE_RAW_PRINT},
	{"APF_AltStream",  SF_CMD_TOGGLE_ALT_PRINT},
	{"APF_Help",       SF_CMD_PRINT_HELP},
};
}

/************************************
 * VARIABLES
 ************************************/

/************************************
 * FUNCTION DECLARATIONS
 ************************************/

/************************************
 * FUNCTION DEFINITIONS
 ************************************/
SensorFrameTask::SensorFrameTask():Task(TASK_SENSORFRAME_DEPTH_OBJS)
{

}

/**
 * @brief Initialize the SensorFrameTask
 *        Do not modify this function aside from adding the task name
 */
void SensorFrameTask::InitTask()
{
    // Make sure the task is not already initialized
    SOAR_ASSERT(rtTaskHandle == nullptr, "Cannot initialize SensorFrame task twice");

    BaseType_t rtValue =
        xTaskCreate((TaskFunction_t)SensorFrameTask::RunTask,
            (const char*)"SensorFrameTask",
            (uint16_t)TASK_SENSORFRAME_DEPTH_WORDS,
            (void*)this,
            (UBaseType_t)TASK_SENSORFRAME_PRIORITY,
            (TaskHandle_t*)&rtTaskHandle);

                SOAR_ASSERT(rtValue == pdPASS, "SensorFrameTask::InitTask() - xTaskCreate() failed");

	DataBroker::Subscribe<IMUData>(this);
	DataBroker::Subscribe<BaroData>(this);
	DataBroker::Subscribe<MagData>(this);
	// DataBroker::Subscribe<GPSData>(this);

}

void SensorFrameTask::Run(void * pvParams){

    while (1) {
        /* Process commands in blocking mode */
        Command cm;
        bool res = qEvtQueue->ReceiveWait(cm);
        if(res){

        	HandleCommand(cm);
        }
    }

}

/**
 * @brief Runs in the DebugTask's context, so it only looks up the name and queues
 *        the command. All printing happens later in this task, which owns the data
 */
bool SensorFrameTask::HandleApfCommand(const char* msg)
{
	for (const ApfCommand& apf : APF_COMMANDS) {
		if (strcmp(msg, apf.name) == 0) {
			Command cmd(TASK_SPECIFIC_COMMAND, static_cast<uint16_t>(apf.command));
			Inst().GetEventQueue()->Send(cmd);
			return true;
		}
	}
	return false;
}

void SensorFrameTask::HandleCommand(Command& cm)
{
	switch (cm.GetCommand())
	{
	case TASK_SPECIFIC_COMMAND:
		HandleTaskCommand(cm.GetTaskCommand());
		break;

	case DATA_BROKER_COMMAND:
		HandleSensorData(cm);
		break;

	default:
		break;
	}

	// Release the command's allocated memory
	cm.Reset();
}

/**
 * @brief Handles commands sent from other tasks (e.g. the DebugTask)
 */
void SensorFrameTask::HandleTaskCommand(uint16_t taskCommand)
{
	switch (taskCommand)
	{
	case SF_CMD_TOGGLE_RAW_PRINT:
		rawPrintEnabled = !rawPrintEnabled;
		SOAR_PRINT("SensorFrame raw print %s\n", rawPrintEnabled ? "ON" : "OFF");
		break;

	case SF_CMD_TOGGLE_ALT_PRINT:
		altPrintEnabled = !altPrintEnabled;
		SOAR_PRINT("SensorFrame altitude print %s (prints once the tare is done, needs land)\n",
				   altPrintEnabled ? "ON" : "OFF");
		break;

	case SF_CMD_RETARE:
	{
		// Re-tare is pad only, in flight it would stop calibrated frames and tare on flight data
		const RocketState state = PollingTask::Inst().GetRocketState();
		if (!IsPadState(state)) {
			SOAR_PRINT("SensorFrame: re-tare refused, rocket state %d is not a pad state\n", static_cast<int>(state));
			break;
		}
		preProcess.Reset();
		launchDetection.Reset();
		SOAR_PRINT("SensorFrame: re-tare started\n");
		break;
	}

	case SF_CMD_PRINT_TARE:
		preProcess.PrintTareStatus();
		break;

	case SF_CMD_PRINT_RAW:
		if (hasRawFrame) {
			PrintFrame(lastRaw, "Raw (sensor axes)");
		} else {
			SOAR_PRINT("APF: no frame yet, is polling running (land/launch)?\n");
		}
		break;

	case SF_CMD_PRINT_CALIBRATED:
		if (hasCalibratedFrame) {
			PrintFrame(lastCalibrated, "Calibrated (body axes)");
		} else {
			SOAR_PRINT("APF: no calibrated frame yet, tare not done\n");
		}
		break;

	case SF_CMD_PRINT_ALIGNMENT:
		if (preProcess.IsTared()) {
			preProcess.PrintAlignment();
		} else {
			SOAR_PRINT("APF: alignment is checked when the tare finishes, tare not done\n");
		}
		break;

	case SF_CMD_PRINT_FILTER_INIT:
		if (preProcess.IsTared()) {
			preProcess.PrintFilterInit();
		} else {
			SOAR_PRINT("APF: UKF init is built when the tare finishes, tare not done\n");
		}
		break;

	case SF_CMD_PRINT_LAUNCH:
		SOAR_PRINT("LaunchDetection:\n");
		launchDetection.PrintStatus();
		break;

	case SF_CMD_PRINT_HELP:
		PrintApfHelp();
		break;

	default:
		SOAR_PRINT("SensorFrameTask - Received Unsupported Command {%d}\n", taskCommand);
		break;
	}
}

/**
 * @brief Stores a DataBroker sensor message and processes the frame once
 *        all sensors of the poll have arrived
 */
void SensorFrameTask::HandleSensorData(Command& cm)
{
    DataBrokerMessageTypes messageType = DataBroker::getMessageType(cm);

    switch (messageType)
    {
        case DataBrokerMessageTypes::IMU_DATA:
        {
            IMUData data = DataBroker::ExtractData<IMUData>(cm);
            if (data.id >= 2) {
            	break;
            }

            // IMU0 is always published first, so it marks the start of a new poll.
            // Anything still waiting from the previous poll is an incomplete frame
            if (data.id == 0) {
            	if (arrivedMask != 0) {
            		droppedFrames++;
            		arrivedMask = 0;
            	}
            	frameTick = xTaskGetTickCount();
            }

            imu[data.id] = data;
            arrivedMask |= (IMU0_BIT << data.id);
            break;
        }

        case DataBrokerMessageTypes::BARO_DATA:
        {
            BaroData data = DataBroker::ExtractData<BaroData>(cm);
            if (data.id >= 2) {
            	break;
            }

            baro[data.id] = data;
            arrivedMask |= (BARO0_BIT << data.id);
            break;
        }

        case DataBrokerMessageTypes::MAG_DATA:
        {
        	mag = DataBroker::ExtractData<MagData>(cm);
        	arrivedMask |= MAG_BIT;
        	break;
        }

        default:
        {
            break;
        }
    }

    if (arrivedMask == ALL_ARRIVED) {
    	ProcessFrame();
    	arrivedMask = 0;
    }
}

/**
 * @brief Converts the stored raw readings to SI units, timestamps them
 *        and passes the SensorFrame to the PreProcess stage
 */
void SensorFrameTask::ProcessFrame()
{
	SensorFrame frame{};

	// Header
	frame.sequence = sequence++;
	frame.timestampMs = frameTick * portTICK_PERIOD_MS;
	if (firstFrame) {
		frame.dt = 0.0f;
		firstFrame = false;
	} else {
		frame.dt = static_cast<float>((frameTick - lastFrameTick) * portTICK_PERIOD_MS) / 1000.0f;
	}
	lastFrameTick = frameTick;

	// IMUs: accel mg -> m/s^2, gyro mdps -> rad/s, temp already degC
	for (uint8_t i = 0; i < 2; i++) {
		frame.imu[i].accel[0] = imu[i].accel.x * MG_TO_MPS2;
		frame.imu[i].accel[1] = imu[i].accel.y * MG_TO_MPS2;
		frame.imu[i].accel[2] = imu[i].accel.z * MG_TO_MPS2;
		frame.imu[i].gyro[0]  = imu[i].gyro.x * MDPS_TO_RADPS;
		frame.imu[i].gyro[1]  = imu[i].gyro.y * MDPS_TO_RADPS;
		frame.imu[i].gyro[2]  = imu[i].gyro.z * MDPS_TO_RADPS;
		frame.imu[i].temp     = static_cast<float>(imu[i].temp);
	}

	// Baros: pressure already Pa, temp 0.01 degC -> degC
	for (uint8_t i = 0; i < 2; i++) {
		frame.baro[i].pressure = static_cast<float>(baro[i].pressure);
		frame.baro[i].temp     = baro[i].temp * CENTI_DEGC_TO_DEGC;
	}

	// Mag: milligauss -> microtesla
	frame.mag[0] = mag.magX * MGAUSS_TO_UT;
	frame.mag[1] = mag.magY * MGAUSS_TO_UT;
	frame.mag[2] = mag.magZ * MGAUSS_TO_UT;

	// Everything arrived; an all-zero mag means the magnetometer is not delivering data
	frame.validMask = ALL_ARRIVED;
	if (mag.magX == 0 && mag.magY == 0 && mag.magZ == 0) {
		frame.validMask &= ~MAG_BIT;
	}

	if (rawPrintEnabled) {
		PrintRawFrame();
	}

	lastRaw = frame;
	hasRawFrame = true;

	// PreProcess tares on the first pad frames, then returns calibrated frames
	const RocketState rocketState = PollingTask::Inst().GetRocketState();
	SensorFrame calibrated{};
	if (preProcess.AddFrame(frame, calibrated, IsPadState(rocketState))) {
		lastCalibrated = calibrated;
		hasCalibratedFrame = true;

		// Detection stage. The UKF (filter task) will also take the calibrated frame here
		launchDetection.Update(calibrated, preProcess.GetTareResult());

		// Rate limited so it is safe to leave on at launch poll rates
		const TickType_t now = xTaskGetTickCount();
		if (altPrintEnabled && (now - lastAltPrintTick) * portTICK_PERIOD_MS >= RAW_PRINT_MIN_PERIOD_MS) {
			lastAltPrintTick = now;
			launchDetection.PrintAltitudes();
		}
	} else {
		// Tare restarted or blocked: the old calibrated frame no longer matches the tare
		hasCalibratedFrame = false;
	}
}

/**
 * @brief Prints one SensorFrame in milli SI units (SOAR_PRINT has no %f)
 */
void SensorFrameTask::PrintFrame(const SensorFrame& frame, const char* label) const
{
	SOAR_PRINT("%s #%lu t=%lums dt=%ldms valid=0x%02X\n",
			   label,
			   static_cast<unsigned long>(frame.sequence),
			   static_cast<unsigned long>(frame.timestampMs),
			   ToMilli(frame.dt),
			   static_cast<unsigned int>(frame.validMask));

	for (uint8_t i = 0; i < 2; i++) {
		const SensorFrameIMU& imuFrame = frame.imu[i];
		SOAR_PRINT("  IMU%u a=(%ld,%ld,%ld) mm/s2 g=(%ld,%ld,%ld) mrad/s T=%ld mC\n",
				   static_cast<unsigned int>(i),
				   ToMilli(imuFrame.accel[0]), ToMilli(imuFrame.accel[1]), ToMilli(imuFrame.accel[2]),
				   ToMilli(imuFrame.gyro[0]), ToMilli(imuFrame.gyro[1]), ToMilli(imuFrame.gyro[2]),
				   ToMilli(imuFrame.temp));
	}

	SOAR_PRINT("  BARO0 p=%ld mPa T=%ld mC | BARO1 p=%ld mPa T=%ld mC\n",
			   ToMilli(frame.baro[0].pressure), ToMilli(frame.baro[0].temp),
			   ToMilli(frame.baro[1].pressure), ToMilli(frame.baro[1].temp));

	SOAR_PRINT("  MAG m=(%ld,%ld,%ld) nT\n",
			   ToMilli(frame.mag[0]), ToMilli(frame.mag[1]), ToMilli(frame.mag[2]));
}

void SensorFrameTask::PrintApfHelp() const
{
	SOAR_PRINT("APF commands: APF_Raw APF_Calibrated APF_Tare APF_Alignment APF_FilterInit\n");
	SOAR_PRINT("  APF_Launch APF_Retare (pad only) APF_RawStream (toggle) APF_AltStream (toggle) APF_Help\n");
}

/**
 * @brief Prints the raw driver values that went into the current frame.
 *        Rate limited so it is safe to leave on at launch poll rates
 */
void SensorFrameTask::PrintRawFrame()
{
	const TickType_t now = xTaskGetTickCount();
	if ((now - lastRawPrintTick) * portTICK_PERIOD_MS < RAW_PRINT_MIN_PERIOD_MS) {
		return;
	}
	lastRawPrintTick = now;

	SOAR_PRINT("SF #%lu t=%lums drop=%lu mag=(%ld,%ld,%ld)mG\n",
			   static_cast<unsigned long>(sequence - 1),
			   static_cast<unsigned long>(frameTick * portTICK_PERIOD_MS),
			   static_cast<unsigned long>(droppedFrames),
			   static_cast<long>(mag.magX),
			   static_cast<long>(mag.magY),
			   static_cast<long>(mag.magZ));

	SOAR_PRINT("  IMU0 a=(%d,%d,%d)mg g=(%d,%d,%d)mdps T=%dC | IMU1 a=(%d,%d,%d)mg g=(%d,%d,%d)mdps T=%dC\n",
			   imu[0].accel.x, imu[0].accel.y, imu[0].accel.z,
			   imu[0].gyro.x, imu[0].gyro.y, imu[0].gyro.z, imu[0].temp,
			   imu[1].accel.x, imu[1].accel.y, imu[1].accel.z,
			   imu[1].gyro.x, imu[1].gyro.y, imu[1].gyro.z, imu[1].temp);

	SOAR_PRINT("  BARO0 p=%luPa T=%d(0.01C) | BARO1 p=%luPa T=%d(0.01C)\n",
			   static_cast<unsigned long>(baro[0].pressure), baro[0].temp,
			   static_cast<unsigned long>(baro[1].pressure), baro[1].temp);
}
