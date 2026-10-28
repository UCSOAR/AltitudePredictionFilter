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
#include "FlightTask.hpp"
#include "Inc/PrintUtils.hpp"
#include <cstring>
#include <cmath>

/************************************
 * PRIVATE MACROS AND DEFINES
 ************************************/
namespace {
/**
 * @brief The tare may only run while the rocket is on the pad, never in flight.
 *        This board has no rocket state, FlightTask's GROUND before any flight is the pad
 */
bool IsPadState()
{
	const FlightTask& flight = FlightTask::Inst();
	return (flight.GetFlightState() == FlightState::GROUND) && !flight.HasFlown();
}

/**
 * @brief Launch is decided by the board (FlightTask GROUND -> FLIGHT), the APF moves DETECTION -> LAUNCH.
 *        Stays true after landing so the APF never falls back to DETECTION
 */
bool IsBoardLaunched()
{
	const FlightTask& flight = FlightTask::Inst();
	return (flight.GetFlightState() == FlightState::FLIGHT) || flight.HasFlown();
}

/**
 * @brief FlightTask is back in GROUND after a flight, the APF moves LAUNCH -> LANDED
 */
bool IsLandedState()
{
	const FlightTask& flight = FlightTask::Inst();
	return (flight.GetFlightState() == FlightState::GROUND) && flight.HasFlown();
}

// Unit conversions from the driver units to SI
constexpr float MG_TO_MPS2         = 9.80665f / 1000.0f;                  // milli-g -> m/s^2
constexpr float DDPS_TO_RADPS      = (3.14159265f / 180.0f) / 10.0f;      // 0.1 deg/s (IMU drivers) -> rad/s
constexpr float CENTI_DEGC_TO_DEGC = 1.0f / 100.0f;                       // 0.01 degC -> degC
constexpr float MGAUSS_TO_UT       = 0.1f;                                // milligauss -> microtesla

// Minimum time between raw debug prints so the UART can keep up at launch poll rates
constexpr uint32_t RAW_PRINT_MIN_PERIOD_MS = 250;

// CSV log stream: max 10 frames/s. One frame is 4 lines (~450 bytes), the debug UART runs at
// 115200 baud (~11.5 KB/s), so this leaves room for the other prints
constexpr uint32_t LOG_STREAM_MIN_PERIOD_MS = 100;

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
	{"APF_State",      SF_CMD_PRINT_STATE},
	{"APF_Dynamite",   SF_CMD_PRINT_DYNAMITE},
	{"APF_UKF",        SF_CMD_PRINT_UKF},
	{"APF_Apogee",     SF_CMD_PRINT_APOGEE},
	{"APF_Predict",    SF_CMD_PRINT_PREDICT},
	{"APF_Output",     SF_CMD_PRINT_OUTPUT},
	{"APF_OutStream",  SF_CMD_TOGGLE_OUT_PRINT},
	{"APF_LogStream",  SF_CMD_TOGGLE_LOG_STREAM},
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

	// Build the scenario KD-trees once at boot (task stack, not the small main stack), not on
	// the first LAUNCH frame. A heap that is too small now fails on the bench, not in flight
	const size_t heapBefore = xPortGetFreeHeapSize();
	const TickType_t buildStart = xTaskGetTickCount();
	ScenarioModel::BuildTrees();
	SOAR_PRINT("APF: KD-trees built in %lu ms, heap used %lu B, free %lu B\n",
			   static_cast<unsigned long>((xTaskGetTickCount() - buildStart) * portTICK_PERIOD_MS),
			   static_cast<unsigned long>(heapBefore - xPortGetFreeHeapSize()),
			   static_cast<unsigned long>(xPortGetFreeHeapSize()));

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

	case SF_CMD_PRINT_STATE:
		stateSelection.PrintStatus(xTaskGetTickCount() * portTICK_PERIOD_MS);
		break;

	case SF_CMD_PRINT_DYNAMITE:
		if (dynamite.IsInitialized()) {
			dynamite.PrintStatus();
		} else {
			SOAR_PRINT("APF: dynamite starts in DETECTION, tare not done\n");
		}
		break;

	case SF_CMD_PRINT_UKF:
		if (ukf.IsInitialized()) {
			ukf.PrintStatus();
		} else {
			SOAR_PRINT("APF: UKF starts in LAUNCH (state %s)\n", StateSelection::Name(stateSelection.GetState()));
		}
		break;

	case SF_CMD_PRINT_APOGEE:
		apogeeDetection.PrintStatus();
		break;

	case SF_CMD_PRINT_PREDICT:
		apogeePredictor.PrintStatus();
		break;

	case SF_CMD_PRINT_OUTPUT:
		filterOutput.PrintLast();
		break;

	case SF_CMD_TOGGLE_OUT_PRINT:
		outPrintEnabled = !outPrintEnabled;
		SOAR_PRINT("SensorFrame output print %s (every published FilterData, max 4 Hz)\n", outPrintEnabled ? "ON" : "OFF");
		break;

	case SF_CMD_TOGGLE_LOG_STREAM:
		logStreamEnabled = !logStreamEnabled;
		if (logStreamEnabled) {
			loggedFrames = 0;
			PrintLogHeader();
		} else {
			SOAR_PRINT("LOGEND,%lu\n", static_cast<unsigned long>(loggedFrames));
		}
		break;

	case SF_CMD_RETARE:
	{
		// Re-tare is pad only, in flight it would stop calibrated frames and tare on flight data
		if (!IsPadState()) {
			SOAR_PRINT("SensorFrame: re-tare refused, not on the pad (flight state %d, flown %d)\n",
					   static_cast<int>(FlightTask::Inst().GetFlightState()), FlightTask::Inst().HasFlown() ? 1 : 0);
			break;
		}
		preProcess.Reset();
		dynamite.Reset();
		launchDetection.Reset();
		ukf.Reset();
		apogeeDetection.Reset();
		apogeePredictor.Reset();
		filterOutput.Reset();
		stateSelection.Reset(xTaskGetTickCount() * portTICK_PERIOD_MS);
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
		SOAR_PRINT("Launch is decided by the board FlightTask: launched=%d (APF LaunchDetection disabled)\n",
				   IsBoardLaunched() ? 1 : 0);
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

	// IMUs: accel mg -> m/s^2, gyro 0.1 deg/s -> rad/s, temp already degC
	for (uint8_t i = 0; i < 2; i++) {
		frame.imu[i].accel[0] = imu[i].accel.x * MG_TO_MPS2;
		frame.imu[i].accel[1] = imu[i].accel.y * MG_TO_MPS2;
		frame.imu[i].accel[2] = imu[i].accel.z * MG_TO_MPS2;
		frame.imu[i].gyro[0]  = imu[i].gyro.x * DDPS_TO_RADPS;
		frame.imu[i].gyro[1]  = imu[i].gyro.y * DDPS_TO_RADPS;
		frame.imu[i].gyro[2]  = imu[i].gyro.z * DDPS_TO_RADPS;
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
	SensorFrame calibrated{};
	if (preProcess.AddFrame(frame, calibrated, IsPadState())) {
		lastCalibrated = calibrated;
		hasCalibratedFrame = true;
		RunStages(calibrated);
	} else {
		// Tare restarted or blocked: the old calibrated frame no longer matches the tare
		hasCalibratedFrame = false;
	}

	// Launch comes from the board's FlightTask, not the APF LaunchDetection (disabled in RunStages)
	stateSelection.Update(preProcess.IsTared(), IsBoardLaunched(),
						  IsLandedState(), frame.timestampMs);

	if (logStreamEnabled) {
		PrintLogFrame(frame);
	}
}

/**
 * @brief Runs the stages of the current APF state on a calibrated frame
 *        (state transitions diagram, one pass per SensorFrame)
 */
void SensorFrameTask::RunStages(const SensorFrame& calibrated)
{
	switch (stateSelection.GetState())
	{
	case ApfState::PRE_LAUNCH:
		// Tare finished on this frame, DETECTION starts from the next one
		break;

	case ApfState::DETECTION:
		// PROCESS keeps the attitude and pad altitudes current until the board detects launch
		dynamite.Update(calibrated, preProcess.GetTareResult(), false);
		// APF launch detection disabled for now, the board's FlightTask decides launch (IsBoardLaunched)
		// launchDetection.Update(dynamite.GetMeasurement(), calibrated.timestampMs);
		break;

	case ApfState::LAUNCH:
	{
		// PROCESS
		const EverestMeasurement& everest = dynamite.Update(calibrated, preProcess.GetTareResult(), true);

		// ESTIMATE: start on the first LAUNCH frame, then predict + update every frame
		if (!ukf.IsInitialized()) {
			ukf.Init(preProcess.GetFilterInit(), everest);
			SOAR_PRINT("UKF: started at alt=%ld mm vel=%ld mm/s\n",
					   ToMilli(everest.altitude), ToMilli(everest.velocity));
		} else {
			ukf.Predict(calibrated.dt);
			ukf.Update(everest);
		}

		// PREDICTION: apogee detection switches the scenario mode once, the predictor projects ahead
		const float sd[3] = {std::sqrt(ukf.Variance(0)), std::sqrt(ukf.Variance(1)), std::sqrt(ukf.Variance(2))};
		if (apogeeDetection.Update(ukf.State(), sd, calibrated.timestampMs)) {
			ukf.SetAfterApogee(true);
			dynamite.DistrustBarosUntil(calibrated.timestampMs + DYNAMITE_EJECTION_DISTRUST_MS);
		}
		apogeePredictor.Update(ukf, calibrated.timestampMs);
		break;
	}

	case ApfState::LANDED:
	default:
		// Hold the last estimate
		break;
	}

	// OUTPUT: publish (every frame in flight, 1 Hz otherwise). The LoggingTask logs it to flash
	if (filterOutput.Due(stateSelection.GetState(), calibrated.timestampMs)) {
		FilterData data = filterOutput.Build(calibrated.timestampMs, stateSelection.GetState(),
											 preProcess.GetFilterInit().groundAltitude, dynamite, launchDetection,
											 ukf, apogeeDetection, apogeePredictor);
		DataBroker::Publish<FilterData>(&data);
		filterOutput.Published(data);

		const TickType_t tick = xTaskGetTickCount();
		if (outPrintEnabled && (tick - lastOutPrintTick) * portTICK_PERIOD_MS >= RAW_PRINT_MIN_PERIOD_MS) {
			lastOutPrintTick = tick;
			FilterOutput::PrintLine(data);
		}
	}

	// Rate limited so it is safe to leave on at launch poll rates
	const TickType_t now = xTaskGetTickCount();
	if (altPrintEnabled && stateSelection.GetState() != ApfState::PRE_LAUNCH
		&& (now - lastAltPrintTick) * portTICK_PERIOD_MS >= RAW_PRINT_MIN_PERIOD_MS) {
		lastAltPrintTick = now;
		SOAR_PRINT("%s launch=%d ", StateSelection::Name(stateSelection.GetState()), IsBoardLaunched() ? 1 : 0);
		dynamite.PrintLine();
		if (ukf.IsInitialized()) {
			ukf.PrintLine();
		}
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
	SOAR_PRINT("  APF_State APF_Dynamite APF_Launch APF_UKF APF_Apogee APF_Predict APF_Output\n");
	SOAR_PRINT("  APF_Retare (pad only) APF_RawStream APF_AltStream APF_OutStream APF_LogStream (toggles) APF_Help\n");
}

/**
 * @brief Column names of the log stream records. Values are integers in milli SI units
 *        (SOAR_PRINT has no %f): mm/s2, mrad/s, mdegC, mPa, nT, mm, mm/s, ms.
 *        Tools/serial_log/apf_log.py joins the records of one frame into one CSV row
 */
void SensorFrameTask::PrintLogHeader() const
{
	SOAR_PRINT("LOGHDR,LR1,seq,t_ms,imu0_ax,imu0_ay,imu0_az,imu0_gx,imu0_gy,imu0_gz,imu0_T,"
			   "imu1_ax,imu1_ay,imu1_az,imu1_gx,imu1_gy,imu1_gz,imu1_T\n");
	SOAR_PRINT("LOGHDR,LR2,seq,baro0_p,baro0_T,baro1_p,baro1_T,mag_x,mag_y,mag_z,valid\n");
	SOAR_PRINT("LOGHDR,LC,seq,c0_ax,c0_ay,c0_az,c0_gx,c0_gy,c0_gz,c1_ax,c1_ay,c1_az,c1_gx,c1_gy,c1_gz\n");
	SOAR_PRINT("LOGHDR,LF,seq,state,flags,alt,vel,acc,alt_sd,vel_sd,apogee,t_event,dyn_alt,dyn_vel,dyn_sigma,ground\n");
}

/**
 * @brief One frame of the log stream, rate limited to LOG_STREAM_MIN_PERIOD_MS.
 *        LR1 / LR2 raw (sensor axes), LC calibrated (body axes, only once the tare is done),
 *        LF filter output (same values FilterOutput publishes) + the dynamite measurement
 */
void SensorFrameTask::PrintLogFrame(const SensorFrame& raw)
{
	const TickType_t now = xTaskGetTickCount();
	if ((now - lastLogTick) * portTICK_PERIOD_MS < LOG_STREAM_MIN_PERIOD_MS) {
		return;
	}
	lastLogTick = now;
	loggedFrames++;

	const unsigned long seq = static_cast<unsigned long>(raw.sequence);
	const SensorFrameIMU& i0 = raw.imu[0];
	const SensorFrameIMU& i1 = raw.imu[1];
	SOAR_PRINT("LR1,%lu,%lu,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld\n",
			   seq, static_cast<unsigned long>(raw.timestampMs),
			   ToMilli(i0.accel[0]), ToMilli(i0.accel[1]), ToMilli(i0.accel[2]),
			   ToMilli(i0.gyro[0]), ToMilli(i0.gyro[1]), ToMilli(i0.gyro[2]), ToMilli(i0.temp),
			   ToMilli(i1.accel[0]), ToMilli(i1.accel[1]), ToMilli(i1.accel[2]),
			   ToMilli(i1.gyro[0]), ToMilli(i1.gyro[1]), ToMilli(i1.gyro[2]), ToMilli(i1.temp));

	SOAR_PRINT("LR2,%lu,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%u\n", seq,
			   ToMilli(raw.baro[0].pressure), ToMilli(raw.baro[0].temp),
			   ToMilli(raw.baro[1].pressure), ToMilli(raw.baro[1].temp),
			   ToMilli(raw.mag[0]), ToMilli(raw.mag[1]), ToMilli(raw.mag[2]),
			   static_cast<unsigned int>(raw.validMask));

	if (hasCalibratedFrame && lastCalibrated.sequence == raw.sequence) {
		const SensorFrameIMU& c0 = lastCalibrated.imu[0];
		const SensorFrameIMU& c1 = lastCalibrated.imu[1];
		SOAR_PRINT("LC,%lu,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld\n", seq,
				   ToMilli(c0.accel[0]), ToMilli(c0.accel[1]), ToMilli(c0.accel[2]),
				   ToMilli(c0.gyro[0]), ToMilli(c0.gyro[1]), ToMilli(c0.gyro[2]),
				   ToMilli(c1.accel[0]), ToMilli(c1.accel[1]), ToMilli(c1.accel[2]),
				   ToMilli(c1.gyro[0]), ToMilli(c1.gyro[1]), ToMilli(c1.gyro[2]));
	}

	// Filter values exist once the tare is done (PRE_LAUNCH still runs the tare)
	if (stateSelection.GetState() != ApfState::PRE_LAUNCH) {
		const FilterData d = filterOutput.Build(raw.timestampMs, stateSelection.GetState(),
												preProcess.GetFilterInit().groundAltitude, dynamite, launchDetection,
												ukf, apogeeDetection, apogeePredictor);
		const EverestMeasurement& ev = dynamite.GetMeasurement();
		SOAR_PRINT("LF,%lu,%u,%u,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld\n", seq,
				   static_cast<unsigned int>(d.apfState), static_cast<unsigned int>(d.flags),
				   ToMilli(d.altitude), ToMilli(d.velocity), ToMilli(d.acceleration),
				   ToMilli(d.altitudeStd), ToMilli(d.velocityStd), ToMilli(d.predictedApogee), ToMilli(d.timeToEvent),
				   ToMilli(ev.altitude), ToMilli(ev.velocity), ToMilli(ev.sigma), ToMilli(d.groundAltitude));
	}
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

	SOAR_PRINT("  IMU0 a=(%d,%d,%d)mg g=(%d,%d,%d)x0.1dps T=%dC | IMU1 a=(%d,%d,%d)mg g=(%d,%d,%d)x0.1dps T=%dC\n",
			   imu[0].accel.x, imu[0].accel.y, imu[0].accel.z,
			   imu[0].gyro.x, imu[0].gyro.y, imu[0].gyro.z, imu[0].temp,
			   imu[1].accel.x, imu[1].accel.y, imu[1].accel.z,
			   imu[1].gyro.x, imu[1].gyro.y, imu[1].gyro.z, imu[1].temp);

	SOAR_PRINT("  BARO0 p=%luPa T=%d(0.01C) | BARO1 p=%luPa T=%d(0.01C)\n",
			   static_cast<unsigned long>(baro[0].pressure), baro[0].temp,
			   static_cast<unsigned long>(baro[1].pressure), baro[1].temp);
}
