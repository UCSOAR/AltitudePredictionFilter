/*
 * SensorFrameTask.hpp
 *
 *  Created on: Oct 14, 2026
 *      Author: Aly Masani
 */

#ifndef SensorFrameTask_SensorFrameTask_HPP_
#define SensorFrameTask_SensorFrameTask_HPP_

#include "Task.hpp"
#include "SensorDataTypes.hpp"
#include "PreProcess.hpp"
#include "LaunchDetection.hpp"
#include "StateSelection.hpp"
#include "../../Process/Inc/Dynamite.hpp"
#include "../../Estimate/Inc/Ukf.hpp"
#include "../../Prediction/Inc/ApogeeDetection.hpp"
#include "../../Prediction/Inc/ApogeePredictor.hpp"
#include "../../Prediction/Inc/AirbrakeCheck.hpp"
#include "../../Output/Inc/FilterOutput.hpp"
#include "../../Output/Inc/FilterLog.hpp"

/************************************
 * MACROS AND DEFINES
 ************************************/


/************************************
 * TYPEDEFS
 ************************************/
/**
 * @brief Task specific commands handled by the SensorFrameTask
 */
enum SENSORFRAME_TASK_COMMANDS
{
	SF_CMD_NONE = 0,
	SF_CMD_TOGGLE_RAW_PRINT,	// Toggles printing the raw driver values of every frame
	SF_CMD_RETARE,				// Restarts the pre-process tare (pad only!)
	SF_CMD_PRINT_TARE,			// Prints the tare progress or the tare results
	SF_CMD_PRINT_RAW,			// Prints the latest uncalibrated SensorFrame (SI units)
	SF_CMD_PRINT_CALIBRATED,	// Prints the latest calibrated SensorFrame (body axes)
	SF_CMD_PRINT_ALIGNMENT,		// Prints the axis alignment / pad tilt check
	SF_CMD_PRINT_FILTER_INIT,	// Prints the UKF initial conditions
	SF_CMD_PRINT_LAUNCH,		// Prints the launch detection status
	SF_CMD_PRINT_HELP,			// Lists the APF_ debug commands
	SF_CMD_TOGGLE_ALT_PRINT,	// Toggles printing every source / fused altitude and the launch flag
	SF_CMD_PRINT_STATE,			// Prints the APF internal state
	SF_CMD_PRINT_DYNAMITE,		// Prints the dynamite sources, trust weights and Madgwick attitude
	SF_CMD_PRINT_UKF,			// Prints the UKF state, covariance and scenario match
	SF_CMD_PRINT_APOGEE,		// Prints the apogee detection status (legacy confidences)
	SF_CMD_PRINT_PREDICT,		// Prints the apogee / landing prediction and trajectory
	SF_CMD_PRINT_OUTPUT,		// Prints the last published FilterData
	SF_CMD_TOGGLE_OUT_PRINT,	// Toggles printing every published FilterData
	SF_CMD_TOGGLE_LOG_STREAM,	// Toggles the CSV log stream (raw + calibrated + filter), saved on the PC by Tools/serial_log
	SF_CMD_PRINT_TIMING			// Prints (and restarts) the frame / UKF / predictor CPU time
};

/************************************
 * CLASS DEFINITIONS
 ************************************/
/**
 * @brief Groups one poll of DataBroker sensor messages into a SensorFrame
 *        (SI units + timestamp) and passes it to the PreProcess stage
 */
class SensorFrameTask: public Task
{
	public:
		static SensorFrameTask& Inst() {
			static SensorFrameTask inst;
			return inst;
		}

		void InitTask();

		/**
		 * @brief Handles an "APF_..." UART debug command (called from the DebugTask).
		 *        Only queues the matching command, the work runs in this task
		 * @return false if msg is not a known APF_ command
		 */
		static bool HandleApfCommand(const char* msg);

	protected:
		static void RunTask(void* pvParams) { SensorFrameTask::Inst().Run(pvParams); } // Static Task Interface, passes control to the instance Run();
		void Run(void * pvParams); // Main run code
		void HandleCommand(Command& cm);

	private:
		SensorFrameTask();        // Private constructor
		SensorFrameTask(const SensorFrameTask&);                        // Prevent copy-construction
		SensorFrameTask& operator=(const SensorFrameTask&);             // Prevent assignment

		void HandleTaskCommand(uint16_t taskCommand);
		void HandleSensorData(Command& cm);

		// Process a complete set of measurements
		void ProcessFrame();
		void PrintRawFrame();
		void PrintFrame(const SensorFrame& frame, const char* label) const;
		void PrintApfHelp() const;
		void PrintLogHeader() const;
		void PrintLogFrame(const SensorFrame& raw);

		// Pre-launch calibration / tare stage, receives every SensorFrame
		PreProcess preProcess;

		// Decides which stages run for every frame
		StateSelection stateSelection;

		// PROCESS stage: Madgwick + baro altitude blended into the Everest measurement
		Dynamite dynamite;

		// Detection stage, reads the Everest measurement while in DETECTION
		LaunchDetection launchDetection;

		// ESTIMATE stage: UKF with the simulated flights as process model, runs in LAUNCH
		Ukf ukf;

		// PREDICTION stage (apogee part): switches the scenario model to post apogee
		ApogeeDetection apogeeDetection;

		// PREDICTION stage: apogee / trajectory / landing predictor
		ApogeePredictor apogeePredictor;

		// PREDICTION stage: airbrake level from the predicted apogee (logged only, nothing actuated)
		AirbrakeCheck airbrakeCheck;

		// OUTPUT stage: FilterData for the DataBroker (logged to flash by the LoggingTask)
		FilterOutput filterOutput;

		// Stages after the tare, dispatched by stateSelection
		void RunStages(const SensorFrame& calibrated);

		// Power loss / false launch (FilterLog::PersistedState, saved by the LoggingTask in its own sector)
		void RestoreSavedState();
		void SaveState(uint8_t phase, uint32_t nowMs);
		uint8_t BoardPhase() const;
		void TrackPhaseAndEvents(uint32_t nowMs);
		void CheckFalseLaunch(uint32_t nowMs);

		void ResetFilter(uint32_t nowMs);			// back to a fresh pad tare (re-tare, false launch)
		void QueueEvent(uint8_t code, int32_t value);

		static constexpr uint32_t RESTORE_SETTLE_MS = 1500;		// Dynamite settles before the UKF restarts
		static constexpr uint32_t RESTORED_FLAG_MS = 2000;		// RESTORED flag after the UKF restarts
		static constexpr uint32_t FLIGHT_SAVE_PERIOD_MS = 1000;	// saved state refresh in flight (max altitude for a missed apogee)
		static constexpr float MISSED_APOGEE_VEL = -5.0f;		// m/s: already falling when power came back
		static constexpr float LAUNCH_CONFIRM_ALT_M = 30.0f;	// a board FLIGHT is a launch once it climbs this high

		static constexpr uint32_t CLIP_EVENT_MIN_MS = 1000;
		static constexpr uint8_t EVENT_QUEUE = 4;

		FilterLog::PersistedState persisted{};
		uint8_t savedPhase = FilterLog::PHASE_PAD;
		uint32_t lastSaveMs = 0;
		bool restoredFlight = false;
		bool restorePending = false;		// waiting for FlightTask to take the restore
		uint32_t restoreMs = 0;
		uint32_t ukfStartMs = 0;
		bool inFlightSeen = false;			// FlightTask in FLIGHT (false launch check)
		bool launchConfirmed = false;		// ... and the altitude climbed past LAUNCH_CONFIRM_ALT_M

		ApfState lastApfState = ApfState::PRE_LAUNCH;
		bool lastBaroDistrust = false;
		bool lastClip = false;
		uint32_t lastClipEventMs = 0;
		uint8_t eventCode[EVENT_QUEUE] = {};
		int32_t eventValue[EVENT_QUEUE] = {};
		uint8_t eventCount = 0;

		// CPU time (DWT cycle counter), APF_Timing
		struct CpuTime {
			uint32_t calls = 0;
			uint64_t totalCycles = 0;
			uint32_t maxCycles = 0;
			void Add(uint32_t cycles) { calls++; totalCycles += cycles; if (cycles > maxCycles) maxCycles = cycles; }
		};
		CpuTime frameTime, ukfTime, predictorTime;
		uint32_t timingStartMs = 0;
		void PrintTiming();



		// Latest raw reading from each sensor, indexed by sensor id
		IMUData imu[2]{};
		BaroData baro[2]{};
		MagData mag{};

		// Track which sensors have sent their measurements this poll
		uint8_t arrivedMask = 0;

		// Each bit represents one sensor, same layout as SensorFrame::validMask
		static constexpr uint8_t IMU0_BIT  = SENSOR_FRAME_VALID_IMU0;
		static constexpr uint8_t BARO0_BIT = SENSOR_FRAME_VALID_BARO0;
		static constexpr uint8_t MAG_BIT   = SENSOR_FRAME_VALID_MAG;

		// All five sensors (a frame is processed with whatever arrived, see HandleDataBrokerCommand)
		static constexpr uint8_t ALL_ARRIVED = SENSOR_FRAME_VALID_ALL;

		// Timing
		TickType_t frameTick = 0;        // Tick when IMU0 (first message of a poll) arrived
		TickType_t lastFrameTick = 0;
		bool firstFrame = true;

		// Frame bookkeeping
		uint32_t sequence = 0;
		uint32_t droppedFrames = 0;

		// Latest frames, kept for the APF_ debug commands
		SensorFrame lastRaw{};
		SensorFrame lastCalibrated{};
		bool hasRawFrame = false;
		bool hasCalibratedFrame = false;

		// Debug
		bool rawPrintEnabled = false;		// Continuous print, toggled with APF_RawStream
		TickType_t lastRawPrintTick = 0;
		bool altPrintEnabled = false;		// Continuous dynamite + launch print, toggled with APF_AltStream
		TickType_t lastAltPrintTick = 0;
		bool outPrintEnabled = false;		// Toggled with APF_OutStream
		TickType_t lastOutPrintTick = 0;
		bool logStreamEnabled = false;		// CSV log stream, toggled with APF_LogStream
		TickType_t lastLogTick = 0;
		uint32_t loggedFrames = 0;
};


#endif /* SensorFrameTask_SensorFrameTask_HPP_ */
