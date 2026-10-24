/*
 * PreProcess.cpp
 *
 *  Created on: Oct 15, 2026
 *      Author: Aly Masani
 */

/**
 ********************************************************************************
 * @file    PreProcess.cpp
 * @author  Aly Masani
 * @date    Oct 15, 2026
 * @brief   Pre-launch stage: tare on the pad, then calibrate every SensorFrame.
 *          The work itself lives in TareStats, Barometric, AxisAlignment and FilterInit
 ********************************************************************************
 */

/************************************
 * INCLUDES
 ************************************/
#include "Inc/PreProcess.hpp"
#include "Inc/Barometric.hpp"
#include "Inc/AxisAlignment.hpp"
#include "Inc/PrintUtils.hpp"
#include "SystemDefines.hpp"

/************************************
 * FUNCTION DEFINITIONS
 ************************************/
/**
 * @brief Feed one SensorFrame. While taring the frame is only accumulated,
 *        once tared the calibrated frame is written to calibrated.
 *        Tare samples are only collected while tareAllowed (rocket on the pad),
 *        so a reboot in flight never tares on flight data
 * @return true if calibrated holds a valid frame
 */
bool PreProcess::AddFrame(const SensorFrame& in, SensorFrame& calibrated, bool tareAllowed)
{
	if (tared) {
		ApplyCalibration(in, calibrated);
		return true;
	}

	if (!tareAllowed) {
		// Throw away any partial window, it must be collected entirely on the pad
		if (tareStats.SampleCount() > 0) {
			Reset();
		}
		if (!tareBlockedReported) {
			SOAR_PRINT("PreProcess: tare blocked, rocket is not in a pad state\n");
			tareBlockedReported = true;
		}
		return false;
	}
	tareBlockedReported = false;

	tareStats.Add(in);

	if (tareStats.WindowComplete()) {
		if (!tareStats.IsStill()) {
			SOAR_PRINT("PreProcess: rocket moving during tare, restarting\n");
			Reset();
		} else if (!FinishTare()) {
			Reset();
		}
	}

	return false;
}

/**
 * @brief Clears the tare so it runs again from the next frame
 */
void PreProcess::Reset()
{
	tareStats.Reset();
	result = TareResult{};
	filterInit = FilterInitState{};
	tared = false;
}

/**
 * @brief Turns the tare window into calibration results, checks the alignment
 *        and builds the filter initial conditions
 * @return false if the window has no usable baro or IMU (tare must restart)
 */
bool PreProcess::FinishTare()
{
	tareStats.FillResult(result);

	if (!(result.baroValid[0] || result.baroValid[1]) || !(result.imuValid[0] || result.imuValid[1])) {
		SOAR_PRINT("PreProcess: tare has no valid baro or IMU (baro %d/%d imu %d/%d), restarting\n",
				   result.baroValid[0], result.baroValid[1], result.imuValid[0], result.imuValid[1]);
		return false;
	}

	// Ground altitude reference and altitude noise for each baro
	for (uint8_t i = 0; i < 2; i++) {
		if (result.baroValid[i]) {
			result.groundAltitude[i] = Barometric::PressureToAltitude(result.groundPressure[i]);
			result.altitudeStd[i] = Barometric::AltitudePerPascal(result.groundPressure[i]) * result.pressureStd[i];
		}
	}

	// Axis alignment: where gravity lands in body axes
	AxisAlignment::CheckGravity(result);

	// UKF initial conditions
	filterInit = FilterInit::Build(result);

	tared = true;

	PrintTareResult();
	PrintAlignment();
	PrintFilterInit();
	return true;
}

/**
 * @brief Calibrates a frame with the tare results.
 *        Gyro: bias removed (sensor axes). Then both IMUs are rotated into body axes.
 *        Accel keeps gravity (needed for attitude); accel scale/offset and mag
 *        hard/soft iron are identity until measured
 */
void PreProcess::ApplyCalibration(const SensorFrame& in, SensorFrame& out) const
{
	out = in;

	for (uint8_t i = 0; i < 2; i++) {
		float gyroUnbiased[3];
		for (uint8_t axis = 0; axis < 3; axis++) {
			gyroUnbiased[axis] = in.imu[i].gyro[axis] - result.gyroBias[i][axis];
		}

		AxisAlignment::ApplyToImu(i, gyroUnbiased, out.imu[i].gyro);
		AxisAlignment::ApplyToImu(i, in.imu[i].accel, out.imu[i].accel);
	}
}

void PreProcess::PrintTareStatus() const
{
	if (tared) {
		PrintTareResult();
		return;
	}

	SOAR_PRINT("PreProcess: taring, n=%lu/%lu t=%lu/%lums\n",
			   static_cast<unsigned long>(tareStats.SampleCount()),
			   static_cast<unsigned long>(TARE_MIN_SAMPLES),
			   static_cast<unsigned long>(tareStats.ElapsedMs()),
			   static_cast<unsigned long>(TARE_MIN_DURATION_MS));
}

void PreProcess::PrintTareResult() const
{
	SOAR_PRINT("PreProcess: tare done, n=%lu t=%lums (values in milli-units)\n",
			   static_cast<unsigned long>(result.sampleCount),
			   static_cast<unsigned long>(result.durationMs));

	for (uint8_t i = 0; i < 2; i++) {
		SOAR_PRINT("  BARO%u valid=%d p0=%ld mPa std=%ld mPa alt0=%ld mm altStd=%ld mm\n",
				   static_cast<unsigned int>(i), result.baroValid[i],
				   ToMilli(result.groundPressure[i]), ToMilli(result.pressureStd[i]),
				   ToMilli(result.groundAltitude[i]), ToMilli(result.altitudeStd[i]));
	}

	for (uint8_t i = 0; i < 2; i++) {
		SOAR_PRINT("  IMU%u valid=%d g=(%ld,%ld,%ld) mm/s2 gyroBias=(%ld,%ld,%ld) mrad/s\n",
				   static_cast<unsigned int>(i), result.imuValid[i],
				   ToMilli(result.gravity[i][0]), ToMilli(result.gravity[i][1]), ToMilli(result.gravity[i][2]),
				   ToMilli(result.gyroBias[i][0]), ToMilli(result.gyroBias[i][1]), ToMilli(result.gyroBias[i][2]));

		SOAR_PRINT("  IMU%u accStd=(%ld,%ld,%ld) mm/s2 gyroStd=(%ld,%ld,%ld) mrad/s\n",
				   static_cast<unsigned int>(i),
				   ToMilli(result.accelStd[i][0]), ToMilli(result.accelStd[i][1]), ToMilli(result.accelStd[i][2]),
				   ToMilli(result.gyroStd[i][0]), ToMilli(result.gyroStd[i][1]), ToMilli(result.gyroStd[i][2]));
	}
}

void PreProcess::PrintAlignment() const
{
	SOAR_PRINT("PreProcess: alignment %s (max pad tilt %ld mdeg)\n",
			   result.alignmentOk ? "OK" : "CHECK MOUNTING",
			   ToMilli(AxisAlignment::PAD_MAX_TILT_DEG));

	for (uint8_t i = 0; i < 2; i++) {
		SOAR_PRINT("  IMU%u gBody=(%ld,%ld,%ld) mm/s2 tilt=%ld mdeg\n",
				   static_cast<unsigned int>(i),
				   ToMilli(result.gravityBody[i][0]), ToMilli(result.gravityBody[i][1]), ToMilli(result.gravityBody[i][2]),
				   ToMilli(result.padTiltDeg[i]));
	}
}

void PreProcess::PrintFilterInit() const
{
	SOAR_PRINT("PreProcess: UKF init valid=%d groundAlt=%ld mm x0=(%ld,%ld,%ld)\n",
			   filterInit.valid, ToMilli(filterInit.groundAltitude),
			   ToMilli(filterInit.x0[0]), ToMilli(filterInit.x0[1]), ToMilli(filterInit.x0[2]));

	SOAR_PRINT("  P0 diag (x1000)=(%ld,%ld,%ld) Q/R = legacy HALO values (review)\n",
			   ToMilli(filterInit.P0.m[0][0]), ToMilli(filterInit.P0.m[1][1]), ToMilli(filterInit.P0.m[2][2]));
}
