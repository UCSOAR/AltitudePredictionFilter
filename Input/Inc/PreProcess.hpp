/*
 * PreProcess.hpp
 *
 *  Created on: Oct 15, 2026
 *      Author: Aly Masani
 */

#ifndef PREPROCESS_PREPROCESS_HPP_
#define PREPROCESS_PREPROCESS_HPP_

#include "SensorFrame.hpp"
#include "TareStats.hpp"
#include "FilterInit.hpp"

/************************************
 * CLASS DEFINITIONS
 ************************************/
/**
 * @brief Pre-launch stage. Order of work:
 *          1. Collect tare samples while on the pad (TareStats)
 *          2. Stationary check, reject moving or invalid windows
 *          3. Calibration results: gyro bias, noise, pressure reference + ground altitude (Barometric)
 *          4. Axis alignment check against gravity (AxisAlignment)
 *          5. UKF initial conditions x0 / P0 / Q / R (FilterInit)
 *        After that every frame is calibrated: gyro bias removed, IMUs rotated into body axes
 */
class PreProcess
{
	public:
		/**
		 * @brief Feed one SensorFrame
		 * @param in          Uncalibrated frame from the SensorFrameTask
		 * @param calibrated  Filled with the calibrated frame once the tare is done
		 * @param tareAllowed true only while the rocket is on the pad. When false no
		 *                    tare samples are collected and a partial tare is discarded
		 * @return true if calibrated holds a valid frame, false while still taring
		 */
		bool AddFrame(const SensorFrame& in, SensorFrame& calibrated, bool tareAllowed);

		// Restart the tare. Only call this on the pad, frames stop until it finishes
		void Reset();

		// Power loss in flight / after landing: use the saved pad tare instead of taring again
		void Restore(const TareResult& saved);

		bool IsTared() const { return tared; }
		const TareResult& GetTareResult() const { return result; }
		const FilterInitState& GetFilterInit() const { return filterInit; }

		// Debug prints (values in milli-units, SOAR_PRINT has no %f)
		void PrintTareStatus() const;
		void PrintTareResult() const;
		void PrintAlignment() const;
		void PrintFilterInit() const;

	private:
		bool FinishTare();
		void ApplyCalibration(const SensorFrame& in, SensorFrame& out) const;

		TareStats tareStats;
		TareResult result{};
		FilterInitState filterInit{};

		bool tared = false;
		bool tareBlockedReported = false;	// Print the "tare blocked" message only once
};

#endif /* PREPROCESS_PREPROCESS_HPP_ */
