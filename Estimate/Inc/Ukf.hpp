/*
 * Ukf.hpp
 *
 *  Created on: Oct 24, 2026
 *      Author: Aly Masani
 */

#ifndef UKF_UKF_HPP_
#define UKF_UKF_HPP_

#include "FilterInit.hpp"
#include "ScenarioModel.hpp"
#include "../../Process/Inc/Dynamite.hpp"

/************************************
 * MACROS AND DEFINES
 ************************************/
constexpr uint8_t UKF_N = 3;					// State: altitude, velocity, acceleration
constexpr uint8_t UKF_SIGMA = 2 * UKF_N + 1;

// Merwe scaled sigma points. alpha = 1, kappa = 0 gives lambda = 0: no negative mean
// weight, which keeps the float maths well behaved for a 3 state filter
constexpr float UKF_ALPHA = 1.0f;
constexpr float UKF_BETA = 2.0f;
constexpr float UKF_KAPPA = 0.0f;

// Legacy scenario model (APF_SCENARIO_MODEL_LEGACY): the state snaps to the scenarios every step, so
// the model error is per STEP, not per second. The legacy Q was measured per 1/3 s step (HALO ran at
// 3 Hz) and FilterInit stores it per second (x3): Q * this is added once per predict at any rate
constexpr float UKF_LEGACY_Q_STEP_S = 1.0f / 3.0f;

/************************************
 * TYPEDEFS
 ************************************/
/**
 * @brief Latest UKF values, kept for debug prints and the Output stage
 */
struct UkfStatus
{
	ScenarioMatch scenario;		// Scenario of the propagated mean
	float innovation[UKF_N];	// Measurement - prediction
	float gain[UKF_N];			// Kalman gain diagonal K[i][i]: 1 = follows the measurement, 0 = follows the model
	uint32_t updates;
	uint32_t resets;			// Covariance resets after a failed Cholesky
	uint8_t gainResets;			// Legacy model: sigma points whose scenario pair changed this predict (blended 50 / 50)
	uint32_t gainResetsTotal;
#ifdef APF_PC_REPLAY
	// Replay diagnostics (Tools/replay/q_validation): last update / predict
	float S[UKF_N][UKF_N];		// Innovation covariance P- + R
	float R[UKF_N][UKF_N];		// Effective R (floor, altitude raised to the Dynamite sigma)
	float K[UKF_N][UKF_N];		// Full Kalman gain
	float nis;					// innovation^T S^-1 innovation
	float qUsed[UKF_N][UKF_N];	// Q spectral density added in the last predict
	uint8_t qPhase;				// 0 boost, 1 coast, 2 descent (SetQFull)
#endif
};

/************************************
 * CLASS DEFINITIONS
 ************************************/
/**
 * @brief ESTIMATE stage: unscented Kalman filter on [altitude, velocity, acceleration],
 *        pad relative, up positive.
 *          Predict: every sigma point is moved with the scenario model acceleration
 *          Update:  the Everest measurement (alt, vel, acc) from Dynamite
 *        Replaces the legacy HALO UKF without Eigen or the KD-tree
 */
class Ukf
{
	public:
		/**
		 * @brief Start at the launch frame. State from the Everest measurement,
		 *        covariance from the tare (FilterInit P0) or the measurement noise, whichever is larger
		 */
		void Init(const FilterInitState& init, const EverestMeasurement& everest);
		void Reset();

		void Predict(float dt);
		void Update(const EverestMeasurement& everest);

		// Scenario mode switch (pre / post apogee)
		void SetAfterApogee(bool after) { scenarioModel.SetAfterApogee(after); }
		bool IsAfterApogee() const { return scenarioModel.IsAfterApogee(); }

		// The scenario model in use (the apogee predictor projects with a copy of it)
		const ScenarioModel& Scenarios() const { return scenarioModel; }

		// Testing only (leave-one-case-out): the scenario model skips this case, -1 = use all
		void ExcludeScenario(int8_t index) { scenarioModel.ExcludeCase(index); }

#ifdef APF_PC_REPLAY
		// Replay tuning only (Tools/apf_gui.py, before Init): scales the Q / R diagonals (off diagonals by
		// the geometric mean, R after the live dynamite sigma)
		void SetNoiseTuning(const float q[UKF_N], const float r[UKF_N]);
		// Altitude R = max(R, k * Dynamite spread^2): R grows when the sources disagree (k = 0 off)
		// pull = true: k * Dynamite pull^2 instead (how far the sources move the fused value)
		void SetSpreadR(float k, bool pull = false) { spreadK = k; spreadPull = pull; }
		// R floor replaced by this symmetric matrix {r00, r01, r02, r11, r12, r22} (before Init)
		void SetRFloor(const float r[6]);
		// Q diagonal per flight phase {boost alt, vel, acc, coast ..., descent ...} as spectral densities,
		// boost until the UKF acceleration first goes negative, descent after apogee (before Init)
		void SetPhaseQ(const float q[3 * UKF_N]);
		// Full Q spectral density per phase {boost 3x3 row major, coast 3x3, descent 3x3} replacing Q / Qafter
		// as given (no floor), boost until the UKF acceleration first goes negative (before Init)
		void SetQFull(const float q[3 * UKF_N * UKF_N]);
#endif

		bool IsInitialized() const { return initialized; }
		const float* State() const { return x; }
		float Variance(uint8_t i) const { return P[i][i]; }
		float Covariance(uint8_t i, uint8_t j) const { return P[i][j]; }
		const UkfStatus& GetStatus() const { return status; }

		void PrintStatus() const;
		void PrintLine() const;		// One line: state + sigma

	private:
		void Propagate(const float in[UKF_N], float dt, float out[UKF_N]) const;
		void MeasurementNoise(const EverestMeasurement& everest, float R[UKF_N][UKF_N]) const;
		void ResetCovariance();

		float x[UKF_N] = {0.0f, 0.0f, 0.0f};
		float P[UKF_N][UKF_N] = {};
		Mat3 Q{};			// Spectral density, scaled by dt, before apogee
		Mat3 Qafter{};		// ... after apogee
		Mat3 Rfloor{};

		float weightMean[UKF_SIGMA] = {};
		float weightCov[UKF_SIGMA] = {};
		float gamma = 0.0f;		// sqrt(n + lambda)

		bool initialized = false;
		ScenarioModel scenarioModel;
		LegacyGainMemory legacyMemory[UKF_SIGMA];	// Legacy model: one per sigma point
		UkfStatus status{};
#ifdef APF_PC_REPLAY
		float qScale[UKF_N] = {1.0f, 1.0f, 1.0f};
		float rScale[UKF_N] = {1.0f, 1.0f, 1.0f};
		float spreadK = 0.0f;		// SetSpreadR
		bool spreadPull = false;
		bool rFloorOn = false;
		float rFloor[6] = {};
		bool phaseQOn = false;
		float phaseQ[3][UKF_N] = {};
		bool burnedOut = false;
		bool qFullOn = false;
		float qFull[3][UKF_N][UKF_N] = {};
#endif
};

#endif /* UKF_UKF_HPP_ */
