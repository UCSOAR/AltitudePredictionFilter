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

// Adaptive covariance inflation: an NIS this large (chi-square, 3 values, 99.9%) means the
// scenario model does not fit this flight right now. P is scaled up so the filter listens to the
// measurement again instead of following a wrong model (e.g. a rocket unlike the sims)
constexpr float UKF_NIS_GATE = 16.27f;
constexpr float UKF_MAX_INFLATION = 1000.0f;

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
	float nis;					// innovation^T S^-1 innovation, averages ~3 when R / P are honest
	float inflation;			// P scale applied this update (1 = model fits)
	uint32_t inflations;		// Updates where the model did not fit and P was inflated
	uint32_t updates;
	uint32_t resets;			// Covariance resets after a failed Cholesky
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
		Mat3 Q{};			// Spectral density, scaled by dt
		Mat3 Rfloor{};

		float weightMean[UKF_SIGMA] = {};
		float weightCov[UKF_SIGMA] = {};
		float gamma = 0.0f;		// sqrt(n + lambda)

		bool initialized = false;
		ScenarioModel scenarioModel;
		UkfStatus status{};
};

#endif /* UKF_UKF_HPP_ */
