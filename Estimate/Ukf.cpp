/*
 * Ukf.cpp
 *
 *  Created on: Oct 24, 2026
 *      Author: Aly Masani
 */

/**
 ********************************************************************************
 * @file    Ukf.cpp
 * @author  Aly Masani
 * @date    Oct 24, 2026
 * @brief   3 state unscented Kalman filter (altitude, velocity, acceleration) with the
 *          simulated flights as the process model. Fixed size 3x3 maths, no heap
 ********************************************************************************
 */

/************************************
 * INCLUDES
 ************************************/
#include "Inc/Ukf.hpp"
#include "PrintUtils.hpp"
#include "SystemDefines.hpp"
#include <cmath>

/************************************
 * PRIVATE MACROS AND DEFINES
 ************************************/
namespace {
typedef float M3[UKF_N][UKF_N];

float MaxF(float a, float b)
{
	return (a > b) ? a : b;
}

/**
 * @brief Lower triangular L with L * L^T = A. false if A is not positive definite
 */
bool Cholesky(const M3 A, M3 L)
{
	for (uint8_t i = 0; i < UKF_N; i++) {
		for (uint8_t j = 0; j < UKF_N; j++) {
			L[i][j] = 0.0f;
		}
	}

	for (uint8_t i = 0; i < UKF_N; i++) {
		for (uint8_t j = 0; j <= i; j++) {
			float sum = A[i][j];
			for (uint8_t k = 0; k < j; k++) {
				sum -= L[i][k] * L[j][k];
			}
			if (i == j) {
				if (sum <= 0.0f) {
					return false;
				}
				L[i][i] = std::sqrt(sum);
			} else {
				L[i][j] = sum / L[j][j];
			}
		}
	}
	return true;
}

/**
 * @brief 3x3 inverse by cofactors. false if singular
 */
bool Invert(const M3 A, M3 inv)
{
	const float c00 = A[1][1] * A[2][2] - A[1][2] * A[2][1];
	const float c01 = A[1][2] * A[2][0] - A[1][0] * A[2][2];
	const float c02 = A[1][0] * A[2][1] - A[1][1] * A[2][0];
	const float det = A[0][0] * c00 + A[0][1] * c01 + A[0][2] * c02;
	if (std::fabs(det) < 1e-12f) {
		return false;
	}
	const float s = 1.0f / det;

	inv[0][0] = c00 * s;
	inv[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) * s;
	inv[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * s;
	inv[1][0] = c01 * s;
	inv[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) * s;
	inv[1][2] = (A[0][2] * A[1][0] - A[0][0] * A[1][2]) * s;
	inv[2][0] = c02 * s;
	inv[2][1] = (A[0][1] * A[2][0] - A[0][0] * A[2][1]) * s;
	inv[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) * s;
	return true;
}

void Multiply(const M3 A, const M3 B, M3 out)
{
	for (uint8_t i = 0; i < UKF_N; i++) {
		for (uint8_t j = 0; j < UKF_N; j++) {
			float sum = 0.0f;
			for (uint8_t k = 0; k < UKF_N; k++) {
				sum += A[i][k] * B[k][j];
			}
			out[i][j] = sum;
		}
	}
}

void Symmetrize(M3 A)
{
	for (uint8_t i = 0; i < UKF_N; i++) {
		for (uint8_t j = i + 1; j < UKF_N; j++) {
			const float avg = 0.5f * (A[i][j] + A[j][i]);
			A[i][j] = avg;
			A[j][i] = avg;
		}
	}
}
}

/************************************
 * FUNCTION DEFINITIONS
 ************************************/
void Ukf::Init(const FilterInitState& init, const EverestMeasurement& everest)
{
	Reset();

	// Sigma point weights
	const float n = static_cast<float>(UKF_N);
	const float lambda = UKF_ALPHA * UKF_ALPHA * (n + UKF_KAPPA) - n;
	gamma = std::sqrt(n + lambda);
	weightMean[0] = lambda / (n + lambda);
	weightCov[0] = weightMean[0] + (1.0f - UKF_ALPHA * UKF_ALPHA + UKF_BETA);
	for (uint8_t i = 1; i < UKF_SIGMA; i++) {
		weightMean[i] = 1.0f / (2.0f * (n + lambda));
		weightCov[i] = weightMean[i];
	}
#if APF_SCENARIO_MODEL == APF_SCENARIO_MODEL_LEGACY
	// Legacy HALO::init: the centre point has weight 0 for the covariance too (no beta term)
	weightCov[0] = weightMean[0];
#endif

	Q = init.Q;
	Rfloor = init.R;
#ifdef APF_PC_REPLAY
	for (uint8_t i = 0; i < UKF_N; i++) {
		for (uint8_t j = 0; j < UKF_N; j++) {
			Q.m[i][j] *= std::sqrt(qScale[i] * qScale[j]);
		}
	}
#endif

	// Start where Dynamite says the rocket is at launch detection
	x[0] = everest.altitude;
	x[1] = everest.velocity;
	x[2] = everest.acceleration;

	M3 R;
	MeasurementNoise(everest, R);
	for (uint8_t i = 0; i < UKF_N; i++) {
		P[i][i] = MaxF(init.P0.m[i][i], R[i][i]);
	}

	initialized = true;
}

void Ukf::Reset()
{
	for (uint8_t i = 0; i < UKF_N; i++) {
		x[i] = 0.0f;
		for (uint8_t j = 0; j < UKF_N; j++) {
			P[i][j] = 0.0f;
		}
	}
	scenarioModel.SetAfterApogee(false);
	for (uint8_t s = 0; s < UKF_SIGMA; s++) {
		legacyMemory[s].Reset();
	}
	initialized = false;
	status = UkfStatus{};
}

/**
 * @brief Process model: the scenario acceleration at the sigma point's altitude / velocity
 *        drives the motion over dt (constant acceleration over the step)
 */
void Ukf::Propagate(const float in[UKF_N], float dt, float out[UKF_N]) const
{
	// Legacy KD-tree search on the sigma point's (altitude, velocity, acceleration)
	const ScenarioMatch match = scenarioModel.Match(in[0], in[1], in[2]);
	const float a = match.acceleration;

	out[0] = in[0] + in[1] * dt + 0.5f * a * dt * dt;
	out[1] = in[1] + a * dt;
	out[2] = a;
}

void Ukf::Predict(float dt)
{
	if (!initialized || dt <= 0.0f) {
		return;
	}

	// Sigma points: x, x +/- gamma * columns of chol(P)
	M3 L;
	if (!Cholesky(P, L)) {
		ResetCovariance();
		Cholesky(P, L);
	}

	float sigma[UKF_SIGMA][UKF_N];
	for (uint8_t j = 0; j < UKF_N; j++) {
		sigma[0][j] = x[j];
	}
	for (uint8_t i = 0; i < UKF_N; i++) {
		for (uint8_t j = 0; j < UKF_N; j++) {
			sigma[1 + i][j] = x[j] + gamma * L[j][i];
			sigma[1 + UKF_N + i][j] = x[j] - gamma * L[j][i];
		}
	}

	// Propagate every sigma point through the scenario model
	float moved[UKF_SIGMA][UKF_N];
#if APF_SCENARIO_MODEL == APF_SCENARIO_MODEL_LEGACY
	// Legacy HALO::calculateSigmaPoints -> dynamicModel: each sigma point becomes the blend of its
	// two closest scenarios dt later, with its own gain memory (reset to 50 / 50 on a pair change)
	status.gainResets = 0;
	for (uint8_t s = 0; s < UKF_SIGMA; s++) {
		const LegacyStep step = scenarioModel.LegacyPredict(sigma[s], dt, legacyMemory[s], true);
		for (uint8_t j = 0; j < UKF_N; j++) {
			moved[s][j] = step.state[j];
		}
		if (step.gainsReset) {
			status.gainResets++;
		}
	}
	status.gainResetsTotal += status.gainResets;
	const float qScale = UKF_LEGACY_Q_STEP_S;	// Once per step, see UKF_LEGACY_Q_STEP_S
#else
	for (uint8_t s = 0; s < UKF_SIGMA; s++) {
		Propagate(sigma[s], dt, moved[s]);
	}
	const float qScale = dt;
#endif

	// Predicted mean
	float mean[UKF_N] = {0.0f, 0.0f, 0.0f};
	for (uint8_t s = 0; s < UKF_SIGMA; s++) {
		for (uint8_t j = 0; j < UKF_N; j++) {
			mean[j] += weightMean[s] * moved[s][j];
		}
	}

	// Predicted covariance + process noise
	M3 Pp = {};
	for (uint8_t s = 0; s < UKF_SIGMA; s++) {
		float d[UKF_N];
		for (uint8_t j = 0; j < UKF_N; j++) {
			d[j] = moved[s][j] - mean[j];
		}
		for (uint8_t i = 0; i < UKF_N; i++) {
			for (uint8_t j = 0; j < UKF_N; j++) {
				Pp[i][j] += weightCov[s] * d[i] * d[j];
			}
		}
	}
	for (uint8_t i = 0; i < UKF_N; i++) {
		for (uint8_t j = 0; j < UKF_N; j++) {
			Pp[i][j] += Q.m[i][j] * qScale;
		}
	}
	Symmetrize(Pp);

	for (uint8_t i = 0; i < UKF_N; i++) {
		x[i] = mean[i];
		for (uint8_t j = 0; j < UKF_N; j++) {
			P[i][j] = Pp[i][j];
		}
	}

	status.scenario = scenarioModel.Match(x[0], x[1], x[2]);
}

/**
 * @brief The measurement is the state itself (h(x) = x), so the unscented update reduces
 *        exactly to the linear one: K = P (P + R)^-1, Joseph form covariance
 */
void Ukf::Update(const EverestMeasurement& everest)
{
	if (!initialized || !everest.valid) {
		return;
	}

	const float z[UKF_N] = {everest.altitude, everest.velocity, everest.acceleration};

	M3 R;
	MeasurementNoise(everest, R);

	for (uint8_t i = 0; i < UKF_N; i++) {
		status.innovation[i] = z[i] - x[i];
	}

	// S = P + R, and the normalised innovation squared (consistency check that needs no truth)
	M3 S, Sinv;
	auto innovationCovariance = [&]() -> bool {
		for (uint8_t i = 0; i < UKF_N; i++) {
			for (uint8_t j = 0; j < UKF_N; j++) {
				S[i][j] = P[i][j] + R[i][j];
			}
		}
		return Invert(S, Sinv);
	};
	auto nis = [&]() {
		float n = 0.0f;
		for (uint8_t i = 0; i < UKF_N; i++) {
			for (uint8_t j = 0; j < UKF_N; j++) {
				n += status.innovation[i] * Sinv[i][j] * status.innovation[j];
			}
		}
		return n;
	};

	if (!innovationCovariance()) {
		return;
	}
	status.nis = nis();

	// The model does not fit this flight right now: widen P so the measurement wins
	status.inflation = 1.0f;
#ifdef APF_PC_REPLAY
	const bool inflateOn = inflate;
#else
	const bool inflateOn = true;
#endif
	if (inflateOn && status.nis > UKF_NIS_GATE) {
		status.inflation = status.nis / static_cast<float>(UKF_N);
		if (status.inflation > UKF_MAX_INFLATION) {
			status.inflation = UKF_MAX_INFLATION;
		}
		for (uint8_t i = 0; i < UKF_N; i++) {
			for (uint8_t j = 0; j < UKF_N; j++) {
				P[i][j] *= status.inflation;
			}
		}
		status.inflations++;
		if (!innovationCovariance()) {
			return;
		}
	}

	M3 K;
	Multiply(P, Sinv, K);
	for (uint8_t i = 0; i < UKF_N; i++) {
		status.gain[i] = K[i][i];
	}

	// State
	for (uint8_t i = 0; i < UKF_N; i++) {
		float correction = 0.0f;
		for (uint8_t j = 0; j < UKF_N; j++) {
			correction += K[i][j] * status.innovation[j];
		}
		x[i] += correction;
	}

	// Joseph form: P = (I - K) P (I - K)^T + K R K^T
	M3 IK;
	for (uint8_t i = 0; i < UKF_N; i++) {
		for (uint8_t j = 0; j < UKF_N; j++) {
			IK[i][j] = ((i == j) ? 1.0f : 0.0f) - K[i][j];
		}
	}
	M3 IKP, IKT, A, KR, KT, B;
	Multiply(IK, P, IKP);
	for (uint8_t i = 0; i < UKF_N; i++) {
		for (uint8_t j = 0; j < UKF_N; j++) {
			IKT[i][j] = IK[j][i];
			KT[i][j] = K[j][i];
		}
	}
	Multiply(IKP, IKT, A);
	Multiply(K, R, KR);
	Multiply(KR, KT, B);
	for (uint8_t i = 0; i < UKF_N; i++) {
		for (uint8_t j = 0; j < UKF_N; j++) {
			P[i][j] = A[i][j] + B[i][j];
		}
	}
	Symmetrize(P);

	status.updates++;
}

/**
 * @brief R = FilterInit floors, altitude raised to the live dynamite sigma when larger
 */
void Ukf::MeasurementNoise(const EverestMeasurement& everest, float R[UKF_N][UKF_N]) const
{
	for (uint8_t i = 0; i < UKF_N; i++) {
		for (uint8_t j = 0; j < UKF_N; j++) {
			R[i][j] = Rfloor.m[i][j];
		}
	}
	R[0][0] = MaxF(R[0][0], everest.sigma * everest.sigma);
#ifdef APF_PC_REPLAY
	for (uint8_t i = 0; i < UKF_N; i++) {
		for (uint8_t j = 0; j < UKF_N; j++) {
			R[i][j] *= std::sqrt(rScale[i] * rScale[j]);
		}
	}
#endif
}

#ifdef APF_PC_REPLAY
void Ukf::SetNoiseTuning(const float q[UKF_N], const float r[UKF_N], bool inflateOn)
{
	for (uint8_t i = 0; i < UKF_N; i++) {
		qScale[i] = q[i];
		rScale[i] = r[i];
	}
	inflate = inflateOn;
}
#endif

/**
 * @brief P lost positive definiteness (float rounding): restart it from the measurement noise
 */
void Ukf::ResetCovariance()
{
	for (uint8_t i = 0; i < UKF_N; i++) {
		for (uint8_t j = 0; j < UKF_N; j++) {
			P[i][j] = (i == j) ? MaxF(Rfloor.m[i][i], 0.01f) : 0.0f;
		}
	}
	status.resets++;
}

void Ukf::PrintLine() const
{
	SOAR_PRINT("UKF %s | alt=%ld mm vel=%ld mm/s acc=%ld mm/s2 | sd=(%ld,%ld,%ld) | case=%u\n",
			   IsAfterApogee() ? "POST" : "PRE",
			   ToMilli(x[0]), ToMilli(x[1]), ToMilli(x[2]),
			   ToMilli(std::sqrt(P[0][0])), ToMilli(std::sqrt(P[1][1])), ToMilli(std::sqrt(P[2][2])),
			   static_cast<unsigned int>(status.scenario.caseA));
}

void Ukf::PrintStatus() const
{
	SOAR_PRINT("UKF: init=%d phase=%s updates=%lu resets=%lu model misfits (P inflated)=%lu nis=%ld/1000\n",
			   initialized ? 1 : 0, IsAfterApogee() ? "POST_APOGEE" : "PRE_APOGEE",
			   static_cast<unsigned long>(status.updates), static_cast<unsigned long>(status.resets),
			   static_cast<unsigned long>(status.inflations), ToMilli(status.nis));
	SOAR_PRINT("  x: alt=%ld mm vel=%ld mm/s acc=%ld mm/s2\n",
			   ToMilli(x[0]), ToMilli(x[1]), ToMilli(x[2]));
	SOAR_PRINT("  sd: alt=%ld mm vel=%ld mm/s acc=%ld mm/s2\n",
			   ToMilli(std::sqrt(P[0][0])), ToMilli(std::sqrt(P[1][1])), ToMilli(std::sqrt(P[2][2])));
	SOAR_PRINT("  innovation: alt=%ld mm vel=%ld mm/s acc=%ld mm/s2\n",
			   ToMilli(status.innovation[0]), ToMilli(status.innovation[1]), ToMilli(status.innovation[2]));

	const ScenarioMatch& m = status.scenario;
	SOAR_PRINT("  scenario: %s (%ld/1000) + %s, acc=%ld mm/s2 velErr=%ld mm/s%s\n",
			   SCENARIOS[m.caseA].name, ToMilli(m.weightA), SCENARIOS[m.caseB].name,
			   ToMilli(m.acceleration), ToMilli(m.velocityError),
			   (IsAfterApogee() && SCENARIOS[m.caseA].afterSynthetic) ? " [synthetic descent]" : "");
}
