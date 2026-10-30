/*
 * FilterProfile.hpp
 *
 *  Created on: Oct 24, 2026
 *      Author: Aly Masani
 */

#ifndef FILTERPROFILE_FILTERPROFILE_HPP_
#define FILTERPROFILE_FILTERPROFILE_HPP_

/************************************
 * MACROS AND DEFINES
 ************************************/
// One switch picks the scenario set AND the covariances the filter is built with.
// Change APF_PROFILE below (or pass -DAPF_PROFILE=n to the compiler, as Tools/taber does)
//
//   APF_PROFILE_STORMFLY      9 Stormfly MAPLEAF cases (Estimate/ScenarioData.cpp), our Q / R (Input/FilterInit.cpp)
//   APF_PROFILE_TABER_LEGACY  legacy HALO sims, Legacy/Data.cpp (2 Taber sims), legacy Q / R
//   APF_PROFILE_TABER_SIX     testSuite/data/beforeSimsF2_Short.csv (6 sims the legacy createScenarios()
//                             TODO meant to load: 5 at a 1404 m pad + 1 Taber), legacy Q / R
//   APF_PROFILE_TABER_ALL     every distinct legacy sim together (the 2 above + the 6 above = 8; the 4
//                             covarianceCalc sims are F2 sims 2-5 again), legacy Q / R
//   APF_PROFILE_TABER_OWN     only the sims made for Taber (983 m pad): Legacy/Data.cpp 1-2 + F2 sim 6,
//                             legacy Q / R. F2 sims 1-5 are the 1404 m pad Quasar (SAC) sims
//   (the Taber sets live in Estimate/ScenarioDataTaber.cpp)
#define APF_PROFILE_STORMFLY		0
#define APF_PROFILE_TABER_LEGACY	1
#define APF_PROFILE_TABER_SIX		2
#define APF_PROFILE_TABER_ALL		3
#define APF_PROFILE_TABER_OWN		4

#ifndef APF_PROFILE
#define APF_PROFILE APF_PROFILE_STORMFLY
#endif

// Covariances follow the profile unless overridden (-DAPF_COVARIANCE=n), so the scenario
// set and the Q / R can be tested separately
#define APF_COVARIANCE_DEFAULT		0	// Input/FilterInit.cpp DEFAULT_Q / DEFAULT_R
#define APF_COVARIANCE_TABER		1	// Legacy HALO::initializeHALO (covarianceCalc.py)
#define APF_COVARIANCE_SCENARIO		2	// DEFAULT_R, Q before / after apogee from the case spread (ScenarioData.hpp)

// The legacy Taber Q / R (covarianceCalc.py) is the default for every profile, the board included.
// APF_COVARIANCE_DEFAULT (Input/FilterInit.cpp DEFAULT_Q / DEFAULT_R) stays selectable with -DAPF_COVARIANCE=0
#ifndef APF_COVARIANCE
#define APF_COVARIANCE APF_COVARIANCE_TABER
#endif

// Q scale (tuning sweeps, Tools/taber/taber_test.py --q-sweep): Q is multiplied by this.
// 1 = the matrices in Input/FilterInit.cpp as they are
#ifndef APF_Q_SCALE
#define APF_Q_SCALE 1.0f
#endif
#ifndef APF_R_SCALE
#define APF_R_SCALE 1.0f		// Same for R
#endif

// UKF process model, also follows the profile unless overridden (-DAPF_SCENARIO_MODEL=n)
//   ACCEL   ours: integrate the blended acceleration of the two closest scenarios over dt
//   LEGACY  legacy HALO::predictNextValues: the predicted state IS the blend of the two closest
//           scenarios' states dt later, gains per axis from the previous step, reset to 50/50
//           when the pair changes (ScenarioModel::LegacyPredict)
#define APF_SCENARIO_MODEL_ACCEL	0
#define APF_SCENARIO_MODEL_LEGACY	1

#ifndef APF_SCENARIO_MODEL
#if APF_PROFILE == APF_PROFILE_STORMFLY
#define APF_SCENARIO_MODEL APF_SCENARIO_MODEL_ACCEL
#else
#define APF_SCENARIO_MODEL APF_SCENARIO_MODEL_LEGACY
#endif
#endif

#endif /* FILTERPROFILE_FILTERPROFILE_HPP_ */
