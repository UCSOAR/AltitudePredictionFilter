#ifndef HALO_HPP
#define HALO_HPP

#include <AirbrakeController.hpp>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <ctime>
#include <chrono>
#include <deque>
#include "CommonDefinesFilter.hpp"
#include "AirbrakeSims.hpp"

// #define TESTING_BUILD

#ifdef TESTING_BUILD
#define SOAR_PRINT(...) printf(__VA_ARGS__)
#endif

// #define printf(...)

#ifndef TESTING_BUILD
#include "SystemDefines.hpp"
#include "UARTDriver.hpp"
#endif

// away
#ifndef HOME
#include "Eigen\Cholesky"
#include "Eigen\Dense"
#endif

using namespace Eigen;
// #define LOGON
// #define LOGMETRICS


/**
 * @brief Kinematics struct to store the kinematics of the rocket
 */
struct kinematicsHalo {
  float altitudeStore;
};

class HALO {
 public:
  void init(VectorXf &X0, MatrixXf &P0, MatrixXf Q_input, MatrixXf &R0);

  void stateUpdate();

  void prediction();

  float fAccel, fVelo, fAlt, GPS_Alt;

  float getFAlt();

  float getFVelo();

  float getFAccel();

  float getGPSAlt();

  bool gpsAvailable = 0;

  void setAlt(float gps_alt);

  VectorXf predictNextValues(std::vector<std::vector<float>> &vectors,
                             VectorXf &X_in, int scenario1Index,
                             int scenario2Index);

  VectorXf predictNextValuesOnce(
      std::vector<std::vector<float>> &vectors, VectorXf &X_in,
      int scenario1Index, int scenario2Index, int firstTimeForPoint,
      std::vector<float> &prevGain1, std::vector<float> &prevGain2,
      std::vector<std::vector<int>> &scenariosGainsList,
      int &counterSigmaPoint);

  void setStateVector(float filteredAcc, float filteredVelo, float filteredAlt,
                      float gpsAlt);

  std::pair<std::vector<int>, std::vector<std::vector<float>>>
  findNearestScenarios(std::vector<Scenario> *scenarios, VectorXf &measurement);

  // Takes Altitude, Velocity, Acceleration
  void calculateSigmaPoints();

  // loop through calculateSigmaOnce n times.
  VectorXf predictNStates(int n);

  // Acceleration, Velocity, Altitude
  VectorXf X;  // state vector

  // Altitude, Velocity, Acceleration
  VectorXf X0;  // current state vector

  MatrixXf observe(MatrixXf sigmaPoints);

  float lambda;

  float N1;

  void init(MatrixXf &X0, MatrixXf &P0, MatrixXf Q_input, VectorXf &Z_input,
            MatrixXf &F);

  VectorXf Z;  // measurement vector

  VectorXf dynamicModel(VectorXf &X);

  VectorXf dynamicModelOnce(VectorXf &X, int firstTimeForPoint,
                            std::vector<float> &prevGain1,
                            std::vector<float> &prevGain2,
                            std::vector<std::vector<int>> &scenariosGainsList,
                            int &counterSigmaPoint,
                            std::vector<Scenario> &scenarios);

  void setScenarios(std::vector<Scenario> &scenarios) {
    this->scenarios = scenarios;
  };

  std::vector<Scenario> *getScenarios() { return &this->scenarios; };

  std::vector<Scenario> scenarios;

  bool isBeforeApogee(float acceleration, float velocity, float altitude,
                      float lastAltitude);

  // baseline rate of change
  float deltaTime = 1.0 / 3.0f;

  void setDeltaTime(float deltaTime) { this->deltaTime = deltaTime; }

  float getDeltaTime() { return this->deltaTime; }

  float time = 0;

  void setTime(float time) { this->time = time; }

  float euclideanDistance(const std::vector<float> &vec1, const VectorXf &vec2);

  void createScenarios(HALO *halo);

  std::vector<std::pair<std::vector<float>, std::vector<float>>>
      listOfGainsSigmaPoints = {{{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}},
                                {{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}},
                                {{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}},
                                {{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}},
                                {{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}},
                                {{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}},
                                {{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}}};

  void overrideStateWithGPS(float GPS);

  std::vector<int> firstTime = {1, 1, 1, 1, 1, 1, 1};

  int firstTimeForPoint = 1;

  FILE *file;

  int scenarioIndex = 0;

  void initializeHALO(float initialAlt, HALO *halo);

  std::vector<float> Halo_Input(HALO *haloPointer, bool isInitialized,
                                float eAccelerationZ, float eVelocity,
                                float eAltitude, float gpsAltitude, float time,
                                float deltaTime);

  // for predictNextValues
  int counterSigmaPoint = 0;
  std::vector<std::vector<int>> scenariosGainsList = {{0, 0}, {0, 0}, {0, 0},
                                                      {0, 0}, {0, 0}, {0, 0}};

  void initializeHALOWithQR(float initialAlt, HALO *halo, MatrixXf &Q,
                            MatrixXf &R0);

 private:
  float Uaccel;
  float Ualt;
  float Uvelo;
  float gpsAlt;

  VectorXf X_in;
  VectorXf X_pred;

  float timeStep = 1.0f / 3.0f;

  std::vector<float> prevGain1 = {0.5, 0.5, 0.5};
  std::vector<float> prevGain2 = {0.5, 0.5, 0.5};

  float altitudeAccumulator = 0;
  float maxAltitude = 0;

  // ---------------------------------------------------------------------------
  // Shared confidence helpers — used by both detectors
  // ---------------------------------------------------------------------------

  float calculateConfidence(float current, float previous, float variance) {
    if (current > maxAltitude) {
      maxAltitude = current;
    }

    float difference = maxAltitude - current;

    if (difference > 0) {
      altitudeAccumulator += difference;
    }

    if (altitudeAccumulator >= variance) {
      return 1;
    }

    return std::max(difference, altitudeAccumulator) / std::abs(variance);
  }

  float calculateVelocityConfidence(float currentVelo, float varianceVelo) {
    if (currentVelo > varianceVelo) {
      return 0;
    }
    return 1 - std::abs((currentVelo) / (2 * std::abs(varianceVelo)));
  }

  float calculateAccelerationConfidence(float currentAcc, float varianceAcc) {
    float targetAcc = -9.81;
    float difference = std::abs(currentAcc - targetAcc);

    float lowerBound = currentAcc - varianceAcc;
    float upperBound = currentAcc + varianceAcc;

    if (lowerBound > targetAcc && upperBound < targetAcc) {
      return 0;
    }

    return 1 - (difference / (2 * varianceAcc));
  }

  // ---------------------------------------------------------------------------
  // WindowStats — computed ONCE per cycle via computeWindowStats(),
  // then passed into both detectors so nothing is recalculated twice.
  // ---------------------------------------------------------------------------
  struct WindowStats {
    float avgAltitude;
    float avgVelocity;
    float avgAcceleration;
    float sqrtP_altitude;
    float sqrtP_velocity;
    float sqrtP_acceleration;
    bool  windowFull;
  };

  // Call once per filter cycle. Owns the single updateBuffer() call.
  WindowStats computeWindowStats(const Measurement &m) {
    updateBuffer(m);

    WindowStats s{};
    s.windowFull = (buffer.size() >= static_cast<size_t>(windowSize));
    if (!s.windowFull) return s;

    float n           = static_cast<float>(buffer.size());
    s.avgAltitude     = altitudeSum     / n;
    s.avgVelocity     = velocitySum     / n;
    s.avgAcceleration = accelerationSum / n;

    s.sqrtP_altitude     = std::sqrt(this->P(0, 0));
    s.sqrtP_velocity     = std::sqrt(this->P(1, 1));
    s.sqrtP_acceleration = std::sqrt(this->P(2, 2));

    return s;
  }

  bool apogeeDetection(const WindowStats &s) {
    if (!s.windowFull) return false;

    float avgAltitude     = s.avgAltitude;
    float avgVelocity     = s.avgVelocity;
    float avgAcceleration = s.avgAcceleration;

    float sqrtP_altitude     = s.sqrtP_altitude;
    float sqrtP_velocity     = s.sqrtP_velocity;
    float sqrtP_acceleration = s.sqrtP_acceleration;

    if (prevAvgAltitude == 0.0) {
      prevAvgAltitude = avgAltitude;
    }

    float altitudeConfidence =
        calculateConfidence(avgAltitude, prevAvgAltitude, sqrtP_altitude);

    float velocityConfidence = 0.0;

    if (avgVelocity < sqrtP_velocity) {
      velocityConfidence =
          calculateVelocityConfidence(avgVelocity, sqrtP_velocity);

      if (std::isinf(velocityConfidence)) {
        velocityConfidence = 0.0;
      }
    }

    float accelerationConfidence =
        calculateAccelerationConfidence(avgAcceleration, sqrtP_acceleration);

    // cap confidence at 1
    if (altitudeConfidence > 1)     altitudeConfidence = 1;
    if (velocityConfidence > 1)     velocityConfidence = 1;
    if (accelerationConfidence > 1) accelerationConfidence = 1;

    // prevent misfires
    if (altitudeConfidence == 0 && velocityConfidence == 0) {
      accelerationConfidence = 0;
    }

    float totalConfidence =
        (altitudeConfidence * 0.5 + velocityConfidence * 0.7 +
         accelerationConfidence * 0.4);

    if (totalConfidence >= 1) {
      hitOne = true;
    }

    if (totalConfidence > maxAvgConfidence) {
      maxAvgConfidence = totalConfidence;
      return false;
    }

    if (hitOne) {
      if (totalConfidence < maxAvgConfidence) {
        return true;
      }
    }

    if (maxAvgConfidence > 0.7 && totalConfidence <= (maxAvgConfidence * 0.8)) {
      return true;
    }

    prevAvgAltitude     = avgAltitude;
    prevAvgVelocity     = avgVelocity;
    prevAvgAcceleration = avgAcceleration;

    return false;
  }

  // ---------------------------------------------------------------------------
  // burnoutDetection — receives pre-computed stats, no buffer/Kalman work.
  //
  // Physics:
  //   During burn  : accel large & positive (thrust >> drag)
  //   At burnout   : accel collapses toward 0 (thrust cuts off)
  //   After burnout: accel goes negative (drag + gravity)
  //   Velocity     : still positive and near peak — NOT yet decelerating
  // ---------------------------------------------------------------------------
  bool burnoutDetection(const WindowStats &s) {
    if (!s.windowFull)     return false;
    if (burnoutDetected_)  return true;   // latch

    // 1. Acceleration — primary signal, collapses toward 0 at burnout
    float accelMag = std::fabs(s.avgAcceleration);
    float accelerationConfidence = (s.sqrtP_acceleration > 0.0f)
        ? std::min(1.0f, 1.0f - accelMag / s.sqrtP_acceleration)
        : 0.0f;
    if (accelerationConfidence < 0.0f) accelerationConfidence = 0.0f;

    // 2. Velocity — still large & positive during burnout
    float velocityConfidence = (s.avgVelocity > 0.0f)
        ? std::min(1.0f, s.avgVelocity / (s.sqrtP_velocity + 1e-6f))
        : 0.0f;

    // 3. Jerk — sharpness of the accel drop distinguishes burnout from drift
    float jerkConfidence = 0.0f;
    if (prevAvgAccelBurnout_ != 0.0f && s.sqrtP_acceleration > 0.0f) {
      float drop = prevAvgAccelBurnout_ - s.avgAcceleration;
      if (drop > 0.0f)
        jerkConfidence = std::min(1.0f, drop / s.sqrtP_acceleration);
    }

    // Misfire guard — same pattern as apogeeDetection
    if (accelerationConfidence == 0.0f && velocityConfidence == 0.0f) {
      jerkConfidence = 0.0f;
    }

    float totalConfidence = accelerationConfidence * 0.6f
                          + velocityConfidence      * 0.3f
                          + jerkConfidence          * 0.3f;

    // Trigger logic — identical pattern to apogeeDetection
    if (totalConfidence >= 1.0f) burnHitPeak_ = true;

    if (totalConfidence > maxAvgConfBurnout_) {
      maxAvgConfBurnout_   = totalConfidence;
      prevAvgAccelBurnout_ = s.avgAcceleration;
      prevAvgVelBurnout_   = s.avgVelocity;
      return false;
    }

    if (burnHitPeak_ && totalConfidence < maxAvgConfBurnout_) {
      burnoutDetected_ = true;
      return true;
    }

    if (maxAvgConfBurnout_ > 0.7f &&
        totalConfidence <= (maxAvgConfBurnout_ * 0.8f)) {
      burnoutDetected_ = true;
      return true;
    }

    prevAvgAccelBurnout_ = s.avgAcceleration;
    prevAvgVelBurnout_   = s.avgVelocity;
    return false;
  }

  // from empirical observations window of 10 is best
  int windowSize = 10;

  void updateBuffer(const Measurement &currentMeasurement) {
    if (buffer.size() == static_cast<size_t>(windowSize)) {
      const Measurement &oldest = buffer.front();
      altitudeSum     -= oldest.altitude;
      velocitySum     -= oldest.velocity;
      accelerationSum -= oldest.acceleration;
      buffer.pop_front();
    }

    buffer.push_back(currentMeasurement);
    altitudeSum     += currentMeasurement.altitude;
    velocitySum     += currentMeasurement.velocity;
    accelerationSum += currentMeasurement.acceleration;
  }

  std::deque<Measurement> buffer;
  float altitudeSum     = 0.0f;
  float velocitySum     = 0.0f;
  float accelerationSum = 0.0f;

  // Apogee detection state
  float prevAvgAltitude     = 0.0f;
  float prevAvgVelocity     = 0.0f;
  float prevAvgAcceleration = 0.0f;
  float maxAvgConfidence    = 0.0f;
  bool  hitOne              = false;
  bool  wait                = false;

  // Burnout detection state
  bool  burnoutDetected_   = false;
  bool  burnHitPeak_       = false;
  float maxAvgConfBurnout_ = 0.0f;
  float prevAvgAccelBurnout_ = 0.0f;
  float prevAvgVelBurnout_   = 0.0f;

  // Cached nearest-scenario vectors from last dynamicModel() call.
  // Populated in dynamicModel() so stateUpdate() can pass them to
  // AirbrakeController without a second findNearestScenarios lookup.
  std::vector<std::vector<float>> lastNearestVectors_;

  // Airbrake state
  AirbrakeController airbrakeController_;
  int airbrakeLevel_ = 0;

 protected:
  MatrixXf sigmaPoints;
  MatrixXf Xprediction;
  MatrixXf Pprediction;
  MatrixXf P;
  MatrixXf Q;
  MatrixXf projectError;

  MatrixXf WeightsUKF;
  VectorXf WeightsForSigmaPoints;

  MatrixXf F;
  MatrixXf H;
  MatrixXf R;
  MatrixXf K;

  kinematicsHalo KinematicsHalo;

  MatrixXf sigPoints;
};

#include "AirbrakeController.hpp"
// This is placed at the BOTTOM of HALO.hpp so AirbrakeController.hpp
// sees Scenario, VectorXf, and all HALO types already defined.
// AirbrakeController.hpp must NOT include HALO.hpp itself (already
// guarded by the comment in that file).

#endif
