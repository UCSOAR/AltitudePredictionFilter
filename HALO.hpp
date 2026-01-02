#ifndef HALO_HPP
#define HALO_HPP

#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <ctime>
#include <chrono>
#include "KDTree.hpp"
#include <deque>

#ifdef HOME
#include "C:\Users\harry\Desktop\soar\eigen-3.4.0\eigen-3.4.0\Eigen\Cholesky"
#include "C:\Users\harry\Desktop\soar\eigen-3.4.0\eigen-3.4.0\Eigen\Dense"
#endif

// away
#ifndef HOME
#include "C:\Users\harry\Desktop\soar\eigen-3.4.0\eigen-3.4.0\Eigen\Cholesky"
#include "C:\Users\harry\Desktop\soar\eigen-3.4.0\eigen-3.4.0\Eigen\Dense"
#endif

using namespace Eigen;
// #define LOGON

/**
 * @brief Measurement struct to store the time, altitude, velocity and
 * acceleration
 */
struct Measurement {
  float altitude;
  float velocity;
  float acceleration;
  float time;
};

/**
 * @brief Scenario struct to store the coefficients of the 3rd degree polynomial
 * for acceleration, velocity and altitude before and after apogee, also
 * evaluates the acceleration, velocity and altitude at a given time
 */
struct Scenario {
  std::vector<float> beforeApogeeAccel;
  std::vector<float> afterApogeeAccel;

  std::vector<float> beforeApogeeVelo;
  std::vector<float> afterApogeeVelo;

  std::vector<float> beforeApogeeAlt;
  std::vector<float> afterApogeeAlt;

  std::vector<std::vector<float>> BeforeList;
  std::vector<std::vector<float>> AfterList;

  KDTree treeBefore;
  KDTree treeAfter;

  int name;

  std::vector<float> measurement;
  bool isBeforeApogeeBool = true;

  
  // conversions between vector and array.
  inline std::array<float, 3> vec2arr(const std::vector<float>& v) {
      assert(v.size() == 3);
      return {v[0], v[1], v[2]};
  }

  inline std::vector<float> arr2vec(const std::array<float, 3>& a) {
      return {a[0], a[1], a[2]};
  }


  Scenario(std::vector<std::vector<float>> beforeList,
           std::vector<std::vector<float>> afterList, int Name)
      : BeforeList(beforeList), AfterList(afterList), name(Name) {}

  // Function to find the vector and split the list
  std::pair<std::vector<std::vector<float>>, std::vector<std::vector<float>>>
  findAndSplitVector(const std::vector<std::vector<float>> &inputList) {
    std::vector<std::vector<float>> firstPart;
    std::vector<std::vector<float>> secondPart;
    bool splitPointFound = false;

    for (const auto &vec : inputList) {
      if (!splitPointFound && vec.size() >= 3 && vec[1] < 1.0f &&
          vec[2] < 1.0f) {
        splitPointFound = true;
      }

      if (splitPointFound) {
        secondPart.push_back(vec);
      } else {
        firstPart.push_back(vec);
      }
    }

    return {firstPart, secondPart};
  }

  /**
   * {Altitude, Velocity, Acceleration}
   */
  void setMeasurement(std::vector<float> measurementVector) {
    this->measurement = measurementVector;
  }

  void setIsBeforeApogee(bool isBeforeApogee) {
    isBeforeApogeeBool = isBeforeApogee;
  }

  /**
   * Returns the nearest vector to the measurement vector
   */
  std::pair<std::vector<float>, size_t> nearestKDTree(std::vector<float> measurement) {
    std::array<float, 3> convertedMeasurement = vec2arr(measurement);
    pointIndex result;
    if (isBeforeApogeeBool) {
      result = treeBefore.nearest_pointIndex(convertedMeasurement);
    } else {
      result = treeAfter.nearest_pointIndex(convertedMeasurement);
    }

    return {arr2vec(result.first), result.second};
  }

  void createTree() {
    std::vector<std::array<float, 3>> beforeArrayOfArrays;
    std::vector<std::array<float, 3>> afterArrayOfArrays;

    for (const auto& pt : BeforeList) {
        // Assuming BeforeList[i] has at least 3 elements
        beforeArrayOfArrays.push_back({pt[0], pt[1], pt[2]});
    }

    for (const auto& pt : AfterList) {
        afterArrayOfArrays.push_back({pt[0], pt[1], pt[2]});
    }

    treeBefore = KDTree(beforeArrayOfArrays);
    treeAfter = KDTree(afterArrayOfArrays);
  }

  /**
   * Returns list of vectors of scenario {Altitude, Velocity, Acceleration}
   * before or after apogee pass index instead
   */
  std::vector<std::vector<float>> *getLists() {
    if (isBeforeApogeeBool) {
      return &BeforeList;
    } else {
      return &AfterList;
    }
  }

  // Binary search function to find the index of the closest time
  int binarySearch(const std::vector<std::vector<float>> &list, float time) {
    int left = 0;
    int right = list.size() - 1;

    while (left <= right) {
      int mid = left + (right - left) / 2;

      if (list[mid][3] == time) {
        return mid;
      } else if (list[mid][3] < time) {
        left = mid + 1;
      } else {
        right = mid - 1;
      }
    }

    // If the exact time is not found, return the closest index
    return (left < list.size()) ? left : right;
  }

  /** finds vector at specified index **/
  std::vector<float> evaluateVectorAt(int index) {
    return (*getLists())[index];
  }

  /** finds vector at specified time **/
  // binary search
  std::vector<float> evaluateVectorAtTime(float time) {
    std::vector<std::vector<float>> *list = getLists();
    int index = 0;
    std::vector<float> vect = {0, 0, 0, 0};

    vect = list->at(binarySearch((*list), time));

    return vect;
  }
};

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
                            int &counterSigmaPoint, std::vector<Scenario>& scenarios);

  void setScenarios(std::vector<Scenario> &scenarios) {
    this->scenarios = scenarios;
  };

  std::vector<Scenario> *getScenarios() { return &this->scenarios; };

  std::vector<Scenario> scenarios;

  bool isBeforeApogee(float acceleration, float velocity, float altitude,
                      float lastAltitude);

  float deltaTime = 1.0 / 3;

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
                                float eAltitude, float gpsAltitude, float time);

  // for predictNextValues
  int counterSigmaPoint = 0;
  std::vector<std::vector<int>> scenariosGainsList = {{0, 0}, {0, 0}, {0, 0},
                                                      {0, 0}, {0, 0}, {0, 0}};

  std::chrono::duration<float> updateTime;
  std::chrono::duration<float> predictTime;
  std::chrono::duration<float> triangulationTime;
  std::chrono::duration<float> dynamicModelTime;
  std::chrono::duration<float> nearestScenariosTime;
  std::chrono::duration<float> KDTreeTime;
  std::chrono::duration<float> euclideanTime;
  std::chrono::duration<float> twoDistancesTime;
  std::chrono::duration<float> vectorsTime;
  std::chrono::duration<float> push_backTime;
  std::chrono::duration<float> loopScenariosTime;
  std::chrono::duration<float> emplaceBackTime;
  std::chrono::duration<float> getListsTime;
  std::chrono::duration<float> othersTime;
  std::chrono::duration<float> PpredictionTime;
  std::chrono::duration<float> projErrorTime;
  std::chrono::duration<float> sPointTime;
  std::chrono::duration<float> preMeanTime;
  std::chrono::duration<float> predictLoopTime;
  std::chrono::duration<float> endPredictLoopTime;
  std::chrono::duration<float> getScenarioTime;
  std::chrono::duration<float> treeCreationTime;

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

  // last time that a prediction was made. done to prevent multiple triggers.
  int lastTriggerTime = -1;

  std::vector<float> prevGain1 = {0.5, 0.5, 0.5};
  std::vector<float> prevGain2 = {0.5, 0.5, 0.5};

  float altitudeAccumulator = 0;

  float maxAltitude = 0;

  float calculateConfidence(float current, float previous, float variance) {
    if (current > maxAltitude) {
      maxAltitude = current;
    }

    float difference = maxAltitude - current;

    if (difference > 0) {
      altitudeAccumulator += difference;
    }

    // override since altitude has consistently been going down in the range of
    // the variance
    if (altitudeAccumulator >= variance) {
      return 1;
    }

    return std::max(difference, altitudeAccumulator) / std::abs(variance);
  }

  float calculateVelocityConfidence(float currentVelo, float varianceVelo) {
    // velocity should be < 0 for apogee
    // range
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

    float confidence = 1 - (difference / (2 * varianceAcc));

    return confidence;
  }

  bool apogeeDetection(const Measurement &currentMeasurement) {
    updateBuffer(currentMeasurement);

    float avgAltitude = altitudeSum / buffer.size();
    float avgVelocity = velocitySum / buffer.size();
    float avgAcceleration = accelerationSum / buffer.size();

    // Square root the P values to get std deviation
    float sqrtP_altitude = std::sqrt(this->P(0, 0));
    float sqrtP_velocity = std::sqrt(this->P(1, 1));
    float sqrtP_acceleration = std::sqrt(this->P(2, 2));

    if (buffer.size() < windowSize) {
      return false;
    }

    if (prevAvgAltitude == 0.0) {
      prevAvgAltitude = avgAltitude;
    }

    float altitudeConfidence =
        calculateConfidence(avgAltitude, prevAvgAltitude, sqrtP_altitude);

    float velocityConfidence = 0.0;

    if (avgVelocity < sqrtP_velocity) {
      velocityConfidence =
          calculateVelocityConfidence(avgVelocity, sqrtP_velocity);

      // if infinity, set to 1
      if (std::isinf(velocityConfidence)) {
        velocityConfidence = 0.0;
      }
    }
    float accelerationConfidence =
        calculateAccelerationConfidence(avgAcceleration, sqrtP_acceleration);

    // cap confidence at 1
    if (altitudeConfidence > 1) {
      altitudeConfidence = 1;
    }

    if (velocityConfidence > 1) {
      velocityConfidence = 1;
    }

    if (accelerationConfidence > 1) {
      accelerationConfidence = 1;
    }

    // prevent misfires
    if (altitudeConfidence == 0 && velocityConfidence == 0) {
      accelerationConfidence = 0;
    }

    float totalConfidence =
        (altitudeConfidence * 0.5 + velocityConfidence * 0.7 +
         accelerationConfidence * 0.4);

#ifdef LOGON
    // write to file confidence values
    FILE *file = fopen("testSuite/results/confidence.txt",
                       "a+");  // Open the file for writing
    if (!file) {
      fprintf(stderr, "Error opening confidence.txt...exiting\n");
      exit(1);
    }

    fprintf(file, "%f,%f,%f,%f\n", altitudeConfidence, velocityConfidence,
            accelerationConfidence, totalConfidence);

    fclose(file);

#endif

    // if (totalConfidence >= 1) {
    //   return true;
    // }

    if (totalConfidence >= 1) {
      hitOne = true;
    }

    // Update the maximum average confidence
    if (totalConfidence > maxAvgConfidence) {
      maxAvgConfidence = totalConfidence;
      return false;
    }

    // if we've hit one and its decreasing then trigger
    // if confidence drops below 1 or is 1 trigger
    if (hitOne) {
      if (totalConfidence < maxAvgConfidence) {
        return true;
      }
    }

    // if we haven't hit one and it dips then trigger
    if (maxAvgConfidence > 0.7 && totalConfidence <= (maxAvgConfidence * 0.8)) {
      return true;
    }

    prevAvgAltitude = avgAltitude;
    prevAvgVelocity = avgVelocity;
    prevAvgAcceleration = avgAcceleration;

    return false;
  }

  // from empirical observations window of 10 is best
  // detected highest apogee with relatively
  // best confidence values
  int windowSize = 10;

 private:
  void updateBuffer(const Measurement &currentMeasurement) {
    if (buffer.size() == windowSize) {
      const Measurement &oldest = buffer.front();
      altitudeSum -= oldest.altitude;
      velocitySum -= oldest.velocity;
      accelerationSum -= oldest.acceleration;
      buffer.pop_front();
    }

    buffer.push_back(currentMeasurement);
    altitudeSum += currentMeasurement.altitude;
    velocitySum += currentMeasurement.velocity;
    accelerationSum += currentMeasurement.acceleration;
  }

  std::deque<Measurement> buffer;
  float altitudeSum;
  float velocitySum;
  float accelerationSum;
  float prevAvgAltitude = 0.0;
  float prevAvgVelocity = 0.0;
  float prevAvgAcceleration = 0.0;
  float maxAvgConfidence = 0.0;
  bool hitOne = false;
  bool wait = false;

 protected:
  MatrixXf sigmaPoints;
  MatrixXf Xprediction;
  MatrixXf Pprediction;
  MatrixXf P;
  MatrixXf Q;
  MatrixXf projectError;

  MatrixXf WeightsUKF;
  VectorXf WeightsForSigmaPoints;

  MatrixXf F;  // state to next state transition matrix
  MatrixXf H;  // state to measurement matrix
  MatrixXf R;  // measurement noise covariance matrix
  MatrixXf K;  // Kalman gain matrix

  kinematicsHalo KinematicsHalo;

  MatrixXf sigPoints;
};

#endif
