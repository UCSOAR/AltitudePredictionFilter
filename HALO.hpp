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
#include "C:\Users\andin\OneDrive\Documents\AllRepos\UnscentedKalmanFilter\eigen-3.4.0\Eigen\Cholesky"
#include "C:\Users\andin\OneDrive\Documents\AllRepos\UnscentedKalmanFilter\eigen-3.4.0\Eigen\Dense"
#endif

// away
#ifndef HOME
#include "C:\Users\Andrey\Documents\UKFRepo\UnscentedKalmanFilter\eigen-3.4.0\eigen-3.4.0\Eigen\Cholesky"
#include "C:\Users\Andrey\Documents\UKFRepo\UnscentedKalmanFilter\eigen-3.4.0\eigen-3.4.0\Eigen\Dense"
#endif

using namespace Eigen;
#define LOGON

/**
 * @brief Measurement struct to store the time, altitude, velocity and
 * acceleration
 */
struct Measurement {
  double altitude;
  double velocity;
  double acceleration;
  double time;
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
  std::pair<std::vector<float>, size_t> nearestKDTree(
      std::vector<float> measurement) {
    if (isBeforeApogeeBool) {
      return treeBefore.nearest_pointIndex(measurement);
    } else {
      return treeAfter.nearest_pointIndex(measurement);
    }
  }

  void createTree() {
    std::vector<std::vector<float>> beforeVectorofVectors;
    std::vector<std::vector<float>> afterVectorofVectors;

    for (int i = 0; i < BeforeList.size(); i++) {
      std::vector<float> vect = {BeforeList[i][0], BeforeList[i][1],
                                 BeforeList[i][2]};
      beforeVectorofVectors.push_back(vect);
    }

    for (int i = 0; i < AfterList.size(); i++) {
      std::vector<float> vect = {AfterList[i][0], AfterList[i][1],
                                 AfterList[i][2]};
      afterVectorofVectors.push_back(vect);
    }

    treeBefore = KDTree(beforeVectorofVectors);
    treeAfter = KDTree(afterVectorofVectors);
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

  void setStateVector(float filteredAcc, float filteredVelo, float filteredAlt);

  std::pair<std::vector<int>, std::vector<std::vector<float>>>
  findNearestScenarios(std::vector<Scenario> *scenarios, VectorXf &measurement);

  // Takes Altitude, Velocity, Acceleration
  void calculateSigmaPoints();

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

  void setScenarios(std::vector<Scenario> &scenarios) {
    this->scenarios = scenarios;
  };

  std::vector<Scenario> *getScenarios() { return &this->scenarios; };

  std::vector<Scenario> scenarios;

  bool isBeforeApogee(float acceleration, float velocity, float altitude,
                      float lastAltitude);

  float deltaTime = 1.0 / 3;

  void setDeltaTime(float deltaTime) { this->deltaTime; }

  float getDeltaTime() { return this->deltaTime; }

  float time = 0;

  void setTime(float time) { this->time; }

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

  std::vector<double> Halo_Input(HALO *haloPointer, bool isInitialized,
                                 double eAccelerationZ, double eVelocity,
                                 double eAltitude, float time);

  // for predictNextValues
  int counterSigmaPoint = 0;
  std::vector<std::vector<int>> scenariosGainsList = {{0, 0}, {0, 0}, {0, 0},
                                                      {0, 0}, {0, 0}, {0, 0}};

  std::chrono::duration<double> updateTime;
  std::chrono::duration<double> predictTime;
  std::chrono::duration<double> triangulationTime;
  std::chrono::duration<double> dynamicModelTime;
  std::chrono::duration<double> nearestScenariosTime;
  std::chrono::duration<double> KDTreeTime;
  std::chrono::duration<double> euclideanTime;
  std::chrono::duration<double> twoDistancesTime;
  std::chrono::duration<double> vectorsTime;
  std::chrono::duration<double> push_backTime;
  std::chrono::duration<double> loopScenariosTime;
  std::chrono::duration<double> emplaceBackTime;
  std::chrono::duration<double> getListsTime;
  std::chrono::duration<double> othersTime;
  std::chrono::duration<double> PpredictionTime;
  std::chrono::duration<double> projErrorTime;
  std::chrono::duration<double> sPointTime;
  std::chrono::duration<double> preMeanTime;
  std::chrono::duration<double> predictLoopTime;
  std::chrono::duration<double> endPredictLoopTime;
  std::chrono::duration<double> getScenarioTime;
  std::chrono::duration<double> treeCreationTime;

 private:
  float Uaccel;
  float Ualt;
  float Uvelo;

  VectorXf X_in;
  VectorXf X_pred;

  float timeStep = 1 / 3;

  std::vector<float> prevGain1 = {0.5, 0.5, 0.5};
  std::vector<float> prevGain2 = {0.5, 0.5, 0.5};

  double altitudeAccumulator = 0;

  double maxAltitude = 0;

  double calculateConfidence(double current, double previous, double variance) {
    if (current > maxAltitude) {
      maxAltitude = current;
    }

    double difference = maxAltitude - current;

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

  double calculateVelocityConfidence(double currentVelo, double varianceVelo) {
    // velocity should be < 0 for apogee
    // range
    double range = 2 * std::abs(varianceVelo);

    return (range - (std::abs(varianceVelo) + currentVelo)) / range;
  }

  double calculateAccelerationConfidence(double currentAcc,
                                         double varianceAcc) {
    double targetAcc = -9.81;
    double difference = std::abs(currentAcc - targetAcc);

    double lowerBound = currentAcc - varianceAcc;
    double upperBound = currentAcc + varianceAcc;

    if (lowerBound > targetAcc && upperBound < targetAcc) {
      return 0;
    }

    double confidence = (difference / (2 * varianceAcc));

    return confidence;
  }

  bool apogeeDetection(const Measurement &currentMeasurement) {
    updateBuffer(currentMeasurement);

    double avgAltitude = altitudeSum / buffer.size();
    double avgVelocity = velocitySum / buffer.size();
    double avgAcceleration = accelerationSum / buffer.size();

    // Square root the P values
    double sqrtP_altitude = std::sqrt(this->P(0, 0));
    double sqrtP_velocity = std::sqrt(this->P(1, 1));
    double sqrtP_acceleration = std::sqrt(this->P(2, 2));

    if (buffer.size() < windowSize) {
      return false;
    }

    if (prevAvgAltitude == 0.0) {
      prevAvgAltitude = avgAltitude;
    }

    double altitudeConfidence =
        calculateConfidence(avgAltitude, prevAvgAltitude, sqrtP_altitude);
    double velocityConfidence = 0.0;

    if (avgVelocity < sqrtP_velocity) {
      velocityConfidence =
          calculateVelocityConfidence(avgVelocity, sqrtP_velocity);

      // if infinity, set to 1
      if (std::isinf(velocityConfidence)) {
        velocityConfidence = 0.0;
      }
    }
    double accelerationConfidence =
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

    // check confidence
    std::cout << "Confidence_alt: " << altitudeConfidence
              << " velo: " << velocityConfidence
              << " acc: " << accelerationConfidence << std::endl;

    double totalConfidence =
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

    std::cout << "Altitude: " << avgAltitude << " Velocity: " << avgVelocity
              << " Acceleration: " << avgAcceleration << std::endl;
    std::cout << "Confidence: " << totalConfidence << std::endl;

    if (totalConfidence >= 1) {
      std::cout << "Apogee detected at time: " << currentMeasurement.time
                << " seconds, Altitude: " << avgAltitude << " meters"
                << std::endl;
      // stop program
      // exit(0);
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
  double altitudeSum;
  double velocitySum;
  double accelerationSum;
  double prevAvgAltitude = 0.0;
  double prevAvgVelocity = 0.0;
  double prevAvgAcceleration = 0.0;

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
