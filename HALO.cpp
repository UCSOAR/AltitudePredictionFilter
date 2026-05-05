#include <StandardSims.hpp>
#include "HALO.hpp"
#include <fstream>
#include <map>

// #define LOGON
// #define LOGMETRICS
// #define TIMERON
// #define TESTING_BUILD

#ifndef TESTING_BUILD
#include "SystemDefines.hpp"
#include "UARTDriver.hpp"
#endif

// home
#ifdef HOME
#include "Eigen\Dense"
#endif

// away
#ifndef HOME
#include "Eigen\Dense"
#endif

#ifndef HALO_CPP
#define HALO_CPP
#define REFRESH_RATE 3

#ifdef TESTING_BUILD
#define SOAR_PRINT(...) printf(__VA_ARGS__)
#endif

#define printf(...) ;
// #define SOAR_PRINT(...) ;

/* Constants for the UKF... do we ever use it?
#define N 6
#define dim 6
#define alpha 0.1
#define beta 2
#define k 3 - dim  // 3 dimensions
*/

#define OBSERVATION_DIMENSIONS 4

std::string directoryPath = "testSuite/results";

bool isInitialized = false;
int counter = 0;
int isAfterApogee = 0;
int predictionCounter = 0;

AirbrakeController airbrakeController_;
bool isAfterBurnout = false;
int  airbrakeLevel_ = 0;

using namespace Eigen;

#ifdef LOGON
FILE* resetGainsFile;
#endif

#if defined(LOGON) || defined(LOGMETRICS)
static FILE* predictnalt = NULL;
static FILE* kalmangains = NULL;

// how many predictions have been made.
#endif

void HALO::init(VectorXf& X0, MatrixXf& P0, MatrixXf Q_input, MatrixXf& R0) {
#ifdef LOGON

  // remove the file if it exists
  std::string filePath = directoryPath + "/resetGainsFile.txt";
  if (std::remove(filePath.c_str()) == 0) {
  } else {
    std::perror("Error deleting file");
  }

  resetGainsFile = fopen((directoryPath + "/resetGainsFile.txt").c_str(),
                         "a+");  // Open the file for writing
  if (!resetGainsFile) {
    fprintf(stderr, "Error opening resetGainsFile.txt...exiting\n");
    exit(1);
  }
#endif

#if defined(LOGON) || defined(LOGMETRICS)
  // deleting P.txt
  std::string filePath2 = directoryPath + "/P.txt";
  if (std::remove(filePath2.c_str()) == 0) {
  } else {
    std::perror("Error deleting file");
  }

  // create P.txt
  FILE* file = fopen((directoryPath + "/P.txt").c_str(),
                     "a+");  // Open the file for writing
  if (!file) {
    fprintf(stderr, "Error opening P.txt...exiting\n");
    exit(1);
  }

  fprintf(file, "time,P00,P01,P02,P10,P11,P12,P20,P21,P22\n");
  fclose(file);
#endif

  // Initial Guess
  this->X0 = X0;
  this->P = P0;
  this->Q = Q_input;
  this->R = R0;

  // Weights for sigma points
  float dime = 3;
  float k1 = 3 - dime;
  float w0_m = k1 / (dime + k1);  // weight for first sPoint when cal covar
  float w_i = 1 / (2 * (dime + k1));

  MatrixXf Weights(7, 7);
  VectorXf W(7, 1);

  W.setConstant(7, w_i);

  for (int i = 1; i < 7; i++) {
    Weights.diagonal()[i] = w_i;
  }

  Weights(0) = w0_m;

  this->WeightsUKF = Weights;

  VectorXf WeightsForSigmaPoints(7, 1);
  WeightsForSigmaPoints.setConstant(7, w_i);
  WeightsForSigmaPoints(0) = w0_m;
  this->WeightsForSigmaPoints = WeightsForSigmaPoints;

  this->KinematicsHalo.altitudeStore = X0(0);
  this->N1 = 3;

  calculateSigmaPoints();
}

// Update Step-------------------------------------
void HALO::stateUpdate() {
  std::chrono::high_resolution_clock::time_point stateUpdateTime;

#ifdef TIMERON
  stateUpdateTime = std::chrono::high_resolution_clock::now();
#endif

  // Observation dimensions is the dimension of our observation vector. At the
  // moment, it is the 3 values from Everest and 1 GPS altitude.
  MatrixXf observedValues(OBSERVATION_DIMENSIONS, 7);
  observedValues.setZero(OBSERVATION_DIMENSIONS, 7);

  for (int i = 0; i < 7; i++) {
    observedValues(0, i) = sigPoints(0, i);
    observedValues(1, i) = sigPoints(1, i);
    observedValues(2, i) = sigPoints(2, i);
    observedValues(3, i) = sigPoints(0, i);
  }

  // calculate the mean of the observed values
  VectorXf zMean(OBSERVATION_DIMENSIONS);
  zMean.setZero(OBSERVATION_DIMENSIONS);
  zMean = observedValues * WeightsForSigmaPoints;
  this->Z = zMean;

  // calculate covariance of Z, find matrix of deviations of observations.
  MatrixXf zCovar(OBSERVATION_DIMENSIONS, 7);
  zCovar.setZero(OBSERVATION_DIMENSIONS, 7);

  for (int i = 0; i < OBSERVATION_DIMENSIONS; i++) {
    zCovar.row(i) =
        (observedValues.row(i).array() - zMean.row(i).value()).matrix();
  }

  // create a new measurement noise matrix that is 4x4, with GPS measurement
  // noise added. GPS is an independent altitude measurement, so it has 0 for
  // its cross values.
  MatrixXf R_GPS(OBSERVATION_DIMENSIONS, OBSERVATION_DIMENSIONS);
  R_GPS.setZero();
  R_GPS.block<3, 3>(0, 0) = this->R;  // original 3x3 R matrix
  if (!gpsAvailable) {
    R_GPS(3, 3) = 1e9;  // high noise when GPS is off or dead.
  } else {
    R_GPS(3, 3) =
        25;  // dummy value for 5m std dev. low R for GPS means higher trust!
  }
  // calculate the innovation covariance, measurement covariance
  MatrixXf Pz(OBSERVATION_DIMENSIONS, OBSERVATION_DIMENSIONS);
  Pz.setZero(OBSERVATION_DIMENSIONS, OBSERVATION_DIMENSIONS);
  Pz = (zCovar * WeightsForSigmaPoints.asDiagonal() * zCovar.transpose()) +
       R_GPS;

  // calculate the cross covariance
  MatrixXf Pxz(3, OBSERVATION_DIMENSIONS);
  Pxz.setZero();

  Pxz = projectError * WeightsForSigmaPoints.asDiagonal() * zCovar.transpose();

  // calculate the Kalman gain
  MatrixXf K(3, OBSERVATION_DIMENSIONS);
  K.setZero();
  K = Pxz * Pz.inverse();

  SOAR_PRINT("  Pxz:\n");
  for (int i = 0; i < Pxz.rows(); i++) {
    for (int j = 0; j < Pxz.cols(); j++) {
      SOAR_PRINT("Pxz(%d,%d) = %f  ", i, j, K(i, j));
    }
    SOAR_PRINT("\n");
  }

  SOAR_PRINT("Kalman Gain K:\n");
  for (int i = 0; i < K.rows(); i++) {
    for (int j = 0; j < K.cols(); j++) {
      SOAR_PRINT("K(%d,%d) = %f  ", i, j, K(i, j));
    }
    SOAR_PRINT("\n");
  }

#if defined(LOGON) || defined(LOGMETRICS)

  if (kalmangains == NULL) {
    kalmangains = fopen((directoryPath + "/kalmangains.txt").c_str(), "a+");
  }

  if (!kalmangains) {
    perror("Error opening kalmangains.txt");
    fprintf(stderr, "Errno: %d\n", errno);
    exit(1);
  }

  for (int rows = 0; rows < 3; rows++) {
    fprintf(kalmangains, "%f,%f,%f,%f,%f\n", time, K(rows, 0), K(rows, 1),
            K(rows, 2), K(rows, 3));
  }
#endif

  bool kZero = false;

  for (int row = 0; row < 3; row++) {
    for (int col = 0; col < OBSERVATION_DIMENSIONS; col++) {
      if (std::isnan(K(row, col))) {
        kZero = true;
        K(row, col) = 0;

#ifdef TESTING_BUILD
        FILE* log = fopen("log.txt", "a+");  // Open the file for appending or
                                             // create it if it doesn't exist

        if (!log) {
          fprintf(stderr, "Error opening log.txt...exiting\n");
          exit(1);
        }

        fprintf(log, "At time %f s NAN detected in the Kalman Gain, ",
                this->time);

        fclose(log);
#endif
      }
    }
  }

#ifdef TESTING_BUILD
  if (kZero) {
    FILE* log = fopen(
        "log.txt",
        "a+");  // Open the file for appending or create it if it doesn't exist

    if (!log) {
      fprintf(stderr, "Error opening log.txt...exiting\n");
      exit(1);
    }

    fprintf(log, "\n");
    fclose(log);
  }
#endif

  VectorXf difference(OBSERVATION_DIMENSIONS, 1);
  difference.setZero();
  // flipped X, order should be Alt, Velo, Accel, thats why
  // the order is 2, 1, 0
  // gpsAlt isn't part of the X vector, so it doesn't get flipped.
  difference << (this->X[2] - zMean(0)), (this->X[1] - zMean(1)),
      (this->X[0] - zMean(2)), (this->gpsAlt - zMean(3));

  float nis = (difference.transpose() * Pz.inverse() * difference)(0, 0);

#if defined(LOGON) || defined(LOGMETRICS)
  // Log NIS to file
  FILE* nisFile = fopen((directoryPath + "/nis.txt").c_str(), "a+");
  if (nisFile) {
    fprintf(nisFile, "%f,%f\n", time, nis);
    fclose(nisFile);
  } else {
    perror("Error opening nis.txt");
  }
#endif

  X0 = this->Xprediction + K * difference;

  if (std::isnan(X0(0)) || std::isnan(X0(1)) || std::isnan(X0(2))) {
#ifdef TESTING_BUILD
    FILE* log = fopen(
        "log.txt",
        "a+");  // Open the file for appending or create it if it doesn't exist

    if (!log) {
      fprintf(stderr, "Error opening log.txt...exiting\n");
      exit(1);
    }

    fprintf(log,
            "At time %f s NAN detected in the state update, defaulting to "
            "Prediction as Estimation\n",
            this->time);

    fclose(log);
#endif
    // default to prediction
    X0 = this->Xprediction;
  }

  // check and update before apogee bool
  if (isAfterApogee == 1) {
#ifdef TIMERON
    std::chrono::high_resolution_clock::time_point getScenario =
        std::chrono::high_resolution_clock::now();
#endif

    std::vector<Scenario>* scenarios = this->getScenarios();

#ifdef TIMERON
    this->getScenarioTime +=
        std::chrono::duration_cast<std::chrono::duration<float>>(
            std::chrono::high_resolution_clock::now() - getScenario);
#endif

    for (int i = 0; i < scenarios->size(); i++) {
      scenarios->at(i).setIsBeforeApogee(false);
    }

  } else {
    // check rocket state with filter values
    WindowStats ws = this->computeWindowStats(
    Measurement{X0[0], X0[1], X0[2], this->time});

    if (!isAfterBurnout) {
      isAfterBurnout = this->burnoutDetection(ws);
      if (isAfterBurnout) {
        // [CANBUS DAQ] request switch from burn → coast state
      }
    }

    if (isAfterBurnout && !isAfterApogee) {

      // ---- Airbrake control ----------------------------------------
      airbrakeLevel_ = airbrakeController_.calculate_level(
        this->prevGain1, this->prevGain2, this->scenario_index_1, this->scenario_index_2, (uint32_t)ws.avgAltitude
      );
      // [CANBUS DAQ] airbrakeLevel_ published to CANBUS for airbrake actuation
      // ---------------------------------------------------------------

      isAfterApogee = this->apogeeDetection(ws);
      if (isAfterApogee) {
        // [CANBUS DAQ] request switch to after-apogee state
      }
    }
  }

  this->KinematicsHalo.altitudeStore = X0(0);

  // update the estimate covariance matrix
  // when variations in the estimate are small then its good
  MatrixXf P1(3, 3);
  P1.setZero();

  P1 = Pprediction - (K * Pz * K.transpose());

  this->P = P1;

#if defined(LOGON) || defined(LOGMETRICS)
  // write diagonal of P to file
  FILE* file = fopen((directoryPath + "/P.txt").c_str(), "a+");
  if (!file) {
    fprintf(stderr, "Error opening P.txt...exiting\n");
    exit(1);
  }
  fprintf(file, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", time, P1(0, 0), P1(0, 1),
          P1(0, 2), P1(1, 0), P1(1, 1), P1(1, 2), P1(2, 0), P1(2, 1), P1(2, 2));

  fclose(file);

#endif

  std::chrono::high_resolution_clock::time_point predictTimer;
  std::chrono::high_resolution_clock::time_point endUpdateTime;

#ifdef TIMERON

  endUpdateTime = std::chrono::high_resolution_clock::now();

  this->updateTime += std::chrono::duration_cast<std::chrono::duration<float>>(
      std::chrono::high_resolution_clock::now() - stateUpdateTime);

  predictTimer = std::chrono::high_resolution_clock::now();

#endif

  calculateSigmaPoints();

#ifdef TIMERON

  this->predictTime += std::chrono::duration_cast<std::chrono::duration<float>>(
      std::chrono::high_resolution_clock::now() - predictTimer);

#endif
}

// ------------------------------------------------

// Prediction--------------------------------------
void HALO::calculateSigmaPoints() {
  float mutliplier = 3;  // N - lambda

  std::chrono::high_resolution_clock::time_point tTime;

#ifdef TIMERON

  tTime = std::chrono::high_resolution_clock::now();

#endif

  MatrixXf L(((mutliplier)*P).llt().matrixL());

#ifdef TIMERON

  this->triangulationTime +=
      std::chrono::duration_cast<std::chrono::duration<float>>(
          std::chrono::high_resolution_clock::now() - tTime);

#endif

#ifdef TIMERON

  std::chrono::high_resolution_clock::time_point startSPoint =
      std::chrono::high_resolution_clock::now();

#endif

  // Initialize sigma points matrix
  MatrixXf sigmaPoints(3, 7);
  sigmaPoints.setZero();

  // Set the first sigma point
  sigmaPoints.col(0) = X0;

  // Set the remaining sigma points
  for (int i = 1; i < this->N1 + 1; i++) {
    sigmaPoints.col(i) = X0 + L.col(i - 1);
  }

  for (int j = this->N1 + 1; j < (2 * this->N1) + 1; j++) {
    sigmaPoints.col(j) = X0 - L.col(j - this->N1 - 1);
  }

#ifdef TIMERON

  this->sPointTime += std::chrono::duration_cast<std::chrono::duration<float>>(
      std::chrono::high_resolution_clock::now() - startSPoint);

#endif

#ifdef LOGON

  FILE* file = fopen((directoryPath + "/gains.txt").c_str(), "a+");
  if (!file) {
    fprintf(stderr, "Error opening gains.txt...exiting\n");
    exit(1);
  }

  FILE* sigmaPointsFile =
      fopen((directoryPath + "/sigmaPoints.txt").c_str(), "a+");
  if (!sigmaPointsFile) {
    fprintf(stderr, "Error opening sigmaPoints.txt...exiting\n");
    exit(1);
  }

  FILE* sigmaPointsFile1 =
      fopen((directoryPath + "/sigmaPoints1.txt").c_str(), "a+");
  if (!sigmaPointsFile1) {
    fprintf(stderr, "Error opening sigmaPoints1.txt...exiting\n");
    exit(1);
  }

  FILE* sigmaPointsFile2 =
      fopen((directoryPath + "/sigmaPoints2.txt").c_str(), "a+");
  if (!sigmaPointsFile2) {
    fprintf(stderr, "Error opening sigmaPoints2.txt...exiting\n");
    exit(1);
  }

  FILE* sigmaPointsFile3 =
      fopen((directoryPath + "/sigmaPoints3.txt").c_str(), "a+");
  if (!sigmaPointsFile3) {
    fprintf(stderr, "Error opening sigmaPoints3.txt...exiting\n");
    exit(1);
  }

  FILE* sigmaPointsFile4 =
      fopen((directoryPath + "/sigmaPoints4.txt").c_str(), "a+");
  if (!sigmaPointsFile4) {
    fprintf(stderr, "Error opening sigmaPoints4.txt...exiting\n");
    exit(1);
  }

  FILE* sigmaPointsFile5 =
      fopen((directoryPath + "/sigmaPoints5.txt").c_str(), "a+");
  if (!sigmaPointsFile5) {
    fprintf(stderr, "Error opening sigmaPoints5.txt...exiting\n");
    exit(1);
  }

  FILE* sigmaPointsFile6 =
      fopen((directoryPath + "/sigmaPoints6.txt").c_str(), "a+");
  if (!sigmaPointsFile6) {
    fprintf(stderr, "Error opening sigmaPoints6.txt...exiting\n");
    exit(1);
  }

#endif

  // propagate sigma points through the dynamic model
  for (int i = 0; i < (2 * this->N1) + 1; i++) {
#ifdef TIMERON

    std::chrono::high_resolution_clock::time_point startPredictLoop =
        std::chrono::high_resolution_clock::now();

#endif

    // load variables
    VectorXf column = sigmaPoints.col(i);
    this->firstTimeForPoint = firstTime[i];
    this->prevGain1 = this->listOfGainsSigmaPoints[i].first;
    this->prevGain2 = this->listOfGainsSigmaPoints[i].second;

#ifdef TIMERON

    this->predictLoopTime +=
        std::chrono::duration_cast<std::chrono::duration<float>>(
            std::chrono::high_resolution_clock::now() - startPredictLoop);

#endif

    std::chrono::high_resolution_clock::time_point dynamicTime;

#ifdef TIMERON

    dynamicTime = std::chrono::high_resolution_clock::now();

#endif

    sigmaPoints.col(i) = dynamicModel(column);

#ifdef TIMERON

    this->dynamicModelTime +=
        std::chrono::duration_cast<std::chrono::duration<float>>(
            std::chrono::high_resolution_clock::now() - dynamicTime);

#endif

#ifdef TIMERON

    std::chrono::high_resolution_clock::time_point endPredictLoop =
        std::chrono::high_resolution_clock::now();

#endif

    this->listOfGainsSigmaPoints[i] = {this->prevGain1, this->prevGain2};
    this->firstTime[i] = this->firstTimeForPoint;
#ifdef LOGON
    fprintf(file, " ");
#endif

#ifdef TIMERON
    this->endPredictLoopTime +=
        std::chrono::duration_cast<std::chrono::duration<float>>(
            std::chrono::high_resolution_clock::now() - endPredictLoop);
#endif
  }

#ifdef LOGON

  fprintf(file, "\n");

  fprintf(sigmaPointsFile, "%f,%f,%f\n", sigmaPoints(0, 0), sigmaPoints(1, 0),
          sigmaPoints(2, 0));
  fprintf(sigmaPointsFile1, "%f,%f,%f\n", sigmaPoints(0, 1), sigmaPoints(1, 1),
          sigmaPoints(2, 1));
  fprintf(sigmaPointsFile2, "%f,%f,%f\n", sigmaPoints(0, 2), sigmaPoints(1, 2),
          sigmaPoints(2, 2));
  fprintf(sigmaPointsFile3, "%f,%f,%f\n", sigmaPoints(0, 3), sigmaPoints(1, 3),
          sigmaPoints(2, 3));
  fprintf(sigmaPointsFile4, "%f,%f,%f\n", sigmaPoints(0, 4), sigmaPoints(1, 4),
          sigmaPoints(2, 4));
  fprintf(sigmaPointsFile5, "%f,%f,%f\n", sigmaPoints(0, 5), sigmaPoints(1, 5),
          sigmaPoints(2, 5));
  fprintf(sigmaPointsFile6, "%f,%f,%f\n", sigmaPoints(0, 6), sigmaPoints(1, 6),
          sigmaPoints(2, 6));

  fclose(file);
  fclose(sigmaPointsFile);
  fclose(sigmaPointsFile1);
  fclose(sigmaPointsFile2);
  fclose(sigmaPointsFile3);
  fclose(sigmaPointsFile4);
  fclose(sigmaPointsFile5);
  fclose(sigmaPointsFile6);

#endif

#ifdef TIMERON

  std::chrono::high_resolution_clock::time_point preMeanStart =
      std::chrono::high_resolution_clock::now();

#endif

  // calculate the mean and covariance of the sigma points
  VectorXf xPreMean(3, 1);
  for (int row = 0; row < this->N1; row++) {
    float sum00 = 0;
    for (int col = 0; col < 2 * this->N1 + 1; col++) {
      sum00 += sigmaPoints(row, col) * WeightsForSigmaPoints(col);
    }
    xPreMean(row) = sum00;
  }

#ifdef TIMERON

  this->preMeanTime += std::chrono::duration_cast<std::chrono::duration<float>>(
      std::chrono::high_resolution_clock::now() - preMeanStart);

#endif

  this->Xprediction = xPreMean;

#ifdef TIMERON

  std::chrono::high_resolution_clock::time_point projErrorStart =
      std::chrono::high_resolution_clock::now();

#endif

  MatrixXf projError(3, 7);
  projError.setZero(3, 7);

  for (int i = 0; i < this->N1; i++) {
    projError.row(i) =
        (sigmaPoints.row(i).array() - (this->Xprediction).row(i).value())
            .matrix();
  }

  this->projectError = projError;

#ifdef TIMERON

  this->projErrorTime +=
      std::chrono::duration_cast<std::chrono::duration<float>>(
          std::chrono::high_resolution_clock::now() - projErrorStart);
#endif

  MatrixXf Pprediction(3, 3);
  Pprediction.setZero(3, 3);

#ifdef TIMERON

  std::chrono::high_resolution_clock::time_point pPredictionStart =
      std::chrono::high_resolution_clock::now();

#endif

  Pprediction =
      projError * WeightsForSigmaPoints.asDiagonal() * projError.transpose() +
      this->Q;

  this->Pprediction = Pprediction;

#ifdef TIMERON

  this->PpredictionTime +=
      std::chrono::duration_cast<std::chrono::duration<float>>(
          std::chrono::high_resolution_clock::now() - pPredictionStart);

#endif

  this->sigPoints = sigmaPoints;
}

VectorXf HALO::predictNStates(int n) {
#ifdef TIMERON
  std::chrono::high_resolution_clock::time_point predictNStates_start =
      std::chrono::high_resolution_clock::now();
#endif

  // values to be referenced
  int firstTimeForPoint_storage = 1;
  std::vector<float> prevGain1_storage = {0.5, 0.5, 0.5};
  std::vector<float> prevGain2_storage = {0.5, 0.5, 0.5};
  std::vector<std::vector<int>> scenariosGainsList_storage = {
      {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}};
  std::vector<std::pair<std::vector<float>, std::vector<float>>>
      listOfGainsSigmaPoints_storage = {{{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}},
                                        {{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}},
                                        {{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}},
                                        {{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}},
                                        {{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}},
                                        {{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}},
                                        {{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}}};
  int counterSigmaPoint_storage = 0;

  // references
  int& firstTimeForPoint = firstTimeForPoint_storage;
  std::vector<float>& prevGain1 = prevGain1_storage;
  std::vector<float>& prevGain2 = prevGain2_storage;
  std::vector<std::vector<int>>& scenariosGainsList =
      scenariosGainsList_storage;
  int& counterSigmaPoint = counterSigmaPoint_storage;

#ifdef TIMERON
  std::chrono::high_resolution_clock::time_point predictNStates_getScenario =
      std::chrono::high_resolution_clock::now();
#endif

  std::vector<Scenario> scenarios = *this->getScenarios();

#ifdef TIMERON
  this->predictNStates_getScenarioTime +=
      std::chrono::duration_cast<std::chrono::duration<float>>(
          std::chrono::high_resolution_clock::now() -
          predictNStates_getScenario);
#endif

#if defined(LOGON) || defined(LOGMETRICS)
  if (predictnalt == nullptr) {
    predictnalt = fopen((directoryPath + "/predictnalt.txt").c_str(), "a+");
    if (!predictnalt) {
      fprintf(stderr, "Error opening predictnalt.txt...exiting\n");
      exit(1);
    }
  }
#endif

  predictionCounter++;

  // prediction at time step 0. Uses state vector and current covariance matrix.
  VectorXf calculation =
      this->dynamicModelOnce(this->X0, firstTimeForPoint, prevGain1, prevGain2,
                             scenariosGainsList, counterSigmaPoint, scenarios);

#if defined(LOGON) || defined(LOGMETRICS)
  fprintf(predictnalt, "%f,%f,%f,%f,%d\n", this->time, calculation(0),
          calculation(1), calculation(2), predictionCounter);
#endif

  bool apogee = false;
  for (int i = 1; i < n; i++) {
    if (apogee) {
      for (int i = 0; i < scenarios.size(); i++) {
        scenarios.at(i).setIsBeforeApogee(false);
      }

    } else {
      WindowStats ws_pred{};
      ws_pred.windowFull         = true;   // forward prediction — treat as full
      ws_pred.avgAltitude        = calculation(0);
      ws_pred.avgVelocity        = calculation(1);
      ws_pred.avgAcceleration    = calculation(2);
      ws_pred.sqrtP_altitude     = std::sqrt(this->P(0, 0));
      ws_pred.sqrtP_velocity     = std::sqrt(this->P(1, 1));
      ws_pred.sqrtP_acceleration = std::sqrt(this->P(2, 2));
      apogee = this->apogeeDetection(ws_pred);
    }
// NOTE: We do NOT call computeWindowStats() here because that would
// push the prediction values into the real sliding buffer. We build
// WindowStats manually with windowFull=true to bypass the buffer check.
    calculation = this->dynamicModelOnce(
        calculation, firstTimeForPoint, prevGain1, prevGain2,
        scenariosGainsList, counterSigmaPoint, scenarios);

#if defined(LOGON) || defined(LOGMETRICS)
    fprintf(predictnalt, "%f,%f,%f,%f,%d\n", this->time + (float)i * timeStep,
            calculation(0), calculation(1), calculation(2), predictionCounter);
#endif
  }

#ifdef TIMERON
  this->predictNStatesTime +=
      std::chrono::duration_cast<std::chrono::duration<float>>(
          std::chrono::high_resolution_clock::now() - predictNStates_start);
#endif

  return calculation;
}

// Function to calculate the Euclidean distance between two 3D vectors
float HALO::euclideanDistance(const std::vector<float>& vec1,
                              const VectorXf& vec2) {
  float x1, y1, z1, x2, y2, z2;

  x1 = vec1[0];
  y1 = vec1[1];
  z1 = vec1[2];

  x2 = vec2(0);
  y2 = vec2(1);
  z2 = vec2(2);

  return std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2) +
                   std::pow(z2 - z1, 2));
}

/**
 * @brief Given a list of scenarios, find the nearest 2 scenarios and returns
 * the vectors of the nearest scenarios
 */
std::pair<std::vector<int>, std::vector<std::vector<float>>>
HALO::findNearestScenarios(std::vector<Scenario>* scenarios,
                           VectorXf& measurement) {
  std::vector<std::pair<float, std::pair<int, int>>> distances;
  distances.reserve(scenarios->size());
  float minDistance = std::numeric_limits<float>::max();
  int i = 0;

#ifdef TIMERON

  std::chrono::high_resolution_clock::time_point loopScenariosStart =
      std::chrono::high_resolution_clock::now();

#endif

  for (size_t s = 0; s < scenarios->size(); s++) {
#ifdef TIMERON

    std::chrono::high_resolution_clock::time_point getListsStart =
        std::chrono::high_resolution_clock::now();

#endif

    i = 0;
    int lowestDistanceIndex = 0;

#ifdef TIMERON

    this->getListsTime +=
        std::chrono::duration_cast<std::chrono::duration<float>>(
            std::chrono::high_resolution_clock::now() - getListsStart);

#endif

#ifdef TIMERON

    std::chrono::high_resolution_clock::time_point othersTimeStart =
        std::chrono::high_resolution_clock::now();

#endif

    minDistance = std::numeric_limits<float>::max();

    std::pair<std::vector<float>, size_t> vect = {{0, 0, 0}, 0};

    std::vector<float> measurementVec = {measurement(0), measurement(1),
                                         measurement(2)};

#ifdef TIMERON

    this->othersTime +=
        std::chrono::duration_cast<std::chrono::duration<float>>(
            std::chrono::high_resolution_clock::now() - othersTimeStart);

#endif

#ifdef TIMERON

    std::chrono::high_resolution_clock::time_point KDTreeTimeStart =
        std::chrono::high_resolution_clock::now();

#endif

    vect = (scenarios->at(s)).nearestKDTree(measurementVec);

    /*std::cout << "Time: " << this->time << " {";
    for (auto& f : vect.first) std::cout << f << " ";
    std::cout << "}, " << vect.second << "\n";*/

#ifdef TIMERON

    this->KDTreeTime +=
        std::chrono::duration_cast<std::chrono::duration<float>>(
            std::chrono::high_resolution_clock::now() - KDTreeTimeStart);

#endif

#ifdef TIMERON

    std::chrono::high_resolution_clock::time_point euclideanStart =
        std::chrono::high_resolution_clock::now();

#endif

    minDistance = euclideanDistance(vect.first, measurement);

    int currentLowestDistanceIndex = static_cast<int>(vect.second);

#ifdef TIMERON

    this->euclideanTime +=
        std::chrono::duration_cast<std::chrono::duration<float>>(
            std::chrono::high_resolution_clock::now() - euclideanStart);

#endif

#ifdef TIMERON

    std::chrono::high_resolution_clock::time_point emplace_BackStart =
        std::chrono::high_resolution_clock::now();

#endif

    // pass index of scenario / struct
    std::pair<float, std::pair<int, int>> vector = {
        minDistance,
        {currentLowestDistanceIndex, (((scenarios->at(s)).name) - 1)}};
    distances.push_back(vector);

#ifdef TIMERON

    this->emplaceBackTime +=
        std::chrono::duration_cast<std::chrono::duration<float>>(
            std::chrono::high_resolution_clock::now() - emplace_BackStart);

#endif
  }

#ifdef TIMERON

  this->loopScenariosTime +=
      std::chrono::duration_cast<std::chrono::duration<float>>(
          std::chrono::high_resolution_clock::now() - loopScenariosStart);

#endif

#ifdef TIMERON

  std::chrono::high_resolution_clock::time_point startTime =
      std::chrono::high_resolution_clock::now();

#endif

  float lowestDistance = std::numeric_limits<float>::max();
  int lowestDistanceIndex = 0;
  float secondLowestDistance = std::numeric_limits<float>::max();
  int secondLowestDistanceIndex = 0;

  for (size_t s = 0; s < scenarios->size(); s++) {
    if (distances[s].first < lowestDistance) {
      secondLowestDistance = lowestDistance;
      secondLowestDistanceIndex = lowestDistanceIndex;

      lowestDistance = distances[s].first;
      lowestDistanceIndex = s;
    } else if (distances[s].first < secondLowestDistance) {
      secondLowestDistance = distances[s].first;
      secondLowestDistanceIndex = s;
    }
  }

#ifdef TIMERON

  this->twoDistancesTime +=
      std::chrono::duration_cast<std::chrono::duration<float>>(
          std::chrono::high_resolution_clock::now() - startTime);

#endif

#ifdef LOGON

  FILE* nearestScenariosFile =
      fopen((directoryPath + "/nearestScenarios.txt").c_str(), "a+");
  if (!nearestScenariosFile) {
    fprintf(stderr, "Error opening nearestScenarios.txt...exiting\n");
    exit(1);
  }

  FILE* nearestScenariosFormattedFile =
      fopen((directoryPath + "/nearestScenariosFormatted.txt").c_str(), "a+");
  if (!nearestScenariosFormattedFile) {
    fprintf(stderr, "Error opening nearestScenariosFormatted.txt...exiting\n");
    fprintf(stderr, "%d\n", errno);
    exit(1);
  }

  fprintf(nearestScenariosFile, "%f,%f,", lowestDistance, secondLowestDistance);
  fprintf(nearestScenariosFile, "%d,%d\n", lowestDistanceIndex,
          secondLowestDistanceIndex);

  fprintf(nearestScenariosFormattedFile,
          "For Meas(%f,%f,%f) lowest(%f),secondL(%f),", measurement[0],
          measurement[1], measurement[2], lowestDistance, secondLowestDistance);
  fprintf(nearestScenariosFormattedFile, "lowestName(%d),secondLN(%d)\n",
          distances[lowestDistanceIndex].second.second,
          distances[secondLowestDistanceIndex].second.second);
  fprintf(nearestScenariosFormattedFile, "list: %f,%f,%f,%f,%f,%f\n",
          distances[0].first, distances[1].first, distances[2].first,
          distances[3].first, distances[4].first, distances[5].first);

  fclose(nearestScenariosFile);
  fclose(nearestScenariosFormattedFile);

#endif

#ifdef TIMERON

  std::chrono::high_resolution_clock::time_point startVectors =
      std::chrono::high_resolution_clock::now();

#endif

  // find current vector (by index) and future vector (by time)
  int indexFirst = distances[lowestDistanceIndex].second.first;

  if (indexFirst < 0 || indexFirst >= 10000) {
    SOAR_PRINT("indexFirst is garbage: %d\n", indexFirst);
    return {};
  }

  Scenario* scenario1 =
      &scenarios->at(distances[lowestDistanceIndex].second.second);
  std::vector<float> currentVector1 = scenario1->evaluateVectorAt(indexFirst);
  float nextTimeStep = currentVector1[3] + deltaTime;

  std::vector<float> futureVector1 =
      scenario1->evaluateVectorAtTime(nextTimeStep);

  int indexSecond = distances[secondLowestDistanceIndex].second.first;
  Scenario* scenario2 =
      &scenarios->at(distances[secondLowestDistanceIndex].second.second);
  std::vector<float> currentVector2 = scenario2->evaluateVectorAt(indexSecond);
  float nextTimeStep2 = currentVector2[3] + deltaTime;
  std::vector<float> futureVector2 =
      scenario2->evaluateVectorAtTime(nextTimeStep2);

#ifdef TIMERON

  this->vectorsTime += std::chrono::duration_cast<std::chrono::duration<float>>(
      std::chrono::high_resolution_clock::now() - startVectors);

#endif

#ifdef TIMERON

  std::chrono::high_resolution_clock::time_point push_backStart =
      std::chrono::high_resolution_clock::now();

#endif

  std::vector<std::vector<float>> nearestVectors;
  nearestVectors.push_back(currentVector1);
  nearestVectors.push_back(futureVector1);
  nearestVectors.push_back(currentVector2);
  nearestVectors.push_back(futureVector2);

#ifdef TIMERON

  this->push_backTime +=
      std::chrono::duration_cast<std::chrono::duration<float>>(
          std::chrono::high_resolution_clock::now() - push_backStart);

#endif

#ifdef LOGON

  FILE* file3 = fopen((directoryPath + "/predictedValues.txt").c_str(), "a+");
  if (!file3) {
    fprintf(stderr, "Error opening predictedValues.txt...exiting\n");
    exit(1);
  }

  fprintf(file3, "%f,%f,%f,%f,", currentVector1[0], currentVector1[1],
          currentVector1[2], currentVector1[3]);
  fprintf(file3, "%f,%f,%f,%f,", futureVector1[0], futureVector1[1],
          futureVector1[2], futureVector1[3]);
  fprintf(file3, "%f,%f,%f,%f,", currentVector2[0], currentVector2[1],
          currentVector2[2], currentVector2[3]);
  fprintf(file3, "%f,%f,%f,%f\n", futureVector2[0], futureVector2[1],
          futureVector2[2], futureVector2[3]);

  fclose(file3);

#endif

  std::pair<std::vector<int>, std::vector<std::vector<float>>>
      nearestVectorsWithIndices = {std::make_pair(
          std::vector<int>{lowestDistanceIndex, secondLowestDistanceIndex},
          nearestVectors)};

  return nearestVectorsWithIndices;
}

/**
 * @brief Predicts the next values based on the interpolated scenarios
 */
VectorXf HALO::predictNextValues(std::vector<std::vector<float>>& vectors,
                                 VectorXf& X_in, int scenario1Index,
                                 int scenario2Index) {
  std::vector<float> gainV1 = {0, 0, 0};
  std::vector<float> gainV2 = {0, 0, 0};
  std::vector<float> vector1 = vectors[0];
  std::vector<float> vector1Future = vectors[1];

  std::vector<float> vector2 = vectors[2];
  std::vector<float> vector2Future = vectors[3];
  bool both = false;
  bool vector1Further = false;

  for (int i = 0; i < 3; i++) {
    float distance = std::abs(vector1[i] - vector2[i]);

    if (distance < 1) {
      // similar so defaulting
      gainV1[i] = 0.5;
      gainV2[i] = 0.5;
    } else {
      if (vector1[i] > X_in(i)) {
        // below vector1
        if (vector2[i] > X_in(i)) {
          // point below both lines
          both = true;
          if (vector1[i] > vector2[i]) {
            // vector1 on top
            vector1Further = true;
          } else {
            vector1Further = false;
          }
        }
      } else {
        if (vector2[i] > X_in(i)) {
          // point between lines
        } else {
          // point above both lines
          if (vector1[i] < vector2[i]) {
            vector1Further = true;
          } else {
            vector1Further = false;
          }
          both = true;
        }
      }

      // point below both lines
      // -----------------------------------------------------------
      if (both) {
        float longestDistance = 0;
        float distance1 = 0;

        if (vector1Further) {
          longestDistance = std::abs(vector1[i] - X_in(i));
          distance1 = std::abs(vector2[i] - X_in(i));
        } else {
          longestDistance = std::abs(vector2[i] - X_in(i));
          distance1 = std::abs(vector1[i] - X_in(i));
        }

        float relativeFactor = longestDistance / distance1;

        if (!vector1Further) {
          gainV1[i] = relativeFactor / (relativeFactor + 1);
        } else {
          gainV1[i] = 1 - (relativeFactor / (relativeFactor + 1));
        }

        gainV2[i] = 1 - gainV1[i];
      } else {
        // point between lines
        // ------------------------------------------------------------
        gainV1[i] =
            1 - (std::abs(vector1[i] - X_in(i)) /
                 distance);  // get distance between vector1 and current state
        gainV2[i] = 1 - gainV1[i];
        //--------------------------------------------------------------------------------
      }
    }
  }

  if (this->firstTimeForPoint == 1) {
#ifdef LOGON
    FILE* file = fopen((directoryPath + "/log.txt").c_str(), "a+");
    if (!file) {
      fprintf(stderr, "Error opening log.txt...exiting\n");
      exit(1);
    }

    fprintf(file, "First time for point (%f, %f, %f)\n", X_in(0), X_in(1),
            X_in(2));

    fclose(file);
#endif

    this->prevGain1 = gainV1;
    this->prevGain2 = gainV2;
    this->firstTimeForPoint = 0;
  }

#ifdef LOGON

  FILE* gainsFile = fopen((directoryPath + "/gains.txt").c_str(), "a+");
  if (!gainsFile) {
    fprintf(stderr, "Error opening gains.txt...exiting\n");
    exit(1);
  }

  fprintf(gainsFile, "%f, %f, %f,", gainV1[0], gainV1[1], gainV1[2]);
  fprintf(gainsFile, "%f, %f, %f\n", gainV2[0], gainV2[1], gainV2[2]);

  fclose(gainsFile);

#endif

  // check if scenario indices are the same
  // if not reset the gains
  if (this->scenariosGainsList[this->counterSigmaPoint][0] != scenario1Index ||
      this->scenariosGainsList[this->counterSigmaPoint][1] != scenario2Index) {
    this->prevGain1 = {0.5, 0.5, 0.5};
    this->prevGain2 = {0.5, 0.5, 0.5};

// write to resetGainsFile.txt, already opened
#ifdef LOGON
    fprintf(resetGainsFile, "%d\n", 6);
#endif
  } else {
// not resetting gains, scenarios are the same from prior iteration
#ifdef LOGON
    // write to resetGainsFile.txt, already opened
    fprintf(resetGainsFile, "%d\n", -1);
#endif
  }

  // interpolate between the two scenarios to get predicted values
  float predicted_interpolated_alt = this->prevGain1[0] * vector1Future[0] +
                                     this->prevGain2[0] * vector2Future[0];
  float predicted_interpolated_velo = this->prevGain1[1] * vector1Future[1] +
                                      this->prevGain2[1] * vector2Future[1];
  float predicted_interpolated_acc = this->prevGain1[2] * vector1Future[2] +
                                     this->prevGain2[2] * vector2Future[2];

  // save scenario indices
  this->scenariosGainsList[this->counterSigmaPoint] = {scenario1Index,
                                                       scenario2Index};

  this->prevGain1 = gainV1;
  this->prevGain2 = gainV2;

  VectorXf X_pred(3, 1);
  X_pred << predicted_interpolated_alt, predicted_interpolated_velo,
      predicted_interpolated_acc;

  // increment counter
  this->counterSigmaPoint = this->counterSigmaPoint + 1;
  // reset counter
  this->counterSigmaPoint = this->counterSigmaPoint % 6;

  return X_pred;
}

/**
 * @brief Predicts the next values based on the interpolated scenarios
 */
VectorXf HALO::predictNextValuesOnce(
    std::vector<std::vector<float>>& vectors, VectorXf& X_in,
    int firstTimeForPoint, int scenario1Index, int scenario2Index,
    std::vector<float>& prevGain1, std::vector<float>& prevGain2,
    std::vector<std::vector<int>>& scenariosGainsList, int& counterSigmaPoint) {
  std::vector<float> gainV1 = {0, 0, 0};
  std::vector<float> gainV2 = {0, 0, 0};
  std::vector<float> vector1 = vectors[0];
  std::vector<float> vector1Future = vectors[1];

  std::vector<float> vector2 = vectors[2];
  std::vector<float> vector2Future = vectors[3];
  bool both = false;
  bool vector1Further = false;

  for (int i = 0; i < 3; i++) {
    float distance = std::abs(vector1[i] - vector2[i]);

    if (distance < 1) {
      // similar so defaulting
      gainV1[i] = 0.5;
      gainV2[i] = 0.5;
    } else {
      if (vector1[i] > X_in(i)) {
        // below vector1
        if (vector2[i] > X_in(i)) {
          // point below both lines
          both = true;
          if (vector1[i] > vector2[i]) {
            // vector1 on top
            vector1Further = true;
          } else {
            vector1Further = false;
          }
        }
      } else {
        if (vector2[i] > X_in(i)) {
          // point between lines
        } else {
          // point above both lines
          if (vector1[i] < vector2[i]) {
            vector1Further = true;
          } else {
            vector1Further = false;
          }
          both = true;
        }
      }

      // point below both lines
      // -----------------------------------------------------------
      if (both) {
        float longestDistance = 0;
        float distance1 = 0;

        if (vector1Further) {
          longestDistance = std::abs(vector1[i] - X_in(i));
          distance1 = std::abs(vector2[i] - X_in(i));
        } else {
          longestDistance = std::abs(vector2[i] - X_in(i));
          distance1 = std::abs(vector1[i] - X_in(i));
        }

        float relativeFactor = longestDistance / distance1;

        if (!vector1Further) {
          gainV1[i] = relativeFactor / (relativeFactor + 1);
        } else {
          gainV1[i] = 1 - (relativeFactor / (relativeFactor + 1));
        }

        gainV2[i] = 1 - gainV1[i];
      } else {
        // point between lines
        // ------------------------------------------------------------
        gainV1[i] =
            1 - (std::abs(vector1[i] - X_in(i)) /
                 distance);  // get distance between vector1 and current state
        gainV2[i] = 1 - gainV1[i];
        //--------------------------------------------------------------------------------
      }
    }
  }

  if (firstTimeForPoint == 1) {
    prevGain1 = gainV1;
    prevGain2 = gainV2;
    firstTimeForPoint = 0;
  }

  /*
  if (scenariosGainsList[counterSigmaPoint][0] != scenario1Index ||
      scenariosGainsList[counterSigmaPoint][1] != scenario2Index) {
    prevGain1 = {0.5, 0.5, 0.5};
    prevGain2 = {0.5, 0.5, 0.5};
  }*/

  float predicted_interpolated_alt =
      prevGain1[0] * vector1Future[0] + prevGain2[0] * vector2Future[0];
  float predicted_interpolated_velo =
      prevGain1[1] * vector1Future[1] + prevGain2[1] * vector2Future[1];
  float predicted_interpolated_acc =
      prevGain1[2] * vector1Future[2] + prevGain2[2] * vector2Future[2];

  scenariosGainsList[counterSigmaPoint] = {scenario1Index, scenario2Index};

  prevGain1 = gainV1;
  prevGain2 = gainV2;

  VectorXf X_pred(3, 1);
  X_pred << predicted_interpolated_alt, predicted_interpolated_velo,
      predicted_interpolated_acc;

  counterSigmaPoint = counterSigmaPoint + 1;
  // reset counter
  counterSigmaPoint = counterSigmaPoint % 6;

  return X_pred;
}
/**
 * @brief Take the filtered values from Everest filter
 */
void HALO::setStateVector(float filteredAcc, float filteredVelo,
                          float filteredAlt, float gpsAlt) {
  this->Uaccel = filteredAcc;
  this->Uvelo = filteredVelo;
  this->Ualt = filteredAlt;
  this->gpsAlt = gpsAlt;

  VectorXf X_in(3);
  X_in << this->Uaccel, this->Uvelo, this->Ualt;

  this->X = X_in;

  this->stateUpdate();
}

// prediction step based on the dynamic model
VectorXf HALO::dynamicModel(VectorXf& X) {
  VectorXf Xprediction(3, 1);

#ifdef TIMERON

  std::chrono::high_resolution_clock::time_point getScenario =
      std::chrono::high_resolution_clock::now();

#endif

  // for every scenario get lists and find nearest 2 vectors to the current
  // state
  std::vector<Scenario>* scenarios = this->getScenarios();

#ifdef TIMERON

  this->getScenarioTime +=
      std::chrono::duration_cast<std::chrono::duration<float>>(
          std::chrono::high_resolution_clock::now() - getScenario);

#endif

  // check if X is nan, if so default to static integration
  if (std::isnan(X(0)) || std::isnan(X(1)) || std::isnan(X(2))) {
#ifdef TESTING_BUILD
    FILE* file = fopen((directoryPath + "/log.txt").c_str(), "a+");
    if (!file) {
      fprintf(stderr, "Error opening log.txt...exiting\n");
      exit(1);
    }
    fprintf(file, "At %f X is nan, defaulting to static integration\n",
            this->time);
    fclose(file);
#endif

    SOAR_PRINT("X is nan, defaulting to static integration\n");

    float finalVelocity = X(1) + X(0) * getDeltaTime();
    float altitude = X(2) + (X(1) + finalVelocity) * getDeltaTime() / 2.0;

    Xprediction(0) = altitude;
    Xprediction(1) = finalVelocity;
    Xprediction(2) = X(0);

    return Xprediction;
  }

#ifdef TIMERON

  std::chrono::high_resolution_clock::time_point nearestVectorsStart =
      std::chrono::high_resolution_clock::now();

#endif

  std::pair<std::vector<int>, std::vector<std::vector<float>>>
      nearestVectorsWithIndex = this->findNearestScenarios(scenarios, X);

  std::vector<std::vector<float>> nearestVectors =
      nearestVectorsWithIndex.second;

  this->lastNearestVectors_ = nearestVectors; // cache for airbrake controller

  int scenario1Index = nearestVectorsWithIndex.first[0];
  int scenario2Index = nearestVectorsWithIndex.first[1];

#ifdef TIMERON

  this->nearestScenariosTime +=
      std::chrono::duration_cast<std::chrono::duration<float>>(
          std::chrono::high_resolution_clock::now() - nearestVectorsStart);

#endif

  Xprediction =
      predictNextValues(nearestVectors, X, scenario1Index, scenario2Index);

  return Xprediction;
}

// prediction step based on the dynamic model
VectorXf HALO::dynamicModelOnce(
    VectorXf& X, int firstTimeForPoint, std::vector<float>& prevGain1,
    std::vector<float>& prevGain2,
    std::vector<std::vector<int>>& scenariosGainsList, int& counterSigmaPoint,
    std::vector<Scenario>& scenarios) {
  VectorXf Xprediction(3, 1);

  // for every scenario get lists and find nearest 2 vectors to the current
  // state

  // check if X is nan, if so default to static integration

  if (std::isnan(X(0)) || std::isnan(X(1)) || std::isnan(X(2))) {
#ifdef TESTING_BUILD
    FILE* file = fopen((directoryPath + "/log.txt").c_str(), "a+");
    if (!file) {
      fprintf(stderr, "Error opening log.txt...exiting\n");
      exit(1);
    }
    fprintf(file, "At %f X is nan, defaulting to static integration\n",
            this->time);
    fclose(file);
#endif
    SOAR_PRINT("X is nan, defaulting to static integration\n");

    float finalVelocity = X(1) + X(0) * getDeltaTime();
    float altitude = X(2) + (X(1) + finalVelocity) * getDeltaTime() / 2.0;

    Xprediction(0) = altitude;
    Xprediction(1) = finalVelocity;
    Xprediction(2) = X(0);

    return Xprediction;
  }


#ifdef TIMERON
  std::chrono::high_resolution_clock::time_point
      predictNStates_nearestVectorsStart =
          std::chrono::high_resolution_clock::now();
#endif

  std::pair<std::vector<int>, std::vector<std::vector<float>>>
      nearestVectorsWithIndex = this->findNearestScenarios(&scenarios, X);

#ifdef TIMERON
  this->predictNStates_nearestScenariosTime +=
      std::chrono::duration_cast<std::chrono::duration<float>>(
          std::chrono::high_resolution_clock::now() -
          predictNStates_nearestVectorsStart);
#endif

  std::vector<std::vector<float>> nearestVectors =
      nearestVectorsWithIndex.second;

  int scenario1Index = nearestVectorsWithIndex.first[0];
  int scenario2Index = nearestVectorsWithIndex.first[1];

  Xprediction = predictNextValuesOnce(
      nearestVectors, X, scenario1Index, scenario2Index, firstTimeForPoint,
      prevGain1, prevGain2, scenariosGainsList, counterSigmaPoint);

  return Xprediction;
}

/**
 * @brief Create scenarios for HALO
 */
void HALO::createScenarios(HALO* halo) {
#ifdef TIMERON

  std::chrono::high_resolution_clock::time_point treeCreation =
      std::chrono::high_resolution_clock::now();

#endif

  initAllSimData();

  Scenario scenario1 = Scenario{beforeApogeeSim1, afterApogeeSim1, 1};
  scenario1.createTree();
  Scenario scenario2 = Scenario{beforeApogeeSim2, afterApogeeSim2, 2};
  scenario2.createTree();
  // TODO: add more as there are more
  // Scenario scenario3 = Scenario{sim3, sim3, 3};
  // scenario3.createTree();

  // TODO: add them here as well {scenario3, scenario4, scenario5, scenario6}
  this->scenarios.clear();
  this->scenarios.push_back(scenario1);
  this->scenarios.push_back(scenario2);

  // lookup for apogee extrapolation
  std::vector<float> apogees;
  apogees.reserve(2);

  apogees.push_back(beforeApogeeSim1[0][-1]);
  apogees.push_back(beforeApogeeSim2[0][-1]);

  this->airbrakeController_.init(apogees);

#ifdef TIMERON

  std::chrono::high_resolution_clock::time_point treeCreationEnd =
      std::chrono::high_resolution_clock::now();

  std::chrono::duration<float> treeCreationTime =
      std::chrono::duration_cast<std::chrono::duration<float>>(treeCreationEnd -
                                                               treeCreation);

  halo->treeCreationTime = treeCreationTime;

#endif

}

/**
 * @brief Initialize HALO filter, created HALO object, scenarios and initial
 * state vector from calibration
 */
void HALO::initializeHALO(float initialAlt, HALO* halo) {
  // set initial state (altitude, velocity, acceleration)
  VectorXf X0(3);
  X0 << initialAlt, 0, 0;

  // process noise Covariance matrix (altitude, velocity, acceleration)
  // Calculated using covarianceCalc.py -> covariance matrix from Altimeter
  // Assuming Altimeter has no process noise (Q = 0)
  // unaccounted for noise in environment (wind, etc)-> residual from sims
  MatrixXf Q(3, 3);
  Q << 74777.41, 4458.13, -2164.91, 4458.13, 2413.02, -4.52, -2164.91, -4.52,
      503.78;

  // Measurement Covariance matrix (altitude, velocity, acceleration)
  // Calculated using covarianceCalc.py -> residual from Everest
  // Below, we expand the matrix with the measurement noise for GPS.
  MatrixXf R0(3, 3);
  R0 << 15438.09, 1528.51, 727.68, 1528.51, 5005.79, -469.20, 727.68, -469.20,
      613.20;

  // Initial state covariance matrix
  MatrixXf P0(3, 3);
  P0 << 50, 0, 0, 0, 0, 0, 0, 0, 0;

  // create scenarios
  createScenarios(halo);

  // Initialize with tare / GPS values
  halo->init(X0, P0, Q, R0);
}

void HALO::initializeHALOWithQR(float initialAlt, HALO* halo, MatrixXf& Q,
                                MatrixXf& R0) {
  // set initial state (altitude, velocity, acceleration)
  VectorXf X0(3);
  X0 << initialAlt, 0, 0;

  // Initial state covariance matrix
  MatrixXf P0(3, 3);
  P0 << 50, 0, 0, 0, 0, 0, 0, 0, 0;

  // create scenarios
  createScenarios(halo);

  // Initialize with tare / GPS values
  halo->init(X0, P0, Q, R0);
}

std::vector<float> HALO::Halo_Input(HALO* haloPointer, bool isInitialized,
                                    float eAccelerationZ, float eVelocity,
                                    float eAltitude, float gpsAltitude,
                                    float time, float deltaTime) {
  setDeltaTime(deltaTime);
  std::vector<float> unitedStates = {0, 0, 0};

  if (isInitialized) {
    haloPointer->setTime(time);
    haloPointer->setStateVector(eAccelerationZ, eVelocity, eAltitude,
                                gpsAltitude);

    // X0 = {eAltitude, eVelocity, eAccelerationZ};
    unitedStates = {haloPointer->X0[0], haloPointer->X0[1], haloPointer->X0[2]};
  }

  if (counter == 524) {
#ifdef TIMERON
    std::cout << "Update time:\t\t\t\t\t\t\t\t\t\t"
              << haloPointer->updateTime.count() << std::endl;
    std::cout << "Predict time:\t\t\t\t\t\t\t\t\t\t"
              << haloPointer->predictTime.count() << std::endl;

    std::cout << "\tTriangulationTime:\t\t\t\t\t\t\t"
              << haloPointer->triangulationTime.count() << std::endl;
    std::cout << "\tdModeltime:\t\t\t\t\t\t\t\t"
              << haloPointer->dynamicModelTime.count() << std::endl;

    std::cout << "\t\tgetScenarioTime:\t\t\t\t\t"
              << haloPointer->getScenarioTime.count() << std::endl;
    std::cout << "\t\tpPredictionTime:\t\t\t\t\t"
              << haloPointer->PpredictionTime.count() << std::endl;
    std::cout << "\t\tprojErrorTime:\t\t\t\t\t\t"
              << haloPointer->projErrorTime.count() << std::endl;
    std::cout << "\t\tpreMeanTime:\t\t\t\t\t\t"
              << haloPointer->preMeanTime.count() << std::endl;
    std::cout << "\t\tsPointTime:\t\t\t\t\t\t"
              << haloPointer->sPointTime.count() << std::endl;
    std::cout << "\t\tpredictLoopTime:\t\t\t\t\t"
              << haloPointer->predictLoopTime.count() << std::endl;
    std::cout << "\t\tendPredictLoopTime:\t\t\t\t\t"
              << haloPointer->endPredictLoopTime.count() << std::endl;
    std::cout << "\t\tnearestScenariosTime:\t\t\t\t\t"
              << haloPointer->nearestScenariosTime.count() << std::endl;

    std::cout << "\t\t\t->loopScenariosTime:\t\t\t"
              << haloPointer->loopScenariosTime.count() << std::endl;
    std::cout << "\t\t\t\t->getListsTime:\t\t"
              << haloPointer->getListsTime.count() << std::endl;
    std::cout << "\t\t\t\t->othersTime:\t\t" << haloPointer->othersTime.count()
              << std::endl;
    std::cout << "\t\t\t\t->KDTreeTime:\t\t" << haloPointer->KDTreeTime.count()
              << std::endl;
    std::cout << "\t\t\t\t->twoDistancesTime:\t"
              << haloPointer->twoDistancesTime.count() << std::endl;
    std::cout << "\t\t\t\t->emplaceBackTime:\t"
              << haloPointer->emplaceBackTime.count() << std::endl;

    std::cout << "\t\t\t->vectorsTime:\t\t\t\t"
              << haloPointer->vectorsTime.count() << std::endl;
    std::cout << "\t\t\t->push_backTime:\t\t\t"
              << haloPointer->push_backTime.count() << std::endl;

    std::cout << "\t\t\t->->treeCreationTime:\t\t\t\t"
              << (haloPointer->treeCreationTime).count() << std::endl;

    // profiling for prediction step
    std::cout << "\t\tpredictNStatesTime:\t\t\t\t\t"
              << haloPointer->predictNStatesTime.count() << std::endl;
    std::cout << "\t\t\t->predictNStates_getScenarioTime:\t\t\t\t"
              << haloPointer->predictNStates_getScenarioTime.count()
              << std::endl;
    std::cout << "\t\t\t->predictNStates_nearestScenariosTime:\t\t\t\t"
              << haloPointer->predictNStates_nearestScenariosTime.count()
              << std::endl;

#endif

#ifdef LOGMETRICS
    fclose(predictnalt);
    fclose(kalmangains);
#endif
  }

  counter++;

  return unitedStates;
}

#endif
