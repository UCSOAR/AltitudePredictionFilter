/*
 * CommonDefines.hpp
 *
 *  Created on: May 3, 2026
 *      Author: Andrey
 */

#ifndef ALTITUDEPREDICTIONFILTER_COMMONDEFINESFILTER_HPP_
#define ALTITUDEPREDICTIONFILTER_COMMONDEFINESFILTER_HPP_

#include "KDTree.hpp"

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

  std::vector<std::vector<float>>* BeforeList;
  std::vector<std::vector<float>>* AfterList;

  KDTree treeBefore;
  KDTree treeAfter;

  int name;

  std::vector<float> measurement;
  bool isBeforeApogeeBool = true;

  Measurement prediction{};

  // conversions between vector and array.
  inline std::array<float, 3> vec2arr(const std::vector<float> &v) {
    assert(v.size() == 3);
    return {v[0], v[1], v[2]};
  }

  inline std::vector<float> arr2vec(const std::array<float, 3> &a) {
    return {a[0], a[1], a[2]};
  }

  Scenario(std::vector<std::vector<float>>* beforeList,
           std::vector<std::vector<float>>* afterList, int Name)
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
    pointIndex result;
    if (isBeforeApogeeBool) {
      return treeBefore.nearest_pointIndex(measurement);
    } else {
      return treeAfter.nearest_pointIndex(measurement);
    }
  }

  void createTree() {
    std::vector<std::vector<float>> beforeVectorofVectors;
    std::vector<std::vector<float>> afterVectorofVectors;

    beforeVectorofVectors.reserve(BeforeList->size());
    afterVectorofVectors.reserve(AfterList->size());

    for (int i = 0; i < BeforeList->size(); i++) {
      std::vector<float> vect = {(*BeforeList)[i][0], (*BeforeList)[i][1],
                                 (*BeforeList)[i][2]};
      beforeVectorofVectors.push_back(vect);
    }

    for (int i = 0; i < AfterList->size(); i++) {
      std::vector<float> vect = {(*AfterList)[i][0], (*AfterList)[i][1],
                                 (*AfterList)[i][2]};
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
      return BeforeList;
    } else {
      return AfterList;
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
    auto *lists = getLists();

    return lists->at(index);
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


#endif /* ALTITUDEPREDICTIONFILTER_COMMONDEFINESFILTER_HPP_ */
