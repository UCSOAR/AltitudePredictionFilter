#ifndef DATA_HPP_
#define DATA_HPP_

#include <vector>

const int num_scenarios_parsed = 2;

extern std::vector<std::vector<float>>* beforeApogeeSim1;
extern std::vector<std::vector<float>>* afterApogeeSim1;
extern std::vector<std::vector<float>>* beforeApogeeSim2;
extern std::vector<std::vector<float>>* afterApogeeSim2;

void initAllSimData();

#endif
