import csv

# Read the IMU data from a CSV file
input_data = []
with open("testSuite/data/Imu_Baro.csv", "r") as file:
    reader = csv.reader(file)
    next(reader)  # Skip header
    for row in reader:
        input_data.append(row)

# Time variables
time_increment = 1 / 3
current_time = 0

# Count valid IMU rows
imu_count = 0
for data in input_data:
    if data[4] == "":
        break
    imu_count += 1

# Count valid baro rows
baro_count = 0
for data in input_data:
    if data[10] == "":
        break
    baro_count += 1


# C++ source
cpp_code = f"""
#include <array>
#include "input_data.hpp"

const std::array<std::array<float, 10>, {imu_count}> taberLaunch = {{
"""

# Generate IMU data
current_time = 0

for i, data in enumerate(input_data[:imu_count]):

    comma = "," if i < imu_count - 1 else ""

    cpp_code += (
        f"    std::array<float, 10>{{"
        f"{current_time}, "
        f"{float(data[3])/1000}, "
        f"{float(data[4])/1000}, "
        f"{float(data[5])/1000}, "
        f"{float(data[0])/1000}, "
        f"{float(data[1])/1000}, "
        f"{float(data[2])/1000}, "
        f"{float(data[6])/1000}, "
        f"{float(data[7])/1000}, "
        f"{float(data[8])/1000}"
        f"}}{comma}\n"
    )

    current_time += time_increment


cpp_code += f"""}};

const std::array<std::array<float, 4>, {baro_count}> baroData = {{
"""


# Generate Baro data
current_time = 0

for i, data in enumerate(input_data[:baro_count]):

    comma = "," if i < baro_count - 1 else ""

    cpp_code += (
        f"    std::array<float, 4>{{"
        f"{current_time}, "
        f"{data[10]}, "
        f"0, 0"
        f"}}{comma}\n"
    )

    current_time += time_increment


cpp_code += "};\n"


# Write cpp
with open("input_data.cpp", "w") as file:
    file.write(cpp_code)


# Header file
hpp_code = f"""
#ifndef INPUT_DATA_HPP_
#define INPUT_DATA_HPP_

#include <array>

extern const std::array<std::array<float, 10>, {imu_count}> taberLaunch;

extern const std::array<std::array<float, 4>, {baro_count}> baroData;

#endif
"""


with open("input_data.hpp", "w") as file:
    file.write(hpp_code)


print("input_data.cpp and input_data.hpp generated.")
