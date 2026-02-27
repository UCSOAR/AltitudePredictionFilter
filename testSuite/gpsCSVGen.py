import pandas as pd

# Read the CSV file
data = pd.read_csv("testSuite/data/altimeter1.csv")

# Extract columns
altitudes = data["altitude"]
times = data["time"]

num_rows = len(data)

# Generate the C++ file
cpp_file_path = "gpsData.cpp"

with open(cpp_file_path, "w") as cpp_file:
    cpp_file.write("// Auto-generated file containing altitude and time data\n")
    cpp_file.write("#include <array>\n")
    cpp_file.write("#include \"gpsData.hpp\"\n\n")

    cpp_file.write(
        f"const std::array<std::array<float, 2>, {num_rows}> gpsData1 = {{\n"
    )

    for i, (time, altitude) in enumerate(zip(times, altitudes)):
        comma = "," if i < num_rows - 1 else ""
        cpp_file.write(f"    std::array<float, 2>{{{time}, {altitude}}}{comma}\n")

    cpp_file.write("};\n")


# Generate header
hpp_file_path = "gpsData.hpp"

with open(hpp_file_path, "w") as hpp_file:
    hpp_code = f"""
#ifndef GPS_DATA_HPP_
#define GPS_DATA_HPP_

#include <array>

extern const std::array<std::array<float, 2>, {num_rows}> gpsData1;

#endif
"""
    hpp_file.write(hpp_code)


print(f"Altitude and time data saved to {cpp_file_path}")
