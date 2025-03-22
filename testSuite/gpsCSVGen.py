import pandas as pd

# Read the CSV file
data = pd.read_csv("testSuite/data/altimeter1.csv")

# Extract the altitude and time columns (assuming the column names are 'altitude' and 'time')
altitudes = data["altitude"]
times = data["time"]

# Generate the C++ file
cpp_file_path = "gpsData.cpp"
with open(cpp_file_path, "w") as cpp_file:
    cpp_file.write("// Auto-generated file containing altitude and time data\n")
    cpp_file.write("#include <vector>\n\n")
    cpp_file.write("std::vector<std::pair<float, float>> gpsData = {\n")

    # Write the altitude and time data as pairs into the C++ vector
    for i, (time, altitude) in enumerate(zip(times, altitudes)):
        if i < len(times) - 1:
            cpp_file.write(f"    {{{time}, {altitude}}},\n")
        else:
            cpp_file.write(f"    {{{time}, {altitude}}}\n")

    cpp_file.write("};\n")

print(f"Altitude and time data have been saved to {cpp_file_path}")
