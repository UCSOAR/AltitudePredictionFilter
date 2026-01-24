import pandas as pd
import math

# avg_alt_1,avg_velo_1,avg_acc_1,avg_time_1
data = pd.read_csv("testSuite/data/final_may_sims_formatted.csv")

# Convert DataFrame to list of lists
data_list = data.values.tolist()

# skip the header row
data_list = data_list[1:]

# Replace invalid strings with NaN and convert to floats
invalid_values = ["", "#DIV/0!", "#REF!"]
data_list = [
    [float(value) if value not in invalid_values else float("nan") for value in row]
    for row in data_list
]

# Generate the C++ code
cpp_code = """
#include <vector>

const int num_scenarios_parsed = {};
""".format((len(data_list[0])) // 4)

# Determine the number of scenarios based on the data
num_scenarios = (len(data_list[0])) // 4  # Assuming each scenario has 4 columns

print(f"Number of scenarios: {num_scenarios}")

# Find max index for the first element of each scenario
max_indices = []
for scenario in range(1, num_scenarios + 1):
    max_value = -math.inf
    max_index = -1
    for i, entry in enumerate(data_list):
        index = (scenario - 1) * 4
        if index >= len(entry):
            break
        try:
            value = float(entry[index])
        except ValueError:
            continue
        if math.isnan(value):
            continue
        if value > max_value:
            max_value = value
            max_index = i
    if max_index == -1:
        raise ValueError(f"Max value not found in scenario {scenario}")
    max_indices.append(max_index)

for scenario in range(1, num_scenarios + 1):
    i = 0
    before_apogee_entries = []
    after_apogee_entries = []
    apogee_reached = False
    for entry in data_list:
        index = (scenario - 1) * 4
        if index + 3 >= len(entry):
            break
        # Convert entries to floats and check for NaN values
        try:
            values = [float(entry[index + j]) for j in range(4)]
        except ValueError:
            continue
        if any(math.isnan(value) for value in values):
            break
        if not apogee_reached:
            before_apogee_entries.append(
                f"    {{ {values[0]}, {values[1]}, {values[2]}, {values[3]}, {max_indices[scenario - 1] - i} }}"
            )
        else:
            after_apogee_entries.append(
                f"    {{ {values[0]}, {values[1]}, {values[2]}, {values[3]}, {max_indices[scenario - 1] - i} }}"
            )
        if max_indices[scenario - 1] - i == 0:
            apogee_reached = True
        i += 1
    cpp_code += (
        f"\nstd::vector<std::vector<float>> beforeApogeeSim{scenario} = {{\n"
        + ",\n".join(before_apogee_entries)
        + "\n};\n"
    )
    cpp_code += (
        f"\nstd::vector<std::vector<float>> afterApogeeSim{scenario} = {{\n"
        + ",\n".join(after_apogee_entries)
        + "\n};\n"
    )

# Write the generated C++ code to a file
with open("Data.cpp", "w") as file:
    file.write(cpp_code)

print("Data.cpp file has been generated.")
