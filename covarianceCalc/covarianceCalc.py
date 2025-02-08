import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re

# import data from csv
Everest_data = pd.read_csv("testSuite/results/HALO.txt")
Quasar_and_sims_data = pd.read_csv("covarianceCalc/covariance_Cal_Data.csv")
HALO_data = pd.read_csv("testSuite/results/HALO.txt")
p_data = pd.read_csv("testSuite/results/P.txt")
confidence = pd.read_csv("testSuite/results/confidence.txt")
altimeter_data = pd.read_csv("testSuite/data/altimeter1.csv")


def calculate_averages(
    time_input, alt_input, velo_input, acc_input, start_time, interval
):
    # for every 0.1 second time interval average velocity, altitude, and acceleration
    avg_Alt = []
    avg_Velo = []
    avg_Acc = []
    time = []

    # for every time thats within 0.1 interval average velocity, altitude, and acceleration
    lower_bound = start_time
    time_interval = start_time + interval
    index = 0

    max_time = time_input.max()
    while time_interval <= max_time:
        sum_Alt = 0
        sum_Velo = 0
        sum_Acc = 0
        counter_Current_Range_Samples = 0
        i = index

        while i < len(time_input):
            if time_input[i] <= time_interval and time_input[i] >= lower_bound:
                if i < len(alt_input):
                    sum_Alt += alt_input.iloc[i]
                if i < len(velo_input):
                    sum_Velo += velo_input.iloc[i]
                if i < len(acc_input):
                    sum_Acc += acc_input.iloc[i]
                counter_Current_Range_Samples += 1
            i += 1

        lower_bound = time_interval
        time_interval += interval

        if counter_Current_Range_Samples == 0:
            continue
        avg_Alt.append(sum_Alt / counter_Current_Range_Samples)
        avg_Velo.append(sum_Velo / counter_Current_Range_Samples)
        avg_Acc.append(sum_Acc / counter_Current_Range_Samples)
        time.append(time_interval)

    return avg_Alt, avg_Velo, avg_Acc, time


start_time = 0
interval = 0.1

residuals_alt = {}
residuals_velo = {}
residuals_acc = {}


# Inputs and graphs
def sims_Residual():
    plt.figure(figsize=(10, 6))

    # Convert DataFrame to list of lists
    data_list = Quasar_and_sims_data.values.tolist()

    # Replace invalid strings with NaN and convert to floats
    invalid_values = ["", "#DIV/0!", "#REF!"]
    data_list = [
        [float(value) if value not in invalid_values else float("nan") for value in row]
        for row in data_list
    ]

    # Determine the number of scenarios using regex
    alt_columns = [
        col for col in Quasar_and_sims_data.columns if re.match(r"alt_\d+", col)
    ]
    num_scenarios = len(alt_columns)

    print(f"Number of scenarios: {num_scenarios}")

    # remove first row
    data_list = data_list[1:]

    # Initialize residuals dictionary
    residuals = {f"sim{scenario}": [] for scenario in range(1, num_scenarios + 1)}

    # Special case for the Altimeter data
    # altimeter altitude, velocity, and acceleration
    alt_Alt = Quasar_and_sims_data["altitude_alt"]
    velo_Alt = Quasar_and_sims_data["velo_alt"]
    acc_Alt = Quasar_and_sims_data["acceleration_alt"]
    time_Alt = Quasar_and_sims_data["new_time_alt"]

    # drop NaN values
    alt_Alt = alt_Alt.dropna()
    velo_Alt = velo_Alt.dropna()
    acc_Alt = acc_Alt.dropna()
    time_Alt = time_Alt.dropna()

    avg_Altimeter_Alt, avg_Altimeter_Velo, avg_Altimeter_Acc, altimeter_Time = (
        calculate_averages(time_Alt, alt_Alt, velo_Alt, acc_Alt, start_time, interval)
    )

    # Plot the Altimeter data for altitude
    plt.plot(altimeter_Time, avg_Altimeter_Alt, label="Altimeter", linestyle="--")

    for scenario in range(1, num_scenarios + 1):
        label_Alt = "alt_" + str(scenario)
        label_Velo = "velo_" + str(scenario)
        label_Acc = "acc_" + str(scenario)
        label_Time = "time_" + str(scenario)

        alt_i = Quasar_and_sims_data[label_Alt].dropna()
        velo_i = Quasar_and_sims_data[label_Velo].dropna()
        acc_i = Quasar_and_sims_data[label_Acc].dropna()
        time_i = Quasar_and_sims_data[label_Time].dropna()

        avg_Alt_i, avg_Velo_i, avg_Acc_i, time_i = calculate_averages(
            time_i, alt_i, velo_i, acc_i, start_time, interval
        )

        residuals_alt[f"sim_time_{scenario}"] = time_i
        residuals_alt[f"sim_avg_alt_{scenario}"] = avg_Alt_i

        # get residuals for each scenario
        residuals = []
        for j in range(len(avg_Alt_i)):
            residuals.append(avg_Altimeter_Alt[j] - avg_Alt_i[j])

        # calculate covariance matrix of residuals
        residuals_alt[f"sim_residual{scenario}"] = residuals

        # Plot each set of averages for altitude
        plt.plot(
            residuals_alt[f"sim_time_{scenario}"],
            residuals_alt[f"sim_avg_alt_{scenario}"],
            label=f"Alt {scenario}",
        )

        # show residuals
        plt.plot(
            residuals_alt[f"sim_time_{scenario}"],
            residuals_alt[f"sim_residual{scenario}"],
            label=f"Residual {scenario}",
        )

    # Add labels, title, legend, and grid for altitude
    plt.xlabel("Time")
    plt.ylabel("Altitude")
    plt.title("Time vs Altitude for Multiple Scenarios")
    plt.legend()
    plt.grid(True)
    plt.show()

    # Plot velocity
    plt.figure(figsize=(10, 6))
    for scenario in range(1, num_scenarios + 1):
        label_Velo = "velo_" + str(scenario)
        label_Time = "time_" + str(scenario)

        velo_i = Quasar_and_sims_data[label_Velo].dropna()
        time_i = Quasar_and_sims_data[label_Time].dropna()

        avg_Alt_i, avg_Velo_i, avg_Acc_i, time_i = calculate_averages(
            time_i, alt_i, velo_i, acc_i, start_time, interval
        )

        # get residuals for each scenario
        residuals_Velo = []
        for j in range(len(avg_Velo_i)):
            # actual - predicted
            residuals_Velo.append(avg_Altimeter_Velo[j] - avg_Velo_i[j])

        residuals_velo[f"sim_time_{scenario}"] = time_i
        residuals_velo[f"sim_avg_velo_{scenario}"] = avg_Velo_i
        residuals_velo[f"sim_residual{scenario}"] = residuals_Velo

        # Plot each set of averages for velocity
        plt.plot(
            residuals_velo[f"sim_time_{scenario}"],
            residuals_velo[f"sim_avg_velo_{scenario}"],
            label=f"Velo {scenario}",
        )

        # plot residuals
        plt.plot(
            residuals_velo[f"sim_time_{scenario}"],
            residuals_velo[f"sim_residual{scenario}"],
            label=f"Residual {scenario}",
        )

    # Plot the Altimeter data for velocity
    plt.plot(altimeter_Time, avg_Altimeter_Velo, label="Altimeter Velo", linestyle="--")

    # Add labels, title, legend, and grid for velocity
    plt.xlabel("Time")
    plt.ylabel("Velocity")
    plt.title("Time vs Velocity for Multiple Altimeters")
    plt.legend()
    plt.grid(True)
    plt.show()

    # Plot acceleration
    plt.figure(figsize=(10, 6))
    for scenario in range(1, num_scenarios + 1):
        label_Acc = "acc_" + str(scenario)
        label_Time = "time_" + str(scenario)

        acc_i = Quasar_and_sims_data[label_Acc].dropna()
        time_i = Quasar_and_sims_data[label_Time].dropna()

        avg_Alt_i, avg_Velo_i, avg_Acc_i, time_i = calculate_averages(
            time_i, alt_i, velo_i, acc_i, start_time, interval
        )

        # get residuals for each scenario
        residuals_Acc = []
        for j in range(len(avg_Acc_i)):
            residuals_Acc.append(avg_Altimeter_Acc[j] - avg_Acc_i[j])

        residuals_acc[f"sim_time_{scenario}"] = time_i
        residuals_acc[f"sim_avg_acc_{scenario}"] = avg_Acc_i
        residuals_acc[f"sim_residual{scenario}"] = residuals_Acc

        # plot residuals
        plt.plot(
            residuals_acc[f"sim_time_{scenario}"],
            residuals_acc[f"sim_residual{scenario}"],
            label=f"Residual {scenario}",
        )

        # Plot each set of averages for acceleration
        plt.plot(
            residuals_acc[f"sim_time_{scenario}"],
            residuals_acc[f"sim_avg_acc_{scenario}"],
            label=f"Acc {scenario}",
        )

    # Plot the Altimeter data for acceleration
    plt.plot(altimeter_Time, avg_Altimeter_Acc, label="Altimeter Acc", linestyle="--")

    # Add labels, title, legend, and grid for acceleration
    plt.xlabel("Time")
    plt.ylabel("Acceleration")
    plt.title("Time vs Acceleration for Multiple Altimeters")
    plt.legend()
    plt.grid(True)
    plt.show()

    scenarios_combined_all = {}
    # Initialize covariance_matrix
    covariance_matrix = np.zeros((3, 3))

    # Combine residuals for all scenarios
    for scenario in range(1, num_scenarios + 1):
        scenarios_combined_all[f"sim_residual{scenario}"] = np.vstack(
            (
                residuals_alt[f"sim_residual{scenario}"],
                residuals_velo[f"sim_residual{scenario}"],
                residuals_acc[f"sim_residual{scenario}"],
            )
        )

    # add all covariance matrices
    for scenario in range(1, num_scenarios + 1):
        covariance_matrix += np.cov(scenarios_combined_all[f"sim_residual{scenario}"])

    # divide by number of scenarios
    covariance_matrix = covariance_matrix / num_scenarios

    print("Covariance Matrix of Altitude, Velocity, and Acceleration Residuals:")
    print(covariance_matrix)

    # plot covariance matrix
    plt.figure(figsize=(8, 6))
    plt.imshow(covariance_matrix, cmap="hot", interpolation="nearest")
    plt.colorbar()
    plt.title("Covariance Matrix of Altitude, Velocity, and Acceleration Residuals")
    plt.xticks([0, 1, 2], ["Altitude", "Velocity", "Acceleration"])
    plt.yticks([0, 1, 2], ["Altitude", "Velocity", "Acceleration"])
    plt.show()


def everest_Residual():
    # HALO data
    alt_HALO = HALO_data["Halo_Alt"][1:]
    velo_HALO = HALO_data["Halo_Velo"][1:]
    acc_HALO = HALO_data["Halo_Accel"][1:]
    time_HALO = HALO_data["Time"][1:]
    halo_alt_np = np.array(alt_HALO)
    halo_velo_np = np.array(velo_HALO)
    halo_acc_np = np.array(acc_HALO)

    # Ensure both arrays have the same shape
    min_length = min(len(halo_alt_np), len(p_data["alt_std"]))
    halo_alt_np = halo_alt_np[:min_length]
    halo_velo_np = halo_velo_np[:min_length]
    halo_acc_np = halo_acc_np[:min_length]

    alt_std_np = np.array(p_data["alt_std"])[:min_length]
    velo_std_np = np.array(p_data["velo_std"])[:min_length]
    acc_std_np = np.array(p_data["acc_std"])[:min_length]

    # get confidence interval
    upper_alt = halo_alt_np + np.sqrt(np.abs(alt_std_np))
    lower_alt = halo_alt_np - np.sqrt(np.abs(alt_std_np))

    upper_velo = halo_velo_np + np.sqrt(np.abs(velo_std_np))
    lower_velo = halo_velo_np - np.sqrt(np.abs(velo_std_np))

    upper_acc = halo_acc_np + np.sqrt(np.abs(acc_std_np))
    lower_acc = halo_acc_np - np.sqrt(np.abs(acc_std_np))

    # substract 2.333 from the time to match the time of the altimeter
    time_HALO = time_HALO - (1 + 2 / 3)
    # time_HALO = time_HALO[1:]

    # Special case of Altimeter averaging (every 0.33333 seconds)
    alt_Alt = Quasar_and_sims_data["alt_Q"]
    velo_Alt = Quasar_and_sims_data["velo_Q"]
    acc_Alt = Quasar_and_sims_data["acc_Q"]
    time_Alt = Quasar_and_sims_data["time_Q"]

    # drop NaN values
    alt_Alt = alt_Alt.dropna()
    velo_Alt = velo_Alt.dropna()
    acc_Alt = acc_Alt.dropna()
    time_Alt = time_Alt.dropna()

    # variables
    start_time = 0
    interval = 1 / 3
    avg_Altimeter_Alt, avg_Altimeter_Velo, avg_Altimeter_Acc, altimeter_Time = (
        calculate_averages(time_Alt, alt_Alt, velo_Alt, acc_Alt, start_time, interval)
    )

    # Everest data
    alt_Everest = Everest_data["Everest_Alt"]
    velo_Everest = Everest_data["Everest_Velo"]
    acc_Everest = Everest_data["Everest_Accel"]
    time_Everest = Everest_data["Time"]

    # remove first row
    alt_Everest = alt_Everest[1:].reset_index(drop=True)
    velo_Everest = velo_Everest[1:].reset_index(drop=True)
    acc_Everest = acc_Everest[1:].reset_index(drop=True)
    time_Everest = time_Everest[1:].reset_index(drop=True)

    # Everest Avg
    avg_Alt_Everest, avg_Velo_Everest, avg_Acc_Everest, time_Everest = (
        calculate_averages(
            time_Everest, alt_Everest, velo_Everest, acc_Everest, start_time, interval
        )
    )

    # drop NaN values
    alt_Everest = np.array(avg_Alt_Everest)
    velo_Everest = np.array(avg_Velo_Everest)
    acc_Everest = np.array(avg_Acc_Everest)
    time_Everest = np.array(time_Everest)

    # cast to np array and drop NaN values
    alt_Everest = alt_Everest[~np.isnan(alt_Everest)]
    velo_Everest = velo_Everest[~np.isnan(velo_Everest)]
    acc_Everest = acc_Everest[~np.isnan(acc_Everest)]
    time_Everest = time_Everest[~np.isnan(time_Everest)]

    # substract 2.333 from the time to match the time of the altimeter
    time_Everest = time_Everest - (1 + 2 / 3)

    # plot
    plt.figure(figsize=(10, 6))
    # plot altimeter
    altimeter_uncut = altimeter_data["altitude"]
    altimeter_Time_uncut = altimeter_data["time"]
    plt.plot(altimeter_Time_uncut, altimeter_uncut, label="Altimeter")
    plt.plot(time_Everest, alt_Everest, label="Everest")
    plt.plot(time_HALO, alt_HALO, label="HALO")
    plt.fill_between(
        time_HALO[: len(upper_alt)],
        lower_alt,
        upper_alt,
        color="gray",
        alpha=0.5,
        label="Confidence Interval",
    )
    plt.xlabel("Time")
    plt.ylabel("Altitude")
    plt.title("Altimeter vs Everest Altitude")
    plt.legend()
    plt.grid(True)

    # get residuals for Everest
    residuals_Alt = []
    for i in range(min(len(avg_Altimeter_Alt), len(alt_Everest))):
        residual = avg_Altimeter_Alt[i] - alt_Everest[i]
        residuals_Alt.append(residual)

    # Ensure both arrays have the same shape for plotting residuals
    min_length_residuals = min(len(time_Everest), len(residuals_Alt))
    time_Everest_residuals_cut = time_Everest[:min_length_residuals]
    residuals_Alt_cut = residuals_Alt[:min_length_residuals]

    # Ensure both arrays have the same shape for velocity
    min_length_velo = min(len(time_Everest), len(velo_Everest))
    time_Everest_velo_cut = time_Everest[:min_length_velo]
    residuals_Velo_cut = velo_Everest[:min_length_velo]

    # plot residuals
    plt.plot(time_Everest_residuals_cut, residuals_Alt_cut, label="Residual")

    # plot velo
    plt.figure(figsize=(10, 6))
    plt.plot(altimeter_Time, avg_Altimeter_Velo, label="Altimeter")
    plt.plot(time_Everest_velo_cut, residuals_Velo_cut, label="Everest")
    plt.plot(time_HALO, velo_HALO, label="HALO")
    plt.fill_between(
        time_HALO[: len(upper_alt)],
        lower_velo,
        upper_velo,
        color="gray",
        alpha=0.5,
        label="Confidence Interval",
    )
    plt.xlabel("Time")
    plt.ylabel("Velocity")
    plt.title("Altimeter vs Everest Velocity")
    plt.legend()
    plt.grid(True)

    # get residuals for Everest
    residuals_Velo = []
    for i in range(min(len(avg_Altimeter_Velo), len(velo_Everest))):
        residuals_Velo.append(avg_Altimeter_Velo[i] - velo_Everest[i])

    # plot residuals
    plt.plot(time_Everest_velo_cut, residuals_Velo_cut, label="Residual")

    # Ensure both arrays have the same shape for acceleration
    min_length_acc = min(len(time_Everest), len(acc_Everest))
    time_Everest_cut = time_Everest[:min_length_acc]
    acc_Everest_cut = acc_Everest[:min_length_acc]

    # plot accel
    plt.figure(figsize=(10, 6))
    plt.plot(altimeter_Time, avg_Altimeter_Acc, label="Altimeter")
    plt.plot(time_Everest_cut, acc_Everest_cut, label="Everest")
    plt.plot(time_HALO, acc_HALO, label="HALO")
    plt.fill_between(
        time_HALO[: len(upper_alt)],
        lower_acc,
        upper_acc,
        color="gray",
        alpha=0.5,
        label="Confidence Interval",
    )
    plt.xlabel("Time")
    plt.ylabel("Acceleration")
    plt.title("Altimeter vs Everest Acceleration")
    plt.legend()
    plt.grid(True)

    # get residuals for Everest
    residuals_Acc = []
    for i in range(min(len(avg_Altimeter_Acc), len(acc_Everest))):
        residuals_Acc.append(avg_Altimeter_Acc[i] - acc_Everest[i])

    residuals_min_length = min(len(time_Everest_cut), len(residuals_Acc))
    time_Everest_cut = time_Everest[:residuals_min_length]
    residuals_Acc_cut = residuals_Acc[:residuals_min_length]

    # plot residuals
    plt.plot(time_Everest_cut, residuals_Acc_cut, label="Residual")

    # plot confidence values
    plt.figure(figsize=(10, 6))
    plt.plot(
        time_HALO[: len(confidence)],
        confidence.iloc[:, 0],
        label="Altitude Confidence",
        color="red",
    )
    plt.plot(
        time_HALO[: len(confidence)],
        confidence.iloc[:, 3],
        label="TotalConfidence",
        color="blue",
    )
    plt.plot(
        time_HALO[: len(confidence)],
        confidence.iloc[:, 1],
        label="Velocity Confidence",
        color="green",
    )
    plt.plot(
        time_HALO[: len(confidence)],
        confidence.iloc[:, 2],
        label="Acceleration Confidence",
        color="purple",
    )
    plt.yticks(np.arange(0, 1.1, step=0.1))
    plt.xlabel("Time")
    plt.ylabel("Confidence")
    plt.title("Confidence vs Time")
    plt.legend()
    plt.grid(True)

    plt.show()

    # Calculate and plot covariance matrix of residuals
    combined_data = np.vstack((residuals_Alt, residuals_Velo, residuals_Acc))
    covariance_matrix = np.cov(combined_data)

    print("Covariance Matrix of Altitude, Velocity, and Acceleration Residuals:")
    print(covariance_matrix)

    plt.figure(figsize=(8, 6))
    plt.imshow(covariance_matrix, cmap="hot", interpolation="nearest")
    plt.colorbar()
    plt.title("Covariance Matrix of Altitude, Velocity, and Acceleration Residuals")
    plt.xticks([0, 1, 2], ["Altitude", "Velocity", "Acceleration"])
    plt.yticks([0, 1, 2], ["Altitude", "Velocity", "Acceleration"])
    plt.show()


# sims_Residual()
everest_Residual()
