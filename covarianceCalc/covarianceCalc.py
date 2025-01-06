import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# import data from csv
data = pd.read_csv("covarianceCalc/covariance_Cal_CSV.csv")
Everest_data = pd.read_csv("testSuite/results/HALO.txt")

# ignore first row of Everest
Everest_data = Everest_data.iloc[1:]


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

        print(time_interval)

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


# Inputs and graphs
def graph():
    plt.figure(figsize=(10, 6))

    # Special case for the Altimeter data
    # altimeter altitude, velocity, and acceleration
    alt_Alt = data["altitude_alt"]
    velo_Alt = data["velo_alt"]
    acc_Alt = data["acceleration_alt"]
    time_Alt = data["new_time_alt"]

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

    for i in range(5, 9):
        label_Alt = "alt_" + str(i)
        label_Velo = "velo_" + str(i)
        label_Acc = "acc_" + str(i)
        label_Time = "time_" + str(i)

        alt_i = data[label_Alt].dropna()
        velo_i = data[label_Velo].dropna()
        acc_i = data[label_Acc].dropna()
        time_i = data[label_Time].dropna()

        avg_Alt_i, avg_Velo_i, avg_Acc_i, time_i = calculate_averages(
            time_i, alt_i, velo_i, acc_i, start_time, interval
        )

        # get residuals for each scenario
        residuals = []
        for j in range(len(avg_Alt_i)):
            residuals.append(abs(avg_Alt_i[j] - avg_Altimeter_Alt[j]))

        # Plot each set of averages for altitude
        plt.plot(time_i, avg_Alt_i, label=f"Alt {i}")

        # show residuals
        plt.plot(time_i, residuals, label=f"Residual {i}")

    # Add labels, title, legend, and grid for altitude
    plt.xlabel("Time")
    plt.ylabel("Altitude")
    plt.title("Time vs Altitude for Multiple Scenarios")
    plt.legend()
    plt.grid(True)
    plt.show()

    # Plot velocity
    plt.figure(figsize=(10, 6))
    for i in range(5, 9):
        label_Velo = "velo_" + str(i)
        label_Time = "time_" + str(i)

        velo_i = data[label_Velo].dropna()
        time_i = data[label_Time].dropna()

        avg_Alt_i, avg_Velo_i, avg_Acc_i, time_i = calculate_averages(
            time_i, alt_i, velo_i, acc_i, start_time, interval
        )

        # get residuals for each scenario
        residuals_Velo = []
        for j in range(len(avg_Velo_i)):
            # actual - predicted
            residuals_Velo.append(avg_Altimeter_Velo[j] - avg_Velo_i[j])

        # Plot each set of averages for velocity
        plt.plot(time_i, avg_Velo_i, label=f"Velo {i}")

        # show residuals
        plt.plot(time_i, residuals, label=f"Residual {i}")

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
    for i in range(5, 9):
        label_Acc = "acc_" + str(i)
        label_Time = "time_" + str(i)

        acc_i = data[label_Acc].dropna()
        time_i = data[label_Time].dropna()

        avg_Alt_i, avg_Velo_i, avg_Acc_i, time_i = calculate_averages(
            time_i, alt_i, velo_i, acc_i, start_time, interval
        )

        # Plot each set of averages for acceleration
        plt.plot(time_i, avg_Acc_i, label=f"Acc {i}")

    # Plot the Altimeter data for acceleration
    plt.plot(altimeter_Time, avg_Altimeter_Acc, label="Altimeter Acc", linestyle="--")

    # Add labels, title, legend, and grid for acceleration
    plt.xlabel("Time")
    plt.ylabel("Acceleration")
    plt.title("Time vs Acceleration for Multiple Altimeters")
    plt.legend()
    plt.grid(True)
    plt.show()

    # Calculate and plot covariance matrix
    combined_data = np.vstack(
        (avg_Altimeter_Alt, avg_Altimeter_Velo, avg_Altimeter_Acc)
    )
    covariance_matrix = np.cov(combined_data)

    print("Covariance Matrix of Altitude, Velocity, and Acceleration:")
    print(covariance_matrix)

    plt.figure(figsize=(8, 6))
    plt.imshow(covariance_matrix, cmap="hot", interpolation="nearest")
    plt.colorbar()
    plt.title("Covariance Matrix of Altitude, Velocity, and Acceleration")
    plt.xticks([0, 1, 2], ["Altitude", "Velocity", "Acceleration"])
    plt.yticks([0, 1, 2], ["Altitude", "Velocity", "Acceleration"])
    plt.show()


def everest_Residual():
    # Special case of Altimeter averaging (every 0.33333 seconds)
    alt_Alt = data["altitude_alt"]
    velo_Alt = data["velo_alt"]
    acc_Alt = data["acceleration_alt"]
    time_Alt = data["new_time_alt"]

    # drop NaN values
    alt_Alt = alt_Alt.dropna()
    velo_Alt = velo_Alt.dropna()
    acc_Alt = acc_Alt.dropna()
    time_Alt = time_Alt.dropna()

    # variables
    start_time = 0
    interval = 1 / 3
    avg_Altimeter_Alt, avg_Altimeter_Velo, avg_Altimeter_Acc, altimeter_Time = (
        calculate_averages(alt_Alt, velo_Alt, acc_Alt, time_Alt, start_time, interval)
    )

    # Everest data
    alt_Everest = Everest_data["Everest_Alt"]
    velo_Everest = Everest_data["Everest_Velo"]
    acc_Everest = Everest_data["Everest_Accel"]
    time_Everest = Everest_data["Time"]

    # drop NaN values
    alt_Everest = alt_Everest.dropna()
    velo_Everest = velo_Everest.dropna()
    acc_Everest = acc_Everest.dropna()
    time_Everest = time_Everest.dropna()

    # substract 2 from the time to match the time of the altimeter
    time_Everest = time_Everest - (2 + 1 / 3)

    # plot
    plt.figure(figsize=(10, 6))
    plt.plot(altimeter_Time, avg_Altimeter_Alt, label="Altimeter")
    plt.plot(time_Everest, alt_Everest, label="Everest")
    plt.xlabel("Time")
    plt.ylabel("Altitude")
    plt.title("Altimeter vs Everest Altitude")
    plt.legend()
    plt.grid(True)

    # get residuals for Everest
    residuals_Alt = []
    i = 0
    for i in range(min(len(avg_Altimeter_Alt), len(alt_Everest))):
        residuals_Alt.append(avg_Altimeter_Alt[i] - alt_Everest.iloc[i])

    # plot residuals
    plt.plot(time_Everest, residuals_Alt, label="Residual")

    # plot velo
    plt.figure(figsize=(10, 6))
    plt.plot(altimeter_Time, avg_Altimeter_Velo, label="Altimeter")
    plt.plot(time_Everest, velo_Everest, label="Everest")
    plt.xlabel("Time")
    plt.ylabel("Velocity")
    plt.title("Altimeter vs Everest Velocity")
    plt.legend()
    plt.grid(True)

    # get residuals for Everest
    residuals_Velo = []
    for i in range(min(len(avg_Altimeter_Velo), len(velo_Everest))):
        residuals_Velo.append(avg_Altimeter_Velo[i] - velo_Everest.iloc[i])

    # plot residuals
    plt.plot(time_Everest, residuals_Velo, label="Residual")

    # plot accel
    plt.figure(figsize=(10, 6))
    plt.plot(altimeter_Time, avg_Altimeter_Acc, label="Altimeter")
    plt.plot(time_Everest, acc_Everest, label="Everest")
    plt.xlabel("Time")
    plt.ylabel("Acceleration")
    plt.title("Altimeter vs Everest Acceleration")
    plt.legend()
    plt.grid(True)

    # get residuals for Everest
    residuals_Acc = []
    for i in range(min(len(avg_Altimeter_Acc), len(acc_Everest))):
        residuals_Acc.append(avg_Altimeter_Acc[i] - acc_Everest.iloc[i])

    # plot residuals
    plt.plot(time_Everest, residuals_Acc, label="Residual")

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


# graph()
everest_Residual()
