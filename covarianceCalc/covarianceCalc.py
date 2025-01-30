import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# import data from csv
data = pd.read_csv("covarianceCalc/covariance_Cal_CSV.csv")
Everest_data = pd.read_csv("testSuite/results/HALO.txt")
Quasar_data = pd.read_csv("covarianceCalc/covariance_Cal_New_CSV.csv")
HALO_data = pd.read_csv("testSuite/results/HALO.txt")
p_data = pd.read_csv("testSuite/results/P.txt")
confidence = pd.read_csv("testSuite/results/confidence.txt")

# ignore first row of Everest
Everest_data = Everest_data.iloc[1:]
HALO_data = HALO_data.iloc[1:]


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

residual_alt_5 = []
residual_velo_5 = []
residual_acc_5 = []

residual_alt_6 = []
residual_velo_6 = []
residual_acc_6 = []

residual_alt_7 = []
residual_velo_7 = []
residual_acc_7 = []

residual_alt_8 = []
residual_velo_8 = []
residual_acc_8 = []


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
            residuals.append(avg_Altimeter_Alt[j] - avg_Alt_i[j])

        # calculate covariance matrix of residuals
        if i == 5:
            residuals_alt_5 = residuals
        elif i == 6:
            residuals_alt_6 = residuals
        elif i == 7:
            residuals_alt_7 = residuals
        elif i == 8:
            residuals_alt_8 = residuals

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

        if i == 5:
            residuals_velo_5 = residuals_Velo
        elif i == 6:
            residuals_velo_6 = residuals_Velo
        elif i == 7:
            residuals_velo_7 = residuals_Velo
        elif i == 8:
            residuals_velo_8 = residuals_Velo

        # Plot each set of averages for velocity
        plt.plot(time_i, avg_Velo_i, label=f"Velo {i}")

        # plot residuals
        plt.plot(time_i, residuals_Velo, label=f"Residual {i}")

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

        # get residuals for each scenario
        residuals_Acc = []
        for j in range(len(avg_Acc_i)):
            residuals_Acc.append(avg_Altimeter_Acc[j] - avg_Acc_i[j])

        # plot residuals
        plt.plot(time_i, residuals_Acc, label=f"Residual {i}")

        if i == 5:
            residuals_acc_5 = residuals_Acc
        elif i == 6:
            residuals_acc_6 = residuals_Acc
        elif i == 7:
            residuals_acc_7 = residuals_Acc
        elif i == 8:
            residuals_acc_8 = residuals_Acc

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

    # Individually calculate covariance matrix of residuals for each scenario
    combined_data_5 = np.vstack((residuals_alt_5, residuals_velo_5, residuals_acc_5))
    covariance_matrix_5 = np.cov(combined_data_5)

    combined_data_6 = np.vstack((residuals_alt_6, residuals_velo_6, residuals_acc_6))
    covariance_matrix_6 = np.cov(combined_data_6)

    combined_data_7 = np.vstack((residuals_alt_7, residuals_velo_7, residuals_acc_7))
    covariance_matrix_7 = np.cov(combined_data_7)

    combined_data_8 = np.vstack((residuals_alt_8, residuals_velo_8, residuals_acc_8))
    covariance_matrix_8 = np.cov(combined_data_8)

    # add all covariance matrices
    covariance_matrix = (
        covariance_matrix_5
        + covariance_matrix_6
        + covariance_matrix_7
        + covariance_matrix_8
    )
    covariance_matrix = covariance_matrix / 4

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
    alt_HALO = HALO_data["Halo_Alt"]
    velo_HALO = HALO_data["Halo_Velo"]
    acc_HALO = HALO_data["Halo_Accel"]
    time_HALO = HALO_data["Time"]

    halo_alt_np = np.array(alt_HALO)
    halo_velo_np = np.array(velo_HALO)
    halo_acc_np = np.array(acc_HALO)

    # get confidence interval
    upper_alt = halo_alt_np[0:] + np.sqrt(np.abs(p_data["alt_std"]))
    lower_alt = halo_alt_np[0:] - np.sqrt(np.abs(p_data["alt_std"]))

    upper_velo = halo_velo_np[0:] + np.sqrt(np.abs(p_data["velo_std"]))
    lower_velo = halo_velo_np[0:] - np.sqrt(np.abs(p_data["velo_std"]))

    upper_acc = halo_acc_np[0:] + np.sqrt(np.abs(p_data["acc_std"]))
    lower_acc = halo_acc_np[0:] - np.sqrt(np.abs(p_data["acc_std"]))

    print("Acc std")
    print(np.sqrt(np.abs(p_data["acc_std"])))

    # cut at 89 for apogee
    # alt_HALO = alt_HALO[1:89]
    # velo_HALO = velo_HALO[1:89]
    # acc_HALO = acc_HALO[1:89]
    # time_HALO = time_HALO[1:89]

    # substract 2.333 from the time to match the time of the altimeter
    time_HALO = time_HALO - (1 + 2 / 3)
    time_HALO = time_HALO[1:]

    # Special case of Altimeter averaging (every 0.33333 seconds)
    alt_Alt = Quasar_data["alt_Q"]
    velo_Alt = Quasar_data["velo_Q"]
    acc_Alt = Quasar_data["acc_Q"]
    time_Alt = Quasar_data["time_Q"]

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

    # cutoff at Time = 27.333334
    alt_Everest = alt_Everest[:79]
    velo_Everest = velo_Everest[:79]
    acc_Everest = acc_Everest[:79]
    time_Everest = time_Everest[:79]

    # drop NaN values
    alt_Everest = alt_Everest.dropna()
    velo_Everest = velo_Everest.dropna()
    acc_Everest = acc_Everest.dropna()
    time_Everest = time_Everest.dropna()

    # substract 2.333 from the time to match the time of the altimeter
    time_Everest = time_Everest - (1 + 2 / 3)

    # plot
    plt.figure(figsize=(10, 6))
    plt.plot(altimeter_Time, avg_Altimeter_Alt, label="Altimeter")
    plt.plot(time_Everest, alt_Everest, label="Everest")
    plt.plot(time_HALO, alt_HALO[1:], label="HALO")
    plt.fill_between(
        time_HALO,
        lower_alt[1:],
        upper_alt[1:],
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
    i = 0
    for i in range(min(len(avg_Altimeter_Alt), len(alt_Everest))):
        residuals_Alt.append(avg_Altimeter_Alt[i] - alt_Everest.iloc[i])

    # plot residuals
    plt.plot(time_Everest, residuals_Alt, label="Residual")

    # plot velo
    plt.figure(figsize=(10, 6))
    plt.plot(altimeter_Time, avg_Altimeter_Velo, label="Altimeter")
    plt.plot(time_Everest, velo_Everest, label="Everest")
    plt.plot(time_HALO, velo_HALO[1:], label="HALO")
    plt.fill_between(
        time_HALO,
        lower_velo[1:],
        upper_velo[1:],
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
        residuals_Velo.append(avg_Altimeter_Velo[i] - velo_Everest.iloc[i])

    # plot residuals
    plt.plot(time_Everest, residuals_Velo, label="Residual")

    # plot accel
    plt.figure(figsize=(10, 6))
    plt.plot(altimeter_Time, avg_Altimeter_Acc, label="Altimeter")
    plt.plot(time_Everest, acc_Everest, label="Everest")
    plt.plot(time_HALO, acc_HALO[1:], label="HALO")
    plt.fill_between(
        time_HALO,
        lower_acc[1:],
        upper_acc[1:],
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
        residuals_Acc.append(avg_Altimeter_Acc[i] - acc_Everest.iloc[i])

    # plot residuals
    plt.plot(time_Everest, residuals_Acc, label="Residual")

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


graph()
# everest_Residual()
