import pandas as pd
import matplotlib.pyplot as plt

# import data from csv
data = pd.read_csv("covarianceCalc/covariance_Cal.csv")


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


# all averages
big_List = []

start_time = 0
interval = 0.1


def get_inputs():
    plt.figure(figsize=(10, 6))

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

        # Plot each set of averages for altitude
        plt.plot(time_i, avg_Alt_i, label=f"Alt {i}")

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

        # Plot each set of averages for velocity
        plt.plot(time_i, avg_Velo_i, label=f"Velo {i}")

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


get_inputs()
