import pandas as pd
import matplotlib.pyplot as plt

# import data from csv
data = pd.read_csv("covarianceCalc/covariance_Cal.csv")

print(data)

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

print(time_Alt)

# for every 0.1 second time interval average velocity, altitude, and acceleration
avg_Alt = []
avg_Velo = []
avg_Acc = []
time = []

# for every time thats within 0.1 interval average velocity, altitude, and acceleration
lower_bound = 0
time_interval = 0.1
index = 0

print(len(time_Alt))

# iterates through every 0.1 second interval
while time_interval <= 34.4:
    sum_Alt = 0
    sum_Velo = 0
    sum_Acc = 0
    counter_Current_Range_Samples = 0
    i = index
    j = index

    while j < len(time_Alt):
        print(j)

        print(len(time_Alt))

        print(time_Alt[j] <= time_interval)
        print(time_Alt[j] >= lower_bound)

        print("Time: %f" % time_Alt[j])
        print("lower_bound: %f" % lower_bound)
        print("time_interval: %f" % time_interval)

        while time_Alt[j] <= time_interval and time_Alt[j] >= lower_bound:
            sum_Alt += alt_Alt[j]
            sum_Velo += velo_Alt[j]
            sum_Acc += acc_Alt[j]
            counter_Current_Range_Samples += 1
            print("Time: %f" % time_Alt[j])
            print("Alt: %f" % alt_Alt[j])
            print("Velo: %f" % velo_Alt[j])
            print("Acc: %f" % acc_Alt[j])

            print("j inside    : %d" % j)

            j += 1
            index += 1

        break

    # print("Sum_Alt: %f" % sum_Alt)
    # print("Sum_Velo: %f" % sum_Velo)
    # print("Sum_Acc: %f" % sum_Acc)

    # print("In time range (%f, %f)" % (lower_bound, time_interval))
    # print("Counter: %d" % counter_Current_Range_Samples)
    lower_bound = time_interval
    time_interval += 0.1

    if counter_Current_Range_Samples == 0:
        continue
    avg_Alt.append(sum_Alt / counter_Current_Range_Samples)
    avg_Velo.append(sum_Velo / counter_Current_Range_Samples)
    avg_Acc.append(sum_Acc / counter_Current_Range_Samples)
    time.append(time_interval)


# plot altimeter altitude, velocity, and acceleration
plt.figure(figsize=(10, 6))
plt.plot(time, avg_Alt, label="Altimeter Alt")
plt.xlabel("Time")
plt.ylabel("Altitude")
plt.title("Time vs Alt")
plt.legend()
plt.grid(True)

# plot velo
plt.figure(figsize=(10, 6))
plt.plot(time, avg_Velo, label="Altimeter speed")
plt.xlabel("Time")
plt.ylabel("speed")
plt.title("Time vs speed")
plt.legend()

plt.grid(True)

# plot accel
plt.figure(figsize=(10, 6))
plt.plot(time, avg_Acc, label="Altimeter Acc")
plt.xlabel("Time")
plt.ylabel("Acc")
plt.title("Time vs Acc")
plt.legend()
plt.grid(True)

plt.show()
