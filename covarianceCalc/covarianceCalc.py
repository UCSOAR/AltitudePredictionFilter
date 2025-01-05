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

# for every 0.1 second time interval average velocity, altitude, and acceleration
avg_Alt = []
avg_Velo = []
avg_Acc = []
time = []

# for every 0.1 second time interval average velocity, altitude, and acceleration
for i in range(0, len(alt_Alt), 10):
    avg_Alt.append(alt_Alt[i : i + 10].mean())
    avg_Velo.append(velo_Alt[i : i + 10].mean())
    avg_Acc.append(acc_Alt[i : i + 10].mean())
    time.append(time_Alt[i])

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
