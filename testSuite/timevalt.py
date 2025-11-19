import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("testSuite/results/HALO.txt")
sims = pd.read_csv("testSuite/data/beforeSimsF2_Short.csv")

# plot altitude (altimeter, baro, imu, everest, halo, gps, sigmaPoints, scenarios)
plt.figure(figsize=(10, 6))
# plt.plot(df['Time'], df['Altimeter_Alt'], label='Altimeter Alt')
# plt.plot(df['Time'], df['Baro_Alt'], label='Baro Alt')
# plt.plot(df['Time'], df['IMU_Alt'], label='IMU Alt')
plt.plot(df["Time"], df["Everest_Alt"], label="Everest Alt")
plt.plot(df["Time"], df["Halo_Alt"], label="HALO Alt")
# plt.plot(df['Time'], df['GPS_Alt'], label='GPS Alt')
# plt.plot(df['Time'], df['Sigma_Alt_Upper'], label='Sigma Alt Upper', marker = "x")
# plt.plot(df['Time'], df['Sigma_Alt_Lower'], label='Sigma Alt Lower', marker = "o")
plt.plot(sims["time_1"], sims["alt_1"], label="Scenarios Alt1")
plt.plot(sims["time_2"], sims["alt_2"], label="Scenarios Alt2")
plt.plot(sims["time_3"], sims["alt_3"], label="Scenarios Alt3")
plt.plot(sims["time_4"], sims["alt_4"], label="Scenarios Alt4")
plt.plot(sims["time_5"], sims["alt_5"], label="Scenarios Alt5")
plt.plot(sims["time_6"], sims["alt_6"], color="red", label="Altimeter")

plt.xlabel("Time")
plt.ylabel("Altitude")
plt.title("Time vs Alt")
plt.legend()
plt.grid(True)

plt.show()