import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("testSuite/results/HALO.txt")
alt = pd.read_csv("testSuite/data/Altimeter(in)(in).csv")

sims = pd.read_csv("testSuite/data/final_may_sims_formatted.csv")

imu = pd.read_csv("testSuite/data/imu_baro.csv")
imu = pd.read_csv("testSuite/data/imu_baro.csv")

print(sims)
sims["avg_alt_1"] = pd.to_numeric(sims["avg_alt_1"], errors="coerce")
sims["avg_velo_1"] = pd.to_numeric(sims["avg_velo_1"], errors="coerce")
sims["avg_acc_1"] = pd.to_numeric(sims["avg_acc_1"], errors="coerce")
sims["avg_time_1"] = pd.to_numeric(sims["avg_time_1"], errors="coerce")

sims["avg_alt_2"] = pd.to_numeric(sims["avg_alt_2"], errors="coerce")
sims["avg_velo_2"] = pd.to_numeric(sims["avg_velo_2"], errors="coerce")
sims["avg_acc_2"] = pd.to_numeric(sims["avg_acc_2"], errors="coerce")
sims["avg_time_2"] = pd.to_numeric(sims["avg_time_2"], errors="coerce")
sims = sims.dropna()

print(sims)


# plot altitude (altimeter, baro, imu, everest, halo, gps, sigmaPoints, scenarios)
plt.figure(figsize=(10, 6))
# plt.plot(df['Time'], df['Altimeter_Alt'], label='Altimeter Alt')
# plt.plot(df['Time'], df['Baro_Alt'], label='Baro Alt')
# plt.plot(df['Time'], df['IMU_Alt'], label='IMU Alt')
plt.plot(df["Time"], df["Everest_Alt"], label="Everest Alt")
plt.plot(df["Time"], df["Halo_Alt"], label="HALO Alt")
plt.plot(alt["time"], alt["altitude"], label="Altimeter Alt")
# plt.plot(df['Time'], df['GPS_Alt'], label='GPS Alt')
# plt.plot(df['Time'], df['Sigma_Alt_Upper'], label='Sigma Alt Upper', marker = "x")
# plt.plot(df['Time'], df['Sigma_Alt_Lower'], label='Sigma Alt Lower', marker = "o")
plt.plot(sims["avg_time_1"], sims["avg_alt_1"], label="Scenarios Alt1")
plt.plot(sims["avg_time_2"], sims["avg_alt_2"], label="Scenarios Alt2")


plt.xlabel("Time")
plt.ylabel("Altitude")
plt.title("Time vs Alt")
plt.legend()
plt.grid(True)


# plot altitude (altimeter, baro, imu, everest, halo, gps, sigmaPoints, scenarios)
plt.figure(figsize=(10, 6))
# plt.plot(df['Time'], df['Altimeter_Alt'], label='Altimeter Alt')
# plt.plot(df['Time'], df['Baro_Alt'], label='Baro Alt')
# plt.plot(df['Time'], df['IMU_Alt'], label='IMU Alt')
plt.plot(df["Time"], df["Everest_Velo"], label="Everest Velo")
plt.plot(df["Time"], df["Halo_Velo"], label="HALO Velo")
# plt.plot(df['Time'], df['GPS_Alt'], label='GPS Alt')
# plt.plot(df['Time'], df['Sigma_Alt_Upper'], label='Sigma Alt Upper', marker = "x")
# plt.plot(df['Time'], df['Sigma_Alt_Lower'], label='Sigma Alt Lower', marker = "o")
plt.plot(sims["avg_time_1"], sims["avg_velo_1"], label="Scenarios Velo1")
plt.plot(sims["avg_time_2"], sims["avg_velo_2"], label="Scenarios Velo2")


plt.xlabel("Time")
plt.ylabel("Velocity")
plt.title("Time vs Velocity")
plt.legend()
plt.grid(True)


# plot altitude (altimeter, baro, imu, everest, halo, gps, sigmaPoints, scenarios)
plt.figure(figsize=(10, 6))
# plt.plot(df['Time'], df['Altimeter_Alt'], label='Altimeter Alt')
# plt.plot(df['Time'], df['Baro_Alt'], label='Baro Alt')
# plt.plot(df['Time'], df['IMU_Alt'], label='IMU Alt')
plt.plot(df["Time"], df["Everest_Accel"], label="Everest Acc")
plt.plot(df["Time"], df["Halo_Accel"], label="HALO Acc")
# plt.plot(df['Time'], df['GPS_Alt'], label='GPS Alt')
# plt.plot(df['Time'], df['Sigma_Alt_Upper'], label='Sigma Alt Upper', marker = "x")
# plt.plot(df['Time'], df['Sigma_Alt_Lower'], label='Sigma Alt Lower', marker = "o")
plt.plot(sims["avg_time_1"], sims["avg_acc_1"], label="Scenarios Acc1")
plt.plot(sims["avg_time_2"], sims["avg_acc_2"], label="Scenarios Acc2")


plt.xlabel("Time")
plt.ylabel("Acceleration")
plt.title("Time vs Acceleration")
plt.legend()
plt.grid(True)


plt.show()
