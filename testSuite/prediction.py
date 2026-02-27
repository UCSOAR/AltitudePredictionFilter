import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("testSuite/results/HALO.txt")
pred = pd.read_csv("testSuite/results/predictnalt.txt")
altimeter = pd.read_csv("testSuite/data/altimeter(in)(in).csv")

sims = pd.read_csv("testSuite/data/final_may_sims_formatted.csv")

sims["avg_alt_1"] = pd.to_numeric(sims["avg_alt_1"], errors="coerce")
sims["avg_velo_1"] = pd.to_numeric(sims["avg_velo_1"], errors="coerce")
sims["avg_acc_1"] = pd.to_numeric(sims["avg_acc_1"], errors="coerce")
sims["avg_time_1"] = pd.to_numeric(sims["avg_time_1"], errors="coerce")

sims["avg_alt_2"] = pd.to_numeric(sims["avg_alt_2"], errors="coerce")
sims["avg_velo_2"] = pd.to_numeric(sims["avg_velo_2"], errors="coerce")
sims["avg_acc_2"] = pd.to_numeric(sims["avg_acc_2"], errors="coerce")
sims["avg_time_2"] = pd.to_numeric(sims["avg_time_2"], errors="coerce")
sims = sims.dropna()


plt.figure(figsize=(10, 6))

plt.plot(df["Time"], df["Halo_Alt"], label="HALO Alt")
plt.plot(altimeter["time"], altimeter["altitude"], label="Altimeter Alt")
plt.plot(sims["avg_time_1"], sims["avg_alt_1"], label="Scenarios Alt1")
plt.plot(sims["avg_time_2"], sims["avg_alt_2"], label="Scenarios Alt2")

for pid, group in pred.groupby("prediction"):
    plt.plot(group["time"], group["predicted_alt"], label=f"Predicted Alt {pid}")

plt.xlabel("Time")
plt.ylabel("Altitude")
plt.title("Predicted Alt")
plt.legend()
plt.grid(True)

plt.figure(figsize=(10, 6))

plt.plot(df["Time"], df["Halo_Velo"], label="HALO Velo")
plt.plot(altimeter["time"], altimeter["speed"], label="Altimeter velo")
plt.plot(sims["avg_time_1"], sims["avg_velo_1"], label="Scenarios Velo1")
plt.plot(sims["avg_time_2"], sims["avg_velo_2"], label="Scenarios Velo2")

for pid, group in pred.groupby("prediction"):
    plt.plot(group["time"], group["predicted_vel"], label=f"Predicted Velo {pid}")

plt.xlabel("Time")
plt.ylabel("Velocity")
plt.title("Predicted velo")
plt.legend()
plt.grid(True)


plt.figure(figsize=(10, 6))

plt.plot(df["Time"], df["Halo_Accel"], label="HALO Accel")
plt.plot(altimeter["time"], altimeter["acceleration"], label="Altimeter accel")
plt.plot(sims["avg_time_1"], sims["avg_acc_1"], label="Scenarios Acc1")
plt.plot(sims["avg_time_2"], sims["avg_acc_2"], label="Scenarios Acc2")

for pid, group in pred.groupby("prediction"):
    plt.plot(group["time"], group["predicted_acc"], label=f"Predicted Acc {pid}")

plt.xlabel("Time")
plt.ylabel("Acceleration")
plt.title("Predicted Acc")
plt.legend()
plt.grid(True)


plt.show()
