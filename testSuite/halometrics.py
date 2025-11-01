import pandas as pd
import matplotlib.pyplot as plt


halo = pd.read_csv("testSuite/results/HALO.txt")
altimeter = pd.read_csv("testSuite/data/altimeter(in)(in).csv")


halo.rename(columns={"Time": "time"}, inplace=True)
halo["time"] = halo["time"].round(2)

merged = halo.merge(altimeter, how="inner", on=["time"])


plt.figure(figsize=(10, 6))
plt.plot(
    merged["time"],
    abs(merged["Halo_Alt"] - merged["altitude"]),
    label="HALO vs Altimeter",
)
plt.xlabel("Time")
plt.ylabel("Altitude")
plt.title("HALO Alt vs Altimeter Alt")
plt.legend()
plt.grid(True)


plt.figure(figsize=(10, 6))
plt.plot(
    merged["time"],
    abs(merged["speed"] - abs(merged["Halo_Velo"])),
    label="HALO vs Altimeter",
)
plt.xlabel("Time")
plt.ylabel("Altitude")
plt.title("Altimeter Speed vs Halo Velocity (absolute) absolute difference")
plt.legend()
plt.grid(True)


plt.figure(figsize=(10, 6))
plt.plot(
    merged["time"],
    abs(merged["Halo_Accel"] - merged["acceleration"]),
    label="HALO vs Altimeter",
)
plt.xlabel("Time")
plt.ylabel("Altitude")
plt.title("Altimeter Acceleration vs Halo acceleration absolute difference")
plt.legend()
plt.grid(True)

error_alt = abs(merged["Halo_Alt"] - merged["altitude"]) / merged["altitude"] * 100
error_speed = abs(merged["speed"] - abs(merged["Halo_Velo"])) / merged["speed"] * 100
error_acc = (
    abs(merged["Halo_Accel"] - merged["acceleration"]) / merged["acceleration"] * 100
)

avg = error_alt * 0.6 + error_acc * 0.4 + error_speed * 0.4

plt.figure(figsize=(10, 6))
plt.plot(
    merged["time"],
    avg,
    label="HALO vs Altimeter",
)
plt.xlabel("Time")
plt.ylabel("Error")
plt.title("Error Metric")
plt.legend()
plt.grid(True)

plt.show()
