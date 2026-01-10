import pandas as pd
import matplotlib.pyplot as plt


halo = pd.read_csv("testSuite/results/HALO.txt")
altimeter = pd.read_csv("testSuite/data/altimeter(in)(in).csv")


halo.rename(columns={"Time": "time"}, inplace=True)
halo["time"] = halo["time"].round(2)

merged = halo.merge(altimeter, how="inner", on=["time"])

# Alt difference
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

# Speed Diff
plt.figure(figsize=(10, 6))
plt.plot(
    merged["time"],
    abs(merged["speed"] - abs(merged["Halo_Velo"])),
    label="HALO vs Altimeter",
)
plt.xlabel("Time")
plt.ylabel("Speed")
plt.title("Altimeter Speed vs Halo Velocity (absolute) absolute difference")
plt.legend()
plt.grid(True)

# Acc dif
plt.figure(figsize=(10, 6))
plt.plot(
    merged["time"],
    abs(merged["Halo_Accel"] - merged["acceleration"]),
    label="HALO vs Altimeter",
)
plt.xlabel("Time")
plt.ylabel("Acceleration")
plt.title("Altimeter Acceleration vs Halo acceleration absolute difference")
plt.legend()
plt.grid(True)


# error

error_alt = (abs(merged["Halo_Alt"] - merged["altitude"]) / merged["altitude"]) * 100
error_speed = (abs(merged["Halo_Velo"] - merged["speed"]) / merged["speed"]) * 100
error_acc = (
    abs(merged["Halo_Accel"] - merged["acceleration"]) / merged["acceleration"]
) * 100


avg = (error_alt + error_acc + error_speed) / 3

plt.figure(figsize=(10, 6))

plt.plot(
    merged["time"],
    error_alt,
    label="Altitude",
)
plt.plot(
    merged["time"],
    error_speed,
    label="Speed",
)
plt.plot(
    merged["time"],
    error_speed,
    label="Acceleration",
)
plt.plot(
    merged["time"],
    avg,
    label="Average",
)
plt.xlabel("Time")
plt.ylabel("Error")
plt.title("Error Metric")
plt.legend()
plt.grid(True)

# Confidence over time

import numpy as np

covariance = pd.read_csv("testSuite/results/P.txt").dropna()

# square root
covariance.loc[:, covariance.columns != "time"] = np.sqrt(
    covariance.loc[:, covariance.columns != "time"]
)

# just for naming
confidence = covariance

confidence_avg = (
    covariance["velo_std"] + covariance["acc_std"] + covariance["alt_std"]
) / 3.0

plt.figure(figsize=(10, 6))

plt.plot(
    confidence["time"],
    confidence["alt_std"],
    label="Altitude",
)
plt.plot(
    confidence["time"],
    confidence["velo_std"],
    label="Speed",
)
plt.plot(
    confidence["time"],
    confidence["acc_std"],
    label="Acceleration",
)
plt.plot(
    confidence["time"],
    confidence_avg,
    label="Average",
)
plt.xlabel("Time")
plt.ylabel("Confidence")
plt.title("Confidence Over Time")
plt.legend()
plt.grid(True)


plt.show()
