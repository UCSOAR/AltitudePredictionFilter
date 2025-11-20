import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("testSuite/results/kalmangains.txt")

df = df.groupby("time").mean().reset_index()
print(df)


plt.figure(figsize=(10, 6))
plt.plot(
    df["time"],
    abs(df["acc"]),
    label="acc",
)
plt.plot(
    df["time"],
    abs(df["velo"]),
    label="velo",
)
plt.plot(
    df["time"],
    abs(df["alt"]),
    label="alt",
)
plt.plot(
    df["time"],
    abs(df["gps_alt"]),
    label="gps alt",
)
plt.xlabel("Time")
plt.ylabel("Gain")
plt.title("Kalman gain over time")
plt.legend()
plt.grid(True)

plt.show()
