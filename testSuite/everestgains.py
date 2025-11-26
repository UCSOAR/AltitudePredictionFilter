import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("testSuite/results/everestgains.txt")
df = df.apply(pd.to_numeric, errors="coerce")

plt.figure(figsize=(10, 6))
plt.plot(df["time"], abs(df["gain_IMU"]), label="gain_IMU")
plt.plot(df["time"], abs(df["gain_Baro1"]), label="gain_Baro1")
plt.plot(df["time"], abs(df["gain_Baro2"]), label="gain_Baro2")

plt.xlabel("Time")
plt.ylabel("Gain")
plt.title("Everest gain over time")
plt.legend()
plt.grid(True)
plt.show()
