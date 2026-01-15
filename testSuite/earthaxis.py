import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("testSuite/results/infusion.txt")

# Alt difference
plt.figure(figsize=(10, 6))
plt.plot(
    df["timestamp"],
    df["earth.axis.z"],
    label="down axis",
)
plt.xlabel("Time")
plt.ylabel("Down Axis")
plt.title("Earth down axis over time")
plt.legend()
plt.grid(True)

plt.show()
