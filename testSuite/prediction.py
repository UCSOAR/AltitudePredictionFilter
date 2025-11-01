import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv("testSuite/results/HALO.txt")
pred = pd.read_csv("testSuite/results/predictnalt.txt")

plt.figure(figsize=(10, 6))
plt.plot(df["Time"], df["Halo_Alt"], label="HALO Alt")
plt.plot(pred["time"], pred["predicted_alt"], label="Predicted Alt")
plt.xlabel("Time")
plt.ylabel("Altitude")
plt.title("Predicted Alt vs HALO Alt")
plt.legend()
plt.grid(True)

plt.show()
