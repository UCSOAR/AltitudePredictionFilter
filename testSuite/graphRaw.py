import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# make numpy array from 0 to 180 by 1/3 step
x = np.arange(0, 185, 1 / 3)


def graphRaw():
    # Load data
    df = pd.read_csv("testSuite/data/Imu_Baro.csv")

    # Plot data
    plt.plot(x, df["accel_x"], label="Ax")
    plt.plot(x, df["accel_y"], label="Ay")
    plt.plot(x, df["accel_z"], label="Az")

    # Add labels
    plt.xlabel("Time (s)")
    plt.ylabel("Acceleration (m/s^2)")

    # Add legend
    plt.legend()
    plt.show()

    return


graphRaw()
