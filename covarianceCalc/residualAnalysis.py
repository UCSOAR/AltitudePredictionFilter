import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# file paths

HALO_FILE = "testSuite/results/HALO.txt"
P_FILE = "testSuite/results/P.txt"
CONFIDENCE_FILE = "testSuite/results/confidence.txt"
ALTIMETER_FILE = "testSuite/data/altimeter(in)(in).csv"
SIMS_FILE = "testSuite/data/final_may_sims_formatted.csv"
QUASAR_SIMS_FILE = "covarianceCalc/covariance_Cal_Data.csv"

# Load data
HALO_data = pd.read_csv(HALO_FILE)
p_data = pd.read_csv(P_FILE)
confidence = pd.read_csv(CONFIDENCE_FILE)
altimeter_data = pd.read_csv(ALTIMETER_FILE)
sims = pd.read_csv(SIMS_FILE)
Quasar_data = pd.read_csv(QUASAR_SIMS_FILE)

# clean sims data
sims = sims.replace(["#DIV/0!", "#REF!"], np.nan).dropna().reset_index(drop=True)


# average data over interval
def calculate_averages(time, alt, vel, acc, start_time=0, interval=0.1):
    """Average altitude, velocity, acceleration over time intervals."""
    time = pd.to_numeric(time, errors="coerce")
    alt = pd.to_numeric(alt, errors="coerce")
    vel = pd.to_numeric(vel, errors="coerce")
    acc = pd.to_numeric(acc, errors="coerce")

    avg_time, avg_alt, avg_vel, avg_acc = [], [], [], []
    lower = start_time
    upper = start_time + interval
    max_time = time.max()

    while lower < max_time:
        mask = (time >= lower) & (time < upper)
        count = mask.sum()
        if count > 0:
            avg_alt.append(alt[mask].mean())
            avg_vel.append(vel[mask].mean())
            avg_acc.append(acc[mask].mean())
            avg_time.append(upper)
        lower = upper
        upper += interval

    return np.array(avg_time), np.array(avg_alt), np.array(avg_vel), np.array(avg_acc)


# residual calculation
def compute_residuals(pred_time, pred_values, meas_time, meas_values):
    """Compute residuals by interpolating predicted values onto measurement times."""
    pred_values_interp = np.interp(meas_time, pred_time, pred_values)
    residuals = meas_values - pred_values_interp
    return residuals


# plot function
def plot_residuals(
    time, predicted, measured, residuals, ylabel="Altitude", title="Residual Analysis"
):
    plt.figure(figsize=(10, 6))
    plt.plot(time, measured, label="Measured")
    plt.plot(time, predicted, label="Predicted", linestyle="--")
    plt.plot(time, residuals, label="Residual", linestyle=":")
    plt.xlabel("Time")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.show()


# main residual analysis
# -----------------------------
def residual_analysis():
    # HALO vs Altimeter

    alt_time, alt_avg, vel_avg, acc_avg = calculate_averages(
        altimeter_data["time"],
        altimeter_data["altitude"],
        altimeter_data["speed"],
        altimeter_data["acceleration"],
        start_time=0,
        interval=0.333,
    )

    e_time, e_alt_avg, e_vel_avg, e_acc_avg = calculate_averages(
        HALO_data["Time"][1:],
        HALO_data["Halo_Alt"][1:],
        HALO_data["Halo_Velo"][1:],
        HALO_data["Halo_Accel"][1:],
        start_time=0,
        interval=0.333,
    )

    # residuals
    alt_residuals = compute_residuals(alt_time, alt_avg, e_time, e_alt_avg)
    vel_residuals = compute_residuals(alt_time, vel_avg, e_time, e_vel_avg)
    acc_residuals = compute_residuals(alt_time, acc_avg, e_time, e_acc_avg)

    # covariance matrix of residuals
    residual_matrix = np.vstack([alt_residuals, vel_residuals, acc_residuals])
    covariance_matrix = np.cov(residual_matrix)
    print("Covariance Matrix of Residuals (Altitude, Velocity, Acceleration):")
    print(covariance_matrix)

    alt_pred_interp = np.interp(e_time, alt_time, alt_avg)
    vel_pred_interp = np.interp(e_time, alt_time, vel_avg)
    acc_pred_interp = np.interp(e_time, alt_time, acc_avg)

    alt_residuals = e_alt_avg - alt_pred_interp
    vel_residuals = e_vel_avg - vel_pred_interp
    acc_residuals = e_acc_avg - acc_pred_interp

    min_len = min(len(e_time), len(alt_residuals))
    e_time_trim = e_time[:min_len]
    alt_residuals = alt_residuals[:min_len]
    vel_residuals = vel_residuals[:min_len]
    acc_residuals = acc_residuals[:min_len]

    plot_residuals(
        e_time_trim,
        alt_pred_interp[:min_len],
        e_alt_avg[:min_len],
        alt_residuals,
        ylabel="Altitude",
        title="Altitude Residuals",
    )
    plot_residuals(
        e_time_trim,
        vel_pred_interp[:min_len],
        e_vel_avg[:min_len],
        vel_residuals,
        ylabel="Velocity",
        title="Velocity Residuals",
    )
    plot_residuals(
        e_time_trim,
        acc_pred_interp[:min_len],
        e_acc_avg[:min_len],
        acc_residuals,
        ylabel="Acceleration",
        title="Acceleration Residuals",
    )


if __name__ == "__main__":
    residual_analysis()
