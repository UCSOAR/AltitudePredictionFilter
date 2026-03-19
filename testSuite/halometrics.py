import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
import os

apogee_time = 30.505

ignore_after_apogee = True

save_graphs = True

SAVE_DIR = "./plots"
os.makedirs(SAVE_DIR, exist_ok=True)


def save_plot(name):
    plt.savefig(os.path.join(SAVE_DIR, name), dpi=300, bbox_inches="tight")


# -----------------------------
# Load data
# -----------------------------
halo = pd.read_csv("testSuite/results/HALO.txt")
altimeter = pd.read_csv("testSuite/data/altimeter(in)(in).csv")
covariance = pd.read_csv("testSuite/results/P.txt")
nis_data = pd.read_csv("testSuite/results/nis.txt")

halo.rename(columns={"Time": "time"}, inplace=True)

merged = pd.merge_asof(halo, altimeter, on="time", direction="nearest", tolerance=0.05)

if ignore_after_apogee:
    merged = merged[merged["time"] <= apogee_time]
    nis_data = nis_data[nis_data["time"] <= apogee_time]

# -----------------------------
# Absolute state errors
# -----------------------------
alt_err = merged["Halo_Alt"] - merged["altitude"]
vel_err = merged["Halo_Velo"] - merged["speed"]
acc_err = merged["Halo_Accel"] - merged["acceleration"]

plt.figure(figsize=(10, 6))
plt.plot(merged["time"], abs(alt_err), label="HALO Altitude")
plt.plot(merged["time"], abs(vel_err), label="HALO Velocity")
plt.plot(merged["time"], abs(acc_err), label="HALO Acceleration")
plt.axvline(x=apogee_time, color="purple", linestyle="--", alpha=0.7, label="Apogee")

# Root Mean Squared Error
rmse_alt = np.sqrt(np.mean(alt_err**2))
rmse_vel = np.sqrt(np.mean(vel_err**2))
rmse_acc = np.sqrt(np.mean(acc_err**2))

plt.text(
    0.02,
    0.95,
    f"Overall RMSE\nAlt: {rmse_alt:.3f}\nVel: {rmse_vel:.3f}\nAcc: {rmse_acc:.3f}",
    transform=plt.gca().transAxes,
    verticalalignment="top",
    bbox=dict(boxstyle="round", alpha=0.3),
)

plt.title("Absolute State Error (HALO)")
plt.xlabel("Time")
plt.ylabel("Error")
plt.legend()
plt.grid(True)
if save_graphs:
    save_plot("Absolute State Error (HALO)")

# Everest Abs Error

e_alt_err = merged["Everest_Alt"] - merged["altitude"]
e_vel_err = merged["Everest_Velo"] - merged["speed"]
e_acc_err = merged["Everest_Accel"] - merged["acceleration"]

plt.figure(figsize=(10, 6))
plt.plot(merged["time"], abs(e_alt_err), label="Everest Altitude")
plt.plot(merged["time"], abs(e_vel_err), label="Everest Velocity")
plt.plot(merged["time"], abs(e_acc_err), label="Everest Acceleration")
plt.axvline(x=apogee_time, color="purple", linestyle="--", alpha=0.7, label="Apogee")

# Root Mean Squared Error
rmse_alt = np.sqrt(np.mean(e_alt_err**2))
rmse_vel = np.sqrt(np.mean(e_vel_err**2))
rmse_acc = np.sqrt(np.mean(e_acc_err**2))

plt.text(
    0.02,
    0.95,
    f"Overall RMSE\nAlt: {rmse_alt:.3f}\nVel: {rmse_vel:.3f}\nAcc: {rmse_acc:.3f}",
    transform=plt.gca().transAxes,
    verticalalignment="top",
    bbox=dict(boxstyle="round", alpha=0.3),
)

plt.title("Absolute State Error (Everest)")
plt.xlabel("Time")
plt.ylabel("Error")
plt.legend()
plt.grid(True)
if save_graphs:
    save_plot("Absolute State Error (Everest)")

# -----------------------------
# Clipped percent error
# -----------------------------
MAX_ERR = 200


def percent_error(pred, truth):
    denom = np.maximum(np.abs(truth), 1e-6)
    err = np.abs(pred - truth) / denom * 100
    return np.clip(err, 0, MAX_ERR)


percent_error_alt = percent_error(merged["Halo_Alt"], merged["altitude"])
percent_error_vel = percent_error(merged["Halo_Velo"], merged["speed"])
percent_error_acc = percent_error(merged["Halo_Accel"], merged["acceleration"])

plt.figure(figsize=(10, 6))
plt.plot(merged["time"], percent_error_alt, label="HALO Altitude")
plt.plot(merged["time"], percent_error_vel, label="HALO Velocity")
plt.plot(merged["time"], percent_error_acc, label="HALO Acceleration")
plt.axvline(x=apogee_time, color="purple", linestyle="--", alpha=0.7, label="Apogee")
plt.title("Clipped Percent Error, HALO (200%)")
plt.xlabel("Time")
plt.ylabel("Error (%)")
plt.legend()
plt.grid(True)
if save_graphs:
    save_plot("Clipped Percent Error, HALO (200%)")


# Everest

MAX_ERR = 200


def percent_error(pred, truth):
    denom = np.maximum(np.abs(truth), 1e-6)
    err = np.abs(pred - truth) / denom * 100
    return np.clip(err, 0, MAX_ERR)


e_percent_error_alt = percent_error(merged["Everest_Alt"], merged["altitude"])
e_percent_error_vel = percent_error(merged["Everest_Velo"], merged["speed"])
e_percent_error_acc = percent_error(merged["Everest_Accel"], merged["acceleration"])

plt.figure(figsize=(10, 6))
plt.plot(merged["time"], e_percent_error_alt, label="Everest Altitude")
plt.plot(merged["time"], e_percent_error_vel, label="Everest Velocity")
plt.plot(merged["time"], e_percent_error_acc, label="Everest Acceleration")
plt.axvline(x=apogee_time, color="purple", linestyle="--", alpha=0.7, label="Apogee")
plt.title("Clipped Percent Error, Everest (200%)")
plt.xlabel("Time")
plt.ylabel("Error (%)")
plt.legend()
plt.grid(True)
if save_graphs:
    save_plot("Clipped Percent Error, Everest (200%)")

# -----------------------------
# NEES calculation
# -----------------------------
nees = []
for i in range(len(merged)):
    e = np.array([alt_err.iloc[i], vel_err.iloc[i], acc_err.iloc[i]])
    P = np.array(
        [
            [
                covariance["P00"].iloc[i],
                covariance["P01"].iloc[i],
                covariance["P02"].iloc[i],
            ],
            [
                covariance["P10"].iloc[i],
                covariance["P11"].iloc[i],
                covariance["P12"].iloc[i],
            ],
            [
                covariance["P20"].iloc[i],
                covariance["P21"].iloc[i],
                covariance["P22"].iloc[i],
            ],
        ]
    )
    try:
        val = e.T @ np.linalg.inv(P) @ e
    except np.linalg.LinAlgError:
        val = np.nan
    nees.append(val)

nees = np.array(nees)

# -----------------------------
# NEES χ² consistency bounds
# -----------------------------
state_dim = 3
alpha = 0.95
nees_lower = chi2.ppf((1 - alpha) / 2, state_dim)
nees_upper = chi2.ppf(1 - (1 - alpha) / 2, state_dim)

# -----------------------------
# Plot NEES
# -----------------------------
plt.figure(figsize=(10, 6))
plt.plot(merged["time"], nees, label="NEES", alpha=0.6)
plt.axhline(
    nees_lower, color="r", linestyle="--", label=f"χ² bounds ({alpha*100:.0f}%)"
)
plt.axhline(nees_upper, color="r", linestyle="--")
plt.axhline(state_dim, color="g", linestyle=":", label=f"Expected mean ({state_dim})")
plt.axvline(x=apogee_time, color="purple", linestyle="--", alpha=0.7, label="Apogee")
plt.title("NEES Consistency Test")
plt.xlabel("Time (s)")
plt.ylabel("NEES")
plt.legend()
plt.grid(True, alpha=0.3)
if save_graphs:
    save_plot("NEES Consistency Test")

# -----------------------------
# Smoothed NEES
# -----------------------------
nees_smooth = pd.Series(nees).rolling(20).mean()

plt.figure(figsize=(10, 6))
plt.plot(merged["time"], nees_smooth, label="NEES (moving avg)")
plt.axhline(state_dim, color="g", linestyle="--", label="Expected value")
plt.axvline(x=apogee_time, color="purple", linestyle="--", alpha=0.7, label="Apogee")
plt.title("Smoothed NEES")
plt.xlabel("Time (s)")
plt.ylabel("NEES")
plt.legend()
plt.grid(True, alpha=0.3)
if save_graphs:
    save_plot("Smoothed NEES")

# -----------------------------
# NIS Analysis
# -----------------------------
measurement_dim = 4  # 4 measurements: alt, vel, acc, gps_alt
alpha_nis = 0.05  # 95% confidence interval
nis_lower = chi2.ppf(alpha_nis / 2, measurement_dim)
nis_upper = chi2.ppf(1 - alpha_nis / 2, measurement_dim)

# Calculate statistics
nis_values = nis_data["nis"].values
valid_nis = nis_values[~np.isnan(nis_values)]

in_bounds = np.sum((valid_nis >= nis_lower) & (valid_nis <= nis_upper))
above_upper = np.sum(valid_nis > nis_upper)
below_lower = np.sum(valid_nis < nis_lower)
total = len(valid_nis)

percentage_in = 100 * in_bounds / total
percentage_above = 100 * above_upper / total
percentage_below = 100 * below_lower / total

mean_nis = np.mean(valid_nis)
std_nis = np.std(valid_nis)
expected_std = np.sqrt(2 * measurement_dim)

# -----------------------------
# Plot NIS
# -----------------------------
plt.figure(figsize=(12, 6))
plt.plot(nis_data["time"], nis_data["nis"], "b-", alpha=0.6, label="NIS")
plt.axhline(
    y=nis_lower, color="r", linestyle="--", label=f"{(1-alpha_nis)*100:.0f}% bounds"
)
plt.axhline(y=nis_upper, color="r", linestyle="--")
plt.axhline(
    y=measurement_dim,
    color="g",
    linestyle=":",
    label=f"Expected mean ({measurement_dim})",
)
plt.axvline(x=apogee_time, color="purple", linestyle="--", alpha=0.7, label="Apogee")
plt.xlabel("Time (s)")
plt.ylabel("NIS")
plt.title("NIS Consistency Check")
plt.legend()
plt.grid(True, alpha=0.3)
if save_graphs:
    save_plot("NIS Consistency Check")

# Add statistics text box
stats_text = f"""NIS Statistics:
In bounds: {percentage_in:.1f}% (expect ~95%)
Above upper: {percentage_above:.1f}%
Below lower: {percentage_below:.1f}%
Mean: {mean_nis:.2f} (expect {measurement_dim})
Std: {std_nis:.2f} (expect ~{expected_std:.2f})"""

plt.text(
    0.02,
    0.98,
    stats_text,
    transform=plt.gca().transAxes,
    verticalalignment="top",
    bbox=dict(boxstyle="round", alpha=0.3, facecolor="white"),
)

# -----------------------------
# Smoothed NIS
# -----------------------------
nis_smooth = pd.Series(nis_data["nis"]).rolling(20).mean()

plt.figure(figsize=(12, 6))
plt.plot(nis_data["time"], nis_smooth, "b-", label="NIS (moving avg)")
plt.axhline(
    y=measurement_dim,
    color="g",
    linestyle="--",
    label=f"Expected mean ({measurement_dim})",
)
plt.axhline(y=nis_lower, color="r", linestyle=":", alpha=0.5, label="Bounds")
plt.axhline(y=nis_upper, color="r", linestyle=":", alpha=0.5)
plt.axvline(x=apogee_time, color="purple", linestyle="--", alpha=0.7, label="Apogee")
plt.xlabel("Time (s)")
plt.ylabel("NIS")
plt.title("Smoothed NIS (20-sample moving average)")
plt.legend()
plt.grid(True, alpha=0.3)
if save_graphs:
    save_plot("Smoothed NIS (20-sample moving average)")

# -----------------------------
# Print Summary Report
# -----------------------------
print("=" * 60)
print("NIS CONSISTENCY REPORT")
print("=" * 60)
print(f"Total samples: {total}")
print(f"Measurement dimension: {measurement_dim}")
print(f"Expected bounds: [{nis_lower:.2f}, {nis_upper:.2f}]")
print(f"\nIn bounds:    {percentage_in:6.2f}% (expect ~95%)")
print(f"Above upper:  {percentage_above:6.2f}%")
print(f"Below lower:  {percentage_below:6.2f}%")
print(f"\nMean NIS:     {mean_nis:6.3f} (expect {measurement_dim:.3f})")
print(f"Std NIS:      {std_nis:6.3f} (expect ~{expected_std:.3f})")
print("=" * 60)

# Interpretation
print("\nINTERPRETATION:")
if percentage_above > 10:
    print("⚠ Filter is OVERCONFIDENT - too many values above upper bound")
    print("  → Consider increasing process noise Q")
elif percentage_below > 10:
    print("⚠ Filter is UNDERCONFIDENT - too many values below lower bound")
    print("  → Consider decreasing process noise Q")
else:
    print("✓ Filter appears well-tuned")

if abs(mean_nis - measurement_dim) > 1:
    print(f"⚠ Mean NIS ({mean_nis:.2f}) deviates from expected ({measurement_dim})")
else:
    print("✓ Mean NIS is close to expected value")

print("=" * 60)

plt.show()
