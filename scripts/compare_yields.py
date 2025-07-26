"""A comparison of Yield outputs from Stark et al.'s Altruistic Yield Optimizer,
specifically for generating 'Yield Difference' plots
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import ipdb

COMPARISON_PTH = "EAC1_XeLiF_550nm_8m"

YIELDS_DIR = Path.home() / "FRIDAY_Outputs"
YIELD1_DIR = YIELDS_DIR / "EAC1_XeLiF_550nm_16m_Scalar"
YIELD2_DIR = YIELDS_DIR / COMPARISON_PTH

yield_1 = 0
yield_2 = 0

obs_1 = pd.read_csv(YIELD1_DIR / "observations.csv", delimiter=",")
obs_2 = pd.read_csv(YIELD2_DIR / "observations.csv", delimiter=",")

targets_1 = pd.read_csv(YIELD1_DIR / "target_list.csv", delimiter=",")
targets_2 = pd.read_csv(YIELD2_DIR / "target_list.csv", delimiter=",")

key = "exoEarth candidate yield"
key = "Exp Time (days)"
key = "Spec char time (days)"

ayo_keys = [
    "exoEarth candidate yield",
    "Exp Time (days)",
    "Spec char time (days)",
]
assert key in ayo_keys

if key == "Exp Time (days)":
    cbar_label = key
    vlim = 3
elif key == "Spec char time (days)":
    cbar_label = key
    vlim = 2
elif key == "exoEarth candidate yield":
    cbar_label = r"$\Delta$ HZ Completeness"
    vlim = 0.15

# Dictionaries keyed by StarID
data_1 = {}
data_2 = {}

# Collapse data from obs_1
for index, row in obs_1.iterrows():

    starID = row["starID"]

    if starID not in data_1.keys():
        # construct a dict that contains xyz

        target_row = targets_1[targets_1["starID"] == starID]

        data = {
            key: row[key],
            "Lstar": target_row["Lstar (Lsun)"],
            "dist (pc)": target_row["dist (pc)"]
        }
        data_1[starID] = data

    else:
        data_1[starID][key] += row[key]

    yield_1 += row["exoEarth candidate yield"]

# Collapse data from obs_2
for index, row in obs_2.iterrows():

    starID = row["starID"]

    if starID not in data_2.keys():
        # construct a dict that contains xyz
        target_row = targets_2[targets_2["starID"] == starID]

        data = {
            key: row[key],
            "Lstar": target_row["Lstar (Lsun)"],
            "dist (pc)": target_row["dist (pc)"]
        }
        data_2[starID] = data
    else:
        data_2[starID][key] += row[key]

    yield_2 += row["exoEarth candidate yield"]

# Brute force search over available keys
# TODO: This could merit from some optimization

plt.figure()

# Loop over starID 1 and look for matches in starID2
targs_missed = 0
new_targs = 0
for sid1 in data_1.keys():
    if sid1 in data_2.keys():
        # Compute the difference of the data
        data = data_1[sid1]
        y = data["Lstar"]
        x = data["dist (pc)"]
        z1 = data[key]
        z2 = data_2[sid1][key]
        z_diff = z1 - z2

        im = plt.scatter(x, y, c=z_diff, marker="o", cmap="RdBu", vmin=-vlim, vmax=vlim, s=100)

    else:
        data = data_1[sid1]
        y = data["Lstar"]
        x = data["dist (pc)"]
        plt.scatter(x, y, color="k", marker="x")
        targs_missed += 1

# Catch extra targets
for sid2 in data_2.keys():
    if sid2 not in data_1.keys():
        data = data_2[sid2]
        y = data["Lstar"]
        x = data["dist (pc)"]
        plt.scatter(x, y, edgecolor="r", facecolor="None", marker="o")
        new_targs += 1


plt.colorbar(im, label=cbar_label)
plt.xlabel("dist (pc)", fontsize=16)
plt.ylabel(r"$L_{\odot}$", fontsize=16, rotation=0)
plt.xlim([0, 25])
plt.ylim([1e-2, 30])
plt.yscale("log")
pct_diff = (1 - (yield_2 / yield_1)) * 100
plt.title(COMPARISON_PTH+f"       Yield = {yield_2:.2f} (-{pct_diff:.2f}%)")
plt.scatter(0, -1, color="k", marker="x", label=f"{targs_missed} Missed Targets")
plt.scatter(0, -1, edgecolor="r", facecolor="None", marker="o", label=f"{new_targs} Additional Targets")
plt.legend(loc="lower right")
plt.show()
