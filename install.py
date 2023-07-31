import os
import sys
import numpy as np

print("Provide the path where your scalefree files are installed.")
print("For example: /home/yourname/code-packages/scalefree/", "\n")
path_sf = input()

os.chdir(path_sf)
sys.path.append(path_sf)
import scalefree

path_input = r"./examp/vlos_profile.txt"
expect_ghermite = [
    {"norm": 0.9998, "mean": -291.659, "dispersion": 8.5193},
    {"norm": 0.9775, "mean": -291.5934, "dispersion": 7.9392},
    {
        "h0": 1.0,
        "h1": 0.0,
        "h2": -0.0,
        "h3": -0.008,
        "h4": 0.056,
        "h5": 0.013,
        "h6": -0.032,
        "h7": -0.013,
        "h8": 0.013,
        "h9": 0.003,
        "h10": 0.009,
    },
]
res = scalefree.hermite(path_input, exec=True)
for i in range(len(res)):
    check = res[i] == expect_ghermite[i]
    if check is False:
        raise ValueError(
            "ERROR: Installation failed at Gauss-Hermite fit. Fit #", i + 1
        )

expect_v = [
    {
        "rho": 1.0,
        "<v_ph>": 0.0,
        "<v_r^2>": 0.1,
        "<v_th^2>": 0.1,
        "<v_ph^2>": 0.1,
        "beta": 0.0,
    },
    {
        "<rho>_p": 1.57079633,
        "<v>_p": 0.0,
        "<v^2>_p": 0.08488264,
        "<v^3>_p": 0.0,
        "<v^4>_p": 0.01875,
    },
    {"norm": 1.0, "mean": 0.0, "dispersion": 0.29134625},
    {"norm": 1.0, "mean": 0.0, "dispersion": 0.29134625},
    {"h0": 1.0, "h1": -0.0, "h2": 0.0, "h3": 0.0, "h4": 0.2439, "h5": 0.0, "h6": 0.0},
]

expect_vp = [
    np.asarray([-0.8, -0.4, 0.0, 0.4, 0.8]),
    np.asarray([0.02103181, 0.57901835, 1.29989968, 0.57901835, 0.02103181]),
]

res3d = scalefree.vprofile(maxmom="4", exec=True, average=True)
for j in range(len(res3d)):
    res = res3d[j]
    for i in range(len(res) - 1):
        check = res[i] == expect_v[i]
        if check is False:
            raise ValueError(
                "ERROR: Installation failed at Gauss-Hermite fit. Fit #",
                j + 1,
                "/",
                i + 1,
            )
    check = res[-1]["x"] == expect_vp[0]
    if check is False:
        raise ValueError(
            "ERROR: Installation failed at Gauss-Hermite fit. Fit #",
            j + 1,
            "/",
            i + 1,
        )
    check = res[-1]["f(x)"] == expect_vp[1]
    if check is False:
        raise ValueError(
            "ERROR: Installation failed at Gauss-Hermite fit. Fit #",
            j + 1,
            "/",
            i + 1,
        )

print("\nInstallation succesfull.")
print("Finished.")
