"""
Created on 2020

@author: Eduardo Vitral
"""

###############################################################################
#
# July 2023, Baltimore
#
# This file contains the main functions to convert the Fortran outputs of
# scalefree into Python readable outputs.
#
# If you have any further questions please email evitral@stsci.edu
#
###############################################################################

import subprocess
import numpy as np

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ---------------------------------------------------------------------------
"Interpreters of scalefree's Fortran outputs"
# ---------------------------------------------------------------------------


def vprofile(
    potential="1",
    gamma="4",
    q="1",
    df="1",
    beta="0",
    s="0.5",
    t="0",
    inclination="90",
    integration="1",
    ngl="0",
    algorithm="1",
    maxmom="0",
    theta="0",
    xi="0",
    dim=None,
    average=False,
    exec=False,
    debug=False,
    usevp=True,
):
    """
    Returns the fits of a Gauss-Hermite adjustment to data.
    Note: all results are at an
    (intrinsic or projected) radius of 1 in
    dimensionless units. Results can be scaled to other
    radii using the scale-free nature of the models.

    Parameters
    ----------
    potential: str
        Kepler (1) or Logarithmic (2) Potential.
    gamma: str
        Power-law slope gamma of the mass density.
    q: str
        Intrinsic axial ratio q of the mass density.
    df: str
        Case I (1) or Case II (2) DF
    beta: str
        Anisotropy parameter beta of the DF
    s: str
        Odd part parameter s for the DF
    t: str
        Odd part parameter t for the DF
    inclination: str
        Viewing inclination i in degrees (90=edge-on)
    integration: str
        Use Romberg (0) or Gauss-Legendre (1) integration
        for line-of-sight projection
    ngl: str
        Number of quadrature points / fractional accuracy
    algorithm: str
        Algorithm to calculate VPs and Gauss-Hermite moments
            1: Solve VanderMonde matrix directly without
            regularization. Resulting VP will be nonsense
            but the GH moments are generally well
            determined.
            2: Use regularization with a fixed regularization
            parameter.
            3: Use regularization. Increase regularization
            parameter until the VP has no more than 3
            significant local maxima. A local maximum is
            significant if it exceeds the value of its
            neighbors on the grid by eps times the
            absolute VP maximum.
    maxmom: str
        Maximum number of projected moments
        to use (should be an even number)
    theta: str
        Angle theta in the meridional plane
        (in degrees) (0 = symmetry axis)
    xi: str
        Angle on the projected plane
        (in degrees) (0 = major axis)
    dim: str
        Specific dimension to be considered.
        ("los", "posr", "post").
    average: boolean
        Weather the moments are averaged over the first quadrand of
        the sky. For even models, that is the same as the average
        over the whole sky.
    exec: boolean
        True, if the user wants to generate new .e files.
    debug: boolean
        True, if the user whises to print the Fortran output.
    usevp: boolean
        True, if the user wishes to have VP information.

    Returns
    -------
    vinfo : objects
        List of fit parameters for each of the three observable
        dimensions: LOS, POSr, POSt (in this order).
        Each dimension is itself a collection of dictionaries.
        In this order:
            intmom: Intrinsic velocity moments:
                - rho
                - <v_ph>
                - <v_r^2>
                - <v_th^2>
                - <v_ph^2>
                - beta (for avareged models only)
            projmom: Projected velocity moments:
                - <rho>_p
                - <v>_p
                - <v^2>_p
                - <v^3>_p
                - <v^4>_p
            gauss_info: real Gaussian model.
                - norm
                - mean
                - dispersion
            gaussh_info: Gauss-Hermite model.
                - norm
                - mean
                - dispersion
            h_moments: First 0-6 moments of the
                       Gauss-Hermite model.
                - hi, with i in [0, 6]
            vinfo: Velocity distribution function
                - x
                - f(x)
    """

    if average is True:
        params = [
            potential,
            gamma,
            q,
            df,
            beta,
            s,
            t,
            inclination,
            integration,
            ngl,
            algorithm,
            maxmom,
            "2",
            "1",
            "3",
            "0",
            "0",
        ]
    else:
        params = [
            potential,
            gamma,
            q,
            df,
            beta,
            s,
            t,
            inclination,
            integration,
            ngl,
            algorithm,
            maxmom,
            "0",
            theta,
            "1",
            "1",
            xi,
            "0",
            "0",
        ]

    for i in range(len(params)):
        if isinstance(params[i], str) is False:
            raise ValueError("ERROR: All inputs should be in 'str' format.")
    inputstr = params[0] + "\n"
    for i in range(1, len(params) - 1):
        inputstr += params[i] + "\n"
    inputstr += params[-1]

    if dim is None:
        dimensions = ["los", "posr", "post"]
    else:
        dimensions = [dim]
    vinfo = list()
    for sufix in dimensions:
        h_moments = np.zeros(7)
        for i in range(len(dimensions)):
            prefix = "./scalefree_" + sufix + "/"
            if exec is True:
                p = subprocess.run(
                    ["rm", prefix + "scalefree.e"],
                    text=True,
                    input="y",
                    capture_output=True,
                )
                p = subprocess.run(
                    ["gfortran", prefix + "scalefree.f", "-o", prefix + "scalefree.e"],
                    text=True,
                    input="y",
                    capture_output=True,
                )
            p = subprocess.run(
                [prefix + "scalefree.e", "scalefree.f"],
                text=True,
                input=inputstr,
                capture_output=True,
            )
            split = str(p).split()
            if debug is True:
                print(p.stdout)
            counter = 0
            for j in range(len(split)):
                if usevp is True:
                    if split[j] == r"ITMAX":
                        if split[j + 1] == r"exceeded":
                            raise ValueError(
                                "ERROR: ITMAX exceeded in amoeba. Choose a smaller value for 'maxmom'."
                            )
                if split[j] == r"<v_ph^2>\n":
                    intmom = {
                        "rho": float(split[j + 1]),
                        "<v_ph>": float(split[j + 2]),
                        "<v_r^2>": float(split[j + 3]),
                        "<v_th^2>": float(split[j + 4]),
                        "<v_ph^2>": float(split[j + 5].replace(r"\n", r"")),
                    }
                if split[j] == r"beta\n":
                    intmom = {
                        "rho": float(split[j + 1]),
                        "<v_ph>": float(split[j + 2]),
                        "<v_r^2>": float(split[j + 3]),
                        "<v_th^2>": float(split[j + 4]),
                        "<v_ph^2>": float(split[j + 5]),
                        "beta": float(split[j + 6].replace(r"\n", r"")),
                    }
                if split[j] == r"<v^4>_p\n":
                    projmom = {
                        "<rho>_p": float(split[j + 1]),
                        "<v>_p": float(split[j + 2]),
                        "<v^2>_p": float(split[j + 3]),
                        "<v^3>_p": float(split[j + 4]),
                        "<v^4>_p": float(split[j + 5].replace(r"\n", r"")),
                    }
                if split[j] == r"<v^4>_p\n************":
                    projmom = {
                        "<rho>_p": np.nan,
                        "<v>_p": float(split[j + 1]),
                        "<v^2>_p": float(split[j + 2]),
                        "<v^3>_p": float(split[j + 3]),
                        "<v^4>_p": float(split[j + 4].replace(r"\n", r"")),
                    }
                if usevp is True:
                    if split[j] == r"(gam,V,sig):" and counter == 0:
                        sig = split[j + 3].replace(r"\n", r"")
                        sig = float(sig.replace("Gauss-fit", r""))

                        gauss_info = {
                            "norm": float(split[j + 1]),
                            "mean": float(split[j + 2]),
                            "dispersion": sig,
                        }
                        counter += 1
                    if split[j] == r"(gam,V,sig):" and counter == 1:
                        sig = split[j + 3].replace(r"\n", r"")
                        sig = float(sig.replace("Gauss-fit", r""))

                        gaussh_info = {
                            "norm": float(split[j + 1]),
                            "mean": float(split[j + 2]),
                            "dispersion": sig,
                        }
                        counter += 1
                    if split[j] == r"...":
                        for k in range(len(h_moments)):
                            h_moments[k] = float(split[j + k + 2].replace(r"\n", r""))
                    if split[j] == r"VP(v)\n":
                        v = np.asarray([float(split[j + 1].replace(r"\n", r""))])
                        vp = np.asarray([float(split[j + 2].replace(r"\n", r""))])
                        k = 1
                        stop = False
                        while stop is False:
                            vx = split[j + 1 + 2 * k]
                            if vx == r"\n":
                                stop = True
                            else:
                                v = np.append(v, float(vx.replace(r"\n", r"")))
                                vp = np.append(
                                    vp, float(split[j + 2 + 2 * k].replace(r"\n", r""))
                                )
                                k += 1

        hi = ["h0", "h1", "h2", "h3", "h4", "h5", "h6"]
        h_moments = {hi[i]: h_moments[i] for i in range(len(hi))}
        if usevp is False:
            v = np.zeros(3)
            vp = np.zeros(3)
            gauss_info = np.zeros(3)
            gaussh_info = np.zeros(3)

        vprof = {"x": v, "f(x)": vp}

        fitinfo = (intmom, projmom, gauss_info, gaussh_info, h_moments, vprof)
        vinfo.append(fitinfo)

    if dim is None:
        vinfo = (vinfo[0], vinfo[1], vinfo[2])
    else:
        vinfo = vinfo[0]

    return vinfo


def hermite(input, exec=False):
    """
    Returns the fits of a Gauss-Hermite adjustment to data.

    Parameters
    ----------
    input: str
        Path where to find the data to be fitted.
        Expected file shape:
        --> First line (header): "	v  VP(v)"
        --> 1st column: x / 2nd column: f(x)
    exec: boolean
        True, if the user wants to generate new .e files.

    Returns
    -------
    fitinfo : dictionaries
        In order:
        gauss_info: real Gaussian fit.
                - norm
                - mean
                - dispersion
        gaussh_info: Gauss-Hermite fit.
            - norm
            - mean
            - dispersion
        h_moments: First 0-10 moments of the
                    Gauss-Hermite fit.
            - hi, with i in [0, 10]
    """
    prefix = "./ghermite/"
    if exec is True:
        p = subprocess.run(
            ["rm", prefix + "fitvp.e"],
            text=True,
            input="y",
            capture_output=True,
        )
        p = subprocess.run(
            ["gfortran", prefix + "fitvp.f", "-o", prefix + "fitvp.e"],
            text=True,
            input="y",
            capture_output=True,
        )
    p = subprocess.run(
        [prefix + "fitvp.e", "fitvp.f"],
        text=True,
        input=input,
        capture_output=True,
        shell=True,
    )
    split = str(p).split()

    counter = 0
    h_moments = np.zeros(11)
    for i in range(len(split)):
        if split[i] == r"dispersion:\n":
            if counter == 0:
                gauss_info = {
                    "norm": float(split[i + 1]),
                    "mean": float(split[i + 2]),
                    "dispersion": float(split[i + 3].replace(r"\n", r"")),
                }
                counter += 1
            elif counter == 1:
                gaussh_info = {
                    "norm": float(split[i + 1]),
                    "mean": float(split[i + 2]),
                    "dispersion": float(split[i + 3].replace(r"\n", r"")),
                }
                counter += 1
        if split[i] == r"0-10\n":
            for j in range(len(h_moments)):
                h_moments[j] = float(split[i + j + 1].replace(r"\n',", r""))

    hi = ["h0", "h1", "h2", "h3", "h4", "h5", "h6", "h7", "h8", "h9", "h10"]
    h_moments = {hi[i]: h_moments[i] for i in range(len(hi))}
    fitinfo = (gauss_info, gaussh_info, h_moments)

    return fitinfo
