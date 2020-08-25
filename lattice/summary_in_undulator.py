import datetime
import time
import numpy as np
import scipy
import pandas as pd
import seaborn as sns
import sys
import os
from config import get_from_config, save_to_config
from acnet_reader.acnet_reader import fetch_data, \
    get_interpolated_df, fetch_interpolated_data
import fur.path_assistant as path_assistant
import lattice.lattice as lattice


def get_summary_in_undulator(lattice_file,
                             camera_sizes_X_um,
                             camera_sizes_Y_um,
                             dpp=None,
                             dpp_err=0
                             ):
    lattice_df = lattice.read_lattice_file(lattice_file)
    cameras_df = lattice.get_cameras_df(
        lattice_df,
        camera_sizes_X_um,
        camera_sizes_Y_um)
    ey_um, ey_err = lattice.get_e_um_Y_scipy_curve_fit(cameras_df)
    popt, perr = lattice.get_e_um_X_scipy_curve_fit(
        cameras_df, dpp, dpp_err)
    ex_um, dpp = popt
    ex_err, dpp_err = perr
    emittance_6D = {
        "ex_um": ex_um,
        "ex_err_um": ex_err,
        "ey_um": ey_um,
        "ey_err_um": ey_err,
        "dp/p": dpp,
        "dp/p_err": dpp_err
    }
    undulator_df = lattice.get_undulator_df(lattice_df, emittance_6D)
    return undulator_df, emittance_6D


def CalcTransverseBeamParams(lattice_df, ex_um, ey_um, dpp):
    emittance_6D = {
        'ex_um': ex_um,
        'ex_err': 0,
        'ey_um': ey_um,
        'ey_err': 0,
        'dp/p': dpp,
        'dp/p_err': 0
    }
    und_summary = lattice.get_undulator_df(lattice_df, emittance_6D)
    #print(und_summary)
    v = und_summary.loc['Middle', ['Beta_cm_X', 'Beta_cm_Y',
                                   'Alpha_X', 'Alpha_Y',
                                   'Angle_spread_rad_X', 'Angle_spread_rad_Y',
                                   'ex_um', 'ey_um',
                                   'Dispersion_cm_X', 'dDx/dS',
                                   'dp/p',
                                   'Sigma_um_X', 'Sigma_um_Y']]
    bx, by, ax, ay, sxp, syp, ex, ey, Dx, Dxp, sp, sx, sy = v
    bx, by, Dx = 1e4*np.array([bx, by, Dx])
    gx, gy = (1+ax**2)/bx, (1+ay**2)/by
    #print("sqrt(gx*ex) = ", np.sqrt(gx*ex))
    #print("sqrt(gy*ey) = ", np.sqrt(gy*ey))
    Sx = np.sqrt(ex/gx+(gx*Dx+Dxp*ax)**2*bx*ex*sp**2/sxp**2)
    #print("sqrt(ex/gx) = ", np.sqrt(ex/gx))
    Sy = np.sqrt(ey/gy)
    dx = (ax*ex-Dx*Dxp*sp**2)/sxp**2
    #print("ax*ex/sxp**2 = ", ax*ex/sxp**2)
    #print("(Dx*Dxp*sp**2)/sxp**2 = ", (Dx*Dxp*sp**2)/sxp**2)
    dy = ay/gy
    return Sx, Sy, dx, dy, sxp, syp
