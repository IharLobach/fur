import numpy as np
import pandas as pd
import os
from scipy.interpolate import RegularGridInterpolator
from wiggler_radiation.wiggler_radiation import get_photon_flux_3D
from wiggler_radiation.Wigrad.wigrad_generator import get_rad_mesh_tuple
from config import get_from_config
import re


def get_My(sigma_y_um, rad_mesh_tuple, photon_flux_3D):
    """Legacy. Returns 'effective' My, to get real M (number of coherent modes)
    multiply this by sigma_x_um and sigma_z_um.
    """
    i_3D = photon_flux_3D
    x_1D, y_1D, l_1D = rad_mesh_tuple
    sigma_y = sigma_y_um
    x_step = (x_1D[-1]-x_1D[0])/(len(x_1D)-1)
    y_step = (y_1D[-1]-y_1D[0])/(len(y_1D)-1)
    l_step = (l_1D[-1]-l_1D[0])/(len(l_1D)-1)
    tot = x_step*y_step*l_step*np.sum(i_3D)
    y1_2D, y2_2D = np.meshgrid(y_1D, y_1D)
    y1m2_3D = np.tile(y1_2D-y2_2D, (len(l_1D), 1, 1))
    k = 2*np.pi/l_1D
    exp = (l_1D**3)[:, None, None] * \
        np.exp(-(sigma_y*k[:, None, None]*y1m2_3D)**2)
    aux = np.einsum('...ji,...mi->...jm', i_3D, i_3D)

    res = x_step*y_step**2*l_step*np.sum(aux*exp)
    My = 4*np.pi*tot**2/res
    return My


def get_M_df(en):
    M_df = pd.read_csv(
        os.path.join(os.path.dirname(__file__),
                     "rcc_midway_precalculated",
                     f"Mxy_Ver_SigmaX_Hor_SigmaY_{en}MeV.csv"),
        index_col=0)
    return M_df


def get_M_interpolator_sxsyEn():
    this_dir = os.path.dirname(os.path.abspath(__file__))
    files = os.listdir(os.path.join(this_dir, 'rcc_midway_precalculated'))
    file_paths = [os.path.join(this_dir, 'rcc_midway_precalculated', fn)
                  for fn in files]

    def get_energy_from_name(fn):
        return int(re.match(
            'Mxy_Ver_SigmaX_Hor_SigmaY_([0-9]*)MeV.csv', fn).groups()[0])

    energies = np.asarray([get_energy_from_name(fn) for fn in files])
    energies = np.sort(energies)
    M_dfs = [get_M_df(en) for en in energies]
    sy_range = M_dfs[0].columns.values.astype(float)
    sx_range = M_dfs[0].index.values.astype(float)
    M_3d = np.asarray([m.values for m in M_dfs])
    Mxy_interpolator = RegularGridInterpolator(
        (energies, sx_range, sy_range), M_3d,
        bounds_error=False, fill_value=None)
    energies_der = (energies[1:]+energies[:-1])/2
    sx_range_der = (sx_range[1:]+sx_range[:-1])/2
    sy_range_der = (sy_range[1:]+sy_range[:-1])/2
    e_step = (energies[-1]-energies[0])/(len(energies)-1)
    sx_step = (sx_range[-1]-sx_range[0])/(len(sx_range)-1)
    sy_step = (sy_range[-1]-sy_range[0])/(len(sy_range)-1)
    Mxy_der_e_values = (M_3d[1:, :, :]-M_3d[:-1, :, :])/e_step
    Mxy_der_x_values = (M_3d[:, 1:, :]-M_3d[:, :-1, :])/sx_step
    Mxy_der_y_values = (M_3d[:, :, 1:]-M_3d[:, :, :-1])/sy_step
    Mxy_der_e_interpolator = RegularGridInterpolator(
        (energies_der, sx_range, sy_range), Mxy_der_e_values,
        bounds_error=False, fill_value=None)
    Mxy_der_x_interpolator = RegularGridInterpolator(
        (energies, sx_range_der, sy_range), Mxy_der_x_values,
        bounds_error=False, fill_value=None)
    Mxy_der_y_interpolator = RegularGridInterpolator(
        (energies, sx_range, sy_range_der), Mxy_der_y_values,
        bounds_error=False, fill_value=None)

    def M(sigma_x_um, sigma_y_um, sigma_z_cm, en_MeV):
        return Mxy_interpolator((en_MeV, sigma_x_um, sigma_y_um))*sigma_z_cm

    def M_der_x(sigma_x_um, sigma_y_um, sigma_z_cm, en_MeV):
        return Mxy_der_x_interpolator((en_MeV, sigma_x_um, sigma_y_um))\
            * sigma_z_cm

    def M_der_y(sigma_x_um, sigma_y_um, sigma_z_cm, en_MeV):
        return Mxy_der_y_interpolator((en_MeV, sigma_x_um, sigma_y_um))\
            * sigma_z_cm

    def M_der_z(sigma_x_um, sigma_y_um, sigma_z_cm, en_MeV):
        return Mxy_interpolator((en_MeV, sigma_x_um, sigma_y_um))

    def M_der_e(sigma_x_um, sigma_y_um, sigma_z_cm, en_MeV):
        return Mxy_der_e_interpolator((en_MeV, sigma_x_um, sigma_y_um))\
            * sigma_z_cm

    return M, M_der_x, M_der_y, M_der_z, M_der_e


def get_M_interpolator_at_fixed_energy(energy_MeV=None):
    """Returns interpolators M(sigma_x_um, sigma_y_um, sigma_z_cm),
    dM/dsigma_x_um, dM/dsigma_y_um, dM/dsigma_z_cm"""
    en = energy_MeV if (energy_MeV is not None) else\
        get_from_config("me_MeV")*get_from_config("gamma")
    M0, Mx, My, Mz, Me = get_M_interpolator_sxsyEn()

    def M(sigma_x_um, sigma_y_um, sigma_z_cm):
        return M0(sigma_x_um, sigma_y_um, sigma_z_cm, en)

    def M_der_x(sigma_x_um, sigma_y_um, sigma_z_cm):
        return Mx(sigma_x_um, sigma_y_um, sigma_z_cm, en)

    def M_der_y(sigma_x_um, sigma_y_um, sigma_z_cm):
        return My(sigma_x_um, sigma_y_um, sigma_z_cm, en)

    def M_der_z(sigma_x_um, sigma_y_um, sigma_z_cm):
        return Mz(sigma_x_um, sigma_y_um, sigma_z_cm, en)

    def M_der_e(sigma_x_um, sigma_y_um, sigma_z_cm):
        return Me(sigma_x_um, sigma_y_um, sigma_z_cm, en)

    return M, M_der_x, M_der_y, M_der_z, M_der_e



def get_M_interpolator(energy_MeV=None):
    """Legacy. Returns interpolators M(sigma_x_um, sigma_y_um, sigma_z_cm),
    dM/dsigma_x_um, dM/dsigma_y_um, dM/dsigma_z_cm"""
    en = energy_MeV if (energy_MeV is not None) else\
        get_from_config("me_MeV")*get_from_config("gamma")
    M_df = pd.read_csv(
        os.path.join(os.path.dirname(__file__),
                     "rcc_midway_precalculated",
                     f"Mxy_Ver_SigmaX_Hor_SigmaY_{en}MeV.csv"),
        index_col=0)
    sy_range = M_df.columns.values.astype(float)
    sx_range = M_df.index.values.astype(float)
    Mxy_values = M_df.values
    Mxy_interpolator = RegularGridInterpolator((sx_range, sy_range),
                                               Mxy_values)
    sx_range_der = (sx_range[1:]+sx_range[:-1])/2
    sy_range_der = (sy_range[1:]+sy_range[:-1])/2
    sx_step = (sx_range[-1]-sx_range[0])/(len(sx_range)-1)
    sy_step = (sy_range[-1]-sy_range[0])/(len(sy_range)-1)
    Mxy_der_x_values = (Mxy_values[1:, :]-Mxy_values[:-1, :])/sx_step
    Mxy_der_y_values = (Mxy_values[:, 1:]-Mxy_values[:, :-1])/sy_step
    Mxy_der_x_interpolator = RegularGridInterpolator(
        (sx_range_der, sy_range), Mxy_der_x_values)
    Mxy_der_y_interpolator = RegularGridInterpolator(
        (sx_range, sy_range_der), Mxy_der_y_values)

    def M(sigma_x_um, sigma_y_um, sigma_z_cm):
        return Mxy_interpolator((sigma_x_um, sigma_y_um))*sigma_z_cm

    def M_der_x(sigma_x_um, sigma_y_um, sigma_z_cm):
        return Mxy_der_x_interpolator((sigma_x_um, sigma_y_um))*sigma_z_cm

    def M_der_y(sigma_x_um, sigma_y_um, sigma_z_cm):
        return Mxy_der_y_interpolator((sigma_x_um, sigma_y_um))*sigma_z_cm

    def M_der_z(sigma_x_um, sigma_y_um, sigma_z_cm):
        return Mxy_interpolator((sigma_x_um, sigma_y_um))

    return M, M_der_x, M_der_y, M_der_z
