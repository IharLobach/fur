import numpy as np
import fur.path_assistant as path_assistant
from config import get_from_config
from wiggler_radiation.transmission_data import transmission_function
from wigrad import Wiggler, WigglerRadiationSimulator


def get_rad_mesh_tuple(config_style_mesh=None):
    if config_style_mesh is None:
        rad_mesh = get_from_config("radiation_mesh")
    else:
        rad_mesh = config_style_mesh
    zobs = get_from_config("z_obs_m")
    theta_xs = rad_mesh[0][0]/zobs
    theta_xf = rad_mesh[0][1]/zobs
    theta_ys = rad_mesh[1][0]/zobs
    theta_yf = rad_mesh[1][1]/zobs
    xbins = rad_mesh[0][2]
    ybins = rad_mesh[1][2]
    ls = rad_mesh[2][0]
    lf = rad_mesh[2][1]
    lbins = rad_mesh[2][2]
    mesh = (np.linspace(theta_xs, theta_xf, xbins),
            np.linspace(theta_ys, theta_yf, ybins),
            np.linspace(ls, lf, lbins))
    return mesh


def generate_wr_sim_for_plotting(config_style_mesh=None):
    wiggler = Wiggler(K_peak=get_from_config("K_peak"))
    mesh = get_rad_mesh_tuple(config_style_mesh)
    spectral_transmission = transmission_function(mesh[2])
    wr_sim = WigglerRadiationSimulator(
        wiggler,
        mesh,
        gamma=get_from_config("gamma"),
        harmonics=[1, 2],
        aperture=None,  # 'ellipse',
        # if False, then both polarizations are calculated separately
        only_calc_sum_of_both_polarizations=False,  # True,
        spectral_transmission=None  # spectral_transmission
    )
    return wr_sim


def generate_wr_sim_with_wigrad_results(config_style_mesh=None):
    wiggler = Wiggler(K_peak=get_from_config("K_peak"))
    mesh0 = get_rad_mesh_tuple(config_style_mesh)
    mesh = (mesh0[0][int(len(mesh0[0])/2):],
            mesh0[1][int(len(mesh0[1])/2):],
            mesh0[2])
    spectral_transmission = transmission_function(mesh[2])
    wr_sim = WigglerRadiationSimulator(
        wiggler,
        mesh,
        gamma=get_from_config("gamma"),
        harmonics=[1, 2],
        aperture='ellipse',
        # if False, then both polarizations are calculated separately
        only_calc_sum_of_both_polarizations=True,
        spectral_transmission=spectral_transmission
    )
    wr_sim.calc_photon_flux_on_meshgrid()
    wr_sim.extend_results_using_symmetries()
    return wr_sim
