import numpy as np
from wiggler_radiation.wiggler_radiation import get_photon_flux_3D
from wiggler_radiation.Wigrad.wigrad_generator import get_rad_mesh_tuple


def get_My(sigma_y_um, rad_mesh_tuple, photon_flux_3D):
    """Returns 'effective' My, to get real M (number of coherent modes)
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
