from wiggler_radiation.Wigrad.wigrad_generator \
    import generate_wr_sim_with_wigrad_results
import numpy as np
import fur.path_assistant as path_assistant
from config import get_from_config


def get_photon_flux_3D(source='wigrad', config_style_mesh=None):
    if source == 'wigrad':
        return generate_wr_sim_with_wigrad_results(config_style_mesh)\
            .get_photon_flux_3D(polarization='sum')
    elif source == 'SRW':
        Ex_3D = np.load(path_assistant.srw_Ex_3D_file_path)
        Ey_3D = np.load(path_assistant.srw_Ey_3D_file_path)
        return np.absolute(Ex_3D)**2+np.absolute(Ey_3D)**2
    else:
        raise ValueError("Unknown source type, choose from "
                         "'wigrad' and 'SRW'")
