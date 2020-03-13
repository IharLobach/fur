import os
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

module_directory = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(module_directory, "spectral_transmission_data")

mirrors_file = os.path.join(data_dir, "mirrors.csv")
lens_file = os.path.join(data_dir, "lens.csv")
photodiode_file = os.path.join(data_dir, "hamamatsu_photodiode_G11193.csv")

mirrors_df = pd.read_csv(mirrors_file)
mirrors_df['reflectance'] = mirrors_df['reflectance_percent']*0.01
mirrors_df['wavelength_um'] = mirrors_df['wavelength_nm']*0.001
lens_df = pd.read_csv(lens_file)
lens_df['transmission'] = 1-lens_df['reflectance_percent']*0.01
lens_df['wavelength_um'] = lens_df['wavelength_nm']*0.001
photodiode_df = pd.read_csv(photodiode_file)
photodiode_df['quantum_efficiency'] =\
    photodiode_df['photo_sensitivity_A_per_W']*1.24\
    / photodiode_df['wavelength_um']

mirror_reflectance = interp1d(mirrors_df['wavelength_um'],
                              mirrors_df['reflectance'],
                              bounds_error=False,
                              fill_value=(0, 0))
lens_transmission = interp1d(lens_df['wavelength_um'],
                             lens_df['transmission'],
                             bounds_error=False,
                             fill_value=(0, 0))
photodiode_quantum_efficiency = interp1d(photodiode_df['wavelength_um'],
                                         photodiode_df['quantum_efficiency'],
                                         bounds_error=False,
                                         fill_value=(0, 0))


def transmission_function(lambda_um):
    return mirror_reflectance(lambda_um)**2\
        * lens_transmission(lambda_um)\
        * photodiode_quantum_efficiency(lambda_um)
