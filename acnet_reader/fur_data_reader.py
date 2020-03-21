from datetime import datetime
import time
import requests
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
from acnet_reader.acnet_reader import fetch_data, get_interpolated_df,\
    fetch_interpolated_data
import fur.path_assistant as path_assistant


def save_acnet_data_for_fur(shift, t1, t2, file_name_to_save):
    acnet_data_dir = shift.get_acnet_data_dir()
    cameras = ["N:ITC1RSV", "N:ITC1RSH",
               "N:ITC2RSV", "N:ITC2RSH",
               "N:ITC3RSV", "N:ITC3RSH",
               "N:ITC4RSV", "N:ITC4RSH",
               "N:ITC4LSV", "N:ITC4LSH",
               "N:ITC3LSV", "N:ITC3LSH",
               "N:ITC2LSV", "N:ITC2LSH",
               "N:ITC1LSV", "N:ITC1LSH"]
    synclight_data = fetch_interpolated_data(t1, t2, cameras)
    wcm_devices = ["N:IWCMBF", "N:IWCMBR", "N:IWCMBP",
                   "N:IWCMI", "N:IRFEPA", "N:IRFEPP", "N:IWCMBE",
                   "N:IWCMBM", "N:IWCMBG", "N:IWCMIG"]
    wcm_data = fetch_interpolated_data(t1, t2, wcm_devices)
    ibeam_data = fetch_interpolated_data(t1, t2, "N:IBEAMA")
    rf_phase_data = fetch_interpolated_data(t1, t2, "N:IRFEPC")
    rf_devices = ["N:IRFEAT", "N:IRFEFP", "N:IRFECG"]
    rf_data = fetch_interpolated_data(t1, t2, rf_devices)
    all_data = pd.concat([synclight_data, wcm_data, ibeam_data,
                                     rf_phase_data, rf_data], axis=1)
    all_data.to_csv(acnet_data_dir.fi(file_name_to_save))
