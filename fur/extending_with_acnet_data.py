import datetime
import time
import requests
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
from config import get_from_config, save_to_config
from acnet_reader.acnet_reader import fetch_data, get_interpolated_df,\
    fetch_interpolated_data
import fur.path_assistant as path_assistant
from lattice.summary_in_undulator import get_summary_in_undulator
import lattice.lattice as lattice

cameras = lattice.acnet_devices_X+lattice.acnet_devices_Y

bunch_parameters = ["Sigma_um_X", "Sigma_um_X_err",
                    "Sigma_um_Y", "Sigma_um_Y_err",
                    "Angle_spread_rad_X", "Angle_spread_rad_Y",
                    "ex_um", "ex_err", "ey_um", "ey_err", "dp/p", "dp/p_err"]

plt.rcParams['figure.figsize'] = [20, 7.5]
plt.rcParams.update({'font.size': 16, 'legend.fontsize': 16})


def extend_fluctuations_df_with_acnet_data(
        fluctuations_df,
        bpm_files_df,
        acnet_data_df,
        show_plot=False):
    fluctuations_df.sort_values("file_datetime")
    bpm_files_df.sort_values("file_datetime")
    if show_plot:
        ax = sns.lineplot(x=acnet_data_df.index, y=acnet_data_df["N:IWCMI"])
        ax.set_ylim(-3.5, -2.0)
        for t in fluctuations_df["file_datetime"]:
            plt.axvline(t, color="brown")
        for t in bpm_files_df["file_datetime"]:
            plt.axvline(t, color="green")
        plt.show()
    fluctuations_df["rounded_bpm_file_datetime"] =\
        bpm_files_df["file_datetime"].apply(lambda t: t.round(freq='S'))
    acnet_addition = acnet_data_df\
        .loc[acnet_data_df.index
             .isin(fluctuations_df["rounded_bpm_file_datetime"])]\
        .loc[:,
             ["N:IWCMI", "N:IBEAMA", "N:IWCMBE",
              "N:IWCMBR", "N:IWCMBF", "N:IWCMBG", "N:IRFEPA"]
             + cameras]
    fluctuations_df = pd.concat(
        [fluctuations_df, acnet_addition.reset_index(drop=True)], axis=1)
    fluctuations_df["N:IWCMI_recalibrated_to_IWCMI_absolute"] =\
        get_from_config("IWCMI_to_WCM_ABSOLUTE")\
        * fluctuations_df["N:IWCMI"]
    fluctuations_df["N:IBEAM_recalibrated_to_IWCMI_absolute"] =\
        get_from_config("IBEAM_to_WCM_ABSOLUTE")\
        * fluctuations_df["N:IBEAMA"]
    return fluctuations_df


def extend_fluctuations_df_with_bunch_size(fluctuations_df, lattice_file,
                                           fit_dpp=False):
    def get_bunch_params_in_middle_of_undulator(row):
        if fit_dpp:
            dpp = None
        else:
            dpp = lattice.get_dpp(row["N:IWCMBR"], row["N:IRFEPA"])
        undulator_df, emittance_6D = get_summary_in_undulator(
            lattice_file,
            row[lattice.acnet_devices_X].values,
            row[lattice.acnet_devices_Y].values,
            dpp,
            dpp_err=0
        )
        return undulator_df.loc["Middle", bunch_parameters]

    bunch_sizes_df = fluctuations_df.apply(
        get_bunch_params_in_middle_of_undulator, axis=1, result_type='expand')
    return pd.concat([fluctuations_df, bunch_sizes_df], axis=1)


def get_fluctuations_df_with_acnet_data(
        shift,
        fluctuations_df_file_name,
        lattice_file_name,
        acnet_data_df_file_name=None,
        fit_dpp=False):
    results_dir = shift.get_results_dir()
    bpm_files_df = shift.get_bpm_files_df()
    if acnet_data_df_file_name is None:
        acnet_fn = "all_acnet_data_for_"+shift.shift_folder_name+".csv"
    else:
        acnet_fn = acnet_data_df_file_name
    acnet_data_df = shift.get_acnet_data_df(acnet_fn)
    fluctuations_df = shift.get_fluctuations_df(
        fluctuations_df_file_name)
    lattice_file = shift.get_6dsim_dir().fi(
        lattice_file_name)

    res_df = extend_fluctuations_df_with_acnet_data(
        fluctuations_df, bpm_files_df, acnet_data_df)
    res_df = extend_fluctuations_df_with_bunch_size(res_df, lattice_file,
                                                    fit_dpp)
    return res_df