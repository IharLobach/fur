import numpy as np
import scipy.signal
from scipy.optimize import minimize
import pandas as pd
from sklearn import linear_model
import matplotlib.pyplot as plt
from datetime import datetime
import seaborn as sns
import sys
import os
import fur.path_assistant as path_assistant
from fur.waveform_reader import read_waveform
from fur.finding_period import get_period
from fur.fluctuations import get_fluctiation_and_noise_var


def analyze_one_dataset(shift_folder, base_name, t1=None, t2=None,
                        period_in=None):
    results_dir = shift_folder.get_results_dir()
    if (t1 is None) and (t2 is None):
        wf_paths = shift_folder.get_waveform_paths()
    else:   
        wf_paths = [p for p in shift_folder.get_waveform_paths() if
                    (t1 < shift_folder.get_datetime(os.path.basename(p)) < t2)]
    n_files = len(wf_paths)
    print("There are {} files in this data set.".format(n_files))
    res_df = pd.DataFrame(columns=["waveform_file",
                                   "ch2_amplitude",
                                   "var_of_ch1_amplitude",
                                   "noise_var"],
                          index=np.arange(n_files))
    res_df["waveform_file"] = [os.path.basename(p) for p in wf_paths]
    for i, p in enumerate(wf_paths):
        status = os.path.basename(p)+" ({}/{})".format(i+1, n_files)
        print("Started working on the file ", status)
        try:
            ch1, ch2 = read_waveform(p)
            if period_in is None:
                period = get_period(ch2)
            else:
                period = period_in
            print("period = {}".format(period))
            res_df.iloc[i, 1:] = get_fluctiation_and_noise_var(ch1, ch2,
                                                               period,
                                                               show_plots=True)
            print("Sum amplitude = {:.3} V".format(res_df.iloc[i, 1]))
        except Exception as e:
            print("Exception happened: ", e)
        print("Finished working on ", status)
    res_df["waveform_file"] = res_df["waveform_file"].astype('str')
    for i in range(1, 4):
        res_df.iloc[:, i] = res_df.iloc[:, i].astype(np.float32)
    if not os.path.exists(shift_folder.shift_results_dir):
        os.mkdir(shift_folder.shift_results_dir)
    dumped_file_name = "res_df_"+base_name+".csv"
    res_df.to_csv(results_dir.fi(dumped_file_name))
    print("Results saved to {}".format(results_dir.fi(dumped_file_name)))
    plt.plot(res_df["ch2_amplitude"], res_df["var_of_ch1_amplitude"], '.')
    plt.xlabel("Sum channel amplitude, V")
    plt.ylabel("Variance of the difference channel amplitude, V$^2$")
    plt.show()
    return res_df
