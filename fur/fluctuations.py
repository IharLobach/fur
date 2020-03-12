import numpy as np
import scipy.signal
from scipy.optimize import minimize
import pandas as pd
from sklearn import linear_model
import matplotlib.pyplot as plt


def get_fluctiation_and_noise_var(ch1, ch2, period, n_bins=1000,
                                  show_plots=False):
    data_len = len(ch1)
    times_within_period = np.arange(data_len) % period
    bin_size = period/n_bins
    df = pd.DataFrame({"time": times_within_period,
                       "ch1": ch1,
                       "ch2": ch2})
    df0 = df.groupby(df["time"]//bin_size)
    df1 = df0['ch1'].apply(np.asarray)
    var_ch1 = np.array([np.var(arr) for arr in df1])
    ch1_mean = df0["ch2"].apply(np.mean)
    ch1_mean = ch1_mean-min(ch1_mean)
    ch1_mean_squared = ch1_mean.values**2
    max_ch1_mean_squared = max(ch1_mean_squared)
    normalized_ch1_mean_squared = ch1_mean_squared/max_ch1_mean_squared
    reg = linear_model.LinearRegression()
    x = normalized_ch1_mean_squared
    y = var_ch1
    reg.fit(x.reshape((len(x), 1)), y)
    fluct_V_f = reg.coef_[0]
    noise_var_f = reg.intercept_
    if show_plots:
        plt.plot(var_ch1)
        plt.plot(noise_var_f+fluct_V_f*normalized_ch1_mean_squared)
        plt.xlabel("Bin number within one revolution period")
        plt.ylabel("Ch1 variance in the bin and fit with Ch2^2, V^2")
        plt.show()
    return np.sqrt(max_ch1_mean_squared), fluct_V_f, noise_var_f

