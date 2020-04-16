import numpy as np
import scipy.signal
from scipy.optimize import minimize, curve_fit
import pandas as pd
from sklearn import linear_model
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [15, 7.5]
plt.rcParams.update({'font.size': 16, 'legend.fontsize': 16})
import statsmodels.api as sm


def get_fluctiation_and_noise_var(ch1, ch2, period, n_bins=1000,
                                  show_plots=False, fit_method='lstsq',
                                  output_dic=None):
    data_len = len(ch1)
    times_within_period = np.arange(data_len) % period
    bin_size = period/n_bins
    df = pd.DataFrame({"time": times_within_period,
                       "ch1": ch1,
                       "ch2": ch2})
    df0 = df.groupby(df["time"]//bin_size)
    df1 = df0['ch1'].apply(np.asarray)
    var_ch1 = np.array([np.var(arr) for arr in df1])
    ch2_mean = df0["ch2"].apply(np.mean)
    ch2_mean = (ch2_mean-min(ch2_mean)).values
    if output_dic is not None:
        output_dic["std_of_ch1_var"] = np.std(var_ch1)


    # least squares:
    ch2_mean_squared = ch2_mean**2
    max_ch1_mean_squared = max(ch2_mean_squared)
    normalized_ch2_mean_squared = ch2_mean_squared/max_ch1_mean_squared
    # reg = linear_model.LinearRegression()
    # x = normalized_ch2_mean_squared
    # y = var_ch1
    # reg.fit(x.reshape((len(x), 1)), y)
    # fluct_V_f = reg.coef_[0]
    # noise_var_f = reg.intercept_
    X = normalized_ch2_mean_squared.reshape(-1, 1)
    X = sm.add_constant(X)
    model = sm.OLS(var_ch1, X)
    results = model.fit()
    noise_var_f, fluct_V_f = results.params
    noise_var_f_err, fluct_V_f_err = results.bse
    if fit_method == 'scipy_curve_fit':
        max_ch2_mean = max(ch2_mean)
        normalized_ch2_mean = ch2_mean/max_ch2_mean

        def func(x, var_noise, fluc):
            return var_noise+fluc*x**2

        popt, pcov = curve_fit(func, normalized_ch2_mean, var_ch1,
                               p0=(fluct_V_f, noise_var_f))
        oise_var_f, fluct_V_f = popt
        noise_var_f_err, fluct_V_f_err = np.sqrt(np.diag(pcov))
    res = np.sqrt(max_ch1_mean_squared), fluct_V_f, noise_var_f,\
        fluct_V_f_err, noise_var_f_err
    df = pd.DataFrame(
        {"Value": res},
        index=["Sum channel amplitude, V", "Difference chanel variance, V^2",
               "Noise variance, V^2", "Difference channel variance error, V^2",
               "Noise variance error, V^2"])

    if show_plots:
        fig, ax = plt.subplots(figsize=[15, 7.5])
        ax.plot(var_ch1)
        ax.plot(noise_var_f+fluct_V_f*normalized_ch2_mean_squared)
        ax.set_xlabel("Bin number within one revolution period")
        ax.set_ylabel("Ch1 variance in the bin and fit with Ch2^2, V^2")
        ax.text(0.9, 0.5, df.__repr__(),
                verticalalignment='bottom', horizontalalignment='right',
                transform=ax.transAxes, fontsize=15)
        plt.show()
    return res


def get_average_signal(ch2, period, n_bins=1000):
    data_len = len(ch2)
    times_within_period = np.arange(data_len) % period
    bin_size = period/n_bins
    df = pd.DataFrame({"time": times_within_period,
                       "ch2": ch2})
    df0 = df.groupby(df["time"]//bin_size)
    ch2_mean = df0["ch2"].apply(np.mean)
    return ch2_mean.values
