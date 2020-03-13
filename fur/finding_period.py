import numpy as np
import scipy.signal
from sklearn import linear_model


def get_fitlered_signal(signal, window_length=101, polyorder=3):
    return scipy.signal.savgol_filter(signal, window_length, polyorder)


def get_trig_times(signal, abolute_trigger_level):
    above_tl = np.where(signal > abolute_trigger_level, 1, 0)
    dif = np.diff(above_tl)
    pos_edges = np.where(dif > 0, 1, 0)
    ts_dif = np.diff(signal)
    timeline = np.arange(len(signal) - 1)
    interpolated_time = timeline - (signal[1:] - abolute_trigger_level)\
        / np.where(ts_dif != 0, ts_dif, 10)
    trig_times = interpolated_time[pos_edges > 0] + 1
    return trig_times


def get_reduced_data(signal, resampling_factor=5):
    return signal[::resampling_factor]


def get_absolute_trig_level_from_relative(signal, relative_trigger_level=0.5):
    top = max(signal)
    bottom = min(signal)
    return bottom+relative_trigger_level*(top-bottom)


def get_period_from_trig_times(trig_times):
    """returns period in units of sampling time of the argument trig_times"""
    reg = linear_model.LinearRegression()
    x = np.arange(len(trig_times))
    y = trig_times
    reg.fit(x.reshape((len(x), 1)), y)
    a = reg.coef_[0]
    b = reg.intercept_
    return a


def get_period(signal, sampling_time=1, relative_trigger_level=0.5,
               filter_window_length=101, filter_polyorder=3,
               resampling_factor=5):
    filtered_signal = get_fitlered_signal(signal, filter_window_length,
                                          filter_polyorder)
    reduced_data = get_reduced_data(filtered_signal, resampling_factor)
    absolute_trig_level = \
        get_absolute_trig_level_from_relative(reduced_data,
                                              relative_trigger_level)
    trig_times = get_trig_times(reduced_data, absolute_trig_level)
    period = get_period_from_trig_times(trig_times)
    return resampling_factor*sampling_time*period