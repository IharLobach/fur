import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
from fur.waveform_reader import read_waveform
from fur.finding_period import get_period
from config import get_from_config
import os


def analyze_one_file_one_channel(file_path, show_plot=False, start_to_trigger=200,
                     show_chosen_period_start=False, ouput_dic=None):
    dt = get_from_config("dt")
    Rf = get_from_config("Rf")
    e = get_from_config("e")
    p = file_path
    ch = read_waveform(p, one_channel=True)
    output = {}
    period = get_period(ch, output_dic=output)
    trig_times = output["trig_times"]
    start_of_first_full_period = int(trig_times[0]-start_to_trigger)
    if show_chosen_period_start:
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(ch[:int(period)])
        ax.axvline(start_of_first_full_period)
        ax.set_title("Chosen period start")
        ax.text(0.0, 0.95, os.path.basename(file_path),
                verticalalignment='bottom', horizontalalignment='left',
                transform=ax.transAxes, fontsize=10)
        plt.show()
    end_of_last_full_period = int(trig_times[-1]-start_to_trigger)
    ch_df = pd.DataFrame({
        "ch": ch[start_of_first_full_period:end_of_last_full_period]})
    ch_df["period_number"] = (ch_df.index.to_numpy()/period).astype(int)
    photoelectrons_df = ch_df.groupby(['period_number']).sum()/Rf/e*dt*1e-9
    photoelectrons_df = photoelectrons_df.rename(
        columns={"ch": "n_photoelectrons"})
    if ouput_dic is not None:
        ouput_dic["photoelectrons_df"] = photoelectrons_df
    desc = photoelectrons_df.describe().n_photoelectrons

    def remove_last_line(s):
        return s[:s.rfind('\n')]

    def gauss(x, a, mu, sigma):
        return a*np.exp(-(x-mu)**2/(2*sigma**2))

    if show_plot:
        fontsize = 16
        fig, ax = plt.subplots(figsize=[15, 7.5])
        distplot = sns.distplot(photoelectrons_df, kde=False, ax=ax)
        ax.set_xlabel("Photoelectron count", fontsize=fontsize)
        ax.set_ylabel("Occurrences (number of revolutions)", fontsize=fontsize)
        bin_centers = [h.get_x()+h.get_width()/2 for h in distplot.patches]
        bars = [h.get_height() for h in distplot.patches]
        success = False
        try:
            popt, pcov = curve_fit(gauss, bin_centers, bars,
                                   p0=(max(bars), desc['mean'], desc['std']))
            success = True
        except Exception as e:
            print("Exception happened while trying to fit"
                  "the distribution with a Gaussian: ", e)
        if success:
            ax.plot(bin_centers, gauss(bin_centers, *popt), 'ro:',
                    label='Gaussian fit')
            ax.legend(fontsize=fontsize)
        ax.text(0.9, 0.5, remove_last_line(desc.__repr__()),
                verticalalignment='bottom', horizontalalignment='right',
                transform=ax.transAxes, fontsize=15)
        ax.text(0.0, 0.95, os.path.basename(file_path),
                verticalalignment='bottom', horizontalalignment='left',
                transform=ax.transAxes, fontsize=15)
        ax.set_title("Distribution of photoelectron count per period")
        plt.show()
    return desc


def analyze_one_file(file_path, show_plot=False, start_to_trigger=200,
                     integration_end_ns=50, ch2_periods=10,
                     max_len=12100,
                     ouput_dic=None):
    dt = get_from_config("dt")
    ch1, ch2 = read_waveform(file_path)
    output = {}
    period = get_period(ch2, output_dic=output)
    trig_times = output["trig_times"]
    int_period = int(period)
    start_of_first_full_period = int(trig_times[0]-start_to_trigger)
    end_of_last_full_period = int(trig_times[-1]-start_to_trigger)
    ch_df = pd.DataFrame({
        "ch1": ch1[start_of_first_full_period:end_of_last_full_period],
        "ch2": ch2[start_of_first_full_period:end_of_last_full_period]})
    ch_df["period_number"] = (ch_df.index.to_numpy()/period).astype(int)
    integration_end = int(integration_end_ns/dt)
    time_arr = dt*np.arange(int_period+1)
    if show_plot:
        fig, ax1 = plt.subplots()

        color = 'tab:red'
        ax1.set_xlabel('time, ns')
        ax1.set_ylabel('Difference channel, V', color=color)
        ax1.plot(time_arr, ch_df.loc[:int_period, "ch1"], color=color)
        ax1.tick_params(axis='y', labelcolor=color)

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        color = 'tab:blue'
        # we already handled the x-label with ax1
        ax2.set_ylabel('sSum channel, V', color=color)
        ax2.plot(time_arr, ch_df.loc[:int_period, "ch2"], color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        ax1.axvline(0, color='red')
        ax1.axvline(integration_end_ns, color='red')

        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.show()
    ch1_int = ch_df.groupby(["period_number"])['ch1'].apply(np.asarray).apply(
        lambda ch1: np.sum(ch1[:integration_end])).to_list()
    ch2_int = ch_df[ch_df["period_number"] < ch2_periods] \
        .groupby(["period_number"])['ch2'].apply(np.asarray) \
        .apply(lambda ch2: np.sum(ch2[:integration_end])).values
    nan_len = max_len-1-len(ch1_int)
    return [np.mean(ch2_int)]+ch1_int+[np.nan for _ in range(nan_len)]
