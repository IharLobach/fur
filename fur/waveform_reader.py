import numpy as np
import pandas as pd


def read_waveform(file_path, one_channel=False, csv=False, nrows=None):
    if csv:
        df = pd.read_csv(file_path, sep=';',
                         header=None, nrows=nrows)
        return df.iloc[:, 1], df.iloc[:, 0]
    else:
        vals = np.fromfile(file_path, dtype=np.float32)
        if one_channel:  # one channel (without the comb filter, directly connected)
            return -vals
        else:  # two channels (with the comb filter)
            ch1 = vals[::2]
            ch2 = vals[1::2]
            return ch1, ch2
