import numpy as np


def read_waveform(file_path, one_channel=False):
    vals = np.fromfile(file_path, dtype=np.float32)
    if one_channel:  # one channel (without the comb filter, directly connected)
        return -vals
    else:  # two channels (with the comb filter)
        ch1 = vals[::2]
        ch2 = vals[1::2]
        return ch1, ch2
