import numpy as np


def read_waveform(file_path):
    vals = np.fromfile(file_path, dtype=np.float32)
    ch1 = vals[::2]
    ch2 = vals[1::2]
    return ch1, ch2
