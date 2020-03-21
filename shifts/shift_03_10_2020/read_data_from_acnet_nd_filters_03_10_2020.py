import os
from datetime import datetime
import time
import requests
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from acnet_reader.acnet_reader import fetch_data, get_interpolated_df,\
    fetch_interpolated_data
import fur.path_assistant as path_assistant
from acnet_reader.fur_data_reader import save_acnet_data_for_fur
shift_03_10_2020 = path_assistant.PathAssistant('shift_03_10_2020')
waveforms_dir = shift_03_10_2020.get_waveforms_dir()
results_dir = shift_03_10_2020.get_results_dir()
acnet_data_dir = shift_03_10_2020.get_acnet_data_dir()
t1 = datetime(2020, 3, 10, 11, 32)
t2 = datetime(2020, 3, 10, 12, 58)
save_acnet_data_for_fur(shift_03_10_2020, t1, t2, "all_acnet_data_03_10_2020.csv")