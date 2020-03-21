from datetime import datetime
import time
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
prefix = 'https://www-bd.fnal.gov/cgi-bin/acl.pl?acl='


def datetime_to_string(t):
    return t.strftime("%d-%b-%Y %H:%M:%S.%f")[:-3]


def string_to_datetime(s):
    return datetime.strptime(s, "%d-%b-%Y %H:%M:%S.%f")


def make_acl_command(t1, t2, device):
    t1_s = datetime_to_string(t1)
    t2_s = datetime_to_string(t2)
    return 'logger_get/start="{}"/end="{}" {}'.format(t1_s, t2_s, device)


def parse_text(text, device="value"):
    lines = text.split(b'\n')
    time_list = []
    value_list = []
    for line in lines[:-1]:
        try:
            time_s, value_s = line.strip().split(b'  ')
            time = string_to_datetime(time_s.decode())
            value = float(value_s)
            time_list.append(time)
            value_list.append(value)
        except Exception as e:
            print("Exception happened while parsing line '{}'"
                  .format(line.decode()))
            print("Exception: ", e)
    return pd.DataFrame({device: value_list}, index=pd.DatetimeIndex(time_list))


def fetch_data_one_device(t1, t2, device):
    acl_command = make_acl_command(t1, t2, device)
    request = prefix+acl_command
    response = requests.get(request)
    if response.status_code != 200:
        raise Exception('Fetching data for {} was unsuccessful'.format(device))
    return parse_text(response.content, device)


def fetch_data(t1, t2, devices):
    if isinstance(devices, str):
        return fetch_data_one_device(t1, t2, devices)
    else:
        return pd.concat([fetch_data_one_device(t1, t2, d) for d in devices],
                         axis=1)


def get_interpolated_df(df, date_range):
    x = date_range.astype('uint64')
    xp = df.index.astype('uint64')
    res_df = pd.DataFrame(index=pd.DatetimeIndex(date_range))
    for col in df.columns:
        res_df[col] = np.interp(x, xp, df[col])
    return res_df


def fetch_interpolated_data(t1, t2, devices, freq='S'):
    date_range = pd.date_range(t1, t2, freq='S')
    df = fetch_data(t1, t2, devices)
    return get_interpolated_df(df, date_range)
