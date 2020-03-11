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
    return datetime.strptime(s,"%d-%b-%Y %H:%M:%S.%f")


def make_acl_command(t1, t2, device):
    t1_s = datetime_to_string(t1)
    t2_s = datetime_to_string(t2)
    return 'logger_get/start="{}"/end="{}" {}'.format(t1_s, t2_s, device)


def parse_text(text, device="value"):
    lines = text.split(b'\n')
    time_list = []
    value_list = []
    for line in lines[:-1]:
        time_s, value_s = line.strip().split(b'  ')
        time = string_to_datetime(time_s.decode())
        value = float(value_s)
        time_list.append(time)
        value_list.append(value)
    return pd.Series(value_list, index=pd.DatetimeIndex(time_list), name=device)


def fetch_data_one_device(t1, t2, device):
    acl_command = make_acl_command(t1, t2, device)
    request = prefix+acl_command
    response = requests.get(request)
    if response.status_code != 200:
        raise Exception('Fetching data was unsuccessful')
    return parse_text(response.content, device)


def fetch_data(t1, t2, devices):
    if isinstance(devices, str):
        return fetch_data_one_device(t1, t2, devices)
    else:
        return pd.concat([fetch_data_one_device(t1, t2, d) for d in devices], axis=1)