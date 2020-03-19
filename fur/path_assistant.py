import os
import sys
from datetime import datetime
from fur.config_requests import get_from_config
data_folder = get_from_config("data_folder")
shifts_folder = os.path.join(data_folder, 'shifts')
shift_folders = os.listdir(shifts_folder)


def show_shift_folders():
    global shift_folders
    for i, f in enumerate(shift_folders):
        print(i, f)


class WorkingDirectory():
    def __init__(self, wd):
        self.wd = wd

    def fi(self, file_name):
        return os.path.join(self.wd, file_name)


class PathAssistant():
    def __init__(self, shift_fodler_name, ignore_files=None):
        self.shift_folder_name = shift_fodler_name
        self.shift_dir = os.path.join(shifts_folder,
                                      self.shift_folder_name)
        self.shift_results_dir = os.path.join(self.shift_dir,
                                              "results")
        self.waveforms_folder_path = \
            os.path.join(self.shift_dir, "waveforms")
        self.acnet_data_dir = os.path.join(self.shift_dir, "acnet_data")
        self.ignore_files = ['desktop.ini']
        if ignore_files:
            self.ignore_files += ignore_files

    def get_waveforms_folder_path(self):
        return os.path.join(shifts_folder, self.shift_folder_name,
                            "waveforms")

    def get_waveform_files(self):
        shift = self.waveforms_folder_path
        files = [f for f in os.listdir(shift)
                 if os.path.isfile(os.path.join(shift, f))]
        waveforms = [f for f in files if ".Wfm." in f]
        return [f for f in waveforms if (f not in self.ignore_files)]

    def show_waveform_file_names(self):
        files = self.get_waveform_files()
        for i, f in enumerate(files):
            print(i, f)

    def get_waveform_paths(self):
        shift = self.waveforms_folder_path
        files = self.get_waveform_files()
        return [os.path.join(shift, f) for f in files]

    def get_waveform_path(self, waveform_file_name):
        shift = self.waveforms_folder_path
        return os.path.join(shift, waveform_file_name)

    def get_waveforms_dir(self):
        return WorkingDirectory(self.waveforms_folder_path)

    def get_results_dir(self):
        return WorkingDirectory(self.shift_results_dir)

    def get_acnet_data_dir(self):
        return WorkingDirectory(self.acnet_data_dir)

    def get_datetime(self, waveform_name):
        try:
            _, day_str, _, time_str = waveform_name.split('_')
            time_str = time_str.split('.')[0]
            datetime_str = day_str+' '+time_str
            return datetime.strptime(datetime_str, "%Y-%m-%d %H%M%S")
        except Exception as e:
            print("Exceptino happened in get_datetime('{}')"
                  .format(waveform_name), e)

