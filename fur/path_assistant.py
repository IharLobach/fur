import os
import sys
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
        self.waveforms_folder_path = \
            os.path.join(shifts_folder,
                         self.shift_folder_name, "waveforms")
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