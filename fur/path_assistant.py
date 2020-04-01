import os
import sys
from datetime import datetime
import pandas as pd
from config import get_from_config, save_to_config
data_folder = get_from_config("data_folder")
shifts_folder = os.path.join(data_folder, 'shifts')
shift_folders = os.listdir(shifts_folder)
additional_data_folder = os.path.join(shifts_folder, 'additional_data')
srw_precalculated_spectrum_folder = \
    os.path.join(data_folder, "SRW_SLAC_undulator_spectrum")
srw_Ex_3D_file_path = os.path.join(
    srw_precalculated_spectrum_folder, "Ex_3D.npy")
srw_Ey_3D_file_path = os.path.join(
    srw_precalculated_spectrum_folder, "Ey_3D.npy")


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
        self.bpm_data_dir = os.path.join(self.shift_dir, "bunch_profile_meas")
        self.lattice_6dsim_dir = os.path.join(self.shift_dir, "6dsim")
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

    def get_bpm_data_dir(self):
        return WorkingDirectory(self.bpm_data_dir)

    def get_additional_data_dir(self):
        return WorkingDirectory(additional_data_folder)

    def get_6dsim_dir(self):
        return WorkingDirectory(self.lattice_6dsim_dir)

    def get_datetime(self, waveform_name):
        try:
            _, day_str, _, time_str = waveform_name.split('_')
            time_str = time_str.split('.')[0]
            datetime_str = day_str+' '+time_str
            return datetime.strptime(datetime_str, "%Y-%m-%d %H%M%S")
        except Exception as e:
            print("Exception happened in get_datetime('{}')"
                  .format(waveform_name), e)

    def get_bpm_files_df(self):
        bpm_wf_files = [f for f in os.listdir(self.bpm_data_dir)
                        if f != 'log.db']
        bpm_data_dir = self.get_bpm_data_dir()

        def string_to_datetime(s):
            return datetime.strptime(s,
                                     "bunch_profile_%m-%d-%Y_%H_%M_%S_%f.csv")

        bpm_data_df = pd.DataFrame(
            {"file_name": bpm_wf_files,
             "file_path": [bpm_data_dir.fi(f) for f in bpm_wf_files],
             "file_datetime": [string_to_datetime(f) for f in bpm_wf_files]})
        bpm_data_df.sort_values("file_datetime")
        return bpm_data_df

    def get_acnet_data_df(self, file_name):
        acnet_data_dir = self.get_acnet_data_dir()
        return pd.read_csv(acnet_data_dir.fi(file_name), index_col=0,
                           parse_dates=True)

    def get_fluctuations_df(self, file_name):
        results_dir = self.get_results_dir()
        res_df = pd.read_csv(
            results_dir.fi(file_name),
            index_col=0)
        res_df["file_datetime"] = res_df["waveform_file"] \
            .apply(self.get_datetime)
        res_df.sort_values("file_datetime")
        return res_df

    def get_fluctuation_waveforms_df(self):
        fluctuation_waveforms_df = pd.DataFrame({
            "file_name": self.get_waveform_files(),
            "file_path": self.get_waveform_paths()})
        fluctuation_waveforms_df["file_datetime"] = \
            fluctuation_waveforms_df["file_name"].apply(self.get_datetime)
        fluctuation_waveforms_df.sort_values("file_datetime")
        return fluctuation_waveforms_df
