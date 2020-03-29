import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import BaggingRegressor
camera_positions = [312.3907724914158,
                    847.5314240270161,
                    1146.7606367662588,
                    1648.5578345089393,
                    2312.048374490958,
                    2846.8364422732257,
                    3146.9275264095036,
                    3646.7941322228266]
camera_names = ['M1R', 'M2R', 'M3R', 'M4R', 'M4L', 'M3L', 'M2L', 'M1L']
acnet_devices_X = ['N:ITC1RSH',
                   'N:ITC2RSH',
                   'N:ITC3RSH',
                   'N:ITC4RSH',
                   'N:ITC4LSH',
                   'N:ITC3LSH',
                   'N:ITC2LSH',
                   'N:ITC1LSH']
acnet_devices_Y = ['N:ITC1RSV',
                   'N:ITC2RSV',
                   'N:ITC3RSV',
                   'N:ITC4RSV',
                   'N:ITC4LSV',
                   'N:ITC3LSV',
                   'N:ITC2LSV',
                   'N:ITC1LSV']


def read_lattice_file(lattice_file_path):
    with open(lattice_file_path) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    content = content[5:]
    content = [x for x in content if len(x)>0]
    # loading the data
    data_chunks = []
    chunk_x = []
    chunk_y = []
    data_names = [content[0]]
    len_cont = len(content)
    i = 2
    while i < len_cont:
        line = content[i]
        if "Data" in line:
            data_names.append(line)
            i += 1
            data_chunks.append([chunk_x, chunk_y])
            chunk_x = []
            chunk_y = []
        else:
            x, sigma_x, y, sigma_y = [float(j) for j in line.split()]
            chunk_x.append(x)
            chunk_y.append(y)
            if i == len_cont-1:
                data_chunks.append([chunk_x, chunk_y])
        i += 1
    data_chunks[0] = [data_chunks[0][0], [0.1*q for q in data_chunks[0][1]]]
    data_names[0] = data_names[0][:-2]+"cm"
    lattice_df = pd.DataFrame({
        "S_cm": data_chunks[0][0],
        "Dispersion_cm_X": data_chunks[0][1],
        "Beta_cm_X": data_chunks[1][1],
        "Beta_cm_Y": data_chunks[2][1],
        "Alpha_X": data_chunks[3][1],
        "Alpha_Y": data_chunks[4][1]
    })
    return lattice_df


def add_vertical_lines_at_camera_positions(ax, color='green'):
    for p in camera_positions:
        ax.axvline(p, color='green')


def annotate_camera_positions(ax):
    y_pos_annotate = np.mean(ax.get_ylim())
    for name, p in zip(camera_names, camera_positions):
        ax.annotate(name, (p, y_pos_annotate))


def plot_lattice(lattice_df):
    fig, (ax0, ax1) = plt.subplots(2, figsize=(20, 15))
    ax0.plot(lattice_df["S_cm"], lattice_df["Beta_cm_X"], label="Beta_cm_X")
    ax0.plot(lattice_df["S_cm"], lattice_df["Beta_cm_Y"], label="Beta_cm_Y")
    fs = 16
    ax0.set_ylabel("Beta_cm_X, Beta_cm_Y", fontsize=fs)
    ax0.legend()
    ax1.plot(lattice_df["S_cm"], lattice_df["Dispersion_cm_X"])
    ax1.set_ylabel("Dispersion_cm_X", fontsize=fs)
    ax1.set_xlabel("S_cm")
    # ax2.plot(lattice_df["S_cm"], lattice_df["Alpha_X"], label="Alpha_X")
    # ax2.plot(lattice_df["S_cm"], lattice_df["Alpha_Y"], label="Alpha_Y")
    # ax2.set_ylabel("Alpha_X, Alpha_Y", fontsize=fs)
    # ax2.set_xlabel("S_cm")
    # ax2.legend()
    for ax in (ax0, ax1):
        add_vertical_lines_at_camera_positions(ax)
    annotate_camera_positions(ax1)
    plt.show()


def __get_sigma_um(beta_cm, e_um, dispersion_cm, dpp):
    return np.sqrt(e_um*1e4*beta_cm+(1e4*dispersion_cm*dpp)**2)


def show_sigma_fit(lattice_df, cameras_df, axis, emittance_um, dpp=None):
    if (axis == "X") and (dpp is None):
        raise ValueError("Chosen axis is 'X', but value"
                         " of 'dpp' was not specified.")
    elif axis == "X":
        dispersion_cm = lattice_df["Dispersion_cm_X"]
        dpp_ = dpp
    elif axis == "Y":
        dispersion_cm = 0
        dpp_ = 0
    fig, ax = plt.subplots(figsize=(20, 7.5))
    ax.plot(lattice_df["S_cm"],
            __get_sigma_um(
                lattice_df["Beta_cm_"+axis],
                emittance_um,
                dispersion_cm,
                dpp_),
            label="sigma_um_"+axis)
    no_m4r = cameras_df[cameras_df["Name"] != "M4R"]
    add_vertical_lines_at_camera_positions(ax)
    annotate_camera_positions(ax)
    ax.plot(no_m4r["S_cm"], no_m4r["Measured_sigma_um_"+axis],
            marker='o', linestyle='none', label="Measured_sigma_um_"+axis)
    fs = 16
    ax.set_ylabel("sigma_um_"+axis, fontsize=fs)
    ax.legend()
    return fig, ax


def get_cameras_df(lattice_df,
                   measured_sigma_um_X,
                   measured_sigma_um_Y):
    cameras_df = pd.DataFrame({
                                "Name": camera_names,
                                "S_cm": camera_positions
                            })
    for col in lattice_df.columns[1:]:
        cameras_df[col] = np.interp(cameras_df["S_cm"],
                                    lattice_df["S_cm"],
                                    lattice_df[col])
    cameras_df["ACNET_device_X"] = acnet_devices_X
    cameras_df["ACNET_device_Y"] = acnet_devices_Y
    cameras_df["Measured_sigma_um_X"] = measured_sigma_um_X
    cameras_df["Measured_sigma_um_Y"] = measured_sigma_um_Y
    return cameras_df


def get_ey_um_least_squares(
        cameras_df, show_plot=False,
        n_estimators=10000,
        n_estimators_to_plot=50):
    no_m4r = cameras_df[cameras_df["Name"] != "M4R"]
    x = no_m4r["Beta_cm_Y"]
    X = x.values.reshape(-1, 1)
    y = no_m4r["Measured_sigma_um_Y"]**2
    model = BaggingRegressor(LinearRegression(fit_intercept=False), 
                             n_estimators=n_estimators,
                             bootstrap=True)
    _ = model.fit(X, y)
    x_model = np.array([0, max(x)])
    coefs = np.zeros(n_estimators)
    for i, m in enumerate(model.estimators_):
        coefs[i] = m.coef_[0]
    coefs = 1e-4*coefs
    ey_um_description = pd.Series(coefs).describe()
    if show_plot:
        fig, (ax0, ax1) = plt.subplots(2, figsize=(15, 15))
        ax0.plot(x, y, 'o')
        ax0.set_ylim(0, ax0.get_ylim()[1])
        ax0.set_ylabel("Measured_sigma_um_Y^2")
        ax0.set_xlabel("Beta_cm_Y")
        for x_, y_, name in zip(x, y, no_m4r["Name"]):
            ax0.annotate(name, (x_, y_))
        for c in np.random.choice(1e4*coefs, size=n_estimators_to_plot):
            ax0.plot(x_model, c*x_model, color='grey', alpha=0.2)
        ax0.set_title("Bootstrapping for vertical emittance")
        sns.boxenplot(coefs, ax=ax1)
        ax1.set_xlabel("Vertical emittance, um")
        s = "{}".format(ey_um_description)
        lines = s.split('\n')[:-1]
        s = ''
        for l in lines:
            s += l+'\n'
        ax1.text(0.9, 0.5, s,
                 verticalalignment='bottom', horizontalalignment='right',
                 transform=ax1.transAxes, fontsize=15)
        ax1.set_xlim(0, ax1.get_xlim()[1])
        plt.show()
    return ey_um_description


def get_e_um_Y_scipy_curve_fit(cameras_df):
    no_m4r = cameras_df[cameras_df["Name"] != "M4R"]
    popt, pcov = scipy.optimize.curve_fit(
        lambda beta_cm, e_um: np.sqrt(e_um*1e4*beta_cm),
        no_m4r["Beta_cm_Y"],
        no_m4r["Measured_sigma_um_Y"])
    perr = np.sqrt(np.diag(pcov))[0]
    return popt[0], perr


def get_e_um_X_scipy_curve_fit(cameras_df):
    no_m4r = cameras_df[cameras_df["Name"] != "M4R"]

    def f(beta_cm_dispersion_cm, e_um, dpp):
        beta_cm, disperison_cm = beta_cm_dispersion_cm
        return __get_sigma_um(beta_cm, e_um,
                              disperison_cm, dpp)

    popt, pcov = scipy.optimize.curve_fit(
        f,
        no_m4r.loc[:, ["Beta_cm_X", "Dispersion_cm_X"]].values.T,
        no_m4r["Measured_sigma_um_X"])
    perr = np.sqrt(np.diag(pcov))
    return popt, pcov, perr
