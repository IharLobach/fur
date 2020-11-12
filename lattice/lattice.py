import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import os
labeb_fs = 30
annot_fs = 20
legend_fs = 20
text_fs = annot_fs
marker_size = 7
params = {'axes.labelsize': labeb_fs}
plt.rcParams.update(params)
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import BaggingRegressor
from config import get_from_config
camera_positions = [312.3907724914158,
                    847.5314240270161,
                    1146.7606367662588,
                    1648.5578345089393,
                    2312.048374490958,
                    2846.8364422732257,
                    3146.9275264095036,
                    3646.7941322228266]
camera_names = ['M1R', 'M2R', 'M3R', 'M4R', 'M4L', 'M3L', 'M2L', 'M1L']
active_cameras = ['M1R', 'M2R', 'M3R', 'M4L', 'M3L', 'M2L', 'M1L']
# these calibrations only make it worse:
# calibrations = [
#     5/0.772,
#     5/0.77,
#     5/0.779,
#     5/0.78,
#     5/0.779,
#     5/0.824,
#     5/0.724,
#     5/0.838
# ]
# cal_av = np.mean(calibrations)
# calibrations = np.array(calibrations)/cal_av

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
undulator_range = (1383.9, 1435.4)
undulator_middle = np.mean(undulator_range)
undulator_name = "Und."




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
    lattice_df = lattice_df.drop_duplicates("S_cm")
    lattice_df = lattice_df.reset_index(drop=True)
    dS = lattice_df["S_cm"].diff().fillna(method='bfill')
    dDx = lattice_df["Dispersion_cm_X"].diff().fillna(method='bfill')
    lattice_df["dS"] = dS
    lattice_df["dDx"] = dDx
    lattice_df["dDx/dS"] = dDx/dS
    lattice_df['Phi_X'] = lattice_df['Dispersion_cm_X']\
        * lattice_df['Alpha_X']/lattice_df['Beta_cm_X'] + lattice_df['dDx/dS']
    lattice_df['H'] = lattice_df['Dispersion_cm_X']**2/lattice_df['Beta_cm_X']\
        + lattice_df['Beta_cm_X']*lattice_df['Phi_X']**2
    # lattice_df["Alpha_X"] = -lattice_df["Beta_cm_X"].diff()/lattice_df["dS"]/2
    # lattice_df["Alpha_Y"] = -lattice_df["Beta_cm_Y"].diff()/lattice_df["dS"]/2
    # lattice_df['Alpha_X'] = lattice_df['Alpha_X'].fillna(method='bfill')
    # lattice_df['Alpha_Y'] = lattice_df['Alpha_Y'].fillna(method='bfill')
    # idx_ac = 793
    # a_ac_cm = 0.38
    # beta_ac_cm = lattice_df.loc[idx_ac, 'Beta_cm_X']
    # dispersion_ac_cm = lattice_df.loc[idx_ac, 'Dispersion_cm_X']
    # lattice_df['daccL'] = a_ac_cm\
    #     / (np.sqrt(lattice_df['H']*beta_ac_cm)\
    #     + np.abs(dispersion_ac_cm))

    return lattice_df


def add_vertical_lines_at_camera_positions(ax, color='green'):
    for p in camera_positions:
        ax.axvline(p, color=color)


def annotate_camera_positions(ax):
    y_pos_annotate = np.mean(ax.get_ylim())
    for name, p in zip(camera_names, camera_positions):
        ax.annotate(name, (p, y_pos_annotate), fontsize=annot_fs)


def add_undulator_shaded_area(ax, color='blue'):
    ax.axvspan(undulator_range[0], undulator_range[1],
               alpha=0.5, color=color)


def annotate_undulator(ax):
    y_pos_annotate = np.mean(ax.get_ylim())
    ax.annotate(undulator_name,
                (undulator_range[1], y_pos_annotate), fontsize=annot_fs)


def plot_lattice(lattice_df):
    fig, (ax0, ax1) = plt.subplots(2, figsize=(20, 15))
    ax0.plot(lattice_df["S_cm"], lattice_df["Beta_cm_X"], label="Beta_cm_X")
    ax0.plot(lattice_df["S_cm"], lattice_df["Beta_cm_Y"], label="Beta_cm_Y")
    ax0.set_ylabel("Beta_cm_X, Beta_cm_Y")
    ax0.legend(fontsize=legend_fs)
    ax1.plot(lattice_df["S_cm"], lattice_df["Dispersion_cm_X"])
    ax1.set_ylabel("Dispersion_cm_X")
    ax1.set_xlabel("S_cm")
    # ax2.plot(lattice_df["S_cm"], lattice_df["Alpha_X"], label="Alpha_X")
    # ax2.plot(lattice_df["S_cm"], lattice_df["Alpha_Y"], label="Alpha_Y")
    # ax2.set_ylabel("Alpha_X, Alpha_Y", fontsize=fs)
    # ax2.set_xlabel("S_cm")
    # ax2.legend()
    for ax in (ax0, ax1):
        add_vertical_lines_at_camera_positions(ax)
    annotate_camera_positions(ax1)
    add_undulator_shaded_area(ax=ax1)
    annotate_undulator(ax1)
    plt.show()


def get_sigma_um(beta_cm, e_um, dispersion_cm, dpp):
    return np.sqrt(e_um*1e4*beta_cm+(1e4*dispersion_cm*dpp)**2)


def show_sigma_fit(lattice_df, cameras_df, axis, emittance_um, dpp=None):
    if (axis == "X") and (dpp is None):
        raise ValueError("Chosen axis is 'X', but value"
                         " of 'dpp' was not specified.")
    elif axis == "X":
        dispersion_cm = lattice_df["Dispersion_cm_X"]
        dpp_ = dpp
        s = "dp/p = {:.3e}".format(dpp)
    elif axis == "Y":
        dispersion_cm = 0
        dpp_ = 0
        s = ''
    fig, ax = plt.subplots(figsize=(20, 7.5))
    ax.plot(lattice_df["S_cm"],
            get_sigma_um(
                lattice_df["Beta_cm_"+axis],
                emittance_um,
                dispersion_cm,
                dpp_),
            label="sigma_um_"+axis)
    no_m4r = cameras_df[cameras_df["Name"].isin(active_cameras)]
    add_vertical_lines_at_camera_positions(ax)
    annotate_camera_positions(ax)
    ax.plot(no_m4r["S_cm"], no_m4r["Measured_sigma_um_"+axis],
            marker='o', linestyle='none', label="Measured_sigma_um_"+axis,
            markersize=marker_size)
    ax.set_ylabel("sigma_um_"+axis)
    ax.set_xlabel("S_cm")
    add_undulator_shaded_area(ax)
    annotate_undulator(ax)
    ax.text(0.9, 0.9, "emittance_"+axis+" = {:3f} um".format(emittance_um),
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes, fontsize=text_fs)
    ax.text(0.9, 0.8, s,
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes, fontsize=text_fs)

    ax.set_xlim(0, ax.get_xlim()[1])
    ax.legend(fontsize=legend_fs)
    plt.show()


def show_angle_spread_X_Y(lattice_df, e_um_x, e_um_y, dpp=0):
    gamma_x_um_m1 = (1+lattice_df["Alpha_X"]**2)/lattice_df["Beta_cm_X"]/1e4
    gamma_y_um_m1 = (1+lattice_df["Alpha_Y"]**2)/lattice_df["Beta_cm_Y"]/1e4
    angle_spread_x = np.sqrt(lattice_df["dDx/dS"]**2*dpp**2
                             + e_um_x*gamma_x_um_m1)
    angle_spread_y = np.sqrt(e_um_y*gamma_y_um_m1)
    fig, ax = plt.subplots(figsize=(20, 7.5))
    ax.plot(lattice_df["S_cm"], angle_spread_x, label="Angle spread X")
    ax.plot(lattice_df["S_cm"], angle_spread_y, label="Angle spread Y")
    ax.axhline(1/200/np.sqrt(10),
               label="Characteristic angular spread of"
                     " the undulator radiation",
               color="blue")
    add_vertical_lines_at_camera_positions(ax)
    annotate_camera_positions(ax)
    ax.set_ylabel("angle spread, rad")
    ax.set_xlabel("S_cm")
    add_undulator_shaded_area(ax)
    annotate_undulator(ax)
    ax.set_xlim(0, ax.get_xlim()[1])
    ax.legend(fontsize=legend_fs)
    plt.show()


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
    no_m4r = cameras_df[cameras_df["Name"].isin(active_cameras)]
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
                 transform=ax1.transAxes, fontsize=text_fs)
        ax1.set_xlim(0, ax1.get_xlim()[1])
        plt.show()
    return ey_um_description


def get_e_um_Y_scipy_curve_fit(cameras_df):
    no_m4r = cameras_df[cameras_df["Name"].isin(active_cameras)]
    popt, pcov = scipy.optimize.curve_fit(
        lambda beta_cm, e_um: np.sqrt(e_um*1e4*beta_cm),
        no_m4r["Beta_cm_Y"],
        no_m4r["Measured_sigma_um_Y"])
    perr = np.sqrt(np.diag(pcov))[0]
    return popt[0], perr


def get_e_um_X_scipy_curve_fit(cameras_df, dpp=None, dpp_err=0):
    no_m4r = cameras_df[cameras_df["Name"].isin(active_cameras)]
    if dpp is None:
        def f(beta_cm_dispersion_cm, e_um, dpp):
            beta_cm, disperison_cm = beta_cm_dispersion_cm
            return get_sigma_um(beta_cm, e_um,
                                  disperison_cm, dpp)

        popt, pcov = scipy.optimize.curve_fit(
            f,
            no_m4r.loc[:, ["Beta_cm_X", "Dispersion_cm_X"]].values.T,
            no_m4r["Measured_sigma_um_X"])
        perr = np.sqrt(np.diag(pcov))
        return popt, perr
    else:
        def f(beta_cm_dispersion_cm, e_um,):
            beta_cm, disperison_cm = beta_cm_dispersion_cm
            return get_sigma_um(beta_cm, e_um,
                                  disperison_cm, dpp)
        
        popt, pcov = scipy.optimize.curve_fit(
            f,
            no_m4r.loc[:, ["Beta_cm_X", "Dispersion_cm_X"]].values.T,
            no_m4r["Measured_sigma_um_X"])
        perr = np.sqrt(np.diag(pcov))
        return (popt[0], dpp), (perr[0], dpp_err)


def get_undulator_df(lattice_df, emittance_6D):
    e_um_x, ex_err, e_um_y, ey_err, dpp, dpp_err = emittance_6D.values()
    undulator_df = pd.DataFrame(index=["Start", "Middle", "End"])
    undulator_df["S_cm"] = [
        undulator_range[0],
        undulator_middle,
        undulator_range[1]]
    for col in lattice_df.columns[1:]:
        undulator_df[col] = np.interp(undulator_df["S_cm"],
                                      lattice_df["S_cm"],
                                      lattice_df[col])
    undulator_df["Sigma_um_X"] = get_sigma_um(
        undulator_df["Beta_cm_X"],
        e_um_x,
        undulator_df["Dispersion_cm_X"],
        dpp)
    undulator_df["Sigma_um_X_err"] = 1/undulator_df["Sigma_um_X"] \
        * (undulator_df["Beta_cm_X"]*1e4*ex_err/2
            + (undulator_df["Dispersion_cm_X"]*1e4)**2*dpp*dpp_err)
    undulator_df["Sigma_um_Y"] = get_sigma_um(
        undulator_df["Beta_cm_Y"],
        e_um_y,
        0,
        0)
    undulator_df["Sigma_um_Y_err"] = 1/undulator_df["Sigma_um_Y"] \
        * (undulator_df["Beta_cm_Y"]*1e4*ey_err/2)
    gamma_x_um_m1 =\
        (1+undulator_df["Alpha_X"]**2)/undulator_df["Beta_cm_X"]/1e4
    gamma_y_um_m1 =\
        (1+undulator_df["Alpha_Y"]**2)/undulator_df["Beta_cm_Y"]/1e4
    undulator_df["Angle_spread_rad_X"]\
        = np.sqrt(undulator_df["dDx/dS"]**2*dpp**2+e_um_x*gamma_x_um_m1)
    undulator_df["Angle_spread_rad_Y"] = np.sqrt(e_um_y*gamma_y_um_m1)
    undulator_df["ex_um"] = np.ones(3)*e_um_x
    undulator_df["ex_err"] = np.ones(3)*ex_err
    undulator_df["ey_um"] = np.ones(3)*e_um_y
    undulator_df["ey_err"] = np.ones(3)*ey_err
    undulator_df["dp/p"] = np.ones(3)*dpp
    undulator_df["dp/p_err"] = np.ones(3)*dpp_err
    return undulator_df


c = get_from_config("c_m/s")
f0 = 1/get_from_config("IOTA_revolution_period")
gamma0 = get_from_config("gamma")
e = get_from_config("e")
me = get_from_config("me_MeV")
alpha = get_from_config("ring_alpha")
q = get_from_config("RF_q")

def get_dpp(sigma_z_cm, Vrf_V, gamma=gamma0):
    E = gamma*me
    eta_s = alpha-1/gamma**2
    beta = np.sqrt(1-1/gamma**2)
    f = q*f0
    return sigma_z_cm*1e-2/c*f*2*np.pi\
        * np.sqrt(Vrf_V/(2*np.pi*1e6*E*beta**2*q*np.abs(eta_s)))


def get_e_um_Y_dict(cameras_df):
    no_m4r = cameras_df[cameras_df["Name"].isin(active_cameras)]
    ey_vals = no_m4r["Measured_sigma_um_Y"]**2/(1e4*no_m4r["Beta_cm_Y"])
    ey_dict = {}
    for cam, val in zip(no_m4r["Name"], ey_vals):
        ey_dict["e_um_Y_"+cam] = val
    return ey_dict


def get_e_um_X_dict(cameras_df, dpp=0):
    no_m4r = cameras_df[cameras_df["Name"].isin(active_cameras)]
    ex_vals = (no_m4r["Measured_sigma_um_X"]**2
               - (1e4*no_m4r["Dispersion_cm_X"]*dpp)**2)\
        / (1e4*no_m4r["Beta_cm_X"])
    ex_dict = {}
    for cam, val in zip(no_m4r["Name"], ex_vals):
        ex_dict["e_um_X_"+cam] = val
    return ex_dict


# Coupled lattice below

U = np.array([
    [0, 1, 0, 0],
    [-1, 0, 0, 0],
    [0, 0, 0, 1],
    [0, 0, -1, 0]
])


def BuildXi(Xi, D, sp):
    res = np.zeros(shape=(5, 5))
    res[:4, :4] = Xi
    res[4, 4] = D@Xi@D+1/sp**2
    res[4, :4] = Xi@D
    res[:4, 4] = res[4, :4]
    return res

class CoupledLattice:
    def __init__(self, file_path, gamma=gamma0):
        self.gamma = gamma
        with open(file_path) as f:
            f.read(1)
            ldf = pd.read_table(f, delim_whitespace=True)
        cameras_df = pd.DataFrame({
            "Name": camera_names,
            "S[cm]": camera_positions
        })
        for col in ldf.columns[3:]:
            cameras_df[col] = np.interp(cameras_df["S[cm]"],
                                        ldf["S[cm]"],
                                        ldf[col])
        cameras_df["ACNET_device_X"] = acnet_devices_X
        cameras_df["ACNET_device_Y"] = acnet_devices_Y
        
        self.cameras_df = cameras_df[cameras_df["Name"].isin(active_cameras)]
        self.ldf = ldf
    
    def get_density_matrices(self, ldf, e1, e2, sp):
        ldf['nu1'] = 2*np.pi*ldf['Teta1/(2*PI)']
        ldf['nu2'] = 2*np.pi*ldf['Teta2/(2*PI)']
        N = len(ldf.index)
        c1 = np.cos(ldf['nu1'])
        c2 = np.cos(ldf['nu2'])
        s1 = np.sin(ldf['nu1'])
        s2 = np.sin(ldf['nu2'])
        V = np.zeros(shape=(N, 4, 4))
        V[:, 0, 0] = np.sqrt(ldf['BetaX1'])
        V[:, 0, 2] = np.sqrt(ldf['BetaX2'])*c2
        V[:, 0, 3] = -np.sqrt(ldf['BetaX2'])*s2
        V[:, 1, 0] = -ldf['AlfaX1']/np.sqrt(ldf['BetaX1'])
        V[:, 1, 1] = (1-ldf['U'])/np.sqrt(ldf['BetaX1'])
        V[:, 1, 2] = (ldf['U']*s2-ldf['AlfaX2']*c2)/np.sqrt(ldf['BetaX2'])
        V[:, 1, 3] = (ldf['U']*c2+ldf['AlfaX2']*s2)/np.sqrt(ldf['BetaX2'])
        V[:, 2, 0] = np.sqrt(ldf['BetaY1'])*c1
        V[:, 2, 1] = -np.sqrt(ldf['BetaY1'])*s1
        V[:, 2, 2] = np.sqrt(ldf['BetaY2'])
        V[:, 3, 0] = (ldf['U']*s1-ldf['AlfaY1']*c1)/np.sqrt(ldf['BetaY1'])
        V[:, 3, 1] = (ldf['U']*c1+ldf['AlfaY1']*s1)/np.sqrt(ldf['BetaY1'])
        V[:, 3, 2] = -ldf['AlfaY2']/np.sqrt(ldf['BetaY2'])
        V[:, 3, 3] = (1-ldf['U'])/np.sqrt(ldf['BetaY2'])
        D = np.zeros(shape=(N, 4))
        D[:, 0] = ldf['DspX']
        D[:, 1] = ldf['DspXp']
        D[:, 2] = ldf['DspY']
        D[:, 3] = ldf['DspYp']
        # Sigma = V@np.diag([e1, e1, e2, e2])@(np.transpose(V, axes=(0, 2, 1)))
        Xi = U@V@np.diag([1/e1, 1/e1, 1/e2, 1/e2]
                         )@np.transpose(V, axes=(0, 2, 1))@np.transpose(U)
        Xitot = np.array([BuildXi(xi, d, sp) for xi, d in zip(Xi, D)])
        Sigmatot = np.array([np.linalg.inv(sl) for sl in Xitot])
        Xitheta = np.array([
            np.array([
                [m[1, 1], m[1, 3], m[1, 4]],
                [m[3, 1], m[3, 3], m[3, 4]],
                [m[4, 1], m[4, 3], m[4, 4]]
            ]) for m in Xitot
        ])
        Sigmax = np.array([
            np.array([
                [m[0, 0], m[0, 2]],
                [m[2, 0], m[2, 2]]
            ]) for m in Sigmatot
        ])
        Gamma = np.diag([1, 1, 1/self.gamma])
        SigmathetaBF = Gamma@np.linalg.inv(Xitheta)@Gamma
        stheta2BF = np.linalg.eigvals(SigmathetaBF)
        return Sigmax, Sigmatot

    def get_cameras_sx_sy(self, e1, e2, sp):
        Sigmax, Sigmatheta = self.get_density_matrices(self.cameras_df, e1, e2, sp)
        return np.sqrt(Sigmax[:,0,0]), np.sqrt(Sigmax[:,1,1])

    def get_s1_s2_th1_th2(self, Sigmax, Sigmatheta):
        part1 = (Sigmax[:, 0, 0]+Sigmax[:, 1, 1])
        part2 = np.sqrt((Sigmax[:, 0, 0]-Sigmax[:, 1, 1])**2
                        + 4*Sigmax[:, 0, 1]**2)
        s2 = 1/np.sqrt(2)*np.sqrt(part1+part2)
        s1 = 1/np.sqrt(2)*np.sqrt(part1-part2)
        
        
    
    
    
    

    




