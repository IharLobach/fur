import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
camera_positions = [
                    313.99411487289865,
                    850.7097357006688,
                    1150.6154142505097,
                    1649.9644633972994,
                    2313.0210097709055,
                    2849.9525071571384,
                    3148.9348982723236,
                    3648.20091797355
                    ]
camera_names = ['M1R', 'M2R', 'M3R', 'M4R', 'M4L', 'M3L', 'M2L', 'M1L']


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
        "Dispersion_X_cm": data_chunks[0][1],
        "Beta_X_cm": data_chunks[1][1],
        "Beta_Y_cm": data_chunks[2][1],
        "Alpha_X": data_chunks[3][1],
        "Alpha_Y": data_chunks[4][1]
    })
    return lattice_df


def plot_lattice(lattice_df):
    fig, (ax0, ax1, ax2) = plt.subplots(3, figsize=(20,15))
    ax0.plot(lattice_df["S_cm"], lattice_df["Beta_X_cm"], label="Beta_X_cm")
    ax0.plot(lattice_df["S_cm"], lattice_df["Beta_Y_cm"], label="Beta_Y_cm")
    fs = 16
    ax0.set_ylabel("Beta_X_cm, Beta_Y_cm", fontsize=fs)
    ax0.legend()
    ax1.plot(lattice_df["S_cm"], lattice_df["Dispersion_X_cm"])
    ax1.set_ylabel("Dispersion_X_cm", fontsize=fs)
    ax2.plot(lattice_df["S_cm"], lattice_df["Alpha_X"], label="Alpha_X")
    ax2.plot(lattice_df["S_cm"], lattice_df["Alpha_Y"], label="Alpha_Y")
    ax2.set_ylabel("Alpha_X, Alpha_Y", fontsize=fs)
    ax2.set_xlabel("S_cm")
    ax2.legend()
    for ax in (ax0, ax1, ax2):
        for p in camera_positions:
            ax.axvline(p, color='green')
    y_pos_annotate = np.mean(ax1.get_ylim())
    for name, p in zip(camera_names, camera_positions):
        ax1.annotate(name, (p, y_pos_annotate))
    plt.show()


def get_cameras_df(lattice_df):
    cameras_df = pd.DataFrame({
                                "Name": camera_names,
                                "S_cm": camera_positions
                            })
    for col in lattice_df.columns[1:]:
        cameras_df[col] = np.interp(cameras_df["S_cm"],
                                    lattice_df["S_cm"],
                                    lattice_df[col])
    return cameras_df
