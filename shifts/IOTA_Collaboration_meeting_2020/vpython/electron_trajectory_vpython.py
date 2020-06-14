from vpython import *
import numpy as np


class ElectronTrajectory():
    def __init__(self, z_start, z_end, period_len, radius, amplitude=1):
        c = curve(color=color.purple, radius=radius)
        z_und_start = z_start
        z_und_end = z_end
        zs = np.linspace(z_und_start, z_und_end, 100)
        for z in zs:
            c.append(vector(-np.cos(2*np.pi*z/period_len), 0, z))
