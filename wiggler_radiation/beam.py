class Beam():
    electron_charge = 1.602e-19
    revolution_period_sec = 1/7.5e6

    def __init__(self, Ibeam_mA, sigma_x_um, sigma_y_um, sigma_z_cm,
                 gamma=100/0.511):
        self.gamma = gamma
        self.Ibeam_mA = Ibeam_mA
        self.sigma_x_um = sigma_x_um
        self.sigma_y_um = sigma_y_um
        self.sigma_z_cm = sigma_z_cm

    @property
    def sigma_z_um(self):
        return 1e4*self.sigma_z_cm

    @property
    def bunch_charge(self):
        return self.Ibeam_mA*0.001*self.revolution_period_sec\
            / self.electron_charge

    @property
    def number_of_electrons(self):
        return self.bunch_charge/self.electron_charge
