class Wiggler():
    def __init__(self, K_peak, N_periods, lambda_wiggler):
        self.K_peak = K_peak
        self.N_periods = N_periods
        self.lambda_wiggler = lambda_wiggler
        self.aux_const = 1+self.K_peak**2/2


class Beam():
    def __init__(self, gamma, Ibeam_mA):
        self.gamma = gamma
        self.Ibeam_mA = Ibeam_mA