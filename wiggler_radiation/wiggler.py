class Wiggler():
    def __init__(self, K_peak=1, N_periods=10, lambda_wiggler_m=0.055):
        self.K_peak = K_peak
        self.N_periods = N_periods
        self.lambda_wiggler_m = lambda_wiggler_m
        self.aux_const = 1+self.K_peak**2/2