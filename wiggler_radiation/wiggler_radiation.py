from scipy.special import sinc, jv
import numpy as np


class WigglerRadiationSimulator():
    alpha = 1/137

    def __init__(self, wiggler, beam, aperture=None, harmonic=1,
                 bessel_cutoff=20):
        self.wiggler = wiggler
        self.beam = beam
        self.harmonic = harmonic
        self.bessel_cutoff = bessel_cutoff
        self.lambda1_um = 1e6*self.wiggler.lambda_wiggler_m\
            / 2/self.beam.gamma**2*self.wiggler.aux_const
        self.aperture = aperture

    def calc_Ih_1el(self, lambda_um, theta_x, theta_y):
        dw = self.lambda1_um/lambda_um-1
        theta_eff_2 = self.beam.gamma**2*(theta_x**2+theta_y**2)
        A = self.wiggler.aux_const+theta_eff_2
        Y = self.harmonic*self.wiggler.K_peak**2/4/A
        X = 2*self.harmonic*self.beam.gamma*self.wiggler.K_peak*theta_x/A
        L = sinc(self.wiggler.N_periods
                 * (theta_eff_2+dw*A)/self.wiggler.aux_const)**2
        sum1 = 0
        sum2 = 0
        sum3 = 0
        for p in range(-self.bessel_cutoff, self.bessel_cutoff+1):
            sum1 += jv(self.harmonic+2*p, X)*jv(p, Y)
            sum2 += jv(self.harmonic+2*p-1, X)*jv(p, Y)
            sum3 += jv(self.harmonic+2*p+1, X)*jv(p, Y)
        return self.alpha*self.harmonic**2*self.beam.gamma**2\
            * self.wiggler.N_periods**2\
            / A**2*L*np.absolute(2*self.beam.gamma*sum1
                                 - self.wiggler.K_peak*(sum2+sum3))**2\
            * self.lambda1_um/lambda_um**2

    def calc_x_y_integral_1el(self, lambda_um):
        if self.aperture is None:
            raise TypeError('Aperture for this instance of'
                            ' WigglerRadiationSimulator is not defined.')
        if self.aperture.mesh_size_1D is None:
            raise TypeError("Aperture's mesh is not defined.")
        Ih_vals = [self.calc_Ih_1el(lambda_um, x, y)
                   if self.aperture.accepted(x, y) else 0 for x, y in
                   zip(self.aperture.mesh_xs, self.aperture.mesh_ys)]
        return self.aperture.step**2*sum(Ih_vals)
