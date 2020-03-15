from scipy.special import sinc, jv
import numpy as np


class WigglerRadiationSimulator():
    alpha = 1/137

    def __init__(self, wiggler, beam, mesh, harmonic=1,
                 bessel_cutoff=20):
        self.wiggler = wiggler
        self.beam = beam
        self.harmonic = harmonic
        self.bessel_cutoff = bessel_cutoff
        self.lambda1_um = 1e6*self.wiggler.lambda_wiggler_m\
            / 2/self.beam.gamma**2*self.wiggler.aux_const
        self.x_range, self.y_range, self.lambda_range = mesh

    @property
    def n_x(self):
        return len(self.x_range)

    @property
    def n_y(self):
        return len(self.y_range)

    @property
    def n_lambda(self):
        return len(self.n_lambda)

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
        p = -self.bessel_cutoff
        jv2pm1 = jv(self.harmonic+2*p-1, X)
        for p in range(-self.bessel_cutoff, self.bessel_cutoff+1):
            jvpY = jv(p, Y)
            sum1 += jv(self.harmonic+2*p, X)*jvpY
            sum2 += jv2pm1*jvpY
            jv2pp1 = jv(self.harmonic+2*p+1, X)
            sum3 += jv2pp1*jvpY
            jv2pm1 = jv2pp1
        return self.alpha*self.harmonic**2*self.beam.gamma**2\
            * self.wiggler.N_periods**2\
            / A**2*L*np.absolute(2*self.beam.gamma*theta_x*sum1
                                 - self.wiggler.K_peak*(sum2+sum3))**2\
            * self.lambda1_um/lambda_um**2

    def calc_Ih_on_mesh(self):
        x_2D = np.tile(self.x_range, (self.n_y, 1))
        y_2D = np.tile(self.y_range.reshape(-1, 1), (1, self.n_x))
        r2_2D = x_2D**2+y_2D**2
        A = self.wiggler.aux_const+self.beam.gamma**2*r2_2D
        Y = self.harmonic*self.wiggler.K_peak**2/4/A
        X = 2*self.harmonic*self.beam.gamma*self.wiggler.K_peak*x_2D/A
        sum1 = 0
        sum2 = 0
        sum3 = 0
        p = -self.bessel_cutoff
        jv2pm1 = jv(self.harmonic+2*p-1, X)
        for p in range(-self.bessel_cutoff, self.bessel_cutoff+1):
            jvpY = jv(p, Y)
            sum1 += jv(self.harmonic+2*p, X)*jvpY
            sum2 += jv2pm1*jvpY
            jv2pp1 = jv(self.harmonic+2*p+1, X)
            sum3 += jv2pp1*jvpY
            jv2pm1 = jv2pp1
        bessel_part = self.alpha*self.harmonic**2*self.beam.gamma**2\
            * self.wiggler.N_periods**2\
            / A**2*np.absolute(2*self.beam.gamma*x_2D*sum1
                               - self.wiggler.K_peak*(sum2+sum3))**2\
            * self.lambda1_um
        dw_arr = self.lambda1_um/self.lambda_range-1
        L = [sinc(self.wiggler.N_periods*(r2_2D+dw*A)/self.wiggler.aux_const)\
             ** 2 / l**2 for dw, l in zip(dw_arr, self.lambda_range)]
        L = np.asarray(L)
        return bessel_part*L



