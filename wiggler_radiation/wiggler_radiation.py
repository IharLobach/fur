from scipy.special import sinc, jv
import numpy as np


class WigglerRadiationSimulator():
    alpha = 1/137

    def __init__(self, wiggler, beam, n0=1, nb=20):
        self.wiggler = wiggler
        self.beam = beam
        self.n0 = n0
        self.nb = nb

    def calc_Ih(self, dw, theta_x, theta_y):
        A = self.wiggler.aux_const+self.beam.gamma**2*(theta_x**2+theta_y**2)
        Y = self.n0*self.wiggler.K_peak**2/4/A
        X = 2*self.n0*self.beam.gamma*self.wiggler.K_peak*theta_x/A
        L = sinc(self.wiggler.K_peak*((1+dw)*A/self.wiggler.aux_const-1))**2
        sum1 = 0
        sum2 = 0
        sum3 = 0
        for p in range(-self.nb, self.nb+1):
            sum1 += jv(self.n0+2*p, X)*jv(p, Y)
            sum2 += jv(self.n0+2*p-1, X)*jv(p, Y)
            sum3 += jv(self.n0+2*p+1, X)*jv(p, Y)
        return self.alpha*self.n0**2*self.beam.gamma**2*self.wiggler.K_peak**2\
            / A**2*L*np.absolute(2*self.beam.gamma*sum1-self.wiggler.K_peak*(sum2+sum3))**2
