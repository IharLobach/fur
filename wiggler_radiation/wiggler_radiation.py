from scipy.special import sinc, jv
import numpy as np


class WigglerRadiationSimulator():
    alpha = 1/137

    def __init__(self, wiggler, beam, mesh, harmonics=[1],
                 bessel_cutoff=10):
        self.wiggler = wiggler
        self.beam = beam
        self.harmonics = harmonics
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
        return len(self.lambda_range)

    @property
    def x_2D(self):
        return np.tile(self.x_range, (self.n_y, 1))

    @property
    def y_2D(self):
        return np.tile(self.y_range.reshape(-1, 1), (1, self.n_x))

    @property
    def x_3D(self):
        return np.tile(self.x_2D, (self.n_lambda, 1, 1))

    @property
    def y_3D(self):
        return np.tile(self.y_2D, (self.n_lambda, 1, 1))

    def calc_Ih_on_mesh_one_harmonic(self, harmonic):
        x_2D = self.x_2D
        y_2D = self.y_2D
        r2_2D = self.beam.gamma**2*(x_2D**2+y_2D**2)
        A = self.wiggler.aux_const+r2_2D
        Y = harmonic*self.wiggler.K_peak**2/4/A
        X = 2*harmonic*self.beam.gamma*self.wiggler.K_peak*x_2D/A
        sum1 = 0
        sum2 = 0
        sum3 = 0
        p = -self.bessel_cutoff
        jv2pm1 = jv(harmonic+2*p-1, X)
        for p in range(-self.bessel_cutoff, self.bessel_cutoff+1):
            jvpY = jv(p, Y)
            sum1 += jv(harmonic+2*p, X)*jvpY
            sum2 += jv2pm1*jvpY
            jv2pp1 = jv(harmonic+2*p+1, X)
            sum3 += jv2pp1*jvpY
            jv2pm1 = jv2pp1
        aux_factor = self.alpha*harmonic**2*self.beam.gamma**2\
            * self.wiggler.N_periods**2\
            / A**2*self.lambda1_um
        bessel_part_x = aux_factor*np.absolute(2*self.beam.gamma*x_2D*sum1
            - self.wiggler.K_peak*(sum2+sum3))**2
        bessel_part_y = aux_factor*np.absolute(2*self.beam.gamma*y_2D*sum1)**2
        dw_arr = self.lambda1_um/self.lambda_range-1
        L = [(sinc(self.wiggler.N_periods*(r2_2D+dw*A)/self.wiggler.aux_const)
              / sinc((r2_2D+dw*A)/self.wiggler.aux_const))
             ** 2 / l**2 for dw, l in zip(dw_arr, self.lambda_range)]
        L = np.asarray(L)
        return bessel_part_x*L, bessel_part_y*L

    def calc_Ih_on_mesh(self):
        res = self.calc_Ih_on_mesh_one_harmonic(self.harmonics[0])
        for h in self.harmonics[1:]:
            res += self.calc_Ih_on_mesh_one_harmonic(h)
        self.res = res
        return res
