from scipy.special import sinc, jv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from enum import Enum


class WigglerRadiationSimulator():
    alpha = 1/137

    def __init__(self, wiggler, beam, mesh, harmonics=[1],
                 bessel_cutoff=10,
                 aperture='ellips',
                 only_calc_sum_of_both_polarizations=True):
        self.wiggler = wiggler
        self.beam = beam
        self.harmonics = harmonics
        self.bessel_cutoff = bessel_cutoff
        self.aperture = aperture
        self.only_calc_sum_of_both_polarizations = \
            only_calc_sum_of_both_polarizations
        self.lambda1_um = 1e6*self.wiggler.lambda_wiggler_m\
            / 2/self.beam.gamma**2*self.wiggler.aux_const
        self.x_range, self.y_range, self.lambda_range = mesh
        self.x_2D = np.tile(self.x_range, (self.n_y, 1))
        self.y_2D = np.tile(self.y_range.reshape(-1, 1), (1, self.n_x))
        self.x_step = (self.x_range[-1]-self.x_range[0])/(self.n_x-1)
        self.y_step = (self.y_range[-1]-self.y_range[0])/(self.n_y-1)
        self.lambda_step = (self.lambda_range[-1]-self.lambda_range[0])\
            / (self.n_lambda-1)

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
    def x_3D(self):
        return np.tile(self.x_2D, (self.n_lambda, 1, 1))

    @property
    def y_3D(self):
        return np.tile(self.y_2D, (self.n_lambda, 1, 1))

    def __calc_photon_flux_on_meshgrid_one_harmonic(self, harmonic):
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
        bessel_part_x = aux_factor \
            * np.absolute(2*self.beam.gamma*x_2D*sum1
                          - self.wiggler.K_peak*(sum2+sum3))**2
        bessel_part_y = aux_factor*np.absolute(2*self.beam.gamma*y_2D*sum1)**2
        dw_arr = self.lambda1_um/self.lambda_range-1
        L = [(sinc(self.wiggler.N_periods*(r2_2D+dw*A)/self.wiggler.aux_const)
              / sinc((r2_2D+dw*A)/self.wiggler.aux_const))
             ** 2 / l**2 for dw, l in zip(dw_arr, self.lambda_range)]
        L = np.asarray(L)
        if self.only_calc_sum_of_both_polarizations:
            return (bessel_part_x + bessel_part_y)*L
        else:
            return bessel_part_x*L, bessel_part_y*L

    def calc_photon_flux_on_meshgrid(self):
        res = \
            self.__calc_photon_flux_on_meshgrid_one_harmonic(self.harmonics[0])
        for h in self.harmonics[1:]:
            res += self.__calc_photon_flux_on_meshgrid_one_harmonic(h)
        if self.only_calc_sum_of_both_polarizations:
            self.__full_photon_flux_3D = res
        else:
            self.photon_flux_3D_polarization_x = res[0]
            self.photon_flux_3D_polarization_y = res[1]
        del res
        if self.aperture == 'ellips':
            x_max = max(self.x_range)
            y_max = max(self.y_range)
            elliptic_aperture = \
                (self.x_3D**2/x_max**2+self.y_3D**2/y_max**2) < 1
            if self.only_calc_sum_of_both_polarizations:
                self.__full_photon_flux_3D = \
                    np.where(elliptic_aperture,
                             self.__full_photon_flux_3D,
                             0)
            else:
                self.photon_flux_3D_polarization_x = \
                    np.where(elliptic_aperture,
                             self.photon_flux_3D_polarization_x,
                             0)
                self.photon_flux_3D_polarization_y = \
                    np.where(elliptic_aperture,
                             self.photon_flux_3D_polarization_y,
                             0)

    @property
    def full_photon_flux_3D(self):
        if self.only_calc_sum_of_both_polarizations:
            return self.__full_photon_flux_3D
        else:
            return self.photon_flux_3D_polarization_x \
                + self.photon_flux_3D_polarization_y

    def __extend_angular_mesh_using_symmetries(self):
        x1 = self.x_2D
        y1 = self.y_2D
        x2 = np.flip(-x1, axis=1)
        y2 = y1
        x21 = np.hstack((x2, x1))
        y21 = np.hstack((y2, y1))
        x34 = x21
        y34 = np.flip(-y21, axis=0)
        self.x_2D = np.vstack((x34, x21))
        self.y_2D = np.vstack((y34, y21))

    def __extend_photon_flux_using_symmetries(self, z):
        z1 = z
        z2 = np.flip(z1, axis=2)
        z21 = np.concatenate((z2, z1), axis=2)
        z34 = np.flip(z21, axis=1)
        return np.concatenate((z34, z21), axis=1)

    def extend_results_using_symmetries(self):
        self.__extend_angular_mesh_using_symmetries()
        if self.only_calc_sum_of_both_polarizations:
            self.__full_photon_flux_3D = \
                self.__extend_photon_flux_using_symmetries(
                    self.full_photon_flux_3D)
        else:
            self.photon_flux_3D_polarization_x = \
                self.__extend_photon_flux_using_symmetries(
                    self.photon_flux_3D_polarization_x)
            self.photon_flux_3D_polarization_y = \
                self.__extend_photon_flux_using_symmetries(
                    self.photon_flux_3D_polarization_y)

    def __show_angular_distribution(self, z):
        fig = plt.figure(figsize=[10, 10])
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(self.x_2D,
                               self.y_2D,
                               z,
                               cmap=cm.coolwarm,
                               linewidth=0,
                               antialiased=False)
        ax.set_xlabel(r"$\theta_x$, rad")
        ax.set_ylabel(r"$\theta_y$, rad")
        return ax

    @property
    def angular_distribution_sum_both_polatizations(self):
        return self.lambda_step *\
            np.apply_over_axes(np.sum, self.full_photon_flux_3D, [0])[0]
    
    @property
    def angular_distribution_x_polatization(self):
        return self.lambda_step *\
            np.apply_over_axes(np.sum, self.photon_flux_3D_polarization_x, [0])[0]

    @property
    def angular_distribution_y_polatization(self):
        return self.lambda_step *\
            np.apply_over_axes(np.sum, self.photon_flux_3D_polarization_y, [0])[0]

    def show_angular_distribution_sum_both_polatizations(self,
                                                         polarization='sum',
                                                         index_of_lambda=None):
        if index_of_lambda is not None:
            if polarization == 'sum':
                z = self.full_photon_flux_3D
            elif polarization == 'x':
                z = self.photon_flux_3D_polarization_x
            elif polarization == 'y':
                z = self.photon_flux_3D_polarization_y
            else:
                raise ValueError("Unknown polarization type. Choose from 'x', 'y' and 'sum'.")
            ax = self.__show_angular_distribution(z[index_of_lambda])
            dim = r" $\frac{\mathrm{Ph}}{\mathrm{rad}^{-2}\mathrm{nm}}$"
            ax.set_zlabel(r"$\frac{dN}{d\theta_x d\theta_y d\lambda}$,"+dim)
            lambda_val = self.lambda_range[index_of_lambda]
            ax.set_title(r"Polarization: "+polarization+", "+
                         r"$\lambda$ = "
                         + "{:.3f} um".format(lambda_val))
            plt.show()
        else:
            if polarization == 'sum':
                z = self.angular_distribution_sum_both_polatizations
            elif polarization == 'x':
                z = self.angular_distribution_x_polatization
            elif polarization == 'y':
                z = self.angular_distribution_y_polatization
            else:
                raise ValueError("Unknown polarization type. Choose from 'x', 'y' and 'sum'.")
            ax = self.__show_angular_distribution(z)
            dim = r" $\frac{\mathrm{Ph}}{\mathrm{rad}^{-2}}$"
            ax.set_zlabel(r"$\frac{dN}{d\theta_x d\theta_y}$,"+dim)
            ax.set_title("Polarization: "+polarization)
            plt.show()
    

