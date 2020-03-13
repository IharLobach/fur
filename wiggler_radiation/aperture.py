import numpy as np

# https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=12767:
# Clear Aperture >90% of Diameter
# cl_ap = 1.0
# lens_diamter = 5.08  # cm
# distance_to_undulator = 350  # cm
# theta_max = 0.007257142857142857


class EllipticAperture():
    def __init__(self, theta_max=0.007257142857142857, mesh_size_1D=None):
        self.theta_x_max = theta_max
        self.theta_y_max = theta_max/np.sqrt(theta_max)
        self.mesh_size_1D = mesh_size_1D
        if self.mesh_size_1D:
            self.step = 2*theta_max/mesh_size_1D
            self.x_range = np.linspace(-self.theta_x_max, self.theta_x_max,
                                       mesh_size_1D)+self.step/2
            self.y_range = np.linspace(-self.theta_y_max, self.theta_y_max,
                                       mesh_size_1D)+self.step/2
            self.mesh_xs =\
                np.hstack((self.x_range for _ in range(mesh_size_1D)))
            self.mesh_ys =\
                np.hstack(np.vstack(
                    (self.y_range for _ in range(mesh_size_1D))).T)

    def accepted(self, theta_x, theta_y):
        return 1-(theta_x/self.theta_x_max)**2\
            - (theta_y/self.theta_y_max)**2 > 0
