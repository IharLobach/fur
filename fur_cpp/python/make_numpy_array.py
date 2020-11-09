import numpy as np
nx, ny, nz = 100, 200, 500
data = np.arange(nx*ny*nz).reshape(nz, ny, nx)
data.dump('data.npy')
