#!/home/ilobach/anaconda3/bin/python

import sys
import numpy as np
import os
from wiggler_radiation.Wigrad.wigrad_generator import get_rad_mesh_tuple
dirpath = os.path.dirname(os.path.abspath(__file__))
foldername = os.path.basename(dirpath)
en = int(foldername[:-3])
i3d = np.load(f'/scratch/midway2/ilobach/field_wigrad_{en}MeV.npy')
rad_mesh_tuple = get_rad_mesh_tuple() 
x1d, y1d, l1d = rad_mesh_tuple
nx, ny, nl = [len(v) for v in rad_mesh_tuple]
def get_step(arr):
    return (arr[-1]-arr[0])/(len(arr)-1)
dx, dy, dl = [get_step(v) for v in get_rad_mesh_tuple()]
tot = dx*dy*dl*np.sum(i3d)
x1_2D, x2_2D = np.meshgrid(x1d, x1d)
x1m2_2d = x1_2D-x2_2D
y1_2D, y2_2D = np.meshgrid(y1d, y1d)
y1m2_2d = y1_2D-y2_2D
def get_Mxy(sx, sy):
    res = 0
    for slice2d, lmda in zip(i3d, l1d):
        Eexpx = np.tensordot(slice2d, lmda*np.exp(-(sx*2*np.pi/lmda*x1m2_2d)**2), (1,0))
        Eexpy = np.tensordot(slice2d, lmda*np.exp(-(sy*2*np.pi/lmda*y1m2_2d)**2), (0,1))
        res += np.tensordot(Eexpx, Eexpy.T)
    return tot**2/(dx*dy*dx*dy*dl*1/(2*np.sqrt(np.pi)*1e4)*res)
sx_range = np.arange(350, 1300, 10)  # um
sy_range = np.arange(30, 370, 5)  # um
sx2d, sy2d = np.meshgrid(sx_range, sy_range)
sxsy_tuples = np.vstack((sx2d.ravel(), sy2d.ravel())).T
idx = int(sys.argv[1])
res = get_Mxy(*sxsy_tuples[idx])
with open(os.path.join(dirpath,"results","{}.out".format(idx)), 'w') as f:
    f.write("{}".format(res))