import numpy as np
from vpython import *
scene.background = color.white
scene.width = 1280
scene.height = 720
from electron_trajectory_vpython import ElectronTrajectory
R = 10/200
#----------------------------------------------
#undulator
from undulator_vpython import Undulator
und = Undulator()
#---------------------------------------------------
# electron trajectory
el_traj = ElectronTrajectory(
    und.z_start, und.z_end, und.period_len, radius=2*R)

# interferometer
from mirror_vpython import Mirror
from beam_splitter_vpython import BeamSplitter
first_bs = vec(0, 0, 12)
interf_side = 5
vz = vec(0, 0, interf_side)
vx = vec(interf_side, 0, 0)
bs1 = BeamSplitter(pos=first_bs)
m1pos = first_bs+vz
m2pos = first_bs+vx

#---------------------------------------
# light ray
kwlr = dict(
    radius=3*R,
    color=color.yellow,
    emissive=True,
    opacity=0.5
)
und_end_pos = vec(0, 0, und.z_end)
ray1 = cylinder(pos=und_end_pos, axis=first_bs+vz-und_end_pos, **kwlr)
ray2 = cylinder(pos=bs1.pos, axis=vx, **kwlr)


#------------------------------------------
# Detectors
from detector_vpython import Detector
d1 = Detector(first_bs+vz)
d1.extrusion.rotate(angle=pi/2, axis=vec(1, 0, 0), origin=d1.pos)
d2 = Detector(first_bs+vx)
d2.extrusion.rotate(angle=-pi/2, axis=vec(0, 0, 1), origin=d2.pos)


#------------------------------------------
# labels
kwl = dict(height=24, border=4, font='sans')
label(pos=bs1.pos,
      text='Beam splitter', xoffset=-40,
      yoffset=-100, space=40, **kwl)
label(pos=d2.pos,
      text='SPAD2', xoffset=-20,
      yoffset=80, space=20, **kwl)
label(pos=d1.pos,
      text='SPAD1', xoffset=80,
      yoffset=-10, space=40, **kwl)
