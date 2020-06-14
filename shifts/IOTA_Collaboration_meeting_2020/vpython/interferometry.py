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
m1 = Mirror(pos=first_bs+vz)
bs2 = BeamSplitter(pos=first_bs+vx+vz)
m2 = Mirror(pos=first_bs+vx)



#---------------------------------------
# light ray
kwlr = dict(
    radius=3*R,
    color=color.yellow,
    emissive=True,
    opacity=0.5
)
und_end_pos = vec(0, 0, und.z_end)
cylinder(pos=und_end_pos, axis=m1.pos-und_end_pos, **kwlr)
ray2 = cylinder(pos=m1.pos, axis=1.8*vx, **kwlr)
cylinder(pos=bs1.pos, axis=vx, **kwlr)
ray1 = cylinder(pos=m2.pos, axis=1.8*vz, **kwlr)


#------------------------------------------
# Detectors
from detector_vpython import Detector
d1 = Detector(ray1.pos+ray1.axis)
d1.extrusion.rotate(angle=pi/2, axis=vec(1, 0, 0), origin=d1.pos)
d2 = Detector(ray2.pos+ray2.axis)
d2.extrusion.rotate(angle=-pi/2, axis=vec(0, 0, 1), origin=d2.pos)


#------------------------------------------
# labels
kwl = dict(height=32, border=4, font='sans')
label(pos=m1.pos,
      text='Mirror', xoffset=20,
      yoffset=-100, space=30, **kwl)
label(pos=bs1.pos,
      text='Beam splitter', xoffset=-40,
      yoffset=-100, space=40, **kwl)
label(pos=d1.pos,
      text='Detector', xoffset=20,
      yoffset=-80, space=20, **kwl)
