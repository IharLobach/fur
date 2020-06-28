import numpy as np
from vpython import *
scene.background = color.white
scene.width = 1280
scene.height = 720

#---------------------------------------------
#  axes
L=10
R = L/200
d = L-2
axis_col = color.black
headlength = 6*R
headwidth = 4*R
zaxis = arrow(pos=vec(0, 0, 0), axis=vec(
    0, 0, 1.25*d), shaftwidth=R, color=axis_col, headlength=headlength, headwidth=headwidth)
xaxis = arrow(pos=vec(0, 0, 0), axis=vec(
    d, 0, 0), shaftwidth=R, color=axis_col, headlength=headlength, headwidth=headwidth)
yaxis = arrow(pos=vec(0, 0, 0), axis=vec(
    0, 0.6*d, 0), shaftwidth=R, color=axis_col, headlength=headlength, headwidth=headwidth)
k = 1.02
h = 0.05*L
text(pos=zaxis.pos+k*zaxis.axis, text='z', height=h,
     align='center', billboard=True, emissive=True,
     color=color.black)
text(pos=xaxis.pos+k*xaxis.axis, text='x', height=h,
     align='center', billboard=True, emissive=True,
     color=color.black)
text(pos=yaxis.pos+k*yaxis.axis, text='y', height=h,
     align='center', billboard=True, emissive=True,
     color=color.black)

#----------------------------------------------
#undulator
from undulator_vpython import Undulator
und = Undulator()


#----------------------------------------------------
# electron trajectory
from electron_trajectory_vpython import ElectronTrajectory
el_traj = ElectronTrajectory(und.z_start, und.z_end, und.period_len, radius=2*R)


#---------------------------------------------------
# radiation
rel_len = 0.9
light_cone = cone(pos=rel_len*(zaxis.pos+zaxis.axis),
                  axis=-rel_len*zaxis.axis,
                  radius=1,
                  color=color.yellow,
                  opacity=0.5)


#---------------------------------------------------
# labels
kwl = dict(height=32, border=4, font='sans')
light_label = label(pos=light_cone.pos,
                    text='Radiation', xoffset=20,
                    yoffset=100, space=60, **kwl)

electron_label = label(pos=vec(0, 0, und.z_start),
                       text='Electron trajectory', xoffset=80,
                       yoffset=20, space=0, **kwl)

magnets_label = label(pos=vec(
                                -und.mag_width/2,
                                und.gap/2+und.mag_height,
                                und.z_end-2*und.mag_len),
                      text='Permanent magnets', xoffset=-20,
                      yoffset=40, space=0, **kwl)




