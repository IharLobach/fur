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
    0, 0.75*d, 0), shaftwidth=R, color=axis_col, headlength=headlength, headwidth=headwidth)
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
col_hot, col_cold, col_neutral = color.red, color.blue, color.green
cols_top = [col_neutral, col_hot, col_neutral, col_cold]
cols_bottom = [col_neutral, col_cold, col_neutral, col_hot]
mag_len = 1
mag_height = 2
mag_width = 4
mag_dimensions = vector(mag_width, mag_height, mag_len)
mag_half_dimensions = vector(mag_width, mag_height, mag_len/2)
gap = 2
period_len = 4*mag_len


def period(z):
    for i, ct, cb in zip(range(4), cols_top, cols_bottom):
        box(pos=vector(0, gap/2+mag_height/2, z+i*mag_len),
            size=mag_dimensions, color=ct)
        box(pos=vector(0, -gap/2-mag_height/2, z+i*mag_len),
            size=mag_dimensions, color=cb)

def beginning_period(z):
    box(pos=vector(0, gap/2+mag_height/2, z),
        size=mag_half_dimensions, color=col_cold)
    box(pos=vector(0, -gap/2-mag_height/2, z),
        size=mag_half_dimensions, color=col_hot)


def ending_period(z):
    l = 2
    for i, ct, cb in zip(range(l), cols_top, cols_bottom):
        if i==l-1:
            dim = mag_half_dimensions
            d = mag_len/4
        else:
            dim = mag_dimensions
            d = 0
        box(pos=vector(0, gap/2+mag_height/2, z+i*mag_len-d),
            size=dim, color=ct)
        box(pos=vector(0, -gap/2-mag_height/2, z+i*mag_len-d),
            size=dim, color=cb)


n_full_per = 2
und_len = (n_full_per+1/2)*period_len
z0 = -period_len*n_full_per/2
beginning_period(z0-3*mag_len/4)
for p in range(n_full_per):
    period(z0+p*period_len)
ending_period(z0+n_full_per*period_len)


#----------------------------------------------------
# electron trajectory
c = curve(color=color.green, radius=R)
z_und_start = -und_len/2
z_und_end = und_len/2
zs = np.linspace(z_und_start, z_und_end, 100)
for z in zs:
    c.append(vector(-np.cos(2*np.pi*z/period_len), 0, z))


#---------------------------------------------------
# radiation
cone(pos=zaxis.pos+0.9*zaxis.axis,
     axis=-zaxis.axis,
     radius=1,
     color=color.yellow,
     opacity=0.5)





