import numpy as np
from vpython import *


class Detector:
    def __init__(self, pos, radius=1):
        cr = shapes.circle(radius=radius, angle1=0, angle2=pi/2)
        rotation = paths.circle(radius=0.05)
        self.pos = pos
        self.extrusion = extrusion(pos=pos,
                        path=rotation,
                        shape=cr,
                        color=color.purple)


if __name__ == "__main__":
    scene.background = color.white
    scene.width = 1280
    scene.height = 720
    Detector(vec(1, 1, 1), 2)
#det.rotate(angle=pi/2, axis=vec(0,0,1), origin=vec(0,0,0))
