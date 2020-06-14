import numpy as np
from vpython import *


class BeamSplitter:
    def __init__(self, pos, a=2, mir_thickness=0.1):
        self.pos = pos
        b = box(pos=pos, axis=a*vec(0, 0, 1), height=a, width=a,
                color=color.cyan, opacity=0.8)
        mir = box(pos=pos,
                  axis=mir_thickness*vec(np.sqrt(0.5), 0, -np.sqrt(0.5)),
                  width=np.sqrt(2)*a,
                  height=1.04*a,
                  color=color.red)


if __name__=="__main__":
    scene.background = color.white
    scene.width = 1280
    scene.height = 720
    BeamSplitter(vec(1, 1, 1))
