import numpy as np
from vpython import *


class Mirror:
    def __init__(self, pos, a=2, mir_thickness=0.3):
        self.pos = pos
        mir = box(pos=pos,
                  axis=mir_thickness*vec(np.sqrt(0.5), 0, -np.sqrt(0.5)),
                  width=np.sqrt(2)*a,
                  height=1.04*a,
                  color=color.blue)


if __name__ == "__main__":
    scene.background = color.white
    scene.width = 1280
    scene.height = 720
    Mirror(vec(1, 1, 1))
