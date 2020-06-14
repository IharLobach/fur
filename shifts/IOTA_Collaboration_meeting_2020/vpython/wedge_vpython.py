import numpy as np
from vpython import *


class Wedge:
    def __init__(self, pos):
        self.pos = pos
        tri = shapes.points(pos=[[-0.25, -1], [0.25, -1], [-0.25, 1]])
        pa = paths.line(vec(0, 0, -1), vec(0, 0, 1))
        self.extrusion = extrusion(pos=pos,
                                   path=pa,
                                   shape=tri,
                                   color=color.magenta,
                                   opacity=0.5)


if __name__ == "__main__":
    scene.background = color.white
    scene.width = 1280
    scene.height = 720
    pos = vec(1, 1, 1)
    Wedge(pos)
