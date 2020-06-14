from vpython import *


class Undulator:
    def __init__(self, mag_len=1, mag_height=2, mag_width=4, gap=2):
        self.mag_len = mag_len
        self.mag_height = mag_height
        self.mag_width = mag_width
        self.gap = gap
        col_hot, col_cold, col_neutral = color.red, color.blue, color.green
        cols_top = [col_neutral, col_hot, col_neutral, col_cold]
        cols_bottom = [col_neutral, col_cold, col_neutral, col_hot]
        mag_dimensions = vector(self.mag_width, self.mag_height, self.mag_len)
        mag_half_dimensions = vector(self.mag_width, self.mag_height,
                                     self.mag_len/2)
        self.period_len = 4*self.mag_len

        def period(z):
            for i, ct, cb in zip(range(4), cols_top, cols_bottom):
                box(pos=vector(0, self.gap/2+self.mag_height/2, z+i*self.mag_len),
                    size=mag_dimensions, color=ct)
                box(pos=vector(0, -self.gap/2-self.mag_height/2, z+i*self.mag_len),
                    size=mag_dimensions, color=cb)

        def beginning_period(z):
            box(pos=vector(0, self.gap/2+self.mag_height/2, z),
                size=mag_half_dimensions, color=col_cold)
            box(pos=vector(0, -self.gap/2-self.mag_height/2, z),
                size=mag_half_dimensions, color=col_hot)

        def ending_period(z):
            l = 2
            for i, ct, cb in zip(range(l), cols_top, cols_bottom):
                if i == l-1:
                    dim = mag_half_dimensions
                    d = self.mag_len/4
                else:
                    dim = mag_dimensions
                    d = 0
                box(pos=vector(0, self.gap/2+self.mag_height/2, z+i*self.mag_len-d),
                    size=dim, color=ct)
                box(pos=vector(0, -self.gap/2-self.mag_height/2, z+i*self.mag_len-d),
                    size=dim, color=cb)

        n_full_per = 2
        und_len = (n_full_per+1/2)*self.period_len
        z0 = -self.period_len*n_full_per/2
        beginning_period(z0-3*self.mag_len/4)

        for p in range(n_full_per):
            period(z0+p*self.period_len)
        ending_period(z0+n_full_per*self.period_len)
        self.z_start = -und_len/2
        self.z_end = und_len/2
