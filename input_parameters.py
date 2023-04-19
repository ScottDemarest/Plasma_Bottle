import numpy as np

# Normalized length scale (wavelength)
l0 = 2*np.pi

Lx = 4*l0
x_min = 0

Ly = 4*l0
y_min = 0

Lz= 10*l0
z_min = 0

B0 = 1
R = 0.7*l0

#Ring 1
x0 = Lx/2
y0 = Ly/2
z0 = 1*l0

#Ring 2
x1 = Lx/2
y1 = Ly/2
z1 = 9*l0

# Resolution
resx         = 4.    # nb of cells in one laser wavelength
resy         = 4.                     
resz         = 4.