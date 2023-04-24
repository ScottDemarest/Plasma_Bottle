import numpy as np

# Normalized length scale (wavelength)
l0 = 2*np.pi

Lx = 6*l0
x_min = 0

Ly = 6*l0
y_min = 0

Lz= 6*l0
z_min = 0

B0 = 1
R = .5*l0

#Ring 1
x0 = Lx/2
y0 = Ly/2
z0 = 0

#Ring 2
x1 = Lx/2
y1 = Ly/2
z1 = Lz

# Resolution
resx         = 100.    # nb of cells in one laser wavelength
resy         = 100.                     
resz         = 100.

# Domain
dx           = l0/resx
dy           = l0/resy
dz           = l0/resz

Nx = Lx/dx
Ny = Ly/dy
Nz = Lz/dz