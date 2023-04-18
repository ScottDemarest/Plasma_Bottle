from scipy.integrate import quad
import numpy as np
import uproot
from math import pi

def get_B(x, y, z, x0, y0, z0, B0, R):
    # Ring lies in xy plane
    def integrand_x(t, x, y, z):
        r_mag = np.sqrt((x-R*np.cos(t))**2 + (y-R*np.sin(t))**2 + z**2)
        r_3 = r_mag**3
        return R*z*np.cos(t)/r_3

    def integrand_y(t, x, y, z):
        r_mag = np.sqrt((x-R*np.cos(t))**2 + (y-R*np.sin(t))**2 + z**2)
        r_3 = r_mag**3
        return R*z*np.sin(t)/r_3

    def integrand_z(t, x, y, z):
        r_mag = np.sqrt((x-R*np.cos(t))**2 + (y-R*np.sin(t))**2 + z**2)
        r_3 = r_mag**3
        return (R - y*np.sin(t) - x*np.cos(t))/r_3
    
    
    a = 0
    b = 2*np.pi
    
    bx = B0*quad(integrand_x, a, b, args=(x-x0, y-y0, z-z0))[0]
    by = B0*quad(integrand_y, a, b, args=(x-x0, y-y0, z-z0))[0]
    bz = B0*quad(integrand_z, a, b, args=(x-x0, y-y0, z-z0))[0]
    
    return (bx, by, bz)

def get_Bx(x, y, z, x0, y0, z0, B0, R):
    # Ring lies in xy plane
    def integrand_x(t, x, y, z):
        r_mag = np.sqrt((x-R*np.cos(t))**2 + (y-R*np.sin(t))**2 + z**2)
        r_3 = r_mag**3
        return R*z*np.cos(t)/r_3
    
    a = 0
    b = 2*np.pi
    
    bx = B0*quad(integrand_x, a, b, args=(x-x0, y-y0, z-z0))[0]
    
    return bx

def get_By(x, y, z, x0, y0, z0, B0, R):
    # Ring lies in xy plane

    def integrand_y(t, x, y, z):
        r_mag = np.sqrt((x-R*np.cos(t))**2 + (y-R*np.sin(t))**2 + z**2)
        r_3 = r_mag**3
        return R*z*np.sin(t)/r_3
    
    
    a = 0
    b = 2*np.pi
    
    by = B0*quad(integrand_y, a, b, args=(x-x0, y-y0, z-z0))[0]
    
    return by

def get_Bz(x, y, z, x0, y0, z0, B0, R):
    # Ring lies in xy plane

    def integrand_z(t, x, y, z):
        r_mag = np.sqrt((x-R*np.cos(t))**2 + (y-R*np.sin(t))**2 + z**2)
        r_3 = r_mag**3
        return (R - y*np.sin(t) - x*np.cos(t))/r_3
    
    a = 0
    b = 2*np.pi
    
    bz = B0*quad(integrand_z, a, b, args=(x-x0, y-y0, z-z0))[0]
    
    return bz

def two_rings_x(x,y,z, x0,y0,z0, x1,y1,z1, B0, R):
    Bx_ring1 = get_Bx(x,y,z, x0,y0,z0, B0,R)
    Bx_ring2 = get_Bx(x,y,z, x1,y1,z1, B0,R)
    return Bx_ring1 + Bx_ring2

def two_rings_y(x,y,z, x0,y0,z0, x1,y1,z1, B0, R):
    By_ring1 = get_By(x,y,z, x0,y0,z0, B0,R)
    By_ring2 = get_By(x,y,z, x1,y1,z1, B0,R)
    return By_ring1 + By_ring2

def two_rings_z(x,y,z, x0,y0,z0, x1,y1,z1, B0, R):
    Bz_ring1 = get_Bz(x,y,z, x0,y0,z0, B0,R)
    Bz_ring2 = get_Bz(x,y,z, x1,y1,z1, B0,R)
    return Bz_ring1 + Bz_ring2


# Set parameters
l0 = 2*pi
B0 = 1
R = 0.25*l0

Lx           = 2.*l0
Ly           = 2.*l0
Lz           = 2.*l0

#Ring 1
x0 = Lx/2
y0 = Ly/2
z0 = 0

#Ring 2
x1 = Lx/2
y1 = Ly/2
z1 = Lz

# Resolution
resx         = 20.      #64               # nb of cells in one laser wavelength
resy         = 20.                     
resz         = 20.

# Domain
dx           = l0/resx
dy           = l0/resy
dz           = l0/resz

Nx = Lx/dx
Ny = Ly/dy
Nz = Lz/dz

print(f"Lx: {Lx}")
print(f"Nx: {int(Nx)}")

x = np.linspace(0, Lx, int(Nx))
y = np.linspace(0, Ly, int(Ny))
z = np.linspace(0, Lz, int(Nz))

print("Calculating B components.")
A = []
num = 0
x_index = 0
for i in x:
    B = []
    num = num + 1
    print(f"{num} out of {len(x)}")
    for j in y:
        C = []
        for k in z:
            C.append(two_rings_x(i,j,k, x0,y0,z0, x1,y1,z1, B0, R))
        B.append(C)
    A.append(B)
    
Bx = np.array(A)

A = []
num = 0
x_index = 0
for i in x:
    B = []
    num = num + 1
    print(f"{num} out of {len(x)}")
    for j in y:
        C = []
        for k in z:
            C.append(two_rings_y(i,j,k, x0,y0,z0, x1,y1,z1, B0, R))
        B.append(C)
    A.append(B)
    
By = np.array(A)

A = []
num = 0
x_index = 0
for i in x:
    B = []
    num = num + 1
    print(f"{num} out of {len(x)}")
    for j in y:
        C = []
        for k in z:
            C.append(two_rings_z(i,j,k, x0,y0,z0, x1,y1,z1, B0, R))
        B.append(C)
    A.append(B)
    
Bz = np.array(A)

# Bx = two_rings_x(x,y,z, x0,y0,z0, x1,y1,z1, B0, R)
# By = np.array(two_rings_y(x,y,z, x0,y0,z0, x1,y1,z1, B0, R))
# Bz = np.array(two_rings_z(x,y,z, x0,y0,z0, x1,y1,z1, B0, R))

print("Finished calculating fields.")

print(f'Bx shape: {np.shape(Bx)}')
print(f'By shape: {np.shape(By)}')
print(f'Bz shape: {np.shape(Bz)}')

with uproot.recreate("B_field_inputs.root") as f:
    f["B_field"] = {"Bx": Bx, "By": By, "Bz": Bz}

print("Fields saved to root file.")