from math import pi
import math
import scipy
import numpy as np
import uproot

# from input_parameters import *  # Does not work with Smilei for some reason.

# =================================================================
# Copied from input_parameters.py
# =================================================================

# Normalized length scale (wavelength)
l0 = 2*np.pi

Lx = 16*l0
x_min = 0

Ly = 16*l0
y_min = 0

Lz= 40*l0
z_min = 0

B0 = 1
R = 3.6*l0

#Ring 1
x0 = Lx/2
y0 = Ly/2
z0 = 4*l0

#Ring 2
x1 = Lx/2
y1 = Ly/2
z1 = 36*l0

# Resolution
resx         = 4.    # nb of cells in one laser wavelength
resy         = 4.                     
resz         = 4.

# Domain
dx           = l0/resx
dy           = l0/resy
dz           = l0/resz

Nx = Lx/dx
Ny = Ly/dy
Nz = Lz/dz

# =================================================================

# with uproot.open("B_field_inputs.root") as f:
#     Bx_f = f["B_field"]["Bx"].array()
#     By_f = f["B_field"]["By"].array()
#     Bz_f = f["B_field"]["Bz"].array()


t0           = 2*pi    # optical cycle duration in normalized units
Tsim         = 1.*t0  # duration of the simulation, orignally 18

from scipy.constants import c, epsilon_0, e, m_e # constants in SI units
wavelength_SI= 10e-2 # laser wavelength, m // 10 cm

####################  Simulated domain and time interval #######################
# Spatial and temporal resolution
rest           = 1.                     # nb of timesteps in one optical cycle. originally 500
# Mesh and integration timestep
dt           = t0/rest


Main(
    geometry = "3Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [dx,dy,dz],  # normalized units
    grid_length  = [Lx,Ly,Lz], # normalized units
    
    number_of_patches = [ 8,8,16 ],
    
    timestep = dt,          # normalized units
    simulation_time = Tsim, # normalized units
     
    EM_boundary_conditions = [
        ['PML'],
        ['PML'],
        ['PML'],
    ],
    
    random_seed = smilei_mpi_rank
)

c           = scipy.constants.c        # speed of light in vacuum
eps0        = scipy.constants.epsilon_0 # Vacuum permittivity, F/m
e           = scipy.constants.e         # Elementary charge, C
me          = scipy.constants.m_e       # Electron mass, kg
mp          = scipy.constants.m_p       # Proton mass, kg
k           = scipy.constants.Boltzmann  # Boltzmann constant

omega_r     = c/wavelength_SI           # reference frequency
B_r         = me * omega_r / e          # reference_magnetic field

# xx = np.linspace(0, Lx, int(Nx))
# yy = np.linspace(0, Ly, int(Ny))
# zz = np.linspace(0, Lz, int(Nz))

# def find_nearest(array, value):
#     array = np.asarray(array)
#     idx = (np.abs(array - value)).argmin()
#     return idx

# def Bx(x,y,z):
#     i = find_nearest(xx, x)
#     j = find_nearest(yy, y)
#     k = find_nearest(zz, z)
#     return Bx_f[i][j][k]

# def By(x,y,z):
#     i = find_nearest(xx, x)
#     j = find_nearest(yy, y)
#     k = find_nearest(zz, z)
#     return By_f[i][j][k]

# def Bz(x,y,z):
#     i = find_nearest(xx, x)
#     j = find_nearest(yy, y)
#     k = find_nearest(zz, z)
#     return Bz_f[i][j][k]

# field_profile = {'Bx': Bx, 'By': By, 'Bz': Bz}

# for field in ['Bx', 'By', 'Bz']:
#     ExternalField(
#         field=field,
#         profile = field_profile[field]
#     )

# plasma parameters
# n0 = 4 * N_r_SI / N_r         # initial plasma density, normalized units
# n0_ion = 4 * N_r_SI_ion / N_r         # initial plasma density, normalized units

n0 = 1

def ne(x,y,z):  # Ball of plasma
    # Radius
    R = Lx/3

    # Center
    x0 = Lx/2
    y0 = Ly/2
    z0 = Lz/2

    if ((x-x0)**2 + (y-y0)**2 + (z-z0)**2) < R**2:
        return 1.
    else:
        return 0.

E_rest = me * c**2

Tk = 1000  # Kelvin temperature
Te = Tk * k / E_rest  

Species(
    name = 'ion',
    position_initialization = 'random',
    momentum_initialization = 'maxwell-juettner',
    particles_per_cell = 5,
    mass = mp/me, 
    charge = 1.,  # normalized units
    number_density = ne,
    temperature = [Te],
    boundary_conditions = [
        ["remove", "remove"],
        ["remove", "remove"],
        ["remove", "remove"],
    ]
)
Species(
    name = 'eon',
    position_initialization = 'random',
    momentum_initialization = 'maxwell-juettner',
    particles_per_cell = 5,
    mass = 1.,    # normalized units
    charge = -1., # normalized units
    number_density = ne,
    temperature = [Te],
    boundary_conditions = [
        ["remove", "remove"],
        ["remove", "remove"],
        ["remove", "remove"],
    ] 
)


globalEvery = 1

DiagScalar(every=globalEvery)

# DiagProbe(
#     every = globalEvery,
#     origin = [0., Main.grid_length[1]/2., 0.],
#     corners = [
#         # [Main.grid_length[0], Main.grid_length[1]/2., 0.],
#         [Main.grid_length[0], Main.grid_length[1]/2., Main.grid_length[2]],
#         # [0., Main.grid_length[1]/2., Main.grid_length[2]],
#     ],
#     number = [int(Lx/dx)],
#     fields = ['Bx','By','Bz']
# )

DiagFields(
    every = globalEvery,
    fields = ['Bx','By','Bz'] #"Env_A_abs" doesn't exist?
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    species = ["eon","ion"],
    axes = [
        ["x", 0., Main.grid_length[0], 400],
        ["y", 0., Main.grid_length[1], 400],
        ["z", 0., Main.grid_length[2], 400],
    ]
)


# DiagParticleBinning(
#     deposited_quantity = "weight",
#     every = globalEvery,
#     species = ["eon","ion"],
#     axes = [
#         ["y", 0., Main.grid_length[1], 400],
#         ["z", 0., Main.grid_length[2], 400],
#     ]
# )

# DiagParticleBinning(
#     deposited_quantity = "weight",
#     every = globalEvery,
#     species = ["ion"],
#     axes = [
#         ["x", 0., Main.grid_length[0], 400],
#         ["y", 0., Main.grid_length[1], 400],
#     ]
# )
