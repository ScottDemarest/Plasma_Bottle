from math import pi
import math
import scipy
import numpy as np
import uproot

with uproot.open("B_field_inputs2.root") as f:
    Bx_f = f["B_field"]["Bx"].array()
    By_f = f["B_field"]["By"].array()
    Bz_f = f["B_field"]["Bz"].array()


# laser wavelength and period
l0           = 2*pi    # laser wavelength in normalized units
t0           = 2*pi    # optical cycle duration in normalized units
# Grid size and total simulated time
Lx           = 2.*l0
Ly           = 2.*l0
Lz           = 2.*l0
Tsim         = 1.*t0  # duration of the simulation, orignally 18

from scipy.constants import c, epsilon_0, e, m_e # constants in SI units
wavelength_SI= 30e-2 # laser wavelength, m // 30 cm
# Normalization quantities


####################  Simulated domain and time interval #######################
# Spatial and temporal resolution
resx         = 32.                     # nb of cells in one laser wavelength
resy         = 32.                     # nb of cells in one laser wavelength
resz         = 32.
rest         = 2.                     # nb of timesteps in one optical cycle. originally 500
# Mesh and integration timestep
dx           = l0/resx
dy           = l0/resy
dz           = l0/resz
dt           = t0/rest


Main(
    geometry = "3Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [dx,dy,dz],  # normalized units
    grid_length  = [Lx,Ly,Lz], # normalized units
    
    number_of_patches = [ 8,8,8],
    
    timestep = dt,          # normalized units
    simulation_time = Tsim, # normalized units
     
    EM_boundary_conditions = [
        ['PML'],
        ['PML'],
        ['PML'],
    ],
    
    random_seed = smilei_mpi_rank
)

sigma = 4e-2/wavelength_SI * l0 # 4cm / (30cm)  #0.133 * l0 #

# laser parameters
# laser_a0    = 150.        # normalized laser peak field
# laser_waist = 2.0*l0      # laser waist
c           = scipy.constants.c        # speed of light in vacuum
eps0        = scipy.constants.epsilon_0 # Vacuum permittivity, F/m
e           = scipy.constants.e         # Elementary charge, C
me          = scipy.constants.m_e       # Electron mass, kg

Nx = Lx/dx
Ny = Ly/dy
Nz = Lz/dz

xx = np.linspace(0, Lx, int(Nx))
yy = np.linspace(0, Ly, int(Ny))
zz = np.linspace(0, Lz, int(Nz))

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def Bx(x,y,z):
    i = find_nearest(xx, x)
    j = find_nearest(yy, y)
    k = find_nearest(zz, z)
    return Bx_f[i][j][k]

def By(x,y,z):
    i = find_nearest(xx, x)
    j = find_nearest(yy, y)
    k = find_nearest(zz, z)
    return By_f[i][j][k]

def Bz(x,y,z):
    i = find_nearest(xx, x)
    j = find_nearest(yy, y)
    k = find_nearest(zz, z)
    return Bz_f[i][j][k]

field_profile = {'Bx': Bx, 'By': By, 'Bz': Bz}

for field in ['Bx', 'By', 'Bz']:
    ExternalField(
        field=field,
        profile = field_profile[field]
    )

# def Bx(x,y,z):  
#     result = x*y*z
#     return result

# def By(x,y,z,t): 
#     result = two_rings_y(x,y,z, x0,y0,z0, x1,y1,z1, B0, R)
#     return result

# def Bz(x,y,z,t):  
#     result = two_rings_z(x,y,z, x0,y0,z0, x1,y1,z1, B0, R)
#     return result

# ExternalField(  
#     field = "Bx",
#     profile = Bx  # f(x,t)
# )

# plasma parameters
# n0 = 4 * N_r_SI / N_r         # initial plasma density, normalized units
# n0_ion = 4 * N_r_SI_ion / N_r         # initial plasma density, normalized units

# n0 = 4

# def ne(x,y):  # Conductor density (electrons)
#     return n0 * math.exp(-x**2 / (2 * sigma**2))

# def ne_ions(x,y):  # Conductor density (ions)
#     return n0 * math.exp(-x**2 / (2 * sigma**2))

# Species(
#     name = 'ion',
#     position_initialization = 'regular',
#     momentum_initialization = 'cold',
#     particles_per_cell = 4,
#     mass = m_i_SI / me, # 52809 = air molecule avg weight = 28.97 g/mol. e mass = 9.109e-28 g
#     # mass=28.97,
#     charge = 1.,  # normalized units
#     number_density = ne,
#     boundary_conditions = [
#         ["reflective", "remove"],
#         ["remove", "remove"],
#     ]
# )
# Species(
#     name = 'eon',
#     position_initialization = 'regular',
#     momentum_initialization = 'cold',
#     particles_per_cell = 4,
#     mass = 1.,    # normalized units
#     charge = -1., # normalized units
#     number_density = ne,
#     boundary_conditions = [
#         ["reflective", "remove"],
#         ["remove", "remove"],
#     ] 
# )


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

# DiagParticleBinning(
#     deposited_quantity = "weight",
#     every = globalEvery,
#     species = ["eon"],
#     axes = [
#         ["x", 0., Main.grid_length[0], 400],
#         ["y", 0., Main.grid_length[1], 400],
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
