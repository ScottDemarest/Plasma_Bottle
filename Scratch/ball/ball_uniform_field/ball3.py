from math import pi

### Choice 1: Physical Inputs in normalized units ########################

# laser wavelength and period
l0           = 2*pi    # laser wavelength in normalized units
t0           = 2*pi    # optical cycle duration in normalized units
# Grid size and total simulated time
Lx           = 6.*l0
Ly           = 6.*l0
Tsim         = 100.*t0  # duration of the simulation

# ### Choice 2: Physical Inputs in normalized units, converted from SI units #######

# from scipy.constants import c, epsilon_0, e, m_e # constants in SI units
# wavelength_SI= 0.8e-6 # laser wavelength, m
# # Normalization quantities
# L_r          = wavelength_SI       # reference length is the wavelength
# T_r          = L_r / c             # reference time
# omega_r      = 2*pi*c/L_r          # reference angular frequency
# # Laser wavelength and period in normalized units
# l0           = wavelength_SI / L_r
# t0           = wavelength_SI/c / T_r
# # Grid size and total simulated time in normalized units
# Lx           = 6. *wavelength_SI / L_r
# Ly           = 10.*wavelength_SI / L_r
# Tsim         = 10.*wavelength_SI/c / T_r   # duration of the simulation


####################  Simulated domain and time interval #######################
# Spatial and temporal resolution
resx         = 100.                     # nb of cells in one laser wavelength
rest         = 150.                     # nb of timesteps in one optical cycle 
# Mesh and integration timestep
dx           = l0/resx
dy           = l0/resx
dt           = t0/rest

Main(
    geometry = "2Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [dx,dy],  # normalized units
    grid_length  = [Lx,Ly], # normalized units
    
    number_of_patches = [ 8, 8 ],
    
    timestep = dt,          # normalized units
    simulation_time = Tsim, # normalized units
     
    EM_boundary_conditions = [
        ['PML'],
        ['PML'],
    ],
    
    random_seed = smilei_mpi_rank
)

# laser parameters
laser_a0    = 150.        # normalized laser peak field
laser_waist = 2.0*l0      # laser waist
# eps0        = scipy.constants.epsilon_0 # Vacuum permittivity, F/m
# e           = scipy.constants.e         # Elementary charge, C
# me          = scipy.constants.m_e       # Electron mass, kg
# E_r         = me*omega_r*c/e            # Reference electric field, V/m

# plasma parameters
# N_r          = eps0*omega0**2*me/e**2 # Reference density, m-3
n0             = 1         # initial plasma density, normalized units
vacuum_length  = l0          # distance between Xmin and the plasma, normalized units
plateau_length = 0.44*l0     # length of plateau of plasma density profile, normalized units

#n0 = 1

amplitude = 1

PrescribedField( # Something is wrong with this field
    field = "Bx_m",
    # profile = cosine(0.5, 1., xvacuum = 0, xphi=0., xnumber=1)
    profile = constant(100.0, xvacuum=Lx/2)
)

def ne(x,y):  # Ball of plasma
    # Radius
    R = Lx/5

    # Center
    x0 = Lx/2
    y0 = Ly/2

    if ((x-x0)**2 + (y-y0)**2) < R**2:
        return 1.
    else:
        return 0.

#E_rest = me * c**2

#Tk = 100000  # Kelvin temperature
Te = 0.1

Species(
    name = 'ion',
    position_initialization = 'random',
    momentum_initialization = 'maxwell-juettner',
    particles_per_cell = 1,
    mass = 1800, 
    charge = 1.,  # normalized units
    number_density = ne,
    temperature = [Te],
    boundary_conditions = [
        ["remove", "remove"],
        ["remove", "remove"],
    ]
)
Species(
    name = 'eon',
    position_initialization = 'random',
    momentum_initialization = 'maxwell-juettner',
    particles_per_cell = 1,
    mass = 1.,    # normalized units
    charge = -1., # normalized units
    number_density = ne,
    temperature = [Te],
    boundary_conditions = [
        ["remove", "remove"],
        ["remove", "remove"],
    ] 
)
globalEvery = 100

DiagFields(
    every = globalEvery,
    fields = ['Bx','Rho_ion','Rho_eon'] #"Env_A_abs" doesn't exist?
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 100,
    species = ["ion","eon"],
    axes = [
        ["x", 0., Main.grid_length[0], 200],
        ["y", 0., Main.grid_length[1], 200],
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 100,
    species = ["eon"],
    axes = [
        ["x", 0., Main.grid_length[0], 200],
        ["y", 0., Main.grid_length[1], 200],
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 100,
    species = ["ion"],
    axes = [
        ["x", 0., Main.grid_length[0], 200],
        ["y", 0., Main.grid_length[1], 200],
    ]
)

DiagProbe(
    every = globalEvery,
    origin = [0., Main.grid_length[1]/2.],
    corners = [
        [Main.grid_length[0], Main.grid_length[1]/2.],
    ],
    number = [int(Lx/dx)],
    fields = ['Bx','By','Bz','Rho_ion','Rho_eon']
)