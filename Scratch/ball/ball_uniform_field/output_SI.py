import happi
import math
import matplotlib.pyplot as plt
import scipy.constants

c = scipy.constants.c


laser_wavelength_SI = 10e-2
omega_r_SI = 2*math.pi * c / laser_wavelength_SI
 
S_SI = happi.Open(".", reference_angular_frequency_SI = omega_r_SI)
f = S_SI.Field.Field0

p = S_SI.Probe.Probe0

Bx = f("Bx", cmap='Spectral',  units=['cm','ns','T'], vmax=3, vmin=-3)

d = S_SI.ParticleBinning(0, cmap='inferno', units=['cm','ns','1/cm^3'])
d2 = S_SI.ParticleBinning(1, cmap='inferno', units=['cm','ns','1/cm^3'])
d3 = S_SI.ParticleBinning(2, cmap='inferno', units=['cm','ns','1/cm^3'])

p('Rho_eon', units=['cm', 'ns', 'e/cm^3']).slide(vmax = 350, vmin=-350)
p('Rho_ion', units=['cm', 'ns', 'e/cm^3']).slide(vmax = 350, vmin=-350)
