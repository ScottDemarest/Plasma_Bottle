"""

import happi
import matplotlib.pyplot as plt


S = happi.Open(".")
Diag = S.ParticleBinning(0, ymax=1.3, vmax=1, figsize=(8,2), cmap='inferno', xticks=x_str)
result = Diag.getData()

#Diag.plot(saveAs="plots/plot.png")
"""
"""

x_values = [0.5, 1, 1.5, 2, 2.5, 3, 3.5]
x_str = ['']*7  # Can't pass string into xticks
i = 0
for x in x_values:
    x = x + 3.5
    x_values[i] = x
    x_str[i] = str(x)
    i = i + 1
print(x_str)

S = happi.Open(".")
Diag = S.ParticleBinning(0, ymax=1.3, vmax=1, figsize=(8,2), cmap='inferno', xticks=x_str)
result = Diag.getData()

Diag.plot(saveAs="plots/plot.png")

# for i in range(0,94):
    
#     for x in x_values:
#         x = x + 3.5
    
#     plt.xticks(x_values)
#     Diag.plot(timestep = i, saveAs="plots/plot.png")

# Diag.animate(movie="movie.gif")
"""

import happi
import math
import matplotlib.pyplot as plt

# Normalized units ----------------------------

S_normalized = happi.Open(".")

#Ex = S_normalized.Field.Field0("Ex")
#Ex.animate(vmax=10, vmin=-10, cmap='seismic', movie='Test.gif', fps=3)  # use either slide or animate function
Diag = S_normalized.ParticleBinning(0, cmap='inferno')
Diag.animate(movie="dens.gif", fps=12, dpi=100)  # When animating have additional arguments "movie" and "fps"