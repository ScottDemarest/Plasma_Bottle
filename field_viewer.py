import happi
import uproot
import numpy as np
import matplotlib.pyplot as plt

print("Accessing Probe...")
S = happi.Open(".")

field = S.Field.Field0
Bx_f = np.array(field("Bx").getData())
By_f = np.array(field("By").getData())
Bz_f = np.array(field("Bz").getData())

dim = 40
dim_half = int(np.round(dim/2))
slice_y = By_f[0][dim_half][:][:]
shape = np.shape(slice_y)
print(f"Before saving shape: {shape}")

slice_z = Bz_f[0][dim_half][:][:]


@np.vectorize
def mapping_y(i, j):
    return slice_y[i][j]

@np.vectorize
def mapping_z(i, j):
    return slice_z[i][j]

Y, Z = np.mgrid[0:shape[0], 0:shape[1]] 
U = mapping_y(Y,Z)
V = mapping_z(Y,Z)
print(f"U shape: {np.shape(U)}")

skip = (slice(None, None, 4), slice(None, None, 4))
plt.quiver(Z[skip], Y[skip], V[skip], U[skip])
plt.savefig("By_Bz.png", bbox_inches="tight", dpi=300, )