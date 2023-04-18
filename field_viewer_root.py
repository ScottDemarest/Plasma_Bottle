import happi
import uproot
import numpy as np
import matplotlib.pyplot as plt

with uproot.open("B_field_inputs.root") as f:
    Bx_f = f["B_field"]["Bx"].array()
    By_f = f["B_field"]["By"].array()
    Bz_f = f["B_field"]["Bz"].array()

dim_half = 20

slice_y = By_f[dim_half][:][:]
shape = np.shape(slice_y)
print(f"Before saving shape: {shape}")

slice_z = Bz_f[dim_half][:][:]


@np.vectorize
def mapping_y(i, j):
    return slice_y[i][j]

@np.vectorize
def mapping_z(i, j):
    return slice_z[i][j]

Y, Z = np.mgrid[0:shape[0], 0:shape[1]] 
# U, V = slice_y[Y][Z], slice_z[Y][Z]
U = mapping_y(Y.T,Z.T)
V = mapping_z(Y.T,Z.T)
print(f"U shape: {np.shape(U)}")

skip = (slice(None, None, 8), slice(None, None, 8))
plt.quiver(Y, Z, U, V)
plt.savefig("By_Bz.png", bbox_inches="tight", dpi=300, )

# with uproot.recreate("vector_fields.root") as f:
#     f["B_field"] = {"Bx" : np.array(Bx), "By" : np.array(By), "Bz" : np.array(Bz)}

# print("Saved to root file.")

# with uproot.open("vector_fields.root") as f:
#     Bx = f["B_field"]["Bx"].array()
    
# print(np.shape(Bx))