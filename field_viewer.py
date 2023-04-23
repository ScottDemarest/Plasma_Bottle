import uproot
import happi
import numpy as np
import matplotlib.pyplot as plt
import argparse

from input_parameters import *


parser = argparse.ArgumentParser()

parser.add_argument("--bool", type=int, required=True,
                help = "To access Smilei, set  = 1; To access Root, set = 0")

args = parser.parse_args() 

if(args.bool == 1):
    print("Accessing Field in Smilei...")
    S = happi.Open(".")

    field = S.Field.Field0
    Bx_f_temp = np.array(field("Bx").getData())
    By_f_temp = np.array(field("By").getData())
    Bz_f_temp = np.array(field("Bz").getData())

    # Select the first time step.
    Bx_f = Bx_f_temp[0, :, :, :]
    By_f = By_f_temp[0, :, :, :]
    Bz_f = Bz_f_temp[0, :, :, :]
else:
    print("Accessing Root file...")
    with uproot.open("B_field_inputs.root") as f:
        Bx_f = f["B_field"]["Bx"].array()
        By_f = f["B_field"]["By"].array()
        Bz_f = f["B_field"]["Bz"].array()


print(f"Shape of each field array: {np.shape(Bx_f)}")

first_dim_half = int(np.shape(Bx_f)[0]/2)
second_dim_half = int(np.shape(Bx_f)[1]/2)
third_dim_half = int(np.shape(Bx_f)[2]/2)
# third_dim_half = 4 # at first ring

slice_y_forz = By_f[first_dim_half, :, :]
slice_z_fory = Bz_f[first_dim_half, :, :]

slice_x_forz = Bx_f[:, second_dim_half, :]
slice_z_forx = Bz_f[:, second_dim_half, :]
print(f"Before saving shape of slice_x: {np.shape(slice_x_forz)}")

slice_y_forx = By_f[:,:,third_dim_half]
slice_x_fory = Bx_f[:,:,third_dim_half]
print(f"Before saving shape y_forx: {np.shape(slice_y_forx)}")
print(f"Before saving shape x_fory: {np.shape(slice_x_fory)}")

# Plot individually without meshgrid

def plot_fields(name_x, x_field_name, x_field, name_y, y_field_name, y_field): 
    fig, ax = plt.subplots(3,3)
    fig.tight_layout()
    plt.subplot(3,1,1)
    plt.imshow(x_field, cmap='seismic', interpolation='nearest')
    plt.gca().set_aspect('equal')
    plt.title(x_field_name)
    plt.colorbar()

    plt.subplot(3,1,2)
    plt.imshow(y_field, cmap='seismic', interpolation='nearest')
    plt.gca().set_aspect('equal')
    plt.title(y_field_name)
    plt.colorbar()

    plt.subplot(3,1,3)
    for i in range(0, np.shape(x_field)[0], 6):
        for j in range(0, np.shape(y_field)[1], 6):
            x_value = x_field[i][j]
            y_value = y_field[i][j]

            r = np.power(np.add(np.power(x_value,2), np.power(y_value,2)),0.5)
            
            plt.quiver(j, i, x_value/r, y_value/r, width=.005, pivot='mid')

    plt.xlabel(name_x)
    plt.ylabel(name_y)
    plt.title(y_field_name + " vs " + x_field_name)
    plt.gca().set_aspect('equal')
    plt.savefig(y_field_name + "_" + x_field_name + "_fields.png", bbox_inches="tight", dpi=600)

plot_fields("z", "Bz", slice_z_forx, "x", "Bx", slice_x_forz)
plot_fields("z", "Bz", slice_z_fory, "y", "By", slice_y_forz)
plot_fields("x", "Bx", slice_x_fory, "y", "By", slice_y_forx)