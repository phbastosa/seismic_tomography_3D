import sys
import numpy as np
import matplotlib.pyplot as plt

from functions import *

filename = sys.argv[1]

x_samples = int(catch_parameter(filename, "x_samples"))
y_samples = int(catch_parameter(filename, "y_samples"))
z_samples = int(catch_parameter(filename, "z_samples"))

x_spacing = float(catch_parameter(filename, "x_spacing"))
y_spacing = float(catch_parameter(filename, "y_spacing"))
z_spacing = float(catch_parameter(filename, "z_spacing"))

vp_location = catch_parameter(filename, "vp_file")[1:-1]

vp = readBinaryVolume(z_samples, x_samples, y_samples, vp_location)

subplots = np.array([1, 1], dtype = int)
slices = np.array([z_samples/2, x_samples/2, y_samples/2], dtype = int) 
dh = np.array([x_spacing, y_spacing, z_spacing])

check_model(vp, dh, slices, subplots)
plt.show()

