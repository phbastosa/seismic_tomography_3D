import numpy as np

n = 256

model = np.ones((n,n,n)) * 2000

model.astype("float32", order = "F").tofile("input_model.bin")