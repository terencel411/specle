import h5py
from PIL import Image
import os

print("Current working directory:", os.getcwd())

f = h5py.File('specle/intensity.h5', 'r')

print("Available keys in the HDF5 file:")
for key in f.keys():
    print(key)

print(f)

dset = f['real'][:]

print(dset.shape)
print(dset.dtype)

# img = Image.fromarray(dset.astype("uint8"), "RGB")
img = Image.fromarray(dset, "RGB")
img.save("specle/test.png")

print("Complete")