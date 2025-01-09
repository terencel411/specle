import h5py
import numpy as np
import os

def compare_h5_files(file1, file2):
    with h5py.File(file1, 'r') as f1, h5py.File(file2, 'r') as f2:
        # Compare keys
        if set(f1.keys()) != set(f2.keys()):
            return False
        
        # Compare datasets
        for key in f1.keys():
            if not np.array_equal(f1[key][:], f2[key][:]):
                return False
    
    return True

print("Current working directory:", os.getcwd())

filename1 = 'specle_data/intensity.h5'
filename2 = 'specle_data/intensity1.h5'

filename1 = 'specle_data/temporal.h5'
filename2 = 'specle_data/temporal1.h5'

# Example usage
file1 = h5py.File(filename1, 'r')
file2 = h5py.File(filename2, 'r')

print("Available keys in the HDF5 file1:")
for key in file1.keys():
    print(key)

if compare_h5_files(filename1, filename2):
    print("The files are the same.")
else:
    print("The files are different.")
