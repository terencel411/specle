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

# print("Current working directory:", os.getcwd())

dir = os.getcwd() + '/proj2/specle_data/'
print("Working directory:", dir)

# filename1 = dir + 'intensity1.h5'
# filename2 = dir + 'intensity2.h5'
# filename2 = dir + 'intensity4.h5'

filename1 = dir + 'temporal1.h5'
# filename2 = dir + 'temporal2.h5'
filename2 = dir + 'temporal4.h5'

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
