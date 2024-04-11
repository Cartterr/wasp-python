import subprocess
import numpy as np
import matplotlib.pyplot as plt
import os
print("Current working directory: ", os.getcwd())

os.chdir('/mnt/c/Users/josec/Desktop/wasp-python/fortran_code/src_dc_f95')

# Step 1: Compile the Fortran program
subprocess.run(["make"], check=True)

# Step 2: Run the compiled binary (assuming it generates an output.txt)
# Adjust the execution command if necessary
subprocess.run(["./green_bank_f95"], check=True)

# Step 3: Load the output data
# Replace 'output.txt' with the actual output file name
# Ensure the file format matches the expected format (e.g., two columns for time and intensity)
data = np.loadtxt('output.txt')

# Assuming the first column is time and the second column is earthquake intensity
time = data[:, 0]
intensity = data[:, 1]

# Step 4: Plot the data
plt.figure(figsize=(10, 6))
plt.plot(time, intensity, label='Earthquake Intensity')
plt.xlabel('Time')
plt.ylabel('Earthquake Intensity')
plt.title('Earthquake Intensity over Time')
plt.legend()
plt.grid(True)
plt.show()
