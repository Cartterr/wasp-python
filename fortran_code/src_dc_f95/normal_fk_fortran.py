import subprocess
import numpy as np
import matplotlib.pyplot as plt
import os

# Initial debug statement to check the current working directory
print("Initial working directory: ", os.getcwd())

# Change directory to where the Fortran code is
os.chdir('/mnt/e/wasp-python/fortran_code/src_dc_f95')
print("Changed to Fortran source directory: ", os.getcwd())

# Step 1: Compile the Fortran program
print("Starting compilation of Fortran program...")
compilation_result = subprocess.run(["make"], check=True)
print("Compilation completed successfully: ", compilation_result)

# Step 2: Run the compiled binary
print("Running the compiled binary...")
execution_result = subprocess.run(["./green_bank_f95"], check=True)
print("Execution completed successfully: ", execution_result)

# Define nt and other parameters correctly, assuming nt is known
nt = 20  # This should match the number of timesteps used in the Fortran code
n_components = 8  # Assuming there are 8 intensity components
print(f"Defined parameters - nt: {nt}, n_components: {n_components}")

# Step 3: Load the output data
output_file_name = 'green_function_bank.bin'
print("Loading data from: ", output_file_name)
dtype = np.dtype([
    ('depth', np.float64),
    ('distance', np.float64),
    ('time', np.float64),
    ('intensity', np.float32, (nt, n_components))
])
data = np.fromfile(output_file_name, dtype=dtype)
print("Data loaded successfully")

# Check the range and statistics of intensity values
print("Intensity statistics:")
print("Min:", np.min(data['intensity']))
print("Max:", np.max(data['intensity']))
print("Mean:", np.mean(data['intensity']))
print("Standard Deviation:", np.std(data['intensity']))

if np.all(data['intensity'] == 0):
    print("All intensity values are zero. Check the Fortran output or calculation method.")


# Printing the shape of the data to understand its structure
# After loading the data
print("Data shape:", data.shape)
print("Sample intensity data:", data['intensity'][0, :, :])  # Display the intensity matrix of the first record

# Adding checks for data consistency
if data.size == 0:
    print("Warning: No data loaded, check binary file and dtype definition.")

print("Intensity data shape example:", data['intensity'][0].shape)

# Access and plot specific fields
time = data['time']
intensity = data['intensity'][:, 0, 0]  # Selecting the first timestep of the first component
print("Extracted time and intensity data for plotting.")

# Step 4: Plot the data
print("Time data (first 10 entries):", time[:10])

# Selecting a different component if the first doesn't show expected results
intensity = data['intensity'][:, 0, 1]  # Change the second index if needed to another component
print("Intensity data sample (first 10 entries):", intensity[:10])

# Normalize time if it's not in expected range, assume time should be in seconds
max_time = np.max(time)
if max_time > 1e9:  # adjust threshold based on your data scale
    time = time / max_time  # normalize to 1
    print("Normalized time for better scaling.")

# Print loaded data dimensions and sample content
print("Data loaded successfully with shape:", data.shape)
print("Sample data from the loaded file:", data['intensity'][0])

# Verify the min and max values in the loaded intensity data
print("Minimum intensity value:", np.min(data['intensity']))
print("Maximum intensity value:", np.max(data['intensity']))

# Check and print specific component data
for i in range(n_components):
    print(f"Component {i+1} intensity values:", data['intensity'][:, :, i])

# Before plotting, debug print to verify the time and intensity arrays
print("Time data for plotting:", time)
print("Intensity data for first component:", intensity)


plt.figure(figsize=(12, 8))
for i in range(n_components):
    plt.plot(time, data['intensity'][:, 0, i], label=f'Component {i+1}')
plt.xlabel('Time (seconds)')
plt.ylabel('Intensity')
plt.title('Earthquake Intensity Profile Across All Components')
plt.legend()
plt.grid(True)

# Adjust the x-axis limits if needed
plt.xlim(0, np.max(time))  # Ensure that x-axis starts at zero and ends at max time

plt.show()
