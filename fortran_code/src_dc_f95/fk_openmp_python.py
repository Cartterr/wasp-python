import numpy as np
import sys
import types
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def check_fortran_functions_in_module(module, func_names):
    for attr_name in dir(module):
        attr = getattr(module, attr_name)
        # Check if the attribute is a Fortran module by checking its type as a simple heuristic
        if type(attr).__name__ == 'fortran':
            print(f"Checking functions in submodule: {attr_name}")
            for func_name in func_names:
                try:
                    func = getattr(attr, func_name)
                    print(f" - {func_name} exists.")
                    pass
                except AttributeError:
                    # print(f" - {func_name} does not exist.")
                    pass
            print("")  # Add a newline for readability

# Assuming fk_openmp_module is compiled and accessible, with necessary submodules
sys.path.append('C:/Users/josec/Desktop/WASP/fortran_code/src_dc_f95/')
import fk_openmp_module

# List of function names you expect to find in the submodules
expected_functions = ["sub_bs_dc", "update_model_for_depth", "load_velocity_model"]

# Example usage
check_fortran_functions_in_module(fk_openmp_module, expected_functions)

def generate_green_function(nx, distances, t0, nt, ndis, disp):
    """
    Generate Green's function using Fortran subroutine.

    Parameters:
    - nx: Number of distance samples.
    - distances: Distances array.
    - t0: Initial times array.
    - nt: Number of time samples.
    - ndis: Number of distance samples.
    - disp: Whether to calculate displacements.
    """
    # Initialize the Green tensor with zeros for simplicity
    green = np.zeros((nt, 8, ndis), dtype=np.complex128, order='F')

    # Generate the Green's function using the Fortran subroutine
    fk_openmp_module.fk_openmp.sub_bs_dc(nx, distances, t0, green, disp, nt, ndis)

    return green

def create_mock_green_function_bank(dist_min, dist_max, d_step, dep_min, dep_max, dep_step, gf_bank_file):
    """
    Creates a mock Green's function bank with simplified data and saves it to a file.
    """
    nx = int((dist_max - dist_min) / d_step) + 1  # Number of distance samples
    nz = int((dep_max - dep_min) / dep_step) + 1  # Number of depth samples
    nt = 1024  # Number of time samples, for simplicity

    # Creating mock Green's function data
    green_function_bank = np.random.rand(nz, nx, nt, 8).astype(np.float32)

    # Save the mock Green's function bank to a file
    np.save(gf_bank_file, green_function_bank)

def visualize_green_function(gf_bank_file, depth_index=0, distance_index=0, component_index=0):
    """
    Visualizes a single component of the Green's function over time for a specific depth and distance.
    """
    # Load the Green's function bank from file
    green_function_bank = np.load(gf_bank_file)

    # Extract the specified component of the Green's function
    green_function = green_function_bank[depth_index, distance_index, :, component_index]

    # Plot the real part of the Green's function over time
    plt.figure(figsize=(10, 6))
    plt.plot(np.real(green_function))
    plt.xlabel('Time Sample')
    plt.ylabel('Amplitude')
    plt.title(f'Green\'s Function Component {component_index+1} for Depth Index {depth_index} and Distance Index {distance_index}')
    plt.grid(True)
    plt.show()

# Example usage
if __name__ == "__main__":
    dist_min = 0
    dist_max = 100
    d_step = 10
    dep_min = 0
    dep_max = 50
    dep_step = 10
    gf_bank_file = "green_function_bank.npy"
    
    # Create and visualize the mock Green's function bank
    create_mock_green_function_bank(dist_min, dist_max, d_step, dep_min, dep_max, dep_step, gf_bank_file)
    print(f"Mock Green's function bank saved to {gf_bank_file}.")
    visualize_green_function(gf_bank_file, depth_index=0, distance_index=0, component_index=0)