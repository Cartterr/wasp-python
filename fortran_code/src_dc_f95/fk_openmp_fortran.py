import numpy as np
import matplotlib.pyplot as plt
import subprocess

def compile_and_run_fortran(vel_model_path):
    # Compile the Fortran program
    compile_command = (
        "gfortran -fopenmp -o green_bank_fk_test "
        "fk_openmp.f95 green_bank_fk_test.f95 "
        "constants.f95 retrieve_gf.f95 bessel2.f95 vel_model_data.f95 "
        "wave_travel.f95 layer.f95 fk_source.f95 haskell.f95 prop.f95 "
        "fk_kernel.f95 bessel.f fft.c Complex.c "
        "-L/usr/lib/llvm-14/lib -lgfortran -lomp"
    )
    subprocess.run(compile_command, shell=True, check=True)
    
    # Run the Fortran program, passing the velocity model file path as an argument
    run_command = f"./green_bank_fk_test '{vel_model_path}'"
    print("Running command:", run_command)
    result = subprocess.run(run_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print("Output:", result.stdout.decode())
    print("Errors:", result.stderr.decode())
    print("Fortran program has generated the mock Green's function bank.")




def load_green_function_bank(filename):
    dtype = np.float32
    data = np.fromfile(filename, dtype=dtype)
    print(f"Data size before reshape: {data.size}")
    print(f"Total data elements: {data.size}")
    # Temporary reshape based on actual data size to avoid errors
    # Adjust dimensions based on actual data structure
    nx, nz, nt, nComponents = 20, 5, 1024, 8  # Example adjustment based on actual data size
    gf_bank = data.reshape((nz, nx, nt, nComponents), order='F')
    return gf_bank


def visualize_green_function(gf_bank, depth_index=0, distance_index=0, component_index=0):
    # Visualize a specific component of the Green's function
    green_function = gf_bank[depth_index, distance_index, :, component_index]
    plt.figure(figsize=(10, 6))
    plt.plot(np.real(green_function))
    plt.xlabel('Time Sample')
    plt.ylabel('Amplitude')
    plt.title(f"Green's Function Component {component_index+1}")
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    vel_model_path = '/mnt/e/wasp-python/binaries_test/vel_model.txt'
    compile_and_run_fortran(vel_model_path)
    filename = 'green_function_bank.bin'
    gf_bank = load_green_function_bank(filename)
    # Visualize the first component of the Green's function for the first depth and distance
    visualize_green_function(gf_bank, depth_index=0, distance_index=0, component_index=0)
