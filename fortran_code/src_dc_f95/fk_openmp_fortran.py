import subprocess
import numpy as np
import time

def run_fortran_code():
    # Compile the Fortran code
    print("Compiling Fortran code...")
    subprocess.run(["make"], check=True)

    # Run the compiled Fortran executable
    print("Running Fortran executable...")
    start_time = time.time()
    subprocess.run(["./green_bank_openmp_f95"], check=True)
    end_time = time.time()
    print(f"Total calculation time: {end_time - start_time} seconds")

    # Assuming the Fortran code saves its output in a format that can be read by numpy
    # Adjust the filename as needed
    output_file = "green_function_bank.npy"
    try:
        data = np.load(output_file)
        print(f"Data shape: {data.shape}")
        return data
    except FileNotFoundError:
        print(f"Output file {output_file} not found.")
        return None

if __name__ == "__main__":
    data = run_fortran_code()
    if data is not None:
        # Optionally, process or visualize the data here
        pass
