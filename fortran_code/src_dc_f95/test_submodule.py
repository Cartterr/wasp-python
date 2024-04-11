import numpy as np
import sys

# Make sure the path is accessible for Python to find the compiled module
sys.path.append('C:/Users/josec/Desktop/WASP/fortran_code/src_dc_f95/')
import fk_openmp_module
print(dir(fk_openmp_module))


def run_fk_openmp_tests():
    # Example inputs for the sub_bs_dc subroutine.
    # Replace these with the actual inputs your subroutine expects.
    nx = 100  # Number of x points
    x = np.linspace(1, 100, nx)  # x values
    t0 = np.zeros(nx)  # Initial times, assuming zeros for simplicity
    nt = 1024  # Number of time samples, example value
    ndis = 100  # Number of distance samples, example value
    disp = True  # Boolean flag, example value
    
    # Adjust the size of the 'green' array based on the subroutine's expectations
    green = np.zeros((nt, 8, ndis), dtype=np.complex64)

    # Call the subroutine
    fk_openmp_module.fk_openmp.sub_bs_dc(nx, x, t0, green, disp, nt, ndis)

    # Assuming 'green' is modified in-place, you can now inspect or use it
    # Here, we're simply printing its shape to confirm the call was successful
    print(green.shape)

if __name__ == "__main__":
    run_fk_openmp_tests()
