import numpy as np
import matplotlib.pyplot as plt
import ipdb
from poke.writing import read_serial_to_rayfront
from eigenpolarization import (
    jones_matrix_eigenpolarizations,
    jones_vector_to_ellipse_params,
    remove_polarization_piston_multiplicative,
    plot_eigenpolarizations
)
from paraxial_aberrations import (
    polarization_defocus,
    polarization_tilt,
    polarization_piston
)
from prysm.x.polarization import quarter_wave_plate, linear_polarizer

# Example usage and test function
def create_test_jones_pupil(nx=20, ny=20):
    """Create a test Jones pupil with interesting polarization patterns."""
    x = np.linspace(-1, 1, nx)
    y = np.linspace(-1, 1, ny)
    X, Y = np.meshgrid(x, y)
    R = np.sqrt(X**2 + Y**2)
    theta = np.arctan2(Y, X)

    # Create aperture mask (circular)
    aperture = R <= 1.0

    # Initialize Jones pupil matrices
    jones_pupil = np.zeros((ny, nx, 2, 2), dtype=complex)

    # Create a test Jones matrix at each point
    # This creates a varying birefringent element
    for i in range(ny):
        for j in range(nx):
            if aperture[i, j]:
                # Create a rotation + retardation matrix
                angle = theta[i, j]  # Rotation angle varies with position
                retardation = np.pi * R[i, j]  # Retardation varies with radius

                # Rotation matrix
                c, s = np.cos(angle), np.sin(angle)
                R_matrix = np.array([[c, -s], [s, c]])

                # Retardation matrix (linear retarder along x)
                retard_matrix = np.array([[np.exp(-1j*retardation/2), 0],
                                        [0, np.exp(1j*retardation/2)]])

                # Combined Jones matrix
                jones_pupil[i, j] = R_matrix @ retard_matrix @ R_matrix.T

    return jones_pupil, X, Y

# Example usage
if __name__ == "__main__":

    wvls = np.arange(500, 625, 5)
    markers = ['o', 'v', '^', '<', '>', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X']
    markers = markers + markers

    # Create test data
    # print("Creating test Jones pupil...")
    jones_pupil, X, Y = create_test_jones_pupil(10, 10)
    jones_pupil = polarization_defocus(diameter=2, n_samples=13)
    jones_tilt = polarization_tilt(diameter=2, n_samples=13)
    qwp = quarter_wave_plate(theta=0, shape=jones_pupil.shape[:-2])
    lpol = linear_polarizer(theta=np.pi/4, shape=jones_pupil.shape[:-2]) 
    jones_pupil = jones_pupil @ qwp
    # Load in Jones Pupils
    for wvl in wvls[:3]:
        rf = read_serial_to_rayfront(f"/home/jarenashcraft/Data/EAC1_Retardance_Dispersion/hwo_to_coro_ep_{wvl}nm_XeLiF.msgpack")
        X, Y = rf.xData[0,0], rf.yData[0, 0]
        X = X.reshape([5, 5])
        Y = Y.reshape([5, 5])
        jones_pupil = rf.jones_pupil[-1][..., :2, :2]
        jones_pupil = jones_pupil.reshape([10, 10, 2, 2])
    jones_pupil, mean_matrix = remove_polarization_piston_multiplicative(jones_pupil)

        # Plot both eigenpolarizations
    fig, axes = plot_eigenpolarizations(jones_pupil, X, Y, which_eigen='both',
                                        title="Eigenpolarizations of Test Jones Pupil")
    plt.show()

    # Plot just the first eigenpolarization
    # fig2, ax2 = plot_eigenpolarizations(jones_pupil, X, Y, which_eigen='first',
    #                                   title="First Eigenpolarization Only")
    # plt.show()

    # Create a simple test case - quarter wave plate
    # print("Creating quarter wave plate test...")
    # simple_jones = np.zeros((10, 10, 2, 2), dtype=complex)
    # # Quarter wave plate at 45 degrees
    # qwp = np.array([[1, -1j], [-1j, 1]]) / np.sqrt(2)
    # simple_jones[:, :] = qwp

    # fig3, ax3 = plot_eigenpolarizations(simple_jones, which_eigen='both',
    #                                   title="Quarter Wave Plate Eigenpolarizations")
    # plt.show()
