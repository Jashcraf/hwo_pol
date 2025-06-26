import numpy as np
from prysm.coordinates import make_xy_grid, cart_to_polar
from prysm.x.polarization import pauli_spin_matrix

def init_spatial_grid(diameter, n_samples):
    """Initialize a spatial grid for polarization aberrations, returns 
    radial and azimuthal coordinates.

    Parameters
    ----------
    diameter : float
       Distance across the array 
    n_samples : int
        Number of samples along one side of the array

    Returns
    -------
    r, t : np.ndarray
        Radial and azimuthal coordinates of the spatial grid
    """
    x, y = make_xy_grid(n_samples, diameter=diameter)
    r, t = cart_to_polar(x, y)

    return r, t


def polarization_piston(diameter, n_samples):
    """Create a polarization piston `polynomial` 

    Parameters
    ----------
    diameter : float
       Distance across the array 
    n_samples : int
        Number of samples along one side of the array

    Returns
    -------
    ndarray
        N x N x 2 x 2 polarization aberration array
    """

    s0 = pauli_spin_matrix(0, shape=[n_samples, n_samples])
    s1 = pauli_spin_matrix(1, shape=[n_samples, n_samples])

    piston = s0 + 0.5 * s1

    return piston


def polarization_tilt(diameter, n_samples):
    """Create a polarization tilt `polynomial` 

    Parameters
    ----------
    diameter : float
       Distance across the array 
    n_samples : int
        Number of samples along one side of the array

    Returns
    -------
    ndarray
        N x N x 2 x 2 polarization aberration array
    """
    r, t = init_spatial_grid(diameter, n_samples)

    r = r[..., None, None]
    t = t[..., None, None]

    s0 = pauli_spin_matrix(0, shape=[n_samples, n_samples])
    s1 = pauli_spin_matrix(1, shape=[n_samples, n_samples])
    s2 = pauli_spin_matrix(2, shape=[n_samples, n_samples])

    tilt_y = -s1 * np.sin(t)
    tilt_x = s2 * np.cos(t)

    tilt = s0 + 0.5 * r * (tilt_x + tilt_y)

    return tilt


def polarization_defocus(diameter, n_samples):
    """Create a polarization defocus `polynomial` 

    Parameters
    ----------
    diameter : float
       Distance across the array 
    n_samples : int
        Number of samples along one side of the array

    Returns
    -------
    ndarray
        N x N x 2 x 2 polarization aberration array
    """
    r, t = init_spatial_grid(diameter, n_samples)

    r = r[..., None, None]
    t = t[..., None, None]

    s0 = pauli_spin_matrix(0, shape=[n_samples, n_samples])
    s1 = pauli_spin_matrix(1, shape=[n_samples, n_samples])
    s2 = pauli_spin_matrix(2, shape=[n_samples, n_samples])

    astig_0 = s1 * np.cos(2*t)
    astig_45 = s2 * np.sin(2*t)
    defocus = s0 + 0.5 * r**2 * (astig_0 + astig_45)

    return defocus
