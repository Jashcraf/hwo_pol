"""
Zernike analysis of Jones pupils
"""
from prysm.mathops import np
from prysm.coordinates import make_xy_grid, cart_to_polar
from prysm.polynomials import noll_to_nm, zernike_nm_sequence, lstsq

def zernike_coefficients_jones(jones_pupil, nmodes=37, mask=None):
    """Generate Zernike coefficients from a Jones pupil

    Parameters
    ----------
    jones_pupil : ndarray
        Jones pupil of shape (ny, nx, 2, 2)
    nmodes : int, optional
        Number of noll-indexed modes to compute, by default 37
    mask : ndarray, optional
        binary mask to apply to the pupil, shape (nx, ny, 2, 2), by default None

    Returns
    -------
    _type_
        _description_
    """

    assert jones_pupil.ndim == 4, "Input must be a 4D array"
    assert jones_pupil.shape[-2:] == (2, 2), "Input contain a 2x2 Jones matrix"

    # construct a basis
    x, y = make_xy_grid(jones_pupil.shape[:-2])
    r, theta = cart_to_polar(x, y)

    # Init zernike indices
    nms = [noll_to_nm(i) for i in range(1, nmodes + 1)]

    # Generate zernike basis
    basis = list(zernike_nm_sequence(nms, r, theta))

    coefficients = np.zeros([nmodes, 2, 2], dtype=np.complex128)

    # lstsq to get coefficients, loop over Jones indices
    for i in range(2):
        for j in range(2):

            jones_select = jones_pupil[..., i, j]

            if mask is not None:
                jones_select = jones_select / mask

            coefficients[:, i, j] = lstsq(basis, jones_select)

    return coefficients
            




