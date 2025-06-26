"""scripts to compute eigenpolarization given a Jones pupil"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

def jones_matrix_eigenpolarizations(jones_matrices):
    """
    Compute eigenpolarizations (eigenvectors) of Jones matrices.

    Parameters
    ----------
    jones_matrices: complex ndarray of shape [..., 2, 2]
        Jones matrices at each spatial location

    Returns
    -------
    dict with keys: 'eigenvectors', 'eigenvalues'
        eigenvectors: shape [..., 2, 2] - columns are eigenvectors
        eigenvalues: shape [..., 2] - corresponding eigenvalues
    """
    # Compute eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eig(jones_matrices)

    # Sort by eigenvalue magnitude (largest first)
    # [..., ::-1] is a bit cryptic, but it flips the last axis
    # start:stop:step is the slicing syntax, so stepping by -1 reverses the order
    idx = np.argsort(np.abs(eigenvalues), axis=-1)[..., ::-1]

    # Apply sorting to eigenvalues
    eigenvalues_sorted = np.take_along_axis(eigenvalues, idx, axis=-1)

    # For eigenvectors, we need to sort the columns (last axis)
    # eigenvectors has shape (..., 2, 2) where columns are eigenvectors
    # idx has shape (..., 2)
    eigenvectors_sorted = np.zeros_like(eigenvectors)

    # Get the spatial shape
    spatial_shape = jones_matrices.shape[:-2]

    # Create indices for all spatial locations
    spatial_indices = np.unravel_index(np.arange(np.prod(spatial_shape)), spatial_shape)

    # Sort eigenvectors for each spatial location
    flat_idx = idx.reshape(-1, 2)
    flat_eigenvectors = eigenvectors.reshape(-1, 2, 2)
    flat_sorted = np.zeros_like(flat_eigenvectors)

    for i in range(flat_eigenvectors.shape[0]):
        for j in range(2):
            flat_sorted[i, :, j] = flat_eigenvectors[i, :, flat_idx[i, j]]

    eigenvectors_sorted = flat_sorted.reshape(eigenvectors.shape)

    return {
        'eigenvectors': eigenvectors_sorted,
        'eigenvalues': eigenvalues_sorted
    }

def jones_vector_to_ellipse_params(jones_vector):
    """
    Convert a Jones vector to polarization ellipse parameters.

    Parameters:
    jones_vector: complex array of shape (..., 2)
        Jones vector [Ex, Ey] where Ex and Ey are complex amplitudes

    Returns:
    dict with keys: 'semi_major', 'semi_minor', 'angle', 'chirality', 'intensity'
    """
    Ex, Ey = jones_vector[..., 0], jones_vector[..., 1]

    # Intensity components
    Ix = np.abs(Ex)**2
    Iy = np.abs(Ey)**2
    I_total = Ix + Iy

    # Avoid division by zero
    mask = I_total > 1e-10

    # Initialize output arrays
    shape = jones_vector.shape[:-1]
    semi_major = np.zeros(shape)
    semi_minor = np.zeros(shape)
    angle = np.zeros(shape)
    chirality = np.zeros(shape)

    # Only compute for non-zero intensity points
    if np.any(mask):
        Ex_m, Ey_m = Ex[mask], Ey[mask]
        Ix_m, Iy_m = Ix[mask], Iy[mask]
        I_total_m = I_total[mask]

        # Stokes parameters (normalized)
        S0 = I_total_m
        S1 = (Ix_m - Iy_m) / S0
        S2 = 2 * np.real(Ex_m * np.conj(Ey_m)) / S0
        S3 = 2 * np.imag(Ex_m * np.conj(Ey_m)) / S0

        # Degree of polarization
        DOP = np.sqrt(S1**2 + S2**2 + S3**2)
        DOLP = np.sqrt(S1**2 + S2**2)

        # Ellipse parameters
        semi_major_m = np.sqrt(S0) * np.sqrt((S0 + DOLP) / 2)
        semi_minor_m = np.sqrt(S0) * np.sqrt((S0 - DOLP) / 2)
        angle_m = 0.5 * np.arctan2(S2, S1)  # Orientation angle
        chirality_m = S3  # Positive for RCP, negative for LCP

        # Assign back to full arrays
        semi_major[mask] = semi_major_m
        semi_minor[mask] = semi_minor_m
        angle[mask] = angle_m
        chirality[mask] = chirality_m

    return {
        'semi_major': semi_major,
        'semi_minor': semi_minor,
        'angle': angle,
        'chirality': chirality,
        'intensity': I_total
    }

def remove_polarization_piston_multiplicative(jones_pupil, method='mean', aperture_mask=None):
    """
    Remove polarization piston using multiplicative correction (more physically meaningful).
    This removes the common Jones matrix by left-multiplication with its inverse.

    Parameters:
    jones_pupil: complex array of shape (N, M, 2, 2)
        Jones pupil matrices
    method: str, 'mean' or 'median'
        Method to compute the piston term
    aperture_mask: boolean array of shape (N, M), optional
        Mask defining the valid aperture

    Returns:
    jones_corrected: complex array of shape (N, M, 2, 2)
        Jones pupil with piston removed
    piston_matrix: complex array of shape (2, 2)
        The removed piston matrix
    """
    # Create aperture mask if not provided
    if aperture_mask is None:
        det = np.linalg.det(jones_pupil)
        aperture_mask = np.abs(det) > 1e-10

    valid_jones = jones_pupil[aperture_mask]

    if len(valid_jones) == 0:
        print("Warning: No valid points found in aperture")
        return jones_pupil, np.eye(2, dtype=complex)

    # Compute piston term
    if method == 'mean':
        piston_matrix = np.mean(valid_jones, axis=0)
    elif method == 'median':
        piston_matrix = (np.median(np.real(valid_jones), axis=0) +
                        1j * np.median(np.imag(valid_jones), axis=0))
    else:
        raise ValueError("Method must be 'mean' or 'median'")

    # Compute inverse of piston matrix
    try:
        piston_inv = np.linalg.inv(piston_matrix)
    except np.linalg.LinAlgError:
        print("Warning: Piston matrix is singular, using additive correction")
        return remove_polarization_piston(jones_pupil, method, aperture_mask)

    # Apply multiplicative correction: J_corrected = J_piston^-1 @ J_original
    jones_corrected = jones_pupil.copy()
    jones_corrected[aperture_mask] = np.einsum('ij,njk->nik', piston_inv, valid_jones)

    return jones_corrected, piston_matrix


def plot_eigenpolarizations(jones_pupil, x_coords=None, y_coords=None,
                           which_eigen='both', scale_factor=0.4, min_ellipse_size=0.01,
                           figsize=(12, 10), title="Eigenpolarizations Across Pupil"):
    """
    Plot eigenpolarization ellipses across a telescope pupil.

    Parameters:
    jones_pupil: complex array of shape (N, M, 2, 2)
        Jones pupil matrices at each spatial location
    x_coords, y_coords: arrays of shape (N, M), optional
        Spatial coordinates. If None, assumes normalized coordinates [-1, 1]
    which_eigen: str, 'both', 'first', 'second'
        Which eigenpolarization(s) to plot
    scale_factor: float
        Scaling factor for ellipse sizes relative to grid spacing
    min_ellipse_size: float
        Minimum ellipse size (prevents tiny ellipses from disappearing)
    figsize: tuple
        Figure size
    title: str
        Plot title

    Returns:
    fig, ax: matplotlib figure and axis objects
    """

    # Get pupil dimensions
    ny, nx = jones_pupil.shape[:2]

    # Create coordinate grids if not provided
    if x_coords is None or y_coords is None:
        x = np.linspace(-1, 1, nx)
        y = np.linspace(-1, 1, ny)
        x_coords, y_coords = np.meshgrid(x, y)

    # Compute eigenpolarizations
    eigen_result = jones_matrix_eigenpolarizations(jones_pupil)
    eigenvectors = eigen_result['eigenvectors']  # Shape: (ny, nx, 2, 2)
    eigenvalues = eigen_result['eigenvalues']    # Shape: (ny, nx, 2)

    # Create figure
    if which_eigen == 'both':
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        axes = [ax1, ax2]
        titles = ['First Eigenpolarization (Largest Eigenvalue)',
                 'Second Eigenpolarization (Smallest Eigenvalue)']
    else:
        fig, ax = plt.subplots(figsize=(figsize[0]//2, figsize[1]))
        if which_eigen == 'first':
            axes = [ax]
            titles = ['First Eigenpolarization (Largest Eigenvalue)']
        else:
            axes = [ax]
            titles = ['Second Eigenpolarization (Smallest Eigenvalue)']

    # Determine grid spacing for scaling
    dx = np.abs(x_coords[0, 1] - x_coords[0, 0]) if nx > 1 else 1
    dy = np.abs(y_coords[1, 0] - y_coords[0, 0]) if ny > 1 else 1
    grid_scale = min(dx, dy) * scale_factor

    # Plot eigenpolarizations
    eigen_indices = []
    if which_eigen == 'both':
        eigen_indices = [0, 1]
    elif which_eigen == 'first':
        eigen_indices = [0]
    else:  # 'second'
        eigen_indices = [1]

    for plot_idx, eigen_idx in enumerate(eigen_indices):
        ax = axes[plot_idx]

        # Extract eigenvectors for this eigenmode
        current_eigenvectors = eigenvectors[:, :, :, eigen_idx]  # Shape: (ny, nx, 2)
        current_eigenvalues = eigenvalues[:, :, eigen_idx]       # Shape: (ny, nx)

        # Convert eigenvectors to ellipse parameters
        ellipse_params = jones_vector_to_ellipse_params(current_eigenvectors)

        # Normalize ellipse sizes by eigenvalue magnitude
        max_eigenvalue = np.max(np.abs(current_eigenvalues))
        if max_eigenvalue > 0:
            eigenvalue_weights = np.abs(current_eigenvalues) / max_eigenvalue
        else:
            eigenvalue_weights = np.ones_like(current_eigenvalues)

        # Normalize ellipse sizes
        max_semi_major = np.max(ellipse_params['semi_major'])
        if max_semi_major > 0:
            norm_factor = grid_scale / max_semi_major
        else:
            norm_factor = 1

        # Plot ellipses
        for i in range(ny):
            for j in range(nx):
                # Skip if eigenvalue is too small
                if eigenvalue_weights[i, j] < min_ellipse_size:
                    continue

                # Skip if ellipse is too small
                if ellipse_params['semi_major'][i, j] < min_ellipse_size * max_semi_major:
                    continue

                # Ellipse parameters
                center = (x_coords[i, j], y_coords[i, j])
                width = 2 * ellipse_params['semi_major'][i, j] * norm_factor * eigenvalue_weights[i, j]
                height = 2 * ellipse_params['semi_minor'][i, j] * norm_factor * eigenvalue_weights[i, j]
                angle_deg = np.degrees(ellipse_params['angle'][i, j])

                # Color based on chirality (handedness)
                chirality = ellipse_params['chirality'][i, j]
                if chirality > 0.1:
                    color = 'red'  # Right circular polarization
                    alpha = 1 #min(0.8, abs(chirality) * eigenvalue_weights[i, j])
                elif chirality < -0.1:
                    color = 'blue'  # Left circular polarization
                    alpha = 1 #min(0.8, abs(chirality) * eigenvalue_weights[i, j])
                else:
                    color = 'black'  # Linear polarization
                    alpha = 1 #0.6 * eigenvalue_weights[i, j]

                # Create and add ellipse
                ellipse = Ellipse(center, width, height, angle=angle_deg,
                                facecolor=color, alpha=alpha, edgecolor=color,
                                linewidth=0.5)
                ax.add_patch(ellipse)

        # Set equal aspect ratio and limits
        ax.set_aspect('equal')
        ax.set_xlim(np.min(x_coords), np.max(x_coords))
        ax.set_ylim(np.min(y_coords), np.max(y_coords))

        # Labels and title
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        ax.set_title(titles[plot_idx])
        ax.grid(True, alpha=0.3)

    # Add colorbar legend to the rightmost plot
    legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor='red',
                             markersize=10, alpha=0.7, label='Right Circular'),
                      Line2D([0], [0], marker='o', color='w', markerfacecolor='blue',
                             markersize=10, alpha=0.7, label='Left Circular'),
                      Line2D([0], [0], marker='o', color='w', markerfacecolor='black',
                             markersize=10, alpha=0.6, label='Linear')]

    rightmost_ax = axes[-1]
    rightmost_ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.05, 0.5))

    plt.suptitle(title, fontsize=14)
    plt.tight_layout()
    return fig, axes
