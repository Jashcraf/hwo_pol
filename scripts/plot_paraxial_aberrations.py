import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from paraxial_aberrations import (
    polarization_piston,
    polarization_tilt,
    polarization_defocus
)
import ipdb

DIAMETER = 2
N_SAMPLES = 32

piston = polarization_piston(DIAMETER, N_SAMPLES)
tilt = polarization_tilt(DIAMETER, N_SAMPLES)
defocus = polarization_defocus(DIAMETER, N_SAMPLES)

piston = piston + 1j*piston
tilt = tilt + 1j*tilt
defocus = defocus + 1j*defocus

def plot_jones_pupil(array):

    fig, ax = plt.subplots(ncols=4, nrows=2, figsize=[12, 6])
    for i in range(2):
        for j in range(2):

            ax_abs = ax[i, j]
            ax_angle = ax[i, j+2]

            # code for the amplitude
            im_abs = ax_abs.imshow(np.abs(array[..., i, j]), cmap='inferno')
            ax_abs.set_title(f'|{i}{j}|')
            ax_abs.set_xticks([])
            ax_abs.set_yticks([])
            div_abs = make_axes_locatable(ax_abs)
            cax_abs = div_abs.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im_abs, cax=cax_abs, orientation='vertical')

            # code for the phase
            im_angle = ax_angle.imshow(np.angle(array[..., i, j]), cmap='twilight')
            ax_angle.set_title(f'∠{i}{j}')
            ax_angle.set_xticks([])
            ax_angle.set_yticks([])
            div_angle = make_axes_locatable(ax_angle)
            cax_angle = div_angle.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im_angle, cax=cax_angle, orientation='vertical')

    plt.show()
    return

if __name__ == '__main__':
    plot_jones_pupil(piston)
    plot_jones_pupil(tilt)
    plot_jones_pupil(defocus)
