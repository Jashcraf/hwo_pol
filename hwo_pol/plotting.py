from poke.poke_math import np
import matplotlib.pyplot as plt
import ipdb
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot_jones_3x3(jones, type='polar', vmin1=None, vmax1=None, vmin2=None, vmax2=None, title=None):
    """Plotting function including the off-diagonals to debug results
    """

    if type == 'polar':
        cm1 = 'inferno'
        cm2 = 'RdBu_r'
        op1 = np.abs
        op2 = np.angle

    elif type == 'cartesian':
        cm1 = 'viridis'
        cm2 = 'plasma'
        op1 = np.real
        op2 = np.imag

    titles = np.array(['XX', 'XY', 'XZ', 'YX', 'YY', 'YZ', 'ZX', 'ZY', 'ZZ'])
    titles = np.reshape(titles, (3, 3))
    jones = np.moveaxis(jones, 0, -1)
    jones_shape = int(np.sqrt(jones.shape[-1]))
    jones = np.reshape(jones, (3, 3, jones_shape, jones_shape))


    fig, axs = plt.subplots(ncols=6, nrows=3, figsize=(14, 6))

    if title is not None:
        plt.suptitle(title)

    for i in range(3):
        for j in range(3):
            ax = axs[i,j]
            div = make_axes_locatable(ax)
            im = ax.imshow(op1(jones[i, j]), cmap=cm1, vmin=vmin1, vmax=vmax1)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f'|{titles[i, j]}|')
            cax = div.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im, cax=cax, orientation='vertical')

    for i in range(3):
        for j in range(3):
            ax = axs[i, j+3]
            div = make_axes_locatable(ax)
            im = ax.imshow(op2(jones[i, j]), cmap=cm2, vmin=vmin2, vmax=vmax2)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f'∠{titles[i, j]}')
            cax = div.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im, cax=cax, orientation='vertical')


    plt.show()
