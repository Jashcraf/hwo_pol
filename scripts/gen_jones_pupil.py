from hwo_pol.writing import read_model
import poke.plotting as plot
from poke.poke_math import np

models = ["EAC1", "EAC4", "EAC5"]
wavelengths = [250, 550, 760, 950, 1500]
bandwidths = [0.1, 0.2, 0.2, 0.2, 0.2]
N_LAMS = 5 # spectral sampling

for model_folder in models:

    eac = read_model(model_folder)

    for wvl, bw in zip(wavelengths, bandwidths):

        wavelengths = np.linspace(wvl * bw / 2, wvl * bw / 2, N_LAMS)

        for w in wavelengths:

            # re-generate surfaces
