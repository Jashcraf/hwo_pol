from hwo_pol.writing import read_model
import poke.plotting as plot

# Define which
model_folder = "EAC1"

# Demonstrate re-hydrating a model from the models folder
eac1 = read_model(model_folder)

for wvl in eac1.wavelengths:
    plot.jones_pupil(eac1.rayfronts[wvl])
