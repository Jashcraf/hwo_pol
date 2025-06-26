"""These data should have near-zero z-components,
with exception of the ZZ term, which should be 1

You may see some "noise", this is machine epsilon
for single-precision floats arising in the computation
of a zero
"""

from hwo_pol.writing import read_model
from hwo_pol.plotting import plot_jones_3x3

# Define which
EACs = [1, 4, 5]

for eac in EACs:
    model_folder = f"EAC{eac}"
    print(f"Reading model from {model_folder}")

    # Demonstrate re-hydrating a model from the models folder
    eac_model = read_model(model_folder)

    for wvl in eac_model.wavelengths:
        plot_jones_3x3(eac_model.rayfronts[wvl].jones_pupil[-1],
            title=f"EAC{eac} Wavelength {wvl} nm")
