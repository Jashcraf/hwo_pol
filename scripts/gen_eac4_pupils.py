from pathlib import Path
import os
from hwo_pol.models import EAC4_150, EAC4_165, EAC4_180
from hwo_pol.plotting import plot_jones_3x3
import numpy as np
from poke.writing import jones_to_fits

N_LAMS = 5
N_RAYS = 32
CENTER_LAMBDA = 550
FOV = -0.1
bw = 0.2
JONES_DIR = Path.home() / "Box/Polarization PSDD/Jones Pupils/EAC4/"

# Set up EAC models
models = [EAC4_150, EAC4_165, EAC4_180]
aois = [150, 165, 180]
coatings_ota = ["XeLiF.json", "XeMgF2.json"]
coating_internal = "ProtectedAg.json"
bw_short = (1 - bw/2) * CENTER_LAMBDA
bw_long = (1 + bw/2) * CENTER_LAMBDA
wavelengths = np.linspace(bw_short, bw_long, N_LAMS)


for model, aoi in zip(models, aois):
    for coating in coatings_ota:

        coating_name, ext = os.path.splitext(coating)

        # Init the EAC model
        eac = model(num_rays=N_RAYS,
                    fov=FOV,
                    wavelengths=wavelengths,
                    coating=coating,
                    coating_internal=coating_internal, # NOTE this is unused
                    ending="OTA")

        # Trace the rays
        eac.trace_rays()

        # Compute the Jones pupil
        eac.compute_jones_pupil_by_wavelength()

        # Save the Jones pupil
        for wvl in eac.wavelengths:
            
            # Grab relevant rayfront
            rf = eac.rayfronts[wvl]

            # Write the rayfront for the vortex sampling
            jones_to_fits(rf, str(JONES_DIR / f"{coating_name}_2048/jones_pupil_{wvl}_{aoi}"), nmodes=37, npix=2048)

            # Write the rayfront for the APLC sampling
            jones_to_fits(rf, str(JONES_DIR / f"{coating_name}_512/jones_pupil_{wvl}_{aoi}"), nmodes=37, npix=512)