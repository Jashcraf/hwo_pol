from pathlib import Path
import os
from hwo_pol.models import EAC1_Barrel_Decrease
from hwo_pol.plotting import plot_jones_3x3
from hwo_pol.conf import eac1_barrel_lens_data, barrel_lens

from poke.writing import jones_to_fits

import numpy as np
import ipdb

N_LAMS = 5
N_RAYS = 32
# center_wavelengths = [250, 550, 760, 950, 1500]
center_wavelengths = [1500]
# bandwidths = [0.1, 0.2, 0.2, 0.2]
bandwidths = [0.2]
FOV = -0.07302
JONES_DIR = Path.home() / "Box/Polarization PSDD/Jones Pupils/EAC1/"

# Set up EAC models
coatings_ota = ["XeLiF.json", "XeMgF2.json"]

# Build up Barrel Length Modesl
barrel_lens = barrel_lens

for i, (CENTER_LAMBDA, bw) in enumerate(zip(center_wavelengths, bandwidths)):

    # Special case for UV coating
    if CENTER_LAMBDA == 250:
        coating_internal = "EnhancedAg.json"
    else:
        coating_internal = "ProtectedAg.json"

    # Set up bandwidth calculation
    bw_short = (1 - bw/2) * CENTER_LAMBDA
    bw_long = (1 + bw/2) * CENTER_LAMBDA

    wavelengths = np.linspace(bw_short, bw_long, N_LAMS)

    for barrel in barrel_lens:

        for coating in coatings_ota:

            coating_name, ext = os.path.splitext(coating)
            
            # Init the EAC model
            eac = EAC1_Barrel_Decrease(
                        barrel_len=barrel,
                        num_rays=N_RAYS,
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
                jones_to_fits(rf, str(JONES_DIR / f"{coating_name}_2048/jones_pupil_{wvl}_{barrel}mm"), nmodes=37, npix=2048)

                # Write the rayfront for the APLC sampling
                jones_to_fits(rf, str(JONES_DIR / f"{coating_name}_512/jones_pupil_{wvl}_{barrel}mm"), nmodes=37, npix=512)
