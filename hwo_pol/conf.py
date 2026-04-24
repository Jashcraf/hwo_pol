"""Defining constants for EAC1 models barrel length study
"""

from poke.poke_math import np
from pathlib import Path

barrel_lens = np.arange(0, 10_000, 1000)
DESIGN_PATH = Path.home() / "Box/Polarization PSDD/Optical Designs/EAC1/"
SURFLIST = [5, 6, 8, 9, 10, 11, 13, 14, 16, 19, 20]
IS_OTA = [True, True, False, False, False, False, False, False, False, False, False]

# Evaluated at FSM plane, surface
ALOCS_BARREL_DECREASE = [
    np.array([0.07893,    -0.06487,    -0.99477]),
    np.array([0.07893,    -0.06198,    -0.99495]),
    np.array([0.07893,    -0.05724,    -0.99523]),
    np.array([0.07893,    -0.05163,    -0.99554]),
    np.array([0.07893,    -0.04599,    -0.99582]),
    np.array([0.07893,    -0.02967,    -0.99644]),
    np.array([0.07893,    -0.01255,    -0.99680]),
    np.array([0.07893,     0.00738,    -0.99685]),
    np.array([0.07893,     0.03261,    -0.99635]),
    np.array([0.07893,     0.07321,    -0.99419])
]

ALOCPS_BARREL_DECREASE = [
    np.array([0.07894,    -0.09788,    -0.99206]),
    np.array([0.07894,    -0.09483,    -0.99236]),
    np.array([0.07894,    -0.09025,    -0.99279]),
    np.array([0.07894,    -0.08480,    -0.99327]),
    np.array([0.07893,    -0.07931,    -0.99372]),
    np.array([0.07893,    -0.06287,    -0.99490]),
    np.array([0.07893,    -0.04560,    -0.99584]),
    np.array([0.07893,    -0.02544,    -0.99656]),
    np.array([0.07893,     0.00019,    -0.99688]),
    np.array([0.07893,     0.04162,    -0.99601])
]

# Construct dictionary to store the above data
eac1_barrel_lens_data = {}

for i, barrel in enumerate(barrel_lens):

    # Build the path
    path = DESIGN_PATH / f"hwo_ota_with_coronagraph_m{barrel}mm.len"

    eac1_barrel_lens_data[barrel] = {
        "aloc": ALOCS_BARREL_DECREASE[i],
        "alocp": ALOCPS_BARREL_DECREASE[i],
        "path": str(path)
    }