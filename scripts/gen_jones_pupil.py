import matplotlib.pyplot as plt
import ipdb

import poke.plotting as plot
from poke.poke_math import np

from hwo_pol.writing import read_model
from hwo_pol.models import EAC1, EAC4, EAC5
from hwo_pol.zernike import zernike_coefficients_jones

from prysm.coordinates import cart_to_polar

N_LAMS = 5 # spectral sampling
N_RAYS = 32
ending = "TCA" # OTA or TCA

# Set up EAC models
models = [EAC1, EAC4, EAC5]

# Set up coatings to load
coating_ota  = "XeLiF.json"

coatings_internal = [
    "EnhancedAg.json",
    "ProtectedAg.json",
    "ProtectedAg.json",
    "ProtectedAg.json",
    "ProtectedAg.json"
]

# Set up spectral data
center_wavelengths = [350, 550, 760, 950, 1500]
bandwidths = [0.1, 0.2, 0.2, 0.2, 0.2]

# Set up Fields of View
fov_1 = -0.07302
fov_4 = -0.025
fov_5 = -0.025
fovs = [fov_1, fov_4, fov_5]

# Init the EAC objects
for EAC, fov in zip(models, fovs):
    for wvl, bw, internal in zip(center_wavelengths, bandwidths, coatings_internal):
        
        # Set up bandpass
        short_wvl = wvl * (1 - bw/2)
        long_wvl = wvl * (1 + bw/2)
        wavelengths = np.linspace(short_wvl, long_wvl, N_LAMS)

        # init the EAC model
        eac = EAC(N_RAYS, fov, wavelengths, coating_ota, internal, ending=ending)
        
        # Trace the rays
        eac.trace_rays()

        # Compute the Jones pupil
        eac.compute_jones_pupil_by_wavelength()

        # Compute the Zernike coefficients
        for wvl in eac.wavelengths:
            
            rf = eac.rayfronts[wvl]
            jones = rf.jones_pupil[-1][..., :2, :2]
            
            # make a circular mask
            x = rf.xData[0,0] / rf.pupil_radius
            y = rf.yData[0, 0] / rf.pupil_radius
            x = x.reshape([N_RAYS, N_RAYS])
            y = y.reshape([N_RAYS, N_RAYS])

            r, t = cart_to_polar(x, y)
            mask = np.zeros_like(r)
            mask[r <= 1] = 1
            jones = jones.reshape([N_RAYS, N_RAYS, 2, 2])

            plot.jones_pupil(eac.rayfronts[wvl], coordinates="cartesian")

            # Compute the Zernike coefficients
            coeffs = zernike_coefficients_jones(jones, mask=mask)
            
            plt.figure()
            plt.title("Complex Zernike Coefficients, "+r"$a + ib$") 

            print(coeffs)
            
            # plot the Zernike coefficients
            plt.plot(coeffs[:, 0, 0].real, label=r"$a_{xx}$", linestyle="None", marker="o")
            plt.plot(coeffs[:, 0, 1].real, label=r"$a_{xy}$", linestyle="None", marker="o")
            plt.plot(coeffs[:, 1, 0].real, label=r"$a_{yx}$", linestyle="None", marker="o")
            plt.plot(coeffs[:, 1, 1].real, label=r"$a_{yy}$", linestyle="None", marker="o")

            plt.plot(coeffs[:, 0, 0].imag, label=r"$\phi_{xx}$", linestyle="None", marker="x")
            plt.plot(coeffs[:, 0, 1].imag, label=r"$\phi_{xy}$", linestyle="None", marker="x")
            plt.plot(coeffs[:, 1, 0].imag, label=r"$\phi_{yx}$", linestyle="None", marker="x")
            plt.plot(coeffs[:, 1, 1].imag, label=r"$\phi_{yy}$", linestyle="None", marker="x")

            plt.legend()
            plt.show()