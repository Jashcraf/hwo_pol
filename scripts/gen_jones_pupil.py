import poke.plotting as plot
from poke.poke_math import np
import ipdb

from hwo_pol.writing import read_model
from hwo_pol.models import EAC1, EAC4, EAC5
from hwo_pol.zernike import zernike_coefficients_jones

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
            jones = eac.rayfronts[wvl].jones_pupil[-1][..., :2, :2]
            jones = jones.reshape([N_RAYS, N_RAYS, 2, 2])

            # Compute the Zernike coefficients
            coeffs = zernike_coefficients_jones(jones)
            
            plt.figure()
            plt.title("Imaginary Zernike Coefficients, "+r"$a + ib$") 
            
            # plot the Zernike coefficients
            plt.plot(coeffs[:, 0, 0].real, label=r"$a_{xx}$")
            plt.plot(coeffs[:, 0, 1].real, label=r"$a_{xy}$")
            plt.plot(coeffs[:, 1, 0].real, label=r"$a_{yx}$")
            plt.plot(coeffs[:, 1, 1].real, label=r"$a_{yy}$")

            plt.plot(coeffs[:, 0, 0].imag, label=r"$\phi_{xx}$")
            plt.plot(coeffs[:, 0, 1].imag, label=r"$\phi_{xy}$")
            plt.plot(coeffs[:, 1, 0].imag, label=r"$\phi_{yx}$")
            plt.plot(coeffs[:, 1, 1].imag, label=r"$\phi_{yy}$")

            plt.legend()
            plt.show()