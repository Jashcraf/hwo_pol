from hwo_pol.models import EAC1, EAC4, EAC5
import poke.plotting as plot
from hwo_pol.writing import write_model

# Init the EAC object
num_rays = 32
fov_1 = -0.07302
fov_4 = -0.025
fov_5 = -0.025
center_wavelengths = [250, 550, 760, 950, 1500]
bandwidths = [0.1, 0.2, 0.2, 0.2, 0.2]
wavelengths = []

for wvl, bw in zip(center_wavelengths, bandwidths):
    short_wvl = wvl * (1 - bw/2)
    long_wvl = wvl * (1 + bw/2)
    bandpass = np.linspace(short_wvl, long_wvl, 32)
    for lam in bandpass:
        wavelengths.append(bandpass)

coating = "XeLiF.json"
coating_internal = "ProtectedAg.json"
ending = "OTA"

eac1 = EAC1(num_rays, fov_1, wavelengths, coating, coating_internal,
            ending=ending)
eac4 = EAC4(num_rays, fov_4, wavelengths, coating, coating_internal,
            ending=ending)
eac5 = EAC5(num_rays, fov_5, wavelengths, coating, coating_internal,
            ending=ending)

eacs = [eac1, eac4, eac5]

for eac, num in zip(eacs, [1, 4, 5]):

    # First trace the rays, performs the trace for each wavelength
    eac.trace_rays()

    # Compute the Jones pupil
    eac.compute_jones_pupil_by_wavelength()

    # # Plot the pupil
    for wavelength in wavelengths:
        plot.jones_pupil(eac.rayfronts[wavelength])

    # Write the rayfront to a serial file
    # write_model(eac, f"EAC{num}")
