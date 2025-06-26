from hwo_pol.models import EAC1, EAC4, EAC5
import poke.plotting as plot

# Init the EAC object
num_rays = 10
fov = -0.07302
wavelengths = [550]
coating = "XeLiF.json"
coating_internal = "ProtectedAg.json"
ending = "OTA"

eac1 = EAC1(num_rays, fov, wavelengths, coating, coating_internal, ending=ending)
eac4 = EAC4(num_rays, fov, wavelengths, coating, coating_internal, ending=ending)
eac5 = EAC5(num_rays, fov, wavelengths, coating, coating_internal, ending=ending)

eacs = [eac1, eac4, eac5]

for eac in eacs:
    
    # First trace the rays
    eac.trace_rays()

    # Compute the Jones pupil
    eac.compute_jones_pupil_by_wavelength()

    # # Plot the pupil
    for wavelength in wavelengths:
        plot.jones_pupil(eac.rayfronts[wavelength])