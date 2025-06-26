from hwo_pol.models import EAC1
import poke.plotting as plot
from hwo_pol.writing import write_model

# Init the EAC object
num_rays = 32
fov_1 = -0.07302
wavelengths = [550]
coating = "XeLiF.json"
coating_internal = "ProtectedAg.json"
ending = "OTA"

eac1 = EAC1(num_rays, fov_1, wavelengths, coating, coating_internal,
            ending=ending)

# First trace the rays with Poke
eac.trace_rays()

# Compute the Jones pupil with Poke
eac.compute_jones_pupil_by_wavelength()

# Plot the pupil
for wavelength in wavelengths:
    
    # Show the Jones Pupil
    plot.jones_pupil(eac.rayfronts[wavelength])
   
    # Write EAC model to models/ directory
    write_model(eac, "EAC1")

