import numpy as np
import matplotlib.pyplot as plt
from hwo_pol.coatings import load_coating_data

from poke.thinfilms import compute_thin_films_broadcasted

linestyles = ['--', '-.', ':', '-']
wavelengths = np.arange(300, 1693, 10)  # Wavelengths in nm
AOI = 18 # Degrees
supported_coatings = ["XeLiF", "XeMgF2",
                      "EnhancedAg", "ProtectedAg"]

supported_withpath = [f"{c}.json" for c in supported_coatings]

def compute_coefficients(coating_data, aoi, wavelengths):

    aoi = np.ones_like(wavelengths) * aoi

    rs, ts = compute_thin_films_broadcasted(coating_data,
                                            aoi,
                                            wavelengths,
                                            polarization='s',
                                            substrate_index=1.5)

    rp, tp = compute_thin_films_broadcasted(coating_data,
                                            aoi,
                                            wavelengths,
                                            polarization='p',
                                            substrate_index=1.5)

    return rs, rp

def compute_reflectivity(coating_data, aoi, wavelengths):

    rs, rp = compute_coefficients(coating_data, aoi, wavelengths)

    reflectivity = np.abs(rs)**2 + np.abs(rp)**2
    reflectivity /= 2

    return reflectivity

def compute_diattenuation(coating_data, aoi, wavelengths):

    rs, rp = compute_coefficients(coating_data, aoi, wavelengths)

    reflectivity = np.abs(rs)**2 + np.abs(rp)**2
    diattenuation = (np.abs(rs)**2 - np.abs(rp)**2) / reflectivity

    return diattenuation

def compute_retardance(coating_data, aoi, wavelengths):

    rs, rp = compute_coefficients(coating_data, aoi, wavelengths)

    retardance = np.angle(rp) - np.angle(rs)

    return retardance

# Nominal reflectivity figure
plt.figure(figsize=[16, 4])
for coating, ls in zip(supported_withpath, linestyles):

    data = load_coating_data(coating, wavelengths)

    # The nominal reflectivity plot
    plt.subplot(131)
    reflectivity = compute_reflectivity(data, np.radians(AOI), wavelengths)
    plt.plot(wavelengths, reflectivity, label=coating, linestyle=ls, linewidth=3)
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Reflectivity")
    plt.title(f"Coating Reflectivity at {AOI}° AOI")
    plt.grid(alpha=0.25)
    plt.ylim([0.6, 1])

    plt.subplot(132)
    diattenuation = compute_diattenuation(data, np.radians(AOI), wavelengths)
    plt.plot(wavelengths, diattenuation, label=coating, linestyle=ls, linewidth=3)
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Diattenuation")
    plt.title(f"Coating Diattenuation at {AOI}° AOI")
    plt.ylim([0, .1])
    plt.grid(alpha=0.25)

    plt.subplot(133)
    retardance = compute_retardance(data, np.radians(AOI), wavelengths)
    retardance = np.degrees(retardance)
    plt.plot(wavelengths, retardance, label=coating, linestyle=ls, linewidth=3)
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Retardance, deg")
    plt.title(f"Coating Retardance at {AOI}° AOI")
    plt.legend()
    plt.grid(alpha=0.25)

plt.show()
