import numpy as np
import matplotlib.pyplot as plt
from hwo_pol.coatings import load_coating_data

from poke.thinfilms import compute_thin_films_broadcasted

linestyles = ['--', '-.', ':', '-']
wavelengths = np.arange(300, 1700, 10)  # Wavelengths in nm
AOI = np.radians(15) # Radians
supported_coatings = ["XeLiF", "XeMgF2",
                      "EnhancedAg", "ProtectedAg"]

supported_withpath = [f"{c}.json" for c in supported_coatings]

def compute_coefficients(coating_data, aoi, wavelengths):
    
    aoi = np.full_like(wavelengths, aoi)

    rs, ts = compute_thin_films_broadcasted(coating_data,
                                            aoi,
                                            wavelengths,
                                            polarization='s')

    rp, tp = compute_thin_films_broadcasted(coating_data,
                                            aoi,
                                            wavelengths,
                                            polarization='p')
    
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
plt.figure()
for coating, ls in zip(supported_withpath, linestyles):

    data = load_coating_data(coating, wavelengths)
    reflectivity = compute_reflectivity(data, AOI, wavelengths)

    plt.plot(wavelengths, reflectivity, label=coating, linestyle=ls, linewidth=3)

plt.xlabel("Wavelength (nm)")
plt.ylabel("Reflectivity")
plt.title(f"Coating Reflectivity at {np.degrees(AOI)}° AOI")
plt.legend()
plt.grid(alpha=0.25)
plt.ylim([0.6, 1])


plt.figure()
for coating, ls in zip(supported_withpath, linestyles):

    data = load_coating_data(coating, wavelengths)
    reflectivity = compute_diattenuation(data, AOI, wavelengths)

    plt.plot(wavelengths, reflectivity, label=coating, linestyle=ls, linewidth=3)

plt.xlabel("Wavelength (nm)")
plt.ylabel("Diattenuation")
plt.title(f"Coating Diattenuation at {np.degrees(AOI)}° AOI")
plt.legend()
plt.grid(alpha=0.25)

plt.figure()
for coating, ls in zip(supported_withpath, linestyles):

    data = load_coating_data(coating, wavelengths)
    retardance = compute_retardance(data, np.degrees(AOI), wavelengths)
    retardance = np.degrees(retardance)
    plt.plot(wavelengths, retardance, label=coating, linestyle=ls, linewidth=3)

plt.xlabel("Wavelength (nm)")
plt.ylabel("Retardance, deg")
plt.title(f"Coating Retardance at {AOI}° AOI")
plt.legend()
plt.grid(alpha=0.25)

plt.show()
