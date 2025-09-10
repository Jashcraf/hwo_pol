import numpy as np
import matplotlib.pyplot as plt
from hwo_pol.coatings import load_coating_data

from poke.thinfilms import compute_thin_films_broadcasted

linestyles = ['--', '-.', ':', '-']
wvl_detect = np.arange(350, 1050, 50) # nm
sc_lambda = [0.727, 0.746, 0.765, 0.783, 0.802, 0.821, 0.839, 0.859, 0.878, 0.896, 0.915, 0.934, 0.953, 0.971, 1.000]
wvl_spec = [l * 1e3 for l in sc_lambda]
wvl_spec = np.asarray(wvl_spec)

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
N_OPTICS = 14
plt.figure()
plt.title("Throughput Detection")
for coating, ls in zip(supported_withpath, linestyles):

    data = load_coating_data(coating, wvl_detect)
    reflectivity = compute_reflectivity(data, AOI, wvl_detect)

    plt.plot(wvl_detect, reflectivity**N_OPTICS, label=coating, linestyle=ls, linewidth=3, marker="o")

plt.xlabel("Wavelength (nm)")
plt.ylabel("Reflectivity")
plt.title(f"Coating Reflectivity at {np.degrees(AOI)}° AOI")
plt.legend()
plt.grid(alpha=0.25)
plt.ylim([0.6, 1])

print("-"*20)
print(f"Toptical = [{reflectivity**N_OPTICS}]")
print("-"*20)

plt.figure()
plt.title("Throughput Characterization")
for coating, ls in zip(supported_withpath, linestyles):

    data = load_coating_data(coating, wvl_spec)
    reflectivity = compute_reflectivity(data, AOI, wvl_spec)

    plt.plot(wvl_spec, reflectivity**N_OPTICS, label=coating, linestyle=ls, linewidth=3, marker="o")

plt.xlabel("Wavelength (nm)")
plt.ylabel("Reflectivity")
plt.title(f"Coating Reflectivity at {np.degrees(AOI)}° AOI")
plt.legend()
plt.grid(alpha=0.25)
plt.ylim([0.6, 1])

print("-"*20)
print(f"sc_Toptical = [{reflectivity**N_OPTICS}]")
print("-"*20)

plt.show()
