import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from pathlib import Path
import ipdb
from pathlib import Path
import json
import warnings

from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

# The Poke things
from poke.poke_core import Rayfront
import poke.plotting as plot
from poke.writing import write_rayfront_to_serial, jones_to_fits
from poke.materials import create_index_model # poke baked in models
from poke.thinfilms import compute_thin_films_broadcasted

# NOTE this returns the project directory, not the coating folder
PROJECT_PATH = Path(__file__).parent.parent 
COATINGS_PATH = PROJECT_PATH / "coating_recipes"
n_BK7 =  1.4801 # Approximation for ULE, though it really shouldn't matter

# NOTE each of these is a material, not a coating stack. The materials denoted
# 'Layer' are protected / not mine to distribute, so their names are not
# explicitly referenced
supported_coatings = ['XeLiF','eLiF','MgF2','F3','Al2O3','Cr','SiO2',
                      'Aluminum',
                      'Layer1', 'Silver', 'Layer3_5', 'Layer4_6', 'Layer7']

# NOTE These are the names of the coating recipes, not the materials
supported_recipes = ['XeLiF', 'XeMgF2', 'EnhancedAg', 'ProtectedAg']


def make_index_model(coating):
    """Generate a refractive index model interpolated from real data. Relies
    on the coating formulae being available, which needs to be populated before
    installation. 

    Parameters
    ----------
    coating : str
        string containing the coating to load, should be in supported_coatings

    Returns
    -------
    callable
        callable function of nanometers that returns the complex index
    """

    # First check that there are files in the coating_recipes folder
    assert COATINGS_PATH.exists(), f"Coatings folder {COATINGS_PATH} not found"
    assert coating in supported_coatings, f"Coating {coating} not supported" 

    # Skips first two rows because it returns a nan. Suspect this has to do
    # with the csv file saving via excel. 
    if coating in ['MgF2', 'Al2O3', 'Cr', 'F3', 'SiO2']:
        coating_data = np.genfromtxt(COATINGS_PATH / f'{coating}.csv',
                                     delimiter=',', skip_header=2)
    else:    
        coating_data = np.genfromtxt(COATINGS_PATH / f'{coating}.csv',
                                     delimiter=',')

    wavelengths = coating_data[:,0]
    n = coating_data[:,1]
    k = coating_data[:,2]

    n_interp = interp1d(wavelengths, n, kind='cubic')
    k_interp = interp1d(wavelengths, k, kind='cubic')

    index_model = lambda x: n_interp(x) + 1j*k_interp(x)

    return index_model


def load_coating_data(name, wavelengths):
    """Load a coating recipe from a json file

    Parameters
    ----------
    name : str or Path
        Name of the coating recipe in the coating_recipes folder
    wavelengths : float or array of floats
        wavelength in nanometers to evaluate the complex refractive index
        at

    Returns
    -------
    list of tuples
        list of tuples of refractive index and thickness, compatible with the
        format Poke likes coatings in. i.e.

        [
            (n1, t1),
            (n2, t2),
            ...
            (nN)
        ]

        The last entry is the substrate, which has no thickness because the
        thin film algorithm assumes it's semi-infinite
    """


    if not isinstance(name, Path):
        name = Path(name)

    assert name.stem in supported_recipes, f"Coating {name} not supported"
    assert name.suffix == ".json", f"Coating extension {name.suffix} not supported"

    with open(name, "r") as f:
        recipe = json.load(f)

    coating_stack = []
    for layer in recipe["Recipe"]:
        index_model = make_index_model(layer["material"])
        eval_index = index_model(wavelengths)
        thickness = np.full_like(wavelengths, layer["thickness_nm"])
        coating_stack.append((eval_index, thickness))

    # append the mirror substrate
    glass = np.full_like(wavelengths, n_BK7)
    coating_stack.append((glass))

    return coating_stack


def make_enhanced_silver(wvl):
    """Helper function to generate the UV enhanced silver recipe

    Parameters
    ----------
    wvl : float
        wavelength in nanometers

    Returns
    -------
    list
        list of tuples of refractive index and thickness,
        compatible with the format Poke likes coatings in
    """
    warnings.warn("This function is deprecated, use load_coating_data instead")

    coating_models = []
    for coating in supported_coatings[5:]:
        coating_models.append(make_index_model(coating))

    # way it is written in excell sheet (backwards)
    thicknesses = [20, 150, 49.1, 33, 53.8, 27.5, 97]
    pag = [
        (coating_models[0](wvl), thicknesses[0]),
        (coating_models[1](wvl), thicknesses[1]),
        (coating_models[2](wvl), thicknesses[2]),
        (coating_models[3](wvl), thicknesses[3]),
        (coating_models[2](wvl), thicknesses[4]),
        (coating_models[3](wvl), thicknesses[5]),
        (coating_models[4](wvl), thicknesses[6])
    ]
    pag.reverse()
    
    # just appending a 1.5 index for substrate, this should have dispersion but it wont matter given the thickness of the reflector
    pag.append(n_BK7)

    return pag
