import numpy as np
from pathlib import Path
import json
from hwo_pol.coatings import (
    make_index_model,
    load_coating_data,
    n_BK7
)

import pytest
import ipdb

COATING_PATH = Path(__file__).parent.parent / 'coating_recipes'
WAVELENGTHS = np.linspace(300, 800, 100)

def test_make_index_model():
    """Test that the index model is created correctly"""

    for file_path in COATING_PATH.iterdir():
        if file_path.suffix == '.csv':
            name = file_path.stem
            index_model = make_index_model(name)

            # also hard-load data
            if name in ['MgF2', 'Al2O3', 'Cr', 'F3', 'SiO2']:
                coating_data = np.genfromtxt(file_path, delimiter=',', skip_header=2)
            else:    
                coating_data = np.genfromtxt(file_path, delimiter=',')

            
            wavelengths = coating_data[:,0]
            n = coating_data[:, 1]
            k = coating_data[:, 2]
            n_data = n + 1j*k
            n_model = index_model(wavelengths)

            np.testing.assert_allclose(n_data, n_model, rtol=1e-3)


def test_load_coating_data():
    """Test that the coating data is loaded correctly
    NOTE that this does not test the recipe itself, just that the data is
    loaded correctly
    """

    for file_path in COATING_PATH.iterdir():
        if file_path.suffix == '.json':
            name = file_path.stem

            # First make sure the data is loaded properly
            try:
                coating_data = load_coating_data(file_path, WAVELENGTHS)

            except Exception as e:
                print(f"Error loading coating {name}: {e}")

            # Next check that the thicknesses are correct, since 
            # we already have `test_make_index_model` to determine the indices

            with open(file_path, "r") as f:
                recipe = json.load(f)
            
            coating_stack = []
            for i, layer in enumerate(recipe["Recipe"]):
                index_model = make_index_model(layer["material"])
                eval_index = index_model(WAVELENGTHS)
                thickness = np.full_like(WAVELENGTHS, layer["thickness_nm"])
                assert eval_index == pytest.approx(coating_data[i][0], abs=1e-3)
                assert thickness == pytest.approx(coating_data[i][1], abs=1e-3)

            # make sure the mirror substrate has no thickness
            np.testing.assert_allclose(coating_data[-1], n_BK7)


if __name__ == '__main__':
    test_load_coating_data()