from poke.writing import serialize, deserialize, write_rayfront_to_serial
import os
from pathlib import Path

def write_model(eac_class, filename):
    """Write the model to a msgpack file"""

    # clear the Pathlib objects
    eac_class.raytrace_path = None

    # make directory to save constituent files
    models_path = Path(__file__).parent.parent / "models" / filename
    for wvl in eac_class.wavelengths:
        rf = eac_class.rayfronts[wvl]
        path = models_path / f"rayfront_{wvl}nm"
        models_path.mkdir(exist_ok=True)
        write_rayfront_to_serial(rf, str(path))

    # clear the rayfronts from eac_class
    eac_class.rayfronts = {}

    # serialize the model
    serdata = serialize(eac_class)

    with open(models_path / "eac_model.msgpack", "wb") as f:
        f.write(serdata)
