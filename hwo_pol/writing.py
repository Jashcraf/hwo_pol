from poke.writing import write_rayfront_to_serial, read_serial_to_rayfront
import os
from pathlib import Path
import msgpack
import msgpack_numpy as m

# Classes that are allowed to be re-hydrated
from hwo_pol.models import (
    EAC1,
    EAC4,
    EAC5
)
m.patch()

# patch msgpack to allow int/float keys
old_unpackb = msgpack.unpackb

def unpackb_allow_non_str_keys(*args, **kwargs):
    kwargs['strict_map_key'] = False
    return old_unpackb(*args, **kwargs)

msgpack.unpackb = unpackb_allow_non_str_keys

def serialize(T):
    """serializes a class using msgpack, taken from Poke

    Parameters
    ----------
    T : class
        class to convert to hex code. Used for rayfronts

    Returns
    -------
    serdat
        serial data corresponding to class T
    """
    glb = globals()
    Tname = T.__class__.__name__
    # assert Tname in glb, 'class must exist in globals in order to be re-hydrateable, with the same constraint'

    # now we make our storage format.  It will be:
    # 1) a header with the class name
    # 2) the content of vars(T)
    serdat = (Tname, vars(T))
    return msgpack.packb(serdat)


class MsgpackTrickerEmpty:
    """dummy class to trick msgpack
    """
    pass


def deserialize(buf):
    """deserializes a class using msgpack, taken from Poke

    Parameters
    ----------
    buf : serdat
        serial data coorresponding to class

    Returns
    -------
    class
        deserialized class, typically a rayfront
    """
    e = MsgpackTrickerEmpty()
    Tname, varzzz = msgpack.unpackb(buf, use_list=True)
    for k, v in varzzz.items():
        setattr(e, k, v)
    e.__class__ = globals()[Tname]
    return e


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


def read_model(folder_name):

    # Read the data
    parent_path = Path(__file__).parent.parent / "models"
    folder_path = parent_path / folder_name

    with open(folder_path / "eac_model.msgpack", "rb") as f:
        eac_class = f.read()

    eac_class = deserialize(eac_class)

    # Put the rayfronts back in the eac model
    for wvl in eac_class.wavelengths:
        load_rf = read_serial_to_rayfront(folder_path / f"rayfront_{wvl}nm.msgpack")
        eac_class.rayfronts[wvl] = load_rf

    return eac_class
