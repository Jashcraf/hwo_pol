"""
Models for saving state of different EACs
"""

from pathlib import Path
from warnings import warn
import ipdb

# The Poke stuff
from poke.poke_core import Rayfront
import poke.plotting as plot
from poke.writing import write_rayfront_to_serial, jones_to_fits
from poke.poke_math import np

# the hwo_pol stuff
from hwo_pol.coatings import load_coating_data
from hwo_pol.conf import eac1_barrel_lens_data, barrel_lens

PROJECT_PATH = Path.home() / "Box/Polarization PSDD"
JONES_PATH = PROJECT_PATH / "jones_pupils"
DESIGN_PATH = PROJECT_PATH / "Optical Designs"
COATINGS_PATH = Path(__file__).parent.parent / "coating_recipes"
n_BK7 = 1.4801

class EAC:
    def __init__(self, num_rays, fov, wavelengths, coating, coating_internal):
        """Initialize EAC Base Class"""
        self.num_rays = num_rays
        self.fov = fov # maximum field of view
        self.wavelengths = wavelengths
        self.coating = coating
        self.coating_internal = coating_internal
        self._init_rayfronts()

    def _init_rayfronts(self):
        # Init the Poke Rayfront objects
        self.rayfronts = {}
        for wavelength in self.wavelengths:
            self.rayfronts[wavelength] = Rayfront(self.num_rays,
                                                  wavelength,
                                                  self.pupil_radius,
                                                  self.fov,
                                                  circle=False,
                                                  fov=[self.fov_x, self.fov_y])

    
    def normalize_aloc(self):
        """Construct the rays that determine the local coordinate system
        for the Jones pupil computation
        """

        assert hasattr(self, "aloc"), "Need to define aloc"
        assert hasattr(self, "aloc_p"), "Need to define aloc_p"

        # Compute the local X-axis from a cross product of the aloc and aloc_p
        self.aloc /= np.linalg.norm(self.aloc)
        self.aloc_p /= np.linalg.norm(self.aloc_p)
        self.exit_x = np.cross(self.aloc_p, self.aloc)
        self.exit_x = np.cross(self.aloc, self.aloc_p)
        self.exit_x /= np.linalg.norm(self.exit_x)
    
    
    def trace_rays(self):
        """Trace rays through the pupil"""
        # convert Winows/Posix path to string
        raytrace_path = str(self.raytrace_path)
        raytrace_path = '"' + raytrace_path + '"'

        try:
            for wavelength in self.wavelengths:
                self.rayfronts[wavelength].trace_rayset(raytrace_path,
                                                        ref_surf=self.ref_surf,
                                                        wave=wavelength)

        except AttributeError as e:
            warn(f"Error in rayfront: {e}, no poke Rayfront initialized")

    def compute_jones_pupil_by_wavelength(self):
        """Compute the Jones pupil from the rayfront"""
        try:
            for wavelength in self.wavelengths:
                self.rayfronts[wavelength].compute_jones_pupil(exit_x=self.exit_x,
                                                               aloc=self.aloc)

        except AttributeError as e:
            warn(f"Error in rayfront: {e}, no poke Rayfront initialized")

    def construct_surflist_by_wavelength(self):
        
        assert hasattr(self, "surfnums"), "Need to define surfnums"
        assert hasattr(self, "is_ota"), "Need to define is_ota"
        
        self.surflist_by_wavelength = {}


        for wavelength in self.wavelengths:
            coating_ota = load_coating_data(COATINGS_PATH / self.coating, wavelength)
            coating_internal = load_coating_data(COATINGS_PATH / self.coating_internal,
                                                wavelength)
            surflist = [] 
            for surfnum, isota in zip(self.surfnums, self.is_ota):
                
                if isota:
                    surf = {
                        "surf": surfnum,
                        "coating": coating_ota,
                        "mode": "reflect",
                    }

                else:
                    surf = {
                        "surf": surfnum,
                        "coating": coating_internal,
                        "mode": "reflect",
                    }

                surflist.append(surf)

            self.rayfronts[wavelength].as_polarized(surflist)
            self.surflist_by_wavelength[wavelength] = surflist 


class EAC1(EAC):
    def __init__(self, num_rays, fov, wavelengths, coating, coating_internal,
                 ending="TCA"):

        assert ending in ["TCA", "OTA"], "Ending must be TCA or OTA"

        # Hard-coded values from CODE V
        self.pupil_radius = 7200 / 2
        self.max_fov = -0.07302
        self.ref_surf = 4

        if ending == "TCA":
            self.aloc = np.array([0.07893,    -0.06487,    -0.99477])
            self.aloc_p = np.array([0.07894,    -0.09788,    -0.99206])
            self.surfnums = [5, 6, 8, 9, 10, 11]
            self.is_ota = [True, True, False, False, False, False]
        
        elif ending == "OTA":
            self.aloc = np.array([0.00000, 0.01504, 0.99989])
            self.aloc_p = np.array([ 0.00000, -0.01049, 0.99995])
            self.surfnums = [5, 6]
            self.is_ota = [True, True]

        # Compute the local X-axis from a cross product of the aloc and aloc_p
        self.normalize_aloc()

        # Joe Howard says this is the coronagraph FOV
        self.fov_x = 0
        self.fov_y = -0.0730

        # Get the design path
        self.raytrace_path = DESIGN_PATH / f"EAC1/EAC1.len"

        # Init methods from the base class
        super().__init__(num_rays, fov, wavelengths, coating, coating_internal)

        # Assemble Surface dictionaries keyed by wavelength
        self.construct_surflist_by_wavelength()


class EAC1_Barrel_Decrease(EAC1):
    def __init__(self, barrel_len, num_rays, fov, wavelengths, coating, coating_internal,
                 ending="OTA"):

        assert barrel_len in barrel_lens, f"Barrel length {barrel_len} not supported"

        super().__init__(num_rays, fov, wavelengths, coating, coating_internal,
                        ending=ending)

        # Modify based on barrel length
        self.aloc = -eac1_barrel_lens_data[barrel_len]["aloc"]
        self.aloc_p = eac1_barrel_lens_data[barrel_len]["alocp"]
        self.raytrace_path = eac1_barrel_lens_data[barrel_len]["path"]
        self.normalize_aloc()

        # Update the surface list
        self.surfnums = [5, 6, 8, 9, 10,
                         11, 13, 14, 16, 19,
                         20]

        self.is_ota = [True, True, False, False, False,
                       False, False, False, False, False,
                       False]

        # Re-compute the surface list
        self.construct_eac1_surflist_by_wavelength()
    
    def construct_eac1_surflist_by_wavelength(self):
        """Alternative constructor to handle the transmission through the dichroic
        """ 
        assert hasattr(self, "surfnums"), "Need to define surfnums"
        assert hasattr(self, "is_ota"), "Need to define is_ota"
        
        self.surflist_by_wavelength = {}

        for wavelength in self.wavelengths:
            coating_ota = load_coating_data(COATINGS_PATH / self.coating, wavelength)
            coating_internal = load_coating_data(COATINGS_PATH / self.coating_internal,
                                                wavelength)
            surflist = [] 
            for surfnum, isota in zip(self.surfnums, self.is_ota):
                
                # Need a special case for the transmission through the dichroic
                if surfnum == 13:
                    surf = {
                        "surf": surfnum,
                        "coating": (1., n_BK7),
                        "mode": "transmit",
                    }

                elif surfnum == 14:
                    surf = {
                        "surf": surfnum,
                        "coating": (n_BK7, 1.),
                        "mode": "transmit",
                    }

                # Remaining surfaces
                else:
                    if isota:
                        surf = {
                            "surf": surfnum,
                            "coating": coating_ota,
                            "mode": "reflect",
                        }

                    else:
                        surf = {
                            "surf": surfnum,
                            "coating": coating_internal,
                            "mode": "reflect",
                        }

                surflist.append(surf)

            self.rayfronts[wavelength].as_polarized(surflist)
            self.surflist_by_wavelength[wavelength] = surflist 


class EAC4(EAC):
    def __init__(self, num_rays, fov, wavelengths, coating, coating_internal,
                 ending="OTA"):

        assert ending in ["TCA", "OTA"], "Ending must be TCA or OTA"

        # Hard-coded values from CODE V
        self.pupil_radius = 6800 / 2
        self.max_fov = -0.025
        self.ref_surf = 3

        if ending == "TCA":
            self.aloc = np.array([0.00000, 0.99994, 0.01072])
            self.aloc_p = np.array([0.00000, 0.99986, -0.01685])
            self.surfnums = [6, 7, 9, 10]
            self.is_ota = [True, True, False, False]
        
        elif ending == "OTA":
            self.aloc = np.array([ 0.00000, 0.03056, 0.99953])
            self.aloc_p = np.array([ 0.00000, -0.00554, 0.99998])
            self.surfnums = [6, 7]
            self.is_ota = [True, True]

        # Compute the local X-axis from a cross product of the aloc and aloc_p
        self.normalize_aloc()
        
        # Joe Howard says this is the coronagraph FOV
        self.fov_x = 0
        self.fov_y = -0.025

        # Get the design path
        self.raytrace_path = DESIGN_PATH / f"EAC4/EAC4.len"

        # Init methods from the base class
        super().__init__(num_rays, fov, wavelengths, coating, coating_internal)

        # Assemble Surface dictionaries keyed by wavelength
        self.construct_surflist_by_wavelength()


# ----------------------------------------------------------------------------
# EAC4 Variants for Barrel Length Study
# ----------------------------------------------------------------------------


class EAC4_150(EAC4):
    def __init__(self, num_rays, fov, wavelengths, coating, coating_internal,
                 ending="OTA"):

        super().__init__(num_rays, fov, wavelengths, coating, coating_internal,
                         ending=ending)

        # Update system data
        self.pupil_radius = 7000 / 2
        self.fov_x = 0
        self.fov_y = -0.1

        # Update aloc, aloc_p, and surface numbers
        self.aloc = np.array([0.00000,     0.01816,     0.99984])
        self.aloc_p = np.array([0.00000,    -0.01182,     0.99993])
        self.surfnums = [7, 8]
        self.is_ota = [True, True]
        self.normalize_aloc()
        self.raytrace_path = DESIGN_PATH / f"EAC4/EAC4_150.len"
        self.ref_surf = 7

        # Re-init rayfronts
        self._init_rayfronts()

        # Re-compute the surface list
        self.construct_surflist_by_wavelength()

class EAC4_165(EAC4):
    def __init__(self, num_rays, fov, wavelengths, coating, coating_internal,
                 ending="OTA"):

        super().__init__(num_rays, fov, wavelengths, coating, coating_internal,
                         ending=ending)

        # Update system data
        self.pupil_radius = 7000 / 2
        self.fov_x = 0
        self.fov_y = -0.1
        
        # Update aloc and aloc_p
        self.aloc = np.array([0.00000,     0.01222,     0.99993])
        self.aloc_p = np.array([0.00000,    -0.02424,     0.99971])
        self.surfnums = [7, 8]
        self.is_ota = [True, True]
        self.normalize_aloc()
        self.raytrace_path = DESIGN_PATH / f"EAC4/EAC4_165.len"
        
        # Re-init rayfronts
        self._init_rayfronts()

        # Re-compute the surface list
        self.construct_surflist_by_wavelength()


class EAC4_180(EAC4):
    def __init__(self, num_rays, fov, wavelengths, coating, coating_internal,
                 ending="OTA"):

        super().__init__(num_rays, fov, wavelengths, coating, coating_internal,
                         ending=ending)
        
        # Update system data
        self.pupil_radius = 7000 / 2
        self.fov_x = 0
        self.fov_y = -0.1

        # Update aloc and aloc_p
        self.aloc = np.array([0.00000,     0.00475,     0.99999])
        self.aloc_p = np.array([0.00000,    -0.03644,     0.99934])
        self.surfnums = [7, 8]
        self.is_ota = [True, True]
        self.normalize_aloc()
        self.raytrace_path = DESIGN_PATH / f"EAC4/EAC4_180.len"
        
        # Re-init rayfronts
        self._init_rayfronts()

        # Re-compute the surface list
        self.construct_surflist_by_wavelength()


class EAC5(EAC):
    def __init__(self, num_rays, fov, wavelengths, coating, coating_internal,
                 ending="OTA"):
        assert ending in ["TCA", "OTA"], "Ending must be TCA or OTA"

        # Hard-coded values from CODE V
        self.pupil_radius = 10_000 / 2
        self.max_fov = -0.025
        self.ref_surf = 3

        if ending == "TCA":
            self.aloc = -np.array([0.00000, 0.99966, -0.02594])
            self.aloc_p = np.array([0.00000, 0.99988, 0.01544])
            self.surfnums = [6, 7, 9, 10, 11]
            self.is_ota = [True, True, False, False, False]
        
        elif ending == "OTA":
            self.aloc = np.array([0.00000, -0.00867, 0.99996])
            self.aloc_p = np.array([0.00000, -0.06192, 0.99808])
            self.surfnums = [6, 7]
            self.is_ota = [True, True]

        # Compute the local X-axis from a cross product of the aloc and aloc_p
        self.normalize_aloc()
        
        # Joe Howard says this is the coronagraph FOV
        self.fov_x = 0
        self.fov_y = -0.025

        # Get the design path
        self.raytrace_path = DESIGN_PATH / f"EAC5/EAC5.len"

        # Init methods from the base class
        super().__init__(num_rays, fov, wavelengths, coating, coating_internal)

        # Assemble Surface dictionaries keyed by wavelength
        self.construct_surflist_by_wavelength()