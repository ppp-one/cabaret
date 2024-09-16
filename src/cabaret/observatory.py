from dataclasses import dataclass


@dataclass
class Camera:
    name: str = "gaia-camera-simulated"
    width: int = 1024  # pixels
    height: int = 1024  # pixels
    bin_x: int = 1  # binning factor in x
    bin_y: int = 1  # binning factor in y
    pitch: float = 13.5  # pixel pitch, microns
    plate_scale: float | None = (
        None  # arcsec/pixel (calculated from pitch+telescope if None)
    )
    max_adu: int = 2**16 - 1  # maximum ADU value
    well_depth: int = 2**16 - 1  # electrons
    bias: int = 300  # ADU
    gain: float = 1.0  # e-/ADU
    read_noise: float = 6.2  # e-
    dark_current: float = 0.2  # e-/s
    average_quantum_efficiency: float = 0.8  # fraction


@dataclass
class Telescope:
    focal_length: float = 8.0  # meters
    diameter: float = 1.0  # meters
    collecting_area: float | None = None  # m2 (calculated from diameter if None)


@dataclass
class Site:
    sky_background: float = 150  # for I+z band in Paranal, e-/m2/arcsec2/s
    seeing: float = 1.3  # arcsec
