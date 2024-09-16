from dataclasses import dataclass


@dataclass
class Camera:
    name: str = "gaia-camera-simulated"
    width: int = 1000
    height: int = 1000
    bin_x: int = 1
    bin_y: int = 1
    pitch: float = 12
    max_adu: int = 2**16 - 1
    well_depth: int = 2**16 - 1
    bias: int = 200
    gain: float = 1.0
    read_noise: float = 4.0
    dark_current: float = 10.0
    average_quantum_efficiency: float = 1.0
    temperature: float = -10


@dataclass
class Telescope:
    focal_length: float = 8.0
    diameter: float = 1.0


@dataclass
class Site:
    sky_background: float = 10.0
    seeing: float = 1.0
