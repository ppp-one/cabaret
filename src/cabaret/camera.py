from dataclasses import dataclass


@dataclass
class Camera:
    device_number: int = 0
    width: int = 1000
    height: int = 1000
    max_adu: int = 2**16
    bin_x: int = 1
    bin_y: int = 1
    gain: float = 1.0
    pitch: float = 12.0
    plate_scale: float = 0.35
    well_depth: int = 2**16
    bias: int = 200
    read_noise: float = 4.0
    dark_current: float = 10.0
    sky_background: float = 10.0
    quantum_efficiency: float = 1.0
    collecting_area: float = 10.0
    temperature: float = -10
    sensor_name: str = "gaia-camera-simulated"
    seeing: float = 1.0
