from dataclasses import dataclass


@dataclass
class Camera:
    device_number: int = 0
    width: int = 1000
    height: int = 1000
    bin_x: int = 1
    bin_y: int = 1
    gain: float = 1.0
    pitch: float = 0.0
    plate_scale: float = 1.0
    well_depth: int = 1.0
    bias: int = 0
    max_adu: int = 1000
    read_noise: float = 0
    dark_current: float = 0
    sky_background: float = 10.0
    quantum_efficiency: float = 1.0
    collecting_area: float = 10.0
    temperature: float = -10
    sensor_name: str = "gaia-camera-simulated"
