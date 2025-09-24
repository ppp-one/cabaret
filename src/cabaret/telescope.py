from dataclasses import dataclass

import numpy as np


@dataclass
class Telescope:
    """
    Telescope configuration.

    Attributes
    ----------
    focal_length : float
        Focal length of the telescope in meters.
    diameter : float
        Diameter of the telescope in meters.
    collecting_area : float or None
        Collecting area of the telescope in square meters. If None, it is calculated
        from the diameter.
    """

    focal_length: float = 8.0  # meters
    diameter: float = 1.0  # meters
    collecting_area: float | None = None  # square meters

    def __post_init__(self):
        if self.collecting_area is None:
            self.collecting_area = np.pi * (self.diameter / 2) ** 2  # square meters
