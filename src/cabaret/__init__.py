import importlib.util

from cabaret.camera import Camera
from cabaret.focuser import Focuser
from cabaret.image import generate_image
from cabaret.observatory import Observatory
from cabaret.queries import Filters, GaiaQuery, GaiaTAPSource
from cabaret.site import Site
from cabaret.sources import Sources
from cabaret.telescope import Telescope

__all__ = [
    "Camera",
    "Filters",
    "Focuser",
    "GaiaTAPSource",
    "generate_image",
    "Observatory",
    "Site",
    "Sources",
    "Telescope",
    "GaiaQuery",
]

if importlib.util.find_spec("matplotlib") is not None:
    from cabaret.plot import plot_image  # noqa: F401

    __all__.append("plot_image")
