import copy
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import numpy.random

from cabaret.camera import Camera
from cabaret.focuser import Focuser
from cabaret.image import generate_image
from cabaret.site import Site
from cabaret.sources import Sources
from cabaret.telescope import Telescope


@dataclass
class Observatory:
    """Observatory configuration.

    Parameters
    ----------
    name : str, optional
        Observatory name.
    camera : Camera, dict, optional
        Camera configuration.
    focuser: Focuser, dict, optional
        Focuser configuration.
    telescope : Telescope, dict, optional
        Telescope configuration.
    site : Site, dict, optional
        Site configuration.

    Examples
    --------
    >>> from datetime import datetime, UTC
    >>> from cabaret.observatory import Observatory
    >>> observatory = Observatory()
    >>> dateobs = datetime.now(UTC)

    Query Gaia for sources and generate an image:
    >>> image = observatory.generate_image(
    ...     ra=12.3323, dec=30.4343, exp_time=10, dateobs=dateobs, seed=0
    ... )

    Or using a set of predefined sources:
    >>> from cabaret.sources import Sources
    >>> sources = Sources.from_arrays(
    ...     ra=[10.64, 10.68], dec=[10.68, 41.22], fluxes=[169435.6, 52203.9]
    ... )
    >>> img = observatory.generate_image(
    ...     ra=sources.ra.deg.mean(),
    ...     dec=sources.dec.deg.mean(),
    ...     exp_time=10,
    ...     seed=0,
    ...     sources=sources,
    ... )
    """

    name: str = "Observatory"
    camera: Camera = field(default_factory=Camera)
    focuser: Focuser = field(default_factory=Focuser)
    telescope: Telescope = field(default_factory=Telescope)
    site: Site = field(default_factory=Site)

    def __post_init__(self):
        if isinstance(self.camera, dict):
            self.camera = Camera(**self.camera)
        if isinstance(self.focuser, dict):
            self.focuser = Focuser(**self.focuser)
        if isinstance(self.telescope, dict):
            self.telescope = Telescope(**self.telescope)
        if isinstance(self.site, dict):
            self.site = Site(**self.site)

        if not isinstance(self.camera, Camera):
            raise ValueError("camera must be an instance of Camera.")
        if not isinstance(self.focuser, Focuser):
            raise ValueError("focuser must be an instance of Focuser.")
        if not isinstance(self.telescope, Telescope):
            raise ValueError("telescope must be an instance of Telescope.")
        if not isinstance(self.site, Site):
            raise ValueError("site must be an instance of Site.")

    def generate_image(
        self,
        ra: float,
        dec: float,
        exp_time: float,
        dateobs: datetime = datetime.now(UTC),
        light: int = 1,
        tmass: bool = True,
        n_star_limit: int = 2000,
        rng: numpy.random.Generator = numpy.random.default_rng(),
        seed: int | None = None,
        timeout: float | None = None,
        sources: Sources | None = None,
        apply_pixel_defects: bool = True,
    ) -> numpy.ndarray:
        """Generate a simulated image of the sky.

        Parameters
        ----------
        ra : float
            Right ascension of the center of the image in degrees.
        dec : float
            Declination of the center of the image in degrees.
        exp_time : float
            Exposure time in seconds.
        dateobs : datetime, optional
            Observation date and time in UTC.
        light : int, optional
            Light pollution level (1-5).
        tmass : bool, optional
            Include 2MASS stars in the image.
        n_star_limit : int, optional
            Maximum number of stars to include in the image.
        rng : numpy.random.Generator, optional
            Random number generator.
        seed : int, optional
            Random number generator seed.
        timeout : float, optional
            The maximum time to wait for the Gaia query to complete, in seconds.
            If None, there is no timeout. By default, it is set to None.
        sources : Sources, optional
            A collection of sources with their sky coordinates and fluxes.
            If provided, these sources will be used instead of querying Gaia.
        apply_pixel_defects : bool, optional
            Whether to apply camera pixel defects (default: True).
        """
        return generate_image(
            ra=ra,
            dec=dec,
            exp_time=exp_time,
            dateobs=dateobs,
            light=light,
            camera=self.camera,
            focuser=self.focuser,
            telescope=self.telescope,
            site=self.site,
            tmass=tmass,
            n_star_limit=n_star_limit,
            rng=rng,
            seed=seed,
            timeout=timeout,
            sources=sources,
            apply_pixel_defects=apply_pixel_defects,
        )

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, config) -> "Observatory":
        return cls(
            name=config.get("name", "Observatory"),
            camera=Camera(**config["camera"]),
            focuser=Focuser(**config.get("focuser", {})),
            telescope=Telescope(**config["telescope"]),
            site=Site(**config["site"]),
        )

    @classmethod
    def load_from_yaml(cls, file_path: str | Path) -> "Observatory":
        """Load Observatory configuration from a YAML file."""
        try:
            import yaml

            with open(file_path) as f:
                config = yaml.safe_load(f)

            return cls.from_dict(config)

        except ImportError:
            raise ImportError(
                "Please install PyYAML to load Observatory configuration from YAML."
            )
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {file_path}")
        except Exception as e:
            raise Exception(f"Error loading Observatory configuration: {e}")

    def save_to_yaml(self, file_path: str | Path):
        try:
            import yaml

            with open(file_path, "w") as f:
                yaml.dump(self.to_dict(), f)

        except ImportError:
            raise ImportError(
                "Please install PyYAML to save Observatory configuration to YAML."
            )
        except Exception as e:
            raise Exception(f"Error saving Observatory configuration: {e}")

    def copy(self) -> "Observatory":
        return copy.deepcopy(self)
