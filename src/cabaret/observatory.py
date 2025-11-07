import copy
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy.random
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

from cabaret.camera import Camera
from cabaret.fits_manager import FITSManager
from cabaret.focuser import Focuser
from cabaret.image import generate_image, generate_image_stack
from cabaret.queries import Filters
from cabaret.site import Site
from cabaret.sources import Sources
from cabaret.telescope import Telescope


@dataclass
class Observatory:
    """
    Observatory configuration.

    Examples
    --------
    >>> from datetime import datetime, UTC
    >>> dateobs = datetime.now(UTC)
    >>> from cabaret.observatory import Observatory
    >>> observatory = Observatory()

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

    If you have matplotlib installed, you can visualize the image using cabaret's plot
    utility:

    >>> import matplotlib.pyplot as plt
    >>> from cabaret.plot import plot_image
    >>> _ = plot_image(image, title="Simulated Image")
    >>> plt.show()
    """

    name: str
    """Observatory name."""

    camera: Camera
    """Camera configuration."""

    focuser: Focuser
    """Focuser configuration."""

    telescope: Telescope
    """Telescope configuration."""

    site: Site
    """Site configuration."""

    def __init__(
        self,
        name: str = "Observatory",
        camera: Camera | dict = Camera(),
        focuser: Focuser | dict = Focuser(),
        telescope: Telescope | dict = Telescope(),
        site: Site | dict = Site(),
    ):
        if isinstance(camera, dict):
            camera = Camera(**camera)
        if isinstance(focuser, dict):
            focuser = Focuser(**focuser)
        if isinstance(telescope, dict):
            telescope = Telescope(**telescope)
        if isinstance(site, dict):
            site = Site(**site)

        self.name = name
        self.camera = camera
        self.focuser = focuser
        self.telescope = telescope
        self.site = site
        self.__post_init__()

    def __post_init__(self):
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
        dateobs: datetime | None = None,
        light: int = 1,
        filter_band: Filters | str = Filters.G,
        airmass: float = 1.5,
        n_star_limit: int = 2000,
        rng: numpy.random.Generator = numpy.random.default_rng(),
        seed: int | None = None,
        timeout: float | None = None,
        sources: Sources | None = None,
        wcs: WCS | None = None,
        fwhm_multiplier: float = 5.0,
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

        filter_band : Filters or str, optional
            Photometric filter to use for the simulation (default: Filters.G).
        airmass : float, optional
            Airmass value for the observation (default: 1.5).
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
        wcs : WCS or None, optional
            World Coordinate System information for the image.
        fwhm_multiplier : float, optional
            Multiplier to determine the rendering radius around each star
            (default: 5.0).
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
            filter_band=filter_band,
            airmass=airmass,
            n_star_limit=n_star_limit,
            rng=rng,
            seed=seed,
            timeout=timeout,
            sources=sources,
            wcs=wcs,
            fwhm_multiplier=fwhm_multiplier,
        )

    def generate_image_stack(
        self,
        ra: float,
        dec: float,
        exp_time: float,
        dateobs: datetime | None = None,
        light: int = 1,
        filter_band: Filters | str = Filters.G,
        airmass: float = 1.5,
        n_star_limit: int = 2000,
        rng: numpy.random.Generator = numpy.random.default_rng(),
        seed: int | None = None,
        timeout: float | None = None,
        sources: Sources | None = None,
        convert_all_to_adu: bool = True,
        wcs: WCS | None = None,
        fwhm_multiplier: float = 5.0,
    ) -> numpy.ndarray:
        """
        Generate a stack of images from different stages in the image simulation
        pipeline.

        From first to last, the images are:
        1. Base image with bias, dark, and flat applied.
        2. Astronomical image with sources, sky background, and noise.
        3. Final image with pixel defects applied.

        Parameters
        ----------
        ra : float
            Right ascension of the image center (degrees).
        dec : float
            Declination of the image center (degrees).
        exp_time : float
            Exposure time in seconds.
        dateobs : datetime, optional
            Observation date and time (default: now, UTC).
        light : int, optional
            If 1, simulate light exposure; if 0, simulate dark exposure.
        camera : Camera, optional
            Camera configuration.
        focuser : Focuser, optional
            Focuser configuration.
        telescope : Telescope, optional
            Telescope configuration.
        site : Site, optional
            Observatory site configuration.
        filter_band : Filters or str, optional
            The filter to use for the flux column. Default is "G".
        airmass : float, optional
            Airmass value for the observation (default: 1.5).
        n_star_limit : int, optional
            Maximum number of stars to simulate.
        rng : numpy.random.Generator, optional
            Random number generator.
        seed : int or None, optional
            Seed for the random number generator.
        timeout : float or None, optional
            Timeout for Gaia query.
        sources : Sources or None, optional
            Precomputed sources to use instead of querying Gaia.
        convert_all_to_adu : bool, optional
            Whether to convert all images to ADU. Default is False.
        wcs : WCS or None, optional
            World Coordinate System information for the image.
        fwhm_multiplier : float, optional
            Multiplier to determine the rendering radius around each star
            (default: 5.0).

        Returns
        -------
        np.ndarray
            Simulated image stack as a 3D array (uint16, shape (3, height, width)).
            The first slice is the base image, the second is the astronomical image,
            and the third is the ADU image with pixel defects applied.


        """
        return generate_image_stack(
            ra=ra,
            dec=dec,
            exp_time=exp_time,
            dateobs=dateobs,
            light=light,
            camera=self.camera,
            focuser=self.focuser,
            telescope=self.telescope,
            site=self.site,
            filter_band=filter_band,
            airmass=airmass,
            n_star_limit=n_star_limit,
            rng=rng,
            seed=seed,
            timeout=timeout,
            sources=sources,
            convert_all_to_adu=convert_all_to_adu,
            wcs=wcs,
            fwhm_multiplier=fwhm_multiplier,
        )

    def generate_fits_image(
        self,
        ra: float,
        dec: float,
        exp_time: float,
        file_path: str | Path | None = None,
        dateobs: datetime | None = None,
        light: int = 1,
        filter_band: Filters | str = Filters.G,
        airmass: float = 1.5,
        n_star_limit: int = 2000,
        rng: numpy.random.Generator = numpy.random.default_rng(),
        seed: int | None = None,
        timeout: float | None = None,
        sources: Sources | None = None,
        wcs: WCS | None = None,
        fwhm_multiplier: float = 5.0,
        user_header: dict[str, Any] | fits.Header | None = None,
        overwrite: bool = True,
    ) -> fits.HDUList:
        """Generate a simulated FITS image of the sky.

        Parameters
        ----------
        ra : float
            Right ascension of the center of the image in degrees.
        dec : float
            Declination of the center of the image in degrees.
        exp_time : float
            Exposure time in seconds.
        file_path : str or Path, optional
            If provided, the path to save the FITS file.
        user_header : dict or fits.Header, optional
            Additional header keywords to add.
        dateobs : datetime, optional
            Observation date and time in UTC.
        light : int, optional
            If 1, simulate light exposure; if 0, simulate dark exposure.
        filter_band : Filters or str, optional
            Photometric filter to use for the simulation (default: Filters.G).
        airmass : float, optional
            Airmass value for the observation (default: 1.5).
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
        wcs : WCS or None, optional
            World Coordinate System information for the image.
        fwhm_multiplier : float, optional
            Multiplier to determine the rendering radius around each star
            (default: 5.0).
        overwrite : bool, optional
            Whether to overwrite existing file (default: True).

        Returns
        -------
        fits.HDUList
            The generated FITS HDU list.
        """
        image = generate_image(
            ra=ra,
            dec=dec,
            exp_time=exp_time,
            dateobs=dateobs,
            light=light,
            camera=self.camera,
            focuser=self.focuser,
            telescope=self.telescope,
            site=self.site,
            filter_band=filter_band,
            airmass=airmass,
            n_star_limit=n_star_limit,
            rng=rng,
            seed=seed,
            timeout=timeout,
            sources=sources,
            wcs=wcs,
            fwhm_multiplier=fwhm_multiplier,
        )

        if wcs is None:
            wcs = self.camera.get_wcs(SkyCoord(ra=ra, dec=dec, unit="deg"))

        hdu_list = FITSManager.to_hdu_list(
            observatory=self,
            image=image,
            user_header=user_header,
            ra=ra,
            dec=dec,
            exp_time=exp_time,
            dateobs=dateobs,
            wcs=wcs,
        )
        if file_path is not None:
            hdu_list = FITSManager.save(
                observatory=self,
                file_path=file_path,
                hdu_list=hdu_list,
                overwrite=overwrite,
            )
        return hdu_list

    def to_dict(self) -> dict:
        """Convert the Observatory configuration to a dictionary."""
        return asdict(self)

    @classmethod
    def from_dict(cls, config) -> "Observatory":
        """Create an Observatory instance from a configuration dictionary."""
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
        """Save Observatory configuration to a YAML file."""
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
        """Create a deep copy of the Observatory instance."""
        return copy.deepcopy(self)
