import logging
from datetime import UTC, datetime

import numpy as np
import numpy.random
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
from astropy.time import Time

from cabaret.camera import Camera
from cabaret.focuser import Focuser
from cabaret.queries import Filters, get_gaia_sources
from cabaret.site import Site
from cabaret.sources import Sources
from cabaret.telescope import Telescope

logger = logging.getLogger("cabaret")


def moffat_profile(
    x: np.ndarray,
    y: np.ndarray,
    x0: float,
    y0: float,
    FWHM: float,
    beta: float = 2.5,
) -> np.ndarray:
    """
    Compute a normalized Moffat profile centered at (x0, y0).

    Parameters
    ----------
    x, y : np.ndarray
        Meshgrid arrays for the pixel coordinates.
    x0, y0 : float
        Center of the profile.
    FWHM : float
        Full width at half maximum of the profile.
    beta : float, optional
        Moffat beta parameter (default: 2.5).

    Returns
    -------
    np.ndarray
        Normalized Moffat profile evaluated on the grid.
    """
    # https://nbviewer.org/github/ysbach/AO_2017/blob/master/04_Ground_Based_Concept.ipynb#1.2.-Moffat
    # FWHM =  2 * R * (2**(1/beta) - 1)**0.5

    R = (FWHM / 2) * (1 / (2 ** (1 / beta) - 1) ** 0.5)
    A = (beta - 1) / (np.pi * R**2)

    r_squared = (x - x0) ** 2 + (y - y0) ** 2

    mp = A * (1 + (r_squared / R**2)) ** (-beta)

    mp_sum = np.sum(mp)

    return mp / mp_sum


def generate_star_image_slow(
    pos: np.ndarray,
    fluxes: list[float],
    FWHM: float,
    frame_size: tuple[int, int],
    rng: numpy.random.Generator,
) -> np.ndarray:
    """
    Render stars onto an image using a slow loop-based approach.

    Parameters
    ----------
    pos : np.ndarray
        Pixel positions of stars (shape: 2, n_stars).
    fluxes : list
        list of fluxes for each star.
    FWHM : float
        Full width at half maximum for the Moffat profile.
    frame_size : tuple
        Size of the output image (width, height).
    rng : numpy.random.Generator
        Random number generator for Poisson noise.

    Returns
    -------
    np.ndarray
        Image with rendered stars.
    """
    x = np.linspace(0, frame_size[0] - 1, frame_size[0])
    y = np.linspace(0, frame_size[1] - 1, frame_size[1])
    xx, yy = np.meshgrid(x, y)

    image = np.zeros(frame_size).T
    for i, flux in enumerate(fluxes):
        x0 = pos[0][i]
        y0 = pos[1][i]
        star = rng.poisson(flux) * moffat_profile(xx, yy, x0, y0, FWHM)
        image += star

    return image


def generate_star_image(
    pos: np.ndarray,
    fluxes: list[float],
    FWHM: float,
    frame_size: tuple[int, int],
    rng: numpy.random.Generator,
) -> np.ndarray:
    """
    Render stars onto an image using a fast, windowed approach.

    Parameters
    ----------
    pos : np.ndarray
        Pixel positions of stars (shape: 2, n_stars).
    fluxes : list
        list of fluxes for each star.
    FWHM : float
        Full width at half maximum for the Moffat profile.
    frame_size : tuple
        Size of the output image (width, height).
    rng : numpy.random.Generator
        Random number generator for Poisson noise.

    Returns
    -------
    np.ndarray
        Image with rendered stars.
    """
    x = np.linspace(0, frame_size[0] - 1, frame_size[0])
    y = np.linspace(0, frame_size[1] - 1, frame_size[1])
    xx, yy = np.meshgrid(x, y)

    render_radius = FWHM * 5  # render 5 FWHM around the star

    image = np.zeros(frame_size).T
    for i, flux in enumerate(fluxes):
        x0 = pos[0][i]
        y0 = pos[1][i]
        if x0 < 0 or x0 >= frame_size[0] or y0 < 0 or y0 >= frame_size[1]:
            # print(f"Star {i} is outside the frame.")
            continue
        x_min, x_max = int(x0 - render_radius), int(x0 + render_radius)
        y_min, y_max = int(y0 - render_radius), int(y0 + render_radius)
        x_min, x_max = max(0, x_min), min(x_max, frame_size[0] - 1)
        y_min, y_max = max(0, y_min), min(y_max, frame_size[1] - 1)

        star = rng.poisson(flux) * moffat_profile(
            xx[y_min : y_max + 1, x_min : x_max + 1],
            yy[y_min : y_max + 1, x_min : x_max + 1],
            x0,
            y0,
            FWHM,
        )
        image[y_min : y_max + 1, x_min : x_max + 1] += star

    return image


def get_sources(
    center: SkyCoord,
    fovx: float,
    fovy: float,
    dateobs: datetime,
    n_star_limit: int,
    filter_band: Filters | str,
    timeout: float | None,
    sources: Sources | None = None,
) -> Sources:
    """Get sources from Gaia or use provided sources."""
    if not isinstance(sources, Sources):
        sources = get_gaia_sources(
            center=center,
            fov=u.Quantity((fovx * 1.5, fovy * 1.5), "deg"),
            dateobs=dateobs,
            limit=n_star_limit,
            filter_band=filter_band,
            timeout=timeout,
        )
    return sources


def add_sun_sky_background(
    image: np.ndarray,
    site: Site,
    telescope: Telescope,
    camera: Camera,
    exp_time: float,
    dateobs: datetime,
    logger: logging.Logger,
) -> np.ndarray:
    """Add sky background and sunlight if location is specified."""
    if site.latitude is not None and site.longitude is not None:
        logger.info(
            "Since location specified, calculating sunlight brightness"
            " based on sun's position"
        )
        location = EarthLocation(
            lat=u.Quantity(site.latitude, "deg"), lon=u.Quantity(site.longitude, "deg")
        )
        obs_time = Time(dateobs, scale="utc")
        sun = get_sun(obs_time)
        altaz_frame = AltAz(obstime=obs_time, location=location)
        sun_altaz = sun.transform_to(altaz_frame)
        sun_altitude: float = sun_altaz.alt.degree  # type: ignore
        logger.info(f"Sun altitude: {sun_altitude:.5f} deg at {dateobs} UTC")

        a, b, c = (
            np.float64(4533508.655833181),
            np.float64(0.3937301435229289),
            np.float64(-0.7907223506021084),
        )  # calibrated for I+z band in Paranal

        sky_brightness = a * b ** (c * sun_altitude)  # e-/m2/arcsec2/s
        logger.info(f"sky_brightness (e-/m2/arcsec2/s): {sky_brightness}")

        sky_e = (
            sky_brightness * telescope.collecting_area * camera.plate_scale**2
        )  # e-/s
        logger.info(f"sky_e (e-/s): {sky_e}")

        image += np.random.poisson(
            np.ones((camera.height, camera.width)).astype(np.float64) * sky_e * exp_time
        ).astype(np.float64)
    return image


def add_stars(
    image: np.ndarray,
    sources: Sources,
    camera: Camera,
    focuser: Focuser,
    telescope: Telescope,
    site: Site,
    exp_time: float,
    rng: numpy.random.Generator,
    ra: float | None,
    dec: float | None,
) -> np.ndarray:
    """Add stars to the image using the Moffat profile and sky background."""
    if len(sources) > 0:
        fluxes = (
            sources.fluxes
            * camera.average_quantum_efficiency
            * telescope.collecting_area
            * exp_time
        )  # [electrons]

        if ra is None or dec is None:
            ra, dec = sources.ra.deg.mean(), sources.dec.deg.mean()

        wcs = camera.get_wcs(SkyCoord(ra=ra, dec=dec, unit="deg"))
        gaias_pixel = sources.to_pixel(wcs)

        stars = generate_star_image(
            gaias_pixel,
            fluxes,
            focuser.seeing_multiplier * site.seeing / camera.plate_scale,
            (camera.width, camera.height),
            rng=rng,
        ).astype(np.float64)

        sky_background = (
            site.sky_background * telescope.collecting_area * camera.plate_scale**2
        )  # [e-/s]

        image = image + rng.poisson(
            np.ones((camera.height, camera.width)).astype(np.float64)
            * sky_background
            * exp_time
        ).astype(np.float64)

        image += stars
    return image


def add_stars_and_sky(
    base: np.ndarray,
    ra: float,
    dec: float,
    exp_time: float,
    dateobs: datetime,
    light: int,
    camera: Camera,
    focuser: Focuser,
    telescope: Telescope,
    site: Site,
    filter_band: Filters | str,
    n_star_limit: int,
    rng: numpy.random.Generator,
    timeout: float | None,
    sources: Sources | None,
) -> np.ndarray:
    """Add stars and sky background to the base image."""
    if light == 1:
        center = SkyCoord(ra=ra, dec=dec, unit="deg")
        fovx = (
            (1 / np.abs(np.cos(center.dec.rad)))  # type: ignore
            * camera.width
            * camera.plate_scale
            / 3600
        )
        fovy = np.sqrt(2) * camera.height * camera.plate_scale / 3600
        logger.info("Querying Gaia for sources...")
        sources = get_sources(
            center=center,
            fovx=fovx,
            fovy=fovy,
            dateobs=dateobs,
            n_star_limit=n_star_limit,
            filter_band=filter_band,
            timeout=timeout,
            sources=sources,
        )
        logger.info(f"Found {len(sources)} sources (user set limit of {n_star_limit}).")
        image = base
        image = add_sun_sky_background(
            image, site, telescope, camera, exp_time, dateobs, logger
        )
        image = add_stars(
            image=image,
            sources=sources,
            camera=camera,
            focuser=focuser,
            telescope=telescope,
            site=site,
            exp_time=exp_time,
            rng=rng,
            ra=ra,
            dec=dec,
        )
    else:
        image = base
    return image


def generate_image(
    ra: float,
    dec: float,
    exp_time: float,
    dateobs: datetime = datetime.now(UTC),
    light: int = 1,
    camera: Camera = Camera(),
    focuser: Focuser = Focuser(),
    telescope: Telescope = Telescope(),
    site: Site = Site(),
    filter_band: Filters | str = Filters.G,
    n_star_limit: int = 2000,
    rng: numpy.random.Generator = numpy.random.default_rng(),
    seed: int | None = None,
    timeout: float | None = None,
    sources: Sources | None = None,
    apply_pixel_defects: bool = True,
) -> np.ndarray:
    """
    Generate a simulated astronomical image.

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
    apply_pixel_defects : bool, optional
        Whether to apply pixel defects to the image.

    Returns
    -------
    np.ndarray
        Simulated image as a 2D numpy array (uint16).
    """
    if seed is not None:
        rng = numpy.random.default_rng(seed)

    if camera.plate_scale is None:
        camera.set_plate_scale_from_focal_length(telescope.focal_length)

    base = camera.make_base_image(exp_time=exp_time, rng=rng)

    if light == 1:
        image = add_stars_and_sky(
            base=base,
            ra=ra,
            dec=dec,
            exp_time=exp_time,
            dateobs=dateobs,
            light=light,
            camera=camera,
            focuser=focuser,
            telescope=telescope,
            site=site,
            filter_band=filter_band,
            n_star_limit=n_star_limit,
            rng=rng,
            timeout=timeout,
            sources=sources,
        )
    else:
        image = base

    if apply_pixel_defects:
        image = camera.apply_pixel_defects(image, exp_time)

    image = camera.to_adu_image(image)

    return image


def generate_image_stack(
    ra: float,
    dec: float,
    exp_time: float,
    dateobs: datetime = datetime.now(UTC),
    light: int = 1,
    camera: Camera = Camera(),
    focuser: Focuser = Focuser(),
    telescope: Telescope = Telescope(),
    site: Site = Site(),
    filter_band: Filters | str = Filters.G,
    n_star_limit: int = 2000,
    rng: numpy.random.Generator = numpy.random.default_rng(),
    seed: int | None = None,
    timeout: float | None = None,
    sources: Sources | None = None,
    apply_pixel_defects: bool = True,
) -> np.ndarray:
    """
    Generate a stack of images from different stages in the image simulation pipeline.

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
    apply_pixel_defects : bool, optional
        Whether to apply pixel defects to the image.

    Returns
    -------
    np.ndarray
        Simulated image stack as a 3D numpy array (uint16, shape (3, height, width)).
        The first slice is the base image, the second is the astronomical image,
        and the third is the ADU image with pixel defects applied.


    """
    if seed is not None:
        rng = numpy.random.default_rng(seed)

    if camera.plate_scale is None:
        camera.set_plate_scale_from_focal_length(telescope.focal_length)

    base = camera.make_base_image(exp_time=exp_time, rng=rng)

    if light == 1:
        image = add_stars_and_sky(
            base=np.zeros_like(base),
            ra=ra,
            dec=dec,
            exp_time=exp_time,
            dateobs=dateobs,
            light=light,
            camera=camera,
            focuser=focuser,
            telescope=telescope,
            site=site,
            filter_band=filter_band,
            n_star_limit=n_star_limit,
            rng=rng,
            timeout=timeout,
            sources=sources,
        )
    else:
        image = base

    if apply_pixel_defects:
        adu_image = camera.apply_pixel_defects(image.copy(), exp_time)

    adu_image = camera.to_adu_image(adu_image)

    return np.stack([base, image, adu_image])


if __name__ == "__main__":
    import importlib.util

    camera = Camera(width=2000, height=2000)
    telescope = Telescope()
    site = Site(seeing=1.3, sky_background=350)
    exp_time = 0.1  # [s]

    logger.info("Generating image...")

    # example usage
    image = generate_image(
        323.36152,
        -0.82325,
        exp_time=exp_time,
        camera=camera,
        telescope=telescope,
        site=site,
    )

    science = image  # - camera.dark_current / camera.gain * exp_time - camera.bias

    if importlib.util.find_spec("matplotlib") is not None:
        import matplotlib.pyplot as plt

        from cabaret.plot import plot_image

        print("Plotting image...")
        med = np.median(science)
        std = np.std(science)
        print(med, std)

        fig, ax = plt.subplots()
        plot_image(image, ax=ax, title="Simulated Image")
        plt.show()
