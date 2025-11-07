import logging
from datetime import UTC, datetime

import numpy as np
import numpy.random
from astropy import units as u
from astropy.coordinates import AltAz, Angle, EarthLocation, SkyCoord, get_sun
from astropy.time import Time
from astropy.wcs import WCS

from cabaret.camera import Camera
from cabaret.focuser import Focuser
from cabaret.queries import Filters, GaiaQuery
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


def scintillation_noise(
    r: float,
    t: float,
    N_star: float,
    h: float = 2440,
    C: float = 1.56,
    airmass: float = 1.5,
) -> float:
    """
    Calculate the scintillation noise for a given set of parameters.

    Parameters
    ----------
    r : float
        Aperture radius in meters.
    t : float
        Exposure time in seconds.
    N_star : float
        Number of stars.
    h : float
        Altitude of the observatory in meters. Default is 2440 for Paranal Observatory.
    C : float
        Empirical coefficient. Default is 1.56, optimized for the
        20-cm NGTS telescopes at Paranal Observatory.
    airmass : float, optional
        Airmass value. Default is 1.5.

    Returns
    -------
    float
        The calculated scintillation noise.

    Reference
    ---------
    https://academic.oup.com/mnras/article/509/4/6111/6442285
    """

    return (
        np.sqrt(
            1e-5
            * C**2
            * pow(2 * r, -4 / 3)
            * t**-1
            * airmass**3
            * np.exp(-2 * h / 8000)
        )
        * N_star
        * t
    )


def generate_star_image(
    pos: np.ndarray,
    fluxes: list[float],
    FWHM: float,
    frame_size: tuple[int, int],
    rng: numpy.random.Generator,
    telescope_aperture: float,
    site_elevation: float,
    exp_time: float,
    airmass: float = 1.5,
    fwhm_multiplier: float = 5.0,
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
    fwhm_multiplier : float, optional
        Multiplier to determine the rendering radius around each star (default: 5.0).

    Returns
    -------
    np.ndarray
        Image with rendered stars.
    """
    x = np.linspace(0, frame_size[0] - 1, frame_size[0])
    y = np.linspace(0, frame_size[1] - 1, frame_size[1])
    xx, yy = np.meshgrid(x, y)

    render_radius = FWHM * fwhm_multiplier  # render n * FWHM around the star

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

        scint_noise = scintillation_noise(
            r=telescope_aperture / 2,
            t=exp_time,
            N_star=flux,
            h=site_elevation,
            C=1.56,
            airmass=airmass,
        )

        star_flux = rng.poisson(flux) + rng.normal(0, scint_noise)

        if star_flux < 0:
            star_flux = 0

        star = star_flux * moffat_profile(
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
    radius: Angle,
    dateobs: datetime,
    n_star_limit: int,
    filter_band: Filters | str,
    timeout: float | None,
    sources: Sources | None = None,
) -> Sources:
    """Get sources from Gaia or use provided sources."""
    if not isinstance(sources, Sources):
        logger.info("Querying Gaia for sources...")
        sources = GaiaQuery.get_sources(
            center=center,
            radius=radius,
            dateobs=dateobs,
            limit=n_star_limit,
            filter_band=filter_band,
            timeout=timeout,
        )
        logger.info(f"Found {len(sources)} sources (user set limit of {n_star_limit}).")
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
    wcs: WCS | None,
    dateobs: datetime | None,
    airmass: float | None,
    fwhm_multiplier: float = 5.0,
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
            ra, dec = sources.center

        if wcs is None:
            wcs = camera.get_wcs(SkyCoord(ra=ra, dec=dec, unit="deg"))
        gaias_pixel = sources.to_pixel(wcs)
        on_camera_mask = camera.on_camera_mask(gaias_pixel)
        logger.debug(
            f"Number of stars on camera: {on_camera_mask.sum()}"
            f" out of {len(sources)} total stars "
            f"({on_camera_mask.sum() / len(sources) * 100})%."
        )

        if airmass is not None:
            logger.debug(f"Using provided airmass: {airmass:.3f}.")
        elif site.latitude and site.longitude:
            location = EarthLocation(
                lat=u.Quantity(site.latitude, "deg"),
                lon=u.Quantity(site.longitude, "deg"),
            )
            altitude_azimuth_frame = AltAz(obstime=Time(dateobs), location=location)
            center_coord = SkyCoord(ra=ra, dec=dec, unit="deg")
            center_altaz = center_coord.transform_to(altitude_azimuth_frame)
            zenith_angle = 90.0 - center_altaz.alt.degree  # degrees
            airmass = 1.0 / np.cos(np.radians(zenith_angle))

            if airmass < 1.0:
                logger.warning(
                    f"Calculated airmass {airmass:.3f} is less than 1.0."
                    " Setting airmass to 1.0."
                )
                airmass = 1.0  # minimum airmass is 1.0
        else:
            airmass = 1.5  # default value

        stars = generate_star_image(
            gaias_pixel[:, on_camera_mask],
            fluxes[on_camera_mask],
            focuser.seeing_multiplier * site.seeing / camera.plate_scale,
            (camera.width, camera.height),
            rng=rng,
            telescope_aperture=telescope.diameter,
            site_elevation=site.elevation if site.elevation is not None else 0,
            exp_time=exp_time,
            airmass=airmass,
            fwhm_multiplier=fwhm_multiplier,
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
    ra: float | None,
    dec: float | None,
    exp_time: float,
    dateobs: datetime | None,
    light: int,
    camera: Camera,
    focuser: Focuser,
    telescope: Telescope,
    site: Site,
    filter_band: Filters | str,
    airmass: float,
    n_star_limit: int,
    rng: numpy.random.Generator,
    timeout: float | None,
    sources: Sources | None,
    wcs: WCS | None,
    fwhm_multiplier: float = 5.0,
) -> np.ndarray:
    """Add stars and sky background to the base image."""
    if light == 1:
        if ra is None and dec is None and not isinstance(sources, Sources):
            raise ValueError("Either ra/dec or sources must be provided for light.")
        if ra is None or dec is None and isinstance(sources, Sources):
            ra, dec = sources.center  # type: ignore
        center = SkyCoord(ra=ra, dec=dec, unit="deg")

        assert camera.plate_scale is not None, "Camera plate scale must be set."
        radius = camera.get_fov_radius()

        if dateobs is None:
            dateobs = datetime.now(UTC)
        sources = get_sources(
            center=center,
            radius=radius,
            dateobs=dateobs,
            n_star_limit=n_star_limit,
            filter_band=filter_band,
            timeout=timeout,
            sources=sources,
        )
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
            wcs=wcs,
            dateobs=dateobs,
            airmass=airmass,
            fwhm_multiplier=fwhm_multiplier,
        )
    else:
        image = base
    return image


def generate_image(
    ra: float,
    dec: float,
    exp_time: float,
    dateobs: datetime | None = None,
    light: int = 1,
    camera: Camera = Camera(),
    focuser: Focuser = Focuser(),
    telescope: Telescope = Telescope(),
    site: Site = Site(),
    filter_band: Filters | str = Filters.G,
    airmass: float = 1.5,
    n_star_limit: int = 2000,
    rng: numpy.random.Generator = numpy.random.default_rng(),
    seed: int | None = None,
    timeout: float | None = None,
    sources: Sources | None = None,
    wcs: WCS | None = None,
    fwhm_multiplier: float = 5.0,
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
    airmass : float, optional
        Airmass at the image center (default: 1.5).
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
    wcs : WCS or None, optional
        World Coordinate System information for the image.
    fwhm_multiplier : float, optional
        Multiplier to determine the rendering radius around each star (default: 5.0).

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
            airmass=airmass,
            n_star_limit=n_star_limit,
            rng=rng,
            timeout=timeout,
            sources=sources,
            wcs=wcs,
            fwhm_multiplier=fwhm_multiplier,
        )
    else:
        image = base

    image = camera.apply_pixel_defects(image, exp_time)

    image = camera.bin_image(image)

    image = camera.to_adu_image(image)

    return image


def generate_image_stack(
    ra: float,
    dec: float,
    exp_time: float,
    dateobs: datetime | None = None,
    light: int = 1,
    camera: Camera = Camera(),
    focuser: Focuser = Focuser(),
    telescope: Telescope = Telescope(),
    site: Site = Site(),
    filter_band: Filters | str = Filters.G,
    airmass: float = 1.5,
    n_star_limit: int = 2000,
    rng: numpy.random.Generator = numpy.random.default_rng(),
    seed: int | None = None,
    timeout: float | None = None,
    sources: Sources | None = None,
    convert_all_to_adu: bool = False,
    wcs: WCS | None = None,
    fwhm_multiplier: float = 5.0,
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
    airmass : float, optional
        Airmass at the image center (default: 1.5).
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
        Multiplier to determine the rendering radius around each star (default: 5.0).

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
            airmass=airmass,
            n_star_limit=n_star_limit,
            rng=rng,
            timeout=timeout,
            sources=sources,
            wcs=wcs,
            fwhm_multiplier=fwhm_multiplier,
        )
    else:
        image = base

    adu_image = camera.apply_pixel_defects(image.copy(), exp_time)

    adu_image = camera.to_adu_image(adu_image)

    if convert_all_to_adu:
        base = camera.to_adu_image(base)
        image = camera.to_adu_image(image)

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
        ra=323.36152,
        dec=-0.82325,
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
