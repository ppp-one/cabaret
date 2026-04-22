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
from cabaret.queries import Filters, GaiaQuery, GaiaSQLiteSource, GaiaTAPSource
from cabaret.site import Site
from cabaret.sources import Sources, _normalize_rate
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


def _sample_trail(
    x0: float,
    y0: float,
    dx_pixel: float,
    dy_pixel: float,
    n_samples: int,
    jitter_sigma_pixels: float,
    rng: numpy.random.Generator,
) -> np.ndarray:
    """Sample n_samples positions uniformly along a linear drift trail.

    The trail runs from (x0, y0) at the start of the exposure to
    (x0 + dx_pixel, y0 + dy_pixel) at the end.

    To implement non-linear (ephemeris-based) motion, provide a callable with
    this same signature as the ``trail_sampler`` argument of
    ``generate_star_image``.

    Parameters
    ----------
    x0, y0 : float
        Pixel position at the start of the exposure (t=0).
    dx_pixel, dy_pixel : float
        Total pixel displacement over the full exposure duration.
    n_samples : int
        Number of sample positions along the trail.
    jitter_sigma_pixels : float
        1-sigma Gaussian jitter applied independently to each sample (pixels).
    rng : numpy.random.Generator
        Random number generator for jitter.

    Returns
    -------
    np.ndarray
        Shape (2, n_samples): row 0 is x positions, row 1 is y positions.
    """
    t = np.linspace(0.0, 1.0, n_samples)
    xs = x0 + t * dx_pixel
    ys = y0 + t * dy_pixel
    if jitter_sigma_pixels > 0.0:
        xs += rng.normal(0.0, jitter_sigma_pixels, n_samples)
        ys += rng.normal(0.0, jitter_sigma_pixels, n_samples)
    return np.stack([xs, ys])


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
    drift_pixels: np.ndarray | None = None,
    n_trail_samples: int = 50,
    jitter_sigma_pixels: float = 0.0,
    trail_sampler=None,
) -> np.ndarray:
    """
    Render stars onto an image using a fast, windowed approach.

    Parameters
    ----------
    pos : np.ndarray
        Pixel positions of stars at the start of the exposure (shape: 2, n_stars).
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
    drift_pixels : np.ndarray, optional
        Total pixel displacement over the exposure for each star (shape: 2, n_stars).
        Row 0 is x-drift, row 1 is y-drift. None means no drift (point PSF).
    n_trail_samples : int, optional
        Maximum number of Moffat samples along the drift trail (default: 50).
        The actual count is adaptive: ``max(1, ceil(trail_length_px / FWHM * 2))``,
        capped at this value.
    jitter_sigma_pixels : float, optional
        1-sigma guiding jitter added independently to each trail sample (pixels).
    trail_sampler : callable, optional
        Function with the same signature as ``_sample_trail`` for non-linear motion.
        Defaults to ``_sample_trail`` (linear interpolation).

    Returns
    -------
    np.ndarray
        Image with rendered stars.
    """
    # Validate inputs
    if not isinstance(n_trail_samples, int) or n_trail_samples < 1:
        raise ValueError(
            f"n_trail_samples must be an integer >= 1, got "
            f"{n_trail_samples} ({type(n_trail_samples).__name__})"
        )

    if jitter_sigma_pixels < 0:
        raise ValueError(f"jitter_sigma_pixels must be >= 0, got {jitter_sigma_pixels}")

    n_stars = len(fluxes)

    if drift_pixels is None:
        _drift = np.zeros((2, n_stars), dtype=np.float64)
    else:
        _drift = np.asarray(drift_pixels, dtype=np.float64)
        if _drift.ndim != 2 or _drift.shape != (2, n_stars):
            raise ValueError(
                "drift_pixels must have shape (2, n_stars) "
                f"where n_stars={n_stars}; got shape {_drift.shape}"
            )
    has_drift = (_drift[0] != 0.0) | (_drift[1] != 0.0)
    _sampler = trail_sampler if trail_sampler is not None else _sample_trail

    x = np.linspace(0, frame_size[0] - 1, frame_size[0])
    y = np.linspace(0, frame_size[1] - 1, frame_size[1])
    xx, yy = np.meshgrid(x, y)

    render_radius = FWHM * fwhm_multiplier  # render n * FWHM around the star

    W, H = frame_size
    image = np.zeros(frame_size).T
    for i, flux in enumerate(fluxes):
        x0 = pos[0][i]
        y0 = pos[1][i]
        dx = _drift[0, i]
        dy = _drift[1, i]
        x_end = x0 + dx
        y_end = y0 + dy

        # Skip if the entire trail bounding box is outside the frame
        if (
            max(x0, x_end) < 0
            or min(x0, x_end) >= W
            or max(y0, y_end) < 0
            or min(y0, y_end) >= H
        ):
            continue

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

        if has_drift[i]:
            trail_length = np.hypot(dx, dy)
            n_samples = max(
                1, min(n_trail_samples, int(np.ceil(trail_length / FWHM * 2)))
            )

            sample_pos = _sampler(x0, y0, dx, dy, n_samples, jitter_sigma_pixels, rng)
            sample_flux = star_flux / sample_pos.shape[1]
            for j in range(n_samples):
                sx, sy = sample_pos[0, j], sample_pos[1, j]
                sx_min = max(0, int(sx - render_radius))
                sx_max = min(int(sx + render_radius), W - 1)
                sy_min = max(0, int(sy - render_radius))
                sy_max = min(int(sy + render_radius), H - 1)
                if sx_max < 0 or sx_min >= W or sy_max < 0 or sy_min >= H:
                    continue
                image[sy_min : sy_max + 1, sx_min : sx_max + 1] += (
                    sample_flux
                    * moffat_profile(
                        xx[sy_min : sy_max + 1, sx_min : sx_max + 1],
                        yy[sy_min : sy_max + 1, sx_min : sx_max + 1],
                        sx,
                        sy,
                        FWHM,
                    )
                )
        else:
            x_min = max(0, int(x0 - render_radius))
            x_max = min(int(x0 + render_radius), W - 1)
            y_min = max(0, int(y0 - render_radius))
            y_max = min(int(y0 + render_radius), H - 1)

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
    center: SkyCoord | None,
    radius: Angle | None,
    dateobs: datetime,
    n_star_limit: int,
    filter_band: Filters | str,
    timeout: float | None,
    sources: Sources | None = None,
    tap_source: GaiaTAPSource | GaiaSQLiteSource | str | None = None,
    bounds: tuple[float, float, float, float] | None = None,
    additional_sources: Sources | None = None,
) -> Sources:
    """Get sources from Gaia or use provided sources."""
    if not isinstance(sources, Sources):
        logger.info("Querying Gaia for sources...")
        sources = GaiaQuery.get_sources(
            center=center,
            radius=radius,
            bounds=bounds,
            dateobs=dateobs,
            limit=n_star_limit,
            filter_band=filter_band,
            timeout=timeout,
            tap_source=tap_source,
        )
        logger.info(f"Found {len(sources)} sources (user set limit of {n_star_limit}).")
    if additional_sources is not None:
        sources = sources + additional_sources
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

    # Sky brightness model due to sunlight scattering
    # calibrated for I+z band in Paranal
    SKY_MODEL_A = np.float64(4533508.655833181)
    SKY_MODEL_B = np.float64(0.3937301435229289)
    SKY_MODEL_C = np.float64(-0.7907223506021084)

    sun_altitude = _get_sun_altitude(site, dateobs, logger)

    if sun_altitude is None:
        return image

    sun_brightness = SKY_MODEL_A * SKY_MODEL_B ** (SKY_MODEL_C * sun_altitude)
    logger.debug(f"sun_brightness (ph/m2/arcsec2/s): {sun_brightness}")

    sun_electrons_per_sec = (
        sun_brightness
        * telescope.collecting_area
        * camera.average_quantum_efficiency
        * camera.plate_scale**2
    )
    logger.debug(f"sun_e (e-/s): {sun_electrons_per_sec}")

    sun_electrons = sun_electrons_per_sec * exp_time
    sun_background = np.random.poisson(
        np.ones((camera.height, camera.width), dtype=np.float64) * sun_electrons
    ).astype(np.float64)

    image += sun_background

    return image


def _get_sun_altitude(
    site: Site,
    dateobs: datetime,
    logger: logging.Logger,
) -> float | None:
    """
    Calculate the sun's altitude angle at the observation time.

    Parameters
    ----------
    site : Site
        Observatory site configuration.
    dateobs : datetime
        Observation date and time.
    logger : logging.Logger
        Logger instance for debug messages.

    Returns
    -------
    float or None
        Sun altitude in degrees, or None if it cannot be calculated.
    """
    if site.sun_altitude is not None:
        return site.sun_altitude

    if site.latitude is None or site.longitude is None:
        return None

    logger.debug("Calculating sunlight brightness based on sun's position")

    location = EarthLocation(
        lat=u.Quantity(site.latitude, "deg"),
        lon=u.Quantity(site.longitude, "deg"),
    )
    obs_time = Time(dateobs, scale="utc")
    sun = get_sun(obs_time)
    altaz_frame = AltAz(obstime=obs_time, location=location)
    sun_altaz = sun.transform_to(altaz_frame)
    sun_altitude: float = sun_altaz.alt.degree  # type: ignore

    logger.debug(f"Sun altitude: {sun_altitude:.5f} deg at {dateobs} UTC")

    return sun_altitude


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
    tracking_ra_rate: float | u.Quantity = 0.0,
    tracking_dec_rate: float | u.Quantity = 0.0,
    n_trail_samples: int = 50,
    jitter_sigma: float = 0.0,
) -> np.ndarray:
    """Add stars to the image using the Moffat profile and sky background.

    Parameters
    ----------
    tracking_ra_rate : float or astropy Quantity, optional
        Telescope RA tracking rate offset from sidereal in on-sky arcsec/s,
        i.e. dα·cos(δ)/dt (default: 0). Accepts any astropy angular velocity
        Quantity (e.g. ``1 * u.arcsec / u.s``).
    tracking_dec_rate : float or astropy Quantity, optional
        Telescope Dec tracking rate offset from sidereal in arcsec/s (default: 0).
        Accepts any astropy angular velocity Quantity.
    n_trail_samples : int, optional
        Number of PSF samples along the drift trail for moving sources (default: 50).
    jitter_sigma : float, optional
        1-sigma guiding jitter in arcsec applied per trail sample (default: 0).
    """
    tracking_ra_rate = _normalize_rate(tracking_ra_rate)
    tracking_dec_rate = _normalize_rate(tracking_dec_rate)
    sources = sources.drop_nan_fluxes()
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

        start_pixel = sources.to_pixel(wcs)

        # rates are on-sky arcsec/s (dα·cosδ/dt), so convert to RA degrees via cosδ
        rel_ra_arcsec = (sources.ra_rates - tracking_ra_rate) * exp_time
        rel_dec_arcsec = (sources.dec_rates - tracking_dec_rate) * exp_time
        cos_dec = np.cos(np.deg2rad(sources.coords.dec.deg))
        end_coords = SkyCoord(
            ra=sources.coords.ra.deg + rel_ra_arcsec / (3600.0 * cos_dec),
            dec=sources.coords.dec.deg + rel_dec_arcsec / 3600.0,
            unit="deg",
        )
        end_pixel = np.array(end_coords.to_pixel(wcs))
        drift_pixels = end_pixel - start_pixel

        on_camera_mask = (
            (np.maximum(start_pixel[0], end_pixel[0]) >= 0)
            & (np.minimum(start_pixel[0], end_pixel[0]) < camera.width)
            & (np.maximum(start_pixel[1], end_pixel[1]) >= 0)
            & (np.minimum(start_pixel[1], end_pixel[1]) < camera.height)
        )
        logger.debug(
            f"Number of stars with trails intersecting camera: {on_camera_mask.sum()}"
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
            start_pixel[:, on_camera_mask],
            fluxes[on_camera_mask],
            focuser.seeing_multiplier * site.seeing / camera.plate_scale,
            (camera.width, camera.height),
            rng=rng,
            telescope_aperture=telescope.diameter,
            site_elevation=site.elevation if site.elevation is not None else 0,
            exp_time=exp_time,
            airmass=airmass,
            fwhm_multiplier=fwhm_multiplier,
            drift_pixels=drift_pixels[:, on_camera_mask],
            n_trail_samples=n_trail_samples,
            jitter_sigma_pixels=jitter_sigma / camera.plate_scale,
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
    tap_source: GaiaTAPSource | GaiaSQLiteSource | str | None = None,
    tracking_ra_rate: float | u.Quantity = 0.0,
    tracking_dec_rate: float | u.Quantity = 0.0,
    n_trail_samples: int = 50,
    jitter_sigma: float = 0.0,
    additional_sources: Sources | None = None,
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
            tap_source=tap_source,
            additional_sources=additional_sources,
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
            tracking_ra_rate=tracking_ra_rate,
            tracking_dec_rate=tracking_dec_rate,
            n_trail_samples=n_trail_samples,
            jitter_sigma=jitter_sigma,
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
    tap_source: GaiaTAPSource | GaiaSQLiteSource | str | None = None,
    tracking_ra_rate: float | u.Quantity = 0.0,
    tracking_dec_rate: float | u.Quantity = 0.0,
    n_trail_samples: int = 50,
    jitter_sigma: float = 0.0,
    additional_sources: Sources | None = None,
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
            tap_source=tap_source,
            tracking_ra_rate=tracking_ra_rate,
            tracking_dec_rate=tracking_dec_rate,
            n_trail_samples=n_trail_samples,
            jitter_sigma=jitter_sigma,
            additional_sources=additional_sources,
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
    tap_source: GaiaTAPSource | GaiaSQLiteSource | str | None = None,
    tracking_ra_rate: float | u.Quantity = 0.0,
    tracking_dec_rate: float | u.Quantity = 0.0,
    n_trail_samples: int = 50,
    jitter_sigma: float = 0.0,
    additional_sources: Sources | None = None,
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
            tap_source=tap_source,
            tracking_ra_rate=tracking_ra_rate,
            tracking_dec_rate=tracking_dec_rate,
            n_trail_samples=n_trail_samples,
            jitter_sigma=jitter_sigma,
            additional_sources=additional_sources,
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
