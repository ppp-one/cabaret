import logging
from datetime import UTC, datetime

import numpy as np
import numpy.random
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
from astropy.time import Time

from cabaret.camera import Camera
from cabaret.focuser import Focuser
from cabaret.queries import get_gaia_sources
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
    fluxes: list,
    FWHM: float,
    frame_size: tuple,
    rng: numpy.random.Generator,
):
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
    fluxes: list,
    FWHM: float,
    frame_size: tuple,
    rng: numpy.random.Generator,
) -> np.ndarray:
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
    tmass: bool = True,
    n_star_limit: int = 2000,
    rng: numpy.random.Generator = numpy.random.default_rng(),
    seed: int | None = None,
    timeout: float | None = None,
    sources: Sources | None = None,
    apply_pixel_defects: bool = True,
) -> np.ndarray:
    if seed is not None:
        rng = numpy.random.default_rng(seed)

    base = np.ones((camera.height, camera.width)).astype(np.float64)

    base += rng.poisson(base * camera.dark_current * exp_time).astype(np.float64)

    base += rng.normal(0, camera.read_noise, (camera.height, camera.width)).astype(
        np.float64
    )

    if camera.plate_scale is None:
        camera.plate_scale = (
            2
            * np.arctan((camera.pitch * 1e-6) / (2 * telescope.focal_length))
            * (180 / np.pi)
            * 3600
        )  # "/pixel

    if telescope.collecting_area is None:
        telescope.collecting_area = np.pi * (telescope.diameter / 2) ** 2  # [m^2]

    if light == 1:
        # call gaia
        center = SkyCoord(ra=ra, dec=dec, unit="deg")

        fovx = (
            (1 / np.abs(np.cos(center.dec.rad)))
            * camera.width
            * camera.plate_scale
            / 3600
        )
        fovy = (
            np.sqrt(2) * camera.height * camera.plate_scale / 3600
        )  # to account for poles, maybe should scale instead

        logger.info("Querying Gaia for sources...")
        # gaias, vals = get_gaia_sources(
        if not isinstance(sources, Sources):
            sources = get_gaia_sources(
                center,
                (fovx * 1.5, fovy * 1.5),
                tmass=tmass,
                dateobs=dateobs,
                limit=n_star_limit,
                timeout=timeout,
            )
        logger.info(f"Found {len(sources)} sources (user set limit of {n_star_limit}).")

        image = base

        if site.latitude is not None and site.longitude is not None:
            logger.info(
                "Since location specified, calculating sunlight brightness"
                " based on sun's position"
            )
            location = EarthLocation(
                lat=site.latitude * u.deg, lon=site.longitude * u.deg
            )
            obs_time = Time(dateobs, scale="utc")

            # Get sun position
            sun = get_sun(obs_time)

            # Transform to altitude/azimuth coordinates for the observatory
            altaz_frame = AltAz(obstime=obs_time, location=location)
            sun_altaz = sun.transform_to(altaz_frame)
            sun_altitude = sun_altaz.alt.degree
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

            image += rng.poisson(
                np.ones((camera.height, camera.width)).astype(np.float64)
                * sky_e
                * exp_time
            ).astype(np.float64)

        if len(sources) > 0:
            fluxes = (
                sources.fluxes
                * camera.average_quantum_efficiency
                * telescope.collecting_area
                * exp_time
            )  # [electrons]

            # convert gaia stars to pixel coordinates
            wcs = camera.get_wcs(center)
            gaias_pixel = sources.to_pixel(wcs)

            # stars within frame and moffat profile
            stars = generate_star_image(
                gaias_pixel,
                fluxes,
                focuser.seeing_multiplier * site.seeing / camera.plate_scale,
                (camera.width, camera.height),
                rng=rng,
            ).astype(np.float64)  # * flat

            sky_background = (
                site.sky_background * telescope.collecting_area * camera.plate_scale**2
            )  # [e-/s]

            # make base image with sky background
            image = base + rng.poisson(
                np.ones((camera.height, camera.width)).astype(np.float64)
                * sky_background
                * exp_time
            ).astype(np.float64)  # * flat

            image += stars

    else:
        # dark exposure
        image = base

    # inject defect pixels
    if apply_pixel_defects:
        for defect in camera.pixel_defects.values():
            image = defect.introduce_pixel_defect(image, camera)

    # convert to adu and add camera's bias
    image = image / camera.gain + camera.bias  # [adu]

    # clip to max adu
    image = np.clip(image, 0, camera.max_adu)

    # make image 16 bit
    image = image.astype(np.uint16)

    return image


if __name__ == "__main__":
    import matplotlib.pyplot as plt

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

    print("Plotting image...")
    med = np.median(science)
    std = np.std(science)
    print(med, std)

    fig, ax = plt.subplots()
    img = ax.imshow(science, cmap="gray", vmin=med - 1 * std, vmax=med + 1 * std)
    cbar = plt.colorbar(img, ax=ax)
    cbar.set_label("ADU")
    plt.show()
