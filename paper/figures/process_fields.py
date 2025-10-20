import logging
from datetime import datetime
from glob import glob

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS, utils
from photutils.detection import DAOStarFinder
from prose import FITSImage, Sequence, blocks

import cabaret

logging.basicConfig(level=logging.INFO)

files = glob("./images/*.fits")


def open_fits_image(filename: str) -> tuple[np.ndarray, fits.Header]:
    """Open a FITS image and return the data and header."""
    with fits.open(filename) as hdul:
        image_data = hdul[0].data
        header = hdul[0].header
    return image_data, header


def get_seeing(filename: str) -> float:
    """Estimate the seeing from the FWHM in arcseconds."""
    ref = FITSImage(filename)

    calibration = Sequence(
        [
            blocks.PointSourceDetection(
                n=21, min_area=6, minor_length=3, min_separation=0, saturation=50000
            ),  # stars detection
            blocks.Cutouts(21),  # stars cutouts
            blocks.MedianEPSF(),  # building EPSF
            blocks.psf.Gaussian2D(),  # modeling EPSF
        ]
    )

    calibration.run(ref, show_progress=False)

    params = ref.epsf.params

    seeing = 0.348 * 2.355 * (params["sigma_x"] + params["sigma_y"]) / 2
    print(f"Seeing: {seeing:.2f} arcsec")
    return seeing


def get_image_wcs(header: fits.Header) -> WCS:
    """Get the WCS from the FITS header."""
    return WCS(header)


def get_plate_scale(wcs: WCS) -> np.ndarray:
    """Get the plate scale in arcseconds per pixel."""
    return utils.proj_plane_pixel_scales(wcs).mean() * 3600.0  # arcsec / pixel


def get_field_center(wcs: WCS, image_shape: tuple[int, int]) -> SkyCoord:
    """Get the field center in sky coordinates."""
    return utils.pixel_to_skycoord(image_shape[1] / 2, image_shape[0] / 2, wcs)


def get_sky_background(
    image_data: np.ndarray,
    gain: float,
    exptime: float,
    plate_scale: float,
    aptarea: float,
) -> float:
    """Estimate the sky background in photons / arcsec^2 / s."""
    return gain * np.percentile(image_data, 50) / (exptime * plate_scale**2 * aptarea)


def set_up_observatory(
    header: fits.Header, image_data: np.ndarray, seeing: float
) -> tuple[cabaret.Observatory, WCS, np.ndarray]:
    """Set up the observatory from the FITS header."""

    wcs = get_image_wcs(header)
    plate_scale = get_plate_scale(wcs)
    gain = header["GAIN"]
    exptime = header["EXPTIME"]
    collecting_area = header["APTAREA"] / 1e6  # in m^2
    diameter = header["APTDIA"] / 1000  # in m
    focal_length = header["FOCALLEN"] / 1000  # in m
    pitch = header["XPIXSZ"]

    # Paranal
    average_quantum_efficiency = 0.6
    sky_background = get_sky_background(
        image_data, gain, exptime, plate_scale, collecting_area
    )
    site = cabaret.Site(sky_background=sky_background, seeing=seeing)

    # Andor iKon-L 936
    camera = cabaret.Camera(
        read_noise=5.8,
        gain=gain,
        pitch=pitch,
        width=header["NAXIS1"],
        height=header["NAXIS2"],
        dark_current=0,
        average_quantum_efficiency=average_quantum_efficiency,
        bias=0,
        plate_scale=plate_scale.mean(),
        max_adu=(2**16 - 1),
    )

    # SPECULOOS
    telescope = cabaret.Telescope(
        diameter=diameter,
        focal_length=focal_length,
        collecting_area=collecting_area,
    )

    # Combine to make an observatory
    observatory = cabaret.Observatory(
        name="Paranal",
        camera=camera,
        telescope=telescope,
        site=site,
    )

    return observatory, wcs, plate_scale


def find_stars_dao(
    data: np.ndarray, threshold: float = 5.0, fwhm: float = 3.0
) -> np.ndarray:
    """
    Find stars using DAOStarFinder algorithm.

    Uses the photutils DAOStarFinder algorithm to detect point sources in
    astronomical images. The function performs background subtraction and
    returns star coordinates sorted by brightness.

    Parameters:
        data (np.ndarray): The 2D image data array.
        threshold (float, optional): Detection threshold in units of background
            standard deviation. Higher values detect fewer, brighter stars.
            Defaults to 5.0.
        fwhm (float, optional): Expected Full Width at Half Maximum of stars
            in pixels. Should match the typical seeing conditions. Defaults to 3.0.

    Returns:
        simulated_image (np.ndarray):
            The simulated image data array.
        real_image (np.ndarray):
            The real image data array.
        real_stars_dao (np.ndarray):
            Array of (x, y) coordinates of detected stars, sorted by brightness.
        simulated_stars_dao (np.ndarray):
            Array of (x, y) coordinates of detected stars, sorted by brightness.
    """
    # Calculate background statistics
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)

    # Use DAOStarFinder for star detection
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold * std, exclude_border=True)
    sources = daofind(data - median)

    if sources is None or len(sources) == 0:
        return np.array([]).reshape(0, 2)

    # Convert to (x, y) coordinates
    coordinates = np.column_stack([sources["xcentroid"], sources["ycentroid"]])

    # Sort by flux (brightness)
    fluxes = sources["flux"]
    return coordinates[np.argsort(fluxes)[::-1]]


def image_comparison(filename: str, threshold: float = 5.0):
    # Open the real image
    real_image, real_image_header = open_fits_image(filename)

    # Estimate the seeing
    seeing = get_seeing(filename)

    # Set up the observatory
    observatory, wcs, plate_scale = set_up_observatory(
        real_image_header, real_image, seeing
    )

    # Extract observation parameters
    real_center = get_field_center(wcs, real_image.shape)
    ra, dec = real_center.ra.deg, real_center.dec.deg  # in degrees
    exposure_time = real_image_header["EXPTIME"]  # in seconds
    dateobs = datetime.strptime(real_image_header["DATE-OBS"], "%Y-%m-%dT%H:%M:%S.%f")
    filter_band = cabaret.Filters.ensure_enum("RP")

    # Query Gaia catalog
    radius = observatory.camera.get_fov_radius()
    # print(f"Querying Gaia catalog around {real_center} with radius {radius:.2f}")
    # print(observatory)

    table = cabaret.GaiaQuery.query(
        center=real_center,
        radius=radius * np.sqrt(2),
        limit=100000,
        # timeout=60,
        filter_band=filter_band,
    )

    # Apply proper motion to the catalog
    table_filt = cabaret.GaiaQuery._apply_proper_motion(table, dateobs)
    fluxes = table_filt[filter_band.value].value.data  # type: ignore
    table_filt.remove_rows(np.isnan(fluxes))
    fluxes = fluxes[~np.isnan(fluxes)]

    sources = cabaret.Sources.from_arrays(
        ra=table_filt["ra"].value.data,  # type: ignore
        dec=table_filt["dec"].value.data,  # type: ignore
        fluxes=fluxes,
    )

    # Simulate the image
    simulated_image = observatory.generate_image(
        ra,
        dec,
        exposure_time,
        dateobs=dateobs,
        sources=sources,
        wcs=wcs,
        fwhm_multiplier=7,
    )

    # Find stars in both images
    real_stars_dao = find_stars_dao(
        real_image, threshold=threshold, fwhm=seeing / plate_scale
    )

    simulated_stars_dao = find_stars_dao(
        simulated_image, threshold=threshold, fwhm=seeing / plate_scale
    )

    # FITSManager.save(
    #     observatory,
    #     file_path="simulated_" + filename.split("/")[-1],
    #     # user_header=real_image_header,
    #     image=simulated_image,
    #     exp_time=exposure_time,
    #     ra=ra,
    #     dec=dec,
    # )

    return real_image, simulated_image, real_stars_dao, simulated_stars_dao


if __name__ == "__main__":
    from datetime import datetime

    import pandas as pd

    df_files = pd.DataFrame(columns=["filename", "real_stars", "simulated_stars"])
    threshold = 7
    for file in files:
        try:
            print(f"Processing {file}...")
            _r, _s, r_stars, s_stars = image_comparison(file, threshold=threshold)

            df_temp = pd.DataFrame(
                {
                    "filename": [file],
                    "real_stars": [len(r_stars)],
                    "simulated_stars": [len(s_stars)],
                }
            )

            df_files = pd.concat([df_files, df_temp], ignore_index=True)

            print(f"Detected {len(r_stars)} stars in real image.")
            print(f"Detected {len(s_stars)} stars in simulated image.")
            print()
        except Exception as e:
            print(f"Error {e}")

    df_files.to_csv(
        f"number_of_stars_comparison_threshold-{threshold}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
    )
