import numpy as np

from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

from datetime import datetime

from cabaret.camera import Camera
from cabaret.queries import gaia_radecs


def moffat_profile(x, y, x0, y0, FWHM, beta=2.5):
    # https://nbviewer.org/github/ysbach/AO_2017/blob/master/04_Ground_Based_Concept.ipynb#1.2.-Moffat
    # FWHM =  2 * R * (2**(1/beta) - 1)**0.5

    R = (FWHM / 2) * (1 / (2 ** (1 / beta) - 1) ** 0.5)
    A = (beta - 1) / (np.pi * R**2)

    r_squared = (x - x0) ** 2 + (y - y0) ** 2

    mp = A * (1 + (r_squared / R**2)) ** (-beta)

    mp_sum = np.sum(mp)

    return mp / mp_sum


def generate_star_image(pos, fluxes, FWHM, frame_size):
    x = np.linspace(0, frame_size[0] - 1, frame_size[0])
    y = np.linspace(0, frame_size[1] - 1, frame_size[1])
    xx, yy = np.meshgrid(x, y)

    image = np.zeros(frame_size).T
    for i, flux in enumerate(fluxes):
        x0 = pos[0][i]
        y0 = pos[1][i]
        star = np.random.poisson(flux) * moffat_profile(xx, yy, x0, y0, FWHM)
        image += star

    return image


def generate_image(
    ra, dec, exp_time, dateobs=datetime.utcnow(), light=1, camera=Camera, tmass=True
):

    base = np.ones((camera.height, camera.width)).astype(np.float64)
    base += camera.bias + np.random.poisson(
        base * camera.dark_current * exp_time
    ).astype(np.float64)
    base += np.random.normal(
        0, camera.read_noise, (camera.height, camera.width)
    ).astype(np.float64)

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

        gaias, mags = gaia_radecs(
            center, (fovx * 1.5, fovy * 1.5), tmass=tmass, dateobs=dateobs
        )

        wcs = WCS(naxis=2)
        wcs.wcs.cdelt = [-camera.plate_scale / 3600, -camera.plate_scale / 3600]
        wcs.wcs.cunit = ["deg", "deg"]
        wcs.wcs.crpix = [int(camera.width / 2), int(camera.height / 2)]
        wcs.wcs.crval = [center.ra.deg, center.dec.deg]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

        if len(gaias) > 0:
            if tmass:
                # convert mags to fluxes. Reference: https://lweb.cfa.harvard.edu/~dfabricant/huchra/ay145/mags.html
                Jy = 1.51e7  # [photons sec^-1 m^-2 (dlambda/lambda)^-1]
                photons = 0.16 * 1600 * Jy  # [photons sec^-1 m^-2] at mag 0
                fluxes = (
                    photons
                    * 10 ** (-0.4 * mags)
                    * camera.quantum_efficiency
                    * camera.collecting_area
                    * exp_time
                )  # [electrons]
            else:
                fluxes = (
                    fluxes
                    * camera.quantum_efficiency
                    * camera.collecting_area
                    * exp_time
                )  # [electrons]

            # convert gaia stars to pixel coordinates
            gaias_pixel = np.array(SkyCoord(gaias, unit="deg").to_pixel(wcs))

            # stars within frame and moffat profile
            stars = generate_star_image(
                gaias_pixel,
                fluxes,
                camera.seeing / camera.plate_scale,
                (camera.width, camera.height),
            ).astype(
                np.float64
            )  # * flat

            # make base image with sky background
            image = base + np.random.poisson(
                base * camera.sky_background * exp_time
            ).astype(
                np.float64
            )  # * flat

            image += stars

    else:
        # dark exposure
        image = base

    # convert to adu and add camera's bias
    image = image / camera.gain + camera.bias  # [adu]

    # clip to max adu
    image = np.clip(image, 0, camera.max_adu)

    # make image 16 bit
    image = image.astype(np.uint16)

    return image


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # example usage
    image = generate_image(323.36152, -0.82325, 1)

    fig = plt.figure(figsize=(10, 10))

    med = np.median(image)
    std = np.std(image)
    plt.imshow(image, cmap="Greys_r", vmax=3 * std + med, vmin=med - 1 * std)
    plt.show()
