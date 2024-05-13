import numpy as np

from astropy.coordinates import SkyCoord
import astropy.units as u

from astropy.wcs import WCS

from datetime import datetime

from camera import Camera

from typing import Optional, Tuple, Union
from astropy.units import Quantity

def moffat_profile(x, y, x0, y0, FWHM, beta = 2.5):
    # https://nbviewer.org/github/ysbach/AO_2017/blob/master/04_Ground_Based_Concept.ipynb#1.2.-Moffat
    # FWHM =  2 * R * (2**(1/beta) - 1)**0.5

    R = (FWHM / 2) * (1 / (2**(1/beta) - 1)**0.5)
    A = (beta - 1) / (np.pi * R**2)

    r_squared = (x - x0)**2 + (y - y0)**2

    mp = A * (1 + (r_squared / R**2))**(-beta)

    mp_sum = np.sum(mp)

    return mp/mp_sum

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

def generate_image(ra, dec, exp_time, dateobs = datetime.utcnow(), light = 1, camera=Camera, tmass=True):

    base = np.ones((camera.numy, camera.numx)).astype(np.float64)
    base += camera.bias + np.random.poisson(base * camera.dark_current * exp_time).astype(np.float64)
    base += np.random.normal(0, camera.read_noise, (camera.numy, camera.numx)).astype(np.float64)
    
    if light == 1:

        # call gaia
        center = SkyCoord(ra=ra, dec=dec, unit="deg")

        fovx = (1/np.abs(np.cos(center.dec.rad)))*camera.numx * camera.plate_scale / 3600
        fovy = np.sqrt(2)*camera.numy * camera.plate_scale / 3600 # to account for poles, maybe should scale instead
        
        gaias, mags = gaia_radecs(center, (fovx*1.5, fovy*1.5), tmass=tmass, dateobs=dateobs)
 
        wcs = WCS(naxis=2)
        wcs.wcs.cdelt = [-camera.plate_scale / 3600, -camera.plate_scale / 3600]
        wcs.wcs.cunit = ["deg", "deg"]
        wcs.wcs.crpix = [int(camera.numx/2), int(camera.numy/2)]
        wcs.wcs.crval = [center.ra.deg, center.dec.deg]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

        if len(gaias) > 0:
            if tmass:
                # convert mags to fluxes. Reference: https://lweb.cfa.harvard.edu/~dfabricant/huchra/ay145/mags.html
                Jy = 1.51e7 # [photons sec^-1 m^-2 (dlambda/lambda)^-1]
                photons = 0.16 * 1600 * Jy # [photons sec^-1 m^-2] at mag 0
                fluxes = photons * 10**(-0.4*mags) * camera.quantum_efficiency * camera.collecting_area * exp_time # [electrons]
            else:
                fluxes = fluxes * camera.quantum_efficiency * camera.collecting_area * exp_time # [electrons]

            # convert gaia stars to pixel coordinates
            gaias_pixel = np.array(SkyCoord(gaias, unit="deg").to_pixel(wcs))

            # stars within frame and moffat profile
            stars = generate_star_image(gaias_pixel, fluxes, camera.seeing/camera.plate_scale, 
                                        (camera.numx, camera.numy)).astype(np.float64) # * flat
            
            # make base image with sky background
            image = base + np.random.poisson(base * camera.sky_background * exp_time).astype(np.float64) # * flat

            image += stars

    else:
        # dark exposure
        image = base


    # convert to adu and add camera's bias
    image = image / camera.gain + camera.bias # [adu]

    # clip to max adu
    image = np.clip(image, 0, camera.maxadu)

    # make image 16 bit
    image = image.astype(np.uint16)
    
    return image


def gaia_radecs(
    center: Union[Tuple[float, float], SkyCoord],
    fov: Union[float, Quantity],
    limit: int = 100000,
    circular: bool = True,
    tmass: bool = False,
    dateobs: Optional[datetime] = None,
) -> np.ndarray:
    """
    Query the Gaia archive to retrieve the RA-DEC coordinates of stars within a given field-of-view (FOV) centered on a given sky position.

    Parameters
    ----------
    center : tuple or astropy.coordinates.SkyCoord
        The sky coordinates of the center of the FOV. If a tuple is given, it should contain the RA and DEC in degrees.
    fov : float or astropy.units.Quantity
        The field-of-view of the FOV in degrees. If a float is given, it is assumed to be in degrees.
    limit : int, optional
        The maximum number of sources to retrieve from the Gaia archive. By default, it is set to 10000.
    circular : bool, optional
        Whether to perform a circular or a rectangular query. By default, it is set to True.
    tmass : bool, optional
        Whether to retrieve the 2MASS J magnitudes catelog. By default, it is set to False.
    dateobs : datetime.datetime, optional
        The date of the observation. If given, the proper motions of the sources will be taken into account. By default, it is set to None.

    Returns
    -------
    np.ndarray
        An array of shape (n, 2) containing the RA-DEC coordinates of the retrieved sources in degrees.

    Raises
    ------
    ImportError
        If the astroquery package is not installed.

    Examples
    --------
    >>> from astropy.coordinates import SkyCoord
    >>> from twirl import gaia_radecs
    >>> center = SkyCoord(ra=10.68458, dec=41.26917, unit='deg')
    >>> fov = 0.1
    >>> radecs = gaia_radecs(center, fov)
    """
    from astroquery.gaia import Gaia

    if isinstance(center, SkyCoord):
        ra = center.ra.deg
        dec = center.dec.deg
    else:
        ra, dec = center

    if not isinstance(fov, u.Quantity):
        fov = fov * u.deg

    if fov.ndim == 1:
        ra_fov, dec_fov = fov.to(u.deg).value
    else:
        ra_fov = dec_fov = fov.to(u.deg).value

    radius = np.max([ra_fov, dec_fov]) / 2

    if circular and not tmass:
        job = Gaia.launch_job(
            f"""
            SELECT top {limit} gaia.ra, gaia.dec, gaia.pmra, gaia.pmdec, gaia.phot_rp_mean_flux
            FROM gaiadr2.gaia_source AS gaia
            WHERE 1=CONTAINS(
                POINT('ICRS', {ra}, {dec}), 
                CIRCLE('ICRS', gaia.ra, gaia.dec, {radius}))
            ORDER BY gaia.phot_rp_mean_flux DESC
            """
        )
    elif circular and tmass:
        job = Gaia.launch_job(
            f"""
            SELECT top {limit} gaia.ra, gaia.dec, gaia.pmra, gaia.pmdec, gaia.phot_rp_mean_flux, tmass.j_m
            FROM gaiadr2.gaia_source AS gaia
            INNER JOIN gaiadr2.tmass_best_neighbour AS tmass_match ON tmass_match.source_id = gaia.source_id
            INNER JOIN gaiadr1.tmass_original_valid AS tmass ON tmass.tmass_oid = tmass_match.tmass_oid
            WHERE 1=CONTAINS(
                POINT('ICRS', {ra}, {dec}), 
                CIRCLE('ICRS', gaia.ra, gaia.dec, {radius}))
            ORDER BY tmass.j_m
            """
        )
    elif not circular and tmass:
        job = Gaia.launch_job(
            f"""
            SELECT top {limit} gaia.ra, gaia.dec, gaia.pmra, gaia.pmdec, gaia.phot_rp_mean_flux, tmass.j_m
            FROM gaiadr2.gaia_source AS gaia
            INNER JOIN gaiadr2.tmass_best_neighbour AS tmass_match ON tmass_match.source_id = gaia.source_id
            INNER JOIN gaiadr1.tmass_original_valid AS tmass ON tmass.tmass_oid = tmass_match.tmass_oid
            WHERE gaia.ra BETWEEN {ra-ra_fov/2} AND {ra+ra_fov/2} AND
            gaia.dec BETWEEN {dec-dec_fov/2} AND {dec+dec_fov/2}
            ORDER BY tmass.j_m
            """
        )
    else:
        job = Gaia.launch_job(
            f"""
            SELECT top {limit} gaia.ra, gaia.dec, gaia.pmra, gaia.pmdec, gaia.phot_rp_mean_flux
            FROM gaiadr2.gaia_source AS gaia
            WHERE gaia.ra BETWEEN {ra-ra_fov/2} AND {ra+ra_fov/2} AND
            gaia.dec BETWEEN {dec-dec_fov/2} AND {dec+dec_fov/2}
            ORDER BY gaia.phot_rp_mean_flux DESC
            """
        )

    table = job.get_results()
    
    # add proper motion to ra and dec
    if dateobs is not None:
        # calculate fractional year
        dateobs = dateobs.year + (dateobs.timetuple().tm_yday - 1) / 365.25 # type: ignore
        
        years = dateobs - 2015.5 # type: ignore
        table["ra"] += years * table["pmra"] / 1000 / 3600
        table["dec"] += years * table["pmdec"] / 1000 / 3600

    
    if tmass:
        table.remove_rows(np.isnan(table['j_m']))
        return np.array([table["ra"].value.data, table["dec"].value.data]).T, table['j_m'].value.data
    else:
        table.remove_rows(np.isnan(table['phot_rp_mean_flux']))
        return np.array([table["ra"].value.data, table["dec"].value.data]).T, table['phot_rp_mean_flux'].value.data
    


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
    # example usage
    image = generate_image(323.36152, -0.82325, 1)

    fig = plt.figure(figsize=(10, 10))

    med = np.median(image)
    std = np.std(image)
    plt.imshow(image, cmap="Greys_r", vmax=3*std + med, vmin=med - 1*std)
    plt.show()
    
