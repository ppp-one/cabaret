from concurrent.futures import ThreadPoolExecutor, TimeoutError
from datetime import datetime

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.units import Quantity
from astroquery.gaia import Gaia

from cabaret.sources import Sources


def tmass_mag_to_photons(mags: np.ndarray) -> np.ndarray:
    """Convert 2MASS J magnitudes to photon fluxes at mag 0.

    Reference: https://lweb.cfa.harvard.edu/~dfabricant/huchra/ay145/mags.html
    Returns photons/sec/m^2 for each magnitude.
    """
    Jy = 1.51e7  # [photons sec^-1 m^-2 (dlambda/lambda)^-1]
    photons = 0.16 * 1600 * Jy  # [photons sec^-1 m^-2] at mag 0
    return photons * 10 ** (-0.4 * mags)


def gaia_launch_job_with_timeout(query, timeout=None, **kwargs) -> Table:
    """
    Launch a Gaia job and return its results, optionally enforcing a timeout.

    Parameters
    ----------
    query : str
        The query string passed to Gaia.launch_job.
    timeout : float or None, optional
        Maximum number of seconds to wait for Gaia.launch_job to complete.
        If None, the job is run on the main thread (no thread overhead).
    **kwargs
        Additional keyword arguments forwarded to Gaia.launch_job.

    Returns
    -------
    object
        The result returned by job.get_results().

    Raises
    ------
    TimeoutError
        If `timeout` is not None and the call does not complete within `timeout`.
    """
    # Run directly on the main thread when no timeout is requested to avoid
    # unnecessary thread creation and to preserve original callstacks/tracebacks.
    if timeout is None:
        job = Gaia.launch_job(query, **kwargs)
        return job.get_results()

    with ThreadPoolExecutor() as executor:
        future = executor.submit(Gaia.launch_job, query, **kwargs)
        try:
            job = future.result(timeout=timeout)
            return job.get_results()
        except TimeoutError:
            raise TimeoutError("Gaia query timed out.")


def gaia_query(
    center: tuple[float, float] | SkyCoord,
    fov: float | Quantity,
    limit: int = 100000,
    circular: bool = True,
    tmass: bool = False,
    timeout: float | None = None,
) -> Table:
    """Query Gaia and return the raw Astropy Table.

    Example
    -------
    >>> from cabaret.queries import gaia_table_query
    >>> from astropy.coordinates import SkyCoord
    >>> center = SkyCoord(ra=10.68458, dec=41.26917, unit='deg')
    >>> table = gaia_table_query(center, fov=0.1, limit=10, timeout=30)
    """
    if isinstance(center, SkyCoord):
        ra = center.ra.deg
        dec = center.dec.deg
    else:
        ra, dec = center

    if not isinstance(fov, u.Quantity):
        fov = u.Quantity(fov, u.deg)

    if fov.ndim == 1:
        ra_fov, dec_fov = fov.to(u.deg).value
    else:
        ra_fov = dec_fov = fov.to(u.deg).value

    radius = np.max([ra_fov, dec_fov]) / 2

    select_cols = [
        "gaia.ra",
        "gaia.dec",
        "gaia.pmra",
        "gaia.pmdec",
        "gaia.phot_rp_mean_flux",
    ]
    joins = []
    where = []
    order_by = "gaia.phot_rp_mean_flux DESC"

    if tmass:
        select_cols.append("tmass.j_m")
        joins.extend(
            [
                "INNER JOIN gaiadr2.tmass_best_neighbour AS tmass_match "
                + "ON tmass_match.source_id = gaia.source_id",
                "INNER JOIN gaiadr1.tmass_original_valid AS tmass "
                + "ON tmass.tmass_oid = tmass_match.tmass_oid",
            ]
        )
        order_by = "tmass.j_m"

    if circular:
        where.append(
            f"1=CONTAINS(POINT('ICRS', {ra}, {dec}), "
            f"CIRCLE('ICRS', gaia.ra, gaia.dec, {radius}))"
        )
    else:
        where.append(
            f"gaia.ra BETWEEN {ra - ra_fov / 2} AND {ra + ra_fov / 2} "
            f"AND gaia.dec BETWEEN {dec - dec_fov / 2} AND {dec + dec_fov / 2}"
        )

    select_clause = ", ".join(select_cols)
    joins_clause = "\n".join(joins)
    where_clause = " AND ".join(where)

    query = f"""
    SELECT TOP {limit} {select_clause}
    FROM gaiadr2.gaia_source AS gaia
    {joins_clause}
    WHERE {where_clause}
    ORDER BY {order_by}
    """

    table = gaia_launch_job_with_timeout(query, timeout=timeout)
    return table


def apply_proper_motion(table: Table, dateobs: datetime):
    """
    Apply proper motion correction to RA and DEC columns for the given observation date.
    """
    dateobs_frac = dateobs.year + (dateobs.timetuple().tm_yday - 1) / 365.25  # type: ignore
    years = dateobs_frac - 2015.5  # type: ignore
    table["ra"] += years * table["pmra"] / 1000 / 3600
    table["dec"] += years * table["pmdec"] / 1000 / 3600
    return table


def get_gaia_sources(
    center: tuple[float, float] | SkyCoord,
    fov: float | Quantity,
    limit: int = 100000,
    circular: bool = True,
    tmass: bool = False,
    dateobs: datetime | None = None,
    timeout: float | None = None,
) -> Sources:
    """
    Query the Gaia archive to retrieve the RA-DEC coordinates of stars
    within a given field-of-view (FOV) centered on a given sky position.

    Parameters
    ----------
    center : tuple or astropy.coordinates.SkyCoord
        The sky coordinates of the center of the FOV.
        If a tuple is given, it should contain the RA and DEC in degrees.
    fov : float or astropy.units.Quantity
        The field-of-view of the FOV in degrees. If a float is given,
        it is assumed to be in degrees.
    limit : int, optional
        The maximum number of sources to retrieve from the Gaia archive.
        By default, it is set to 10000.
    circular : bool, optional
        Whether to perform a circular or a rectangular query.
        By default, it is set to True.
    tmass : bool, optional
        Whether to retrieve the 2MASS J magnitudes catelog.
        By default, it is set to False.
    dateobs : datetime.datetime, optional
        The date of the observation. If given, the proper motions of the sources
        will be taken into account. By default, it is set to None.
    timeout : float, optional
        The maximum time to wait for the Gaia query to complete, in seconds.
        If None, there is no timeout. By default, it is set to None.

    Returns
    -------
    Sources
        A Sources instance containing the coordinates and fluxes of the retrieved
        sources.


    Notes
    -----
    If `tmass` is True, the fluxes are calculated from the 2MASS J-band magnitudes, see
    `cabaret.queries.tmass_mag_to_photons`.

    Raises
    ------
    ImportError
        If the astroquery package is not installed.

    Examples
    --------
    >>> from cabaret.queries import get_gaia_sources
    >>> from astropy.coordinates import SkyCoord
    >>> center = SkyCoord(ra=10.68458, dec=41.26917, unit='deg')
    >>> sources = get_gaia_sources(center, fov=0.1, timeout=30)
    Sources(coords=<SkyCoord (ICRS): (ra, dec) in deg
        [(10.63950247, 41.26393165), (10.6880729 , 41.22524785),
        (10.70349581, 41.25357386), (10.70022208, 41.26019689),
        (10.71333998, 41.29943347), (10.73974676, 41.2942209 ),
        (10.71181048, 41.29130279), (10.68780207, 41.31717482),
        (10.63804045, 41.27468757), (10.64397532, 41.25237352)]>, fluxes=array(
            [169435.62814443,  52203.9396396 ,  41716.18126449,  29035.89106422,
            22990.85994301,  17672.53437883,  15953.21022642,  15077.12262318,
            14004.42013396,  12271.11779953]))

    """
    table = gaia_query(
        center=center,
        fov=fov,
        limit=limit,
        circular=circular,
        tmass=tmass,
        timeout=timeout,
    )

    if dateobs is not None:
        table = apply_proper_motion(table, dateobs)

    if tmass:
        fluxes = tmass_mag_to_photons(table["j_m"].value.data)
    else:
        fluxes = table["phot_rp_mean_flux"].value.data
    table.remove_rows(np.isnan(fluxes))
    fluxes = fluxes[~np.isnan(fluxes)]

    return Sources.from_arrays(
        ra=table["ra"].value.data,
        dec=table["dec"].value.data,
        fluxes=fluxes,
    )
