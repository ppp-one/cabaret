from concurrent.futures import ThreadPoolExecutor, TimeoutError
from datetime import datetime
from enum import Enum

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.units import Quantity
from astroquery.gaia import Gaia

from cabaret.sources import Sources

__all__ = [
    "Filters",
    "tmass_mag_to_photons",
    "gaia_query",
    "get_gaia_sources",
]


class Filters(Enum):
    """Allowed Gaia and 2MASS flux filter_band strings.

    Examples
    --------
    >>> from cabaret.queries import Filters
    >>> Filters.G
    <Filters.G: 'phot_g_mean_flux'>
    >>> Filters.from_string('RP')
    <Filters.RP: 'phot_rp_mean_flux'>
    >>> Filters.is_tmass('J')
    True
    >>> Filters.options()
    ('G', 'BP', 'RP', 'J', 'H', 'KS')
    """

    G = "phot_g_mean_flux"
    """ Gaia G band flux in [e-/s] """
    BP = "phot_bp_mean_flux"
    """ Gaia BP band flux in [e-/s] """
    RP = "phot_rp_mean_flux"
    """ Gaia RP band flux in [e-/s] """
    J = "j_m"
    """ 2MASS J band magnitude """
    H = "h_m"
    """ 2MASS H band magnitude """
    KS = "ks_m"
    """ 2MASS KS band magnitude """

    @classmethod
    def options(cls) -> tuple[str, ...]:
        """Return all valid filter_band options."""
        return tuple(cls.__members__.keys())

    @classmethod
    def from_string(cls, value: str) -> "Filters":
        """Return the Filters enum member for a given string."""
        try:
            return cls[value.upper()]
        except KeyError:
            raise ValueError(
                f"Invalid filter_band string: {value}. "
                f"Valid options are: {cls.options()}"
            )

    @classmethod
    def is_tmass(cls, value: "Filters | str") -> bool:
        """Check if the filter_band string is a 2MASS filter_band."""
        if isinstance(value, cls):
            return value.name in ("J", "H", "KS")
        elif isinstance(value, str):
            return value.upper() in ("J", "H", "KS")
        else:
            raise ValueError(
                f"Value must be an Filters enum or string, got {type(value)}"
            )

    @classmethod
    def ensure_enum(cls, value: "Filters | str") -> "Filters":
        """Convert a string or Filters to Filters enum."""
        if isinstance(value, cls):
            return value
        elif isinstance(value, str):
            return cls.from_string(value)
        else:
            raise ValueError(
                f"Value must be a Filters enum or string, got {type(value)}"
            )

    @classmethod
    def is_valid(cls, value: str) -> bool:
        """Check if the filter_band string is valid."""
        return value.upper() in cls.__members__


def tmass_mag_to_photons(mags: np.ndarray, filter_band: Filters) -> np.ndarray:
    """Convert 2MASS J magnitudes to photon fluxes at mag 0.

    Reference: https://lweb.cfa.harvard.edu/~dfabricant/huchra/ay145/mags.html
    Returns photons/sec/m^2 for each magnitude.
    """
    # Use a dict for 2MASS filter properties
    dlam_lam_map = {"J": 0.16, "H": 0.23, "KS": 0.23}
    flux_m0_map = {"J": 1600, "H": 1080, "KS": 670}

    try:
        dlam_lam = dlam_lam_map[filter_band.name]
        flux_m0 = flux_m0_map[filter_band.name]
    except KeyError:
        raise ValueError(
            f"tmass_mag_to_photons expects a 2MASS filter (J,H,KS), got {filter_band}"
        )

    Jy = 1.51e7  # [photons sec^-1 m^-2 (dlambda/lambda)^-1]
    photons = dlam_lam * flux_m0 * Jy  # [photons sec^-1 m^-2] at mag 0
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
        return job.get_results()  # type: ignore

    with ThreadPoolExecutor() as executor:
        future = executor.submit(Gaia.launch_job, query, **kwargs)
        try:
            job = future.result(timeout=timeout)
            return job.get_results()  # type: ignore
        except TimeoutError:
            raise TimeoutError(
                "Gaia query timed out."
                " You may want to increase the timeout or reduce the query size."
                f"Query was: {query}"
            )


def gaia_query(
    center: tuple[float, float] | SkyCoord,
    radius: float | Quantity,
    limit: int = 100000,
    timeout: float | None = None,
    filter_band: Filters = Filters.G,
) -> Table:
    """Query Gaia and return the raw Astropy Table.

    Example
    -------
    >>> from cabaret.queries import gaia_query
    >>> from astropy.coordinates import SkyCoord
    >>> center = SkyCoord(ra=10.68458, dec=41.26917, unit='deg')
    >>> table = gaia_query(center, fov=0.1, limit=10, timeout=30)
    """
    filter_band = Filters.ensure_enum(filter_band)

    if isinstance(center, SkyCoord):
        ra = center.ra.deg  # type: ignore
        dec = center.dec.deg  # type: ignore
    else:
        ra, dec = center

    gaia_flux_column = (
        filter_band.value if not Filters.is_tmass(filter_band.name) else Filters.G.value
    )
    select_cols = [
        "gaia.ra",
        "gaia.dec",
        "gaia.pmra",
        "gaia.pmdec",
        f"{gaia_flux_column}",
    ]
    joins = []
    where = []
    order_by = f"{filter_band.value} DESC"  # Prefer brighter stars

    if Filters.is_tmass(filter_band.name):
        select_cols.append(filter_band.value)
        joins.extend(
            [
                "INNER JOIN gaiadr2.tmass_best_neighbour AS tmass_match "
                + "ON tmass_match.source_id = gaia.source_id",
                "INNER JOIN gaiadr1.tmass_original_valid AS tmass "
                + "ON tmass.tmass_oid = tmass_match.tmass_oid",
            ]
        )
        order_by = f"{filter_band.value} ASC"  # <-- fix here
        where.append(f"{filter_band.value} IS NOT NULL")  # <-- fix here
    else:
        where.append(f"{filter_band.value} IS NOT NULL")

    radius_val = radius.value if isinstance(radius, Quantity) else float(radius)
    where.append(
        f"1=CONTAINS(POINT('ICRS', {ra}, {dec}), "
        f"CIRCLE('ICRS', gaia.ra, gaia.dec, {radius_val}))"
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
    table["ra"] += years * table["pmra"] / 1000 / 3600  # type: ignore
    table["dec"] += years * table["pmdec"] / 1000 / 3600  # type: ignore
    return table


def get_gaia_sources(
    center: tuple[float, float] | SkyCoord,
    radius: float | Quantity,
    limit: int = 100000,
    dateobs: datetime | None = None,
    timeout: float | None = None,
    filter_band: Filters | str = Filters.G,
) -> Sources:
    """
    Query the Gaia archive to retrieve the RA-DEC coordinates of stars
    within a given field-of-view (FOV) centered on a given sky position.

    Parameters
    ----------
    center : tuple or astropy.coordinates.SkyCoord
        The sky coordinates of the center of the FOV.
        If a tuple is given, it should contain the RA and DEC in degrees.
    radius : float or astropy.units.Quantity
        The radius of the FOV in degrees. If a Quantity is given, it must be
        convertible to degrees.
    limit : int, optional
        The maximum number of sources to retrieve from the Gaia archive.
        By default, it is set to 10000.
    dateobs : datetime.datetime, optional
        The date of the observation. If given, the proper motions of the sources
        will be taken into account. By default, it is set to None.
    filter_band : Filters or str, optional
        The filter to use for the flux column. Default is Filters.G.
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
    If `filter_band` is a 2MASS filter (J, H, KS), the fluxes are calculated
    from the 2MASS magnitudes using `cabaret.queries.tmass_mag_to_photons`.

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
    filter_band_instance = Filters.ensure_enum(filter_band)

    table = gaia_query(
        center=center,
        radius=radius,
        limit=limit,
        timeout=timeout,
        filter_band=filter_band_instance,
    )

    if dateobs is not None:
        table = apply_proper_motion(table, dateobs)

    if Filters.is_tmass(filter_band_instance.name):
        fluxes = tmass_mag_to_photons(
            table[filter_band_instance.value].value.data,  # type: ignore
            filter_band_instance,
        )
    else:
        fluxes = table[filter_band_instance.value].value.data  # type: ignore
    table.remove_rows(np.isnan(fluxes))
    fluxes = fluxes[~np.isnan(fluxes)]

    return Sources.from_arrays(
        ra=table["ra"].value.data,  # type: ignore
        dec=table["dec"].value.data,  # type: ignore
        fluxes=fluxes,
    )
