from concurrent.futures import ThreadPoolExecutor, TimeoutError
from datetime import datetime
from enum import Enum

import numpy as np
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table
from astropy.time import Time
from astropy.units import Quantity

from cabaret.sources import Sources

__all__ = [
    "Filters",
    "GaiaTAPSource",
    "GaiaQuery",
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


class GaiaTAPSource(Enum):
    """TAP service endpoints for Gaia DR3 data.

    Examples
    --------
    >>> from cabaret.queries import GaiaTAPSource
    >>> GaiaTAPSource.VIZIER
    <GaiaTAPSource.VIZIER: 'https://tapvizier.cds.unistra.fr/TAPVizieR/tap'>
    >>> GaiaTAPSource.ensure_enum("GAIA")
    <GaiaTAPSource.GAIA: 'https://gea.esac.esa.int/tap-server/tap'>
    """

    GAIA = "https://gea.esac.esa.int/tap-server/tap"
    """ESA Gaia Archive TAP service."""
    VIZIER = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap"
    """CDS VizieR TAP service (hosts a copy of Gaia DR3)."""

    @classmethod
    def ensure_enum(cls, value: "GaiaTAPSource | str") -> "GaiaTAPSource":
        """Convert a string or GaiaTAPSource to GaiaTAPSource enum."""
        if isinstance(value, cls):
            return value
        elif isinstance(value, str):
            try:
                return cls[value.upper()]
            except KeyError:
                raise ValueError(
                    f"Invalid TAP source {value!r}. "
                    f"Valid options are: {[m.name for m in cls]}"
                )
        else:
            raise ValueError(
                f"Value must be a GaiaTAPSource enum or string, got {type(value)}"
            )


# Per-source ADQL building blocks. All sources expose the same normalised
# column names via AS aliases so the rest of the code needs no changes.
_TAP_CONFIG: dict[GaiaTAPSource, dict] = {
    GaiaTAPSource.GAIA: {
        "from": "gaiadr3.gaia_source AS gaia",
        "ra": "gaia.ra",
        "dec": "gaia.dec",
        "pmra": "gaia.pmra",
        "pmdec": "gaia.pmdec",
        "g_flux": "phot_g_mean_flux",
        "bp_flux": "phot_bp_mean_flux",
        "rp_flux": "phot_rp_mean_flux",
        "tmass_joins": [
            "INNER JOIN gaiadr3.tmass_psc_xsc_best_neighbour AS tmass_match"
            " ON tmass_match.source_id = gaia.source_id",
            "INNER JOIN external.tmass_psc AS tmass"
            " ON tmass.designation = tmass_match.original_ext_source_id",
        ],
        "j_col": "tmass.j_m",
        "h_col": "tmass.h_m",
        "ks_col": "tmass.ks_m",
    },
    GaiaTAPSource.VIZIER: {
        "from": '"I/355/gaiadr3" AS g',
        "ra": "g.RA_ICRS",
        "dec": "g.DE_ICRS",
        "pmra": "g.pmRA",
        "pmdec": "g.pmDE",
        "g_flux": "g.FG",
        "bp_flux": "g.FBP",
        "rp_flux": "g.FRP",
        "tmass_joins": [
            'INNER JOIN "II/246/out" AS t ON g."2MASS" = t."2MASS"',
        ],
        "j_col": "t.Jmag",
        "h_col": "t.Hmag",
        "ks_col": "t.Kmag",
    },
}

# Map Filters enum names to the per-source 2MASS column key in _TAP_CONFIG.
_TMASS_COL_KEY = {"J": "j_col", "H": "h_col", "KS": "ks_col"}


class GaiaQuery:
    """Class to query Gaia DR3 data and retrieve sources.

    The class provides methods to query a configurable TAP service and return
    either the raw Astropy Table or a Sources instance with RA-DEC coordinates
    and fluxes.

    The TAP service is selected via ``GaiaQuery.DEFAULT_TAP_SOURCE`` (class
    level) or the per-call ``tap_source`` argument.  The default is
    ``GaiaTAPSource.VIZIER`` (CDS VizieR), which hosts a copy of Gaia DR3 and
    is a reliable fallback when the ESA Gaia Archive is unavailable.

    Examples
    --------
    >>> from cabaret.queries import GaiaQuery
    >>> from astropy.coordinates import SkyCoord
    >>> center = SkyCoord(ra=10.68458, dec=41.269, unit="deg")

    The Astropy Table from Gaia can be obtained with:

    >>> table = GaiaQuery.query(center, radius=0.05, limit=10, timeout=30)

    Whereas a Sources instance carrying coordinates and fluxes can be queried for with:

    >>> sources = GaiaQuery.get_sources(center, radius=0.05, limit=10, timeout=30)
    """

    DEFAULT_TAP_SOURCE: GaiaTAPSource = GaiaTAPSource.VIZIER
    """Default TAP service used when ``tap_source=None`` is passed to query methods."""

    @staticmethod
    def query(
        center: tuple[float, float] | SkyCoord,
        radius: float | Angle,
        filter_band: Filters = Filters.G,
        limit: int = 100000,
        timeout: float | None = None,
        tap_source: GaiaTAPSource | str | None = None,
    ) -> Table:
        """Query a Gaia DR3 TAP service within a given radius around the center.

        Parameters
        ----------
        center : tuple or astropy.coordinates.SkyCoord
            The sky coordinates of the center of the FOV.
            If a tuple is given, it should contain the RA and DEC in degrees.
        radius : float or astropy.units.Quantity
            The radius of the FOV in degrees. If a Quantity is given, it must be
            convertible to degrees.
        filter_band : Filters or str, optional
            The filter to use for the flux column. Default is Filters.G.
        limit : int, optional
            The maximum number of sources to retrieve from the Gaia archive.
            By default, it is set to 100000.
        timeout : float, optional
            The maximum time to wait for the Gaia query to complete, in seconds.
            If None, there is no timeout. By default, it is set to None.
        tap_source : GaiaTAPSource, str, or None, optional
            TAP service to query. Accepts a ``GaiaTAPSource`` member or its name
            as a string (e.g. ``"GAIA"`` or ``"VIZIER"``). If None, uses
            ``GaiaQuery.DEFAULT_TAP_SOURCE`` (default: ``GaiaTAPSource.VIZIER``).

        Returns
        -------
        astropy.table.Table
            The raw Astropy Table returned by the TAP service, with columns
            normalised to ``ra``, ``dec``, ``pmra``, ``pmdec``, and the flux
            column named after the filter (e.g. ``phot_g_mean_flux``).

        Example
        -------
        >>> from cabaret.queries import GaiaQuery
        >>> from astropy.coordinates import SkyCoord
        >>> center = SkyCoord(ra=10.68458, dec=41.26917, unit='deg')
        >>> table = GaiaQuery.query(center, radius=0.1, limit=10, timeout=30)
        """
        if tap_source is None:
            tap_source = GaiaQuery.DEFAULT_TAP_SOURCE
        tap_source = GaiaTAPSource.ensure_enum(tap_source)

        filter_band = Filters.ensure_enum(filter_band)

        if isinstance(center, SkyCoord):
            ra = center.ra.deg  # type: ignore
            dec = center.dec.deg  # type: ignore
        else:
            ra, dec = center

        cfg = _TAP_CONFIG[tap_source]

        # Flux column expression (source-dependent) aliased to the standard name.
        if Filters.is_tmass(filter_band.name):
            gaia_flux_expr = cfg["g_flux"]
            gaia_flux_alias = Filters.G.value
            tmass_col_expr = cfg[_TMASS_COL_KEY[filter_band.name]]
            tmass_col_alias = filter_band.value
        else:
            flux_key = filter_band.name.lower() + "_flux"  # g_flux, bp_flux, rp_flux
            gaia_flux_expr = cfg[flux_key]
            gaia_flux_alias = filter_band.value
            tmass_col_expr = None
            tmass_col_alias = None

        select_cols = [
            f"{cfg['ra']} AS ra",
            f"{cfg['dec']} AS dec",
            f"{cfg['pmra']} AS pmra",
            f"{cfg['pmdec']} AS pmdec",
            f"{gaia_flux_expr} AS {gaia_flux_alias}",
        ]
        joins = []
        where = []

        if Filters.is_tmass(filter_band.name):
            select_cols.append(f"{tmass_col_expr} AS {tmass_col_alias}")
            joins.extend(cfg["tmass_joins"])
            order_by = f"{tmass_col_alias} ASC"
            where.append(f"{tmass_col_expr} IS NOT NULL")
        else:
            order_by = f"{gaia_flux_alias} DESC"
            where.append(f"{gaia_flux_expr} IS NOT NULL")

        radius = radius.value if isinstance(radius, Quantity) else float(radius)
        where.append(
            f"1=CONTAINS("
            f"POINT('ICRS', {cfg['ra']}, {cfg['dec']}), "
            f"CIRCLE('ICRS', {ra}, {dec}, {radius}))"
        )

        select_clause = ", ".join(select_cols)
        joins_clause = "\n".join(joins)
        where_clause = " AND ".join(where)

        query = f"""
        SELECT TOP {limit} {select_clause}
        FROM {cfg["from"]}
        {joins_clause}
        WHERE {where_clause}
        ORDER BY {order_by}
        """

        table = GaiaQuery._launch_job_with_timeout(
            query, tap_source=tap_source, timeout=timeout
        )
        return table

    @staticmethod
    def get_sources(
        center: tuple[float, float] | SkyCoord,
        radius: float | Angle,
        filter_band: Filters | str = Filters.G,
        dateobs: datetime | None = None,
        limit: int = 100000,
        timeout: float | None = None,
        tap_source: GaiaTAPSource | str | None = None,
    ) -> Sources:
        """
        Query a Gaia DR3 TAP service to retrieve the RA-DEC coordinates of stars
        within a given radius of a center position, along with their fluxes.

        Parameters
        ----------
        center : tuple or astropy.coordinates.SkyCoord
            The sky coordinates of the center of the FOV.
            If a tuple is given, it should contain the RA and DEC in degrees.
        radius : float or astropy.units.Quantity
            The radius of the FOV in degrees. If a Quantity is given, it must be
            convertible to degrees.
        filter_band : Filters or str, optional
            The filter to use for the flux column. Default is Filters.G.
        dateobs : datetime.datetime, optional
            The date of the observation. If given, the proper motions of the sources
            will be taken into account. By default, it is set to None.
        limit : int, optional
            The maximum number of sources to retrieve from the Gaia archive.
            By default, it is set to 10000.
        timeout : float, optional
            The maximum time to wait for the Gaia query to complete, in seconds.
            If None, there is no timeout. By default, it is set to None.
        tap_source : GaiaTAPSource, str, or None, optional
            TAP service to query. Accepts a ``GaiaTAPSource`` member or its name
            as a string (e.g. ``"GAIA"`` or ``"VIZIER"``). If None, uses
            ``GaiaQuery.DEFAULT_TAP_SOURCE`` (default: ``GaiaTAPSource.VIZIER``).

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
        >>> from cabaret.queries import GaiaQuery
        >>> from astropy.coordinates import SkyCoord
        >>> center = SkyCoord(ra=10.68458, dec=41.26917, unit='deg')
        >>> sources = GaiaQuery.get_sources(
        ...     center, radius=0.1, timeout=30
        ... )  # doctest: +SKIP
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
        filter_band = Filters.ensure_enum(filter_band)

        table = GaiaQuery.query(
            center=center,
            radius=radius,
            limit=limit,
            timeout=timeout,
            filter_band=filter_band,
            tap_source=tap_source,
        )

        if dateobs is not None:
            table = GaiaQuery._apply_proper_motion(table, dateobs)

        if Filters.is_tmass(filter_band.name):
            fluxes = GaiaQuery._tmass_mag_to_photons(
                table[filter_band.value].value.data,  # type: ignore
                filter_band,
            )
        else:
            fluxes = table[filter_band.value].value.data  # type: ignore
        table.remove_rows(np.isnan(fluxes))
        fluxes = fluxes[~np.isnan(fluxes)]

        return Sources.from_arrays(
            ra=table["ra"].value.data,  # type: ignore
            dec=table["dec"].value.data,  # type: ignore
            fluxes=fluxes,
        )

    @staticmethod
    def _launch_job_with_timeout(
        query: str,
        tap_source: GaiaTAPSource,
        timeout: float | None = None,
        **kwargs,
    ) -> Table:
        """
        Launch a TAP job and return its results, optionally enforcing a timeout.

        Parameters
        ----------
        query : str
            The ADQL query string.
        tap_source : GaiaTAPSource
            The TAP service to use.
        timeout : float or None, optional
            Maximum number of seconds to wait for the job to complete.
            If None, the job is run on the main thread (no thread overhead).
        **kwargs
            Additional keyword arguments forwarded to ``TapPlus.launch_job``.

        Returns
        -------
        object
            The result returned by job.get_results().

        Raises
        ------
        TimeoutError
            If `timeout` is not None and the call does not complete within `timeout`.
        """
        from astroquery.utils.tap.core import TapPlus

        def _run():
            tap = TapPlus(url=tap_source.value)
            job = tap.launch_job(query, **kwargs)
            return job.get_results()  # type: ignore

        if timeout is None:
            return _run()

        with ThreadPoolExecutor() as executor:
            future = executor.submit(_run)
            try:
                return future.result(timeout=timeout)
            except TimeoutError:
                raise TimeoutError(
                    "Gaia query timed out."
                    " You may want to increase the timeout or reduce the query size."
                    f" Query was: {query}"
                )

    @staticmethod
    def _tmass_mag_to_photons(mags: np.ndarray, filter_band: Filters) -> np.ndarray:
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
                f"tmass_mag_to_photons expects a 2MASS filter (J,H,KS), "
                f"got {filter_band}."
            )

        Jy = 1.51e7  # [photons sec^-1 m^-2 (dlambda/lambda)^-1]
        photons = dlam_lam * flux_m0 * Jy  # [photons sec^-1 m^-2] at mag 0
        return photons * 10 ** (-0.4 * mags)

    # Gaia DR3 reference epoch J2016.0 expressed as a decimal year.
    _GAIA_DR3_EPOCH = float(Time(2016.0, format="jyear").decimalyear)

    @staticmethod
    def _apply_proper_motion(table: Table, dateobs: datetime | Time) -> Table:
        """
        Apply proper motion correction to RA and DEC columns
        for the given observation date.
        """
        if isinstance(dateobs, datetime):
            dateobs = Time(dateobs)

        if not isinstance(dateobs, Time):
            raise ValueError(
                f"dateobs must be an astropy.time.Time or datetime, got {type(dateobs)}"
            )

        years = float(dateobs.decimalyear) - GaiaQuery._GAIA_DR3_EPOCH
        table["ra"] += years * table["pmra"] / 1000 / 3600  # type: ignore
        table["dec"] += years * table["pmdec"] / 1000 / 3600  # type: ignore

        return table
