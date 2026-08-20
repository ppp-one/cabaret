import math
import re
import sqlite3
import threading
import warnings
from collections.abc import Callable, Sequence
from contextlib import closing
from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from typing import TypeVar

import numpy as np
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table
from astropy.time import Time
from astropy.units import Quantity

from cabaret.sources import Sources

__all__ = [
    "Filters",
    "GaiaSQLiteSource",
    "GaiaTAPSource",
    "GaiaQuery",
]

_T = TypeVar("_T")


@dataclass(frozen=True)
class GaiaSQLiteSource:
    """Configuration for querying a local SQLite catalog.

    Parameters
    ----------
    database : str
        Path to the SQLite database file.
    table : str, optional
        Catalog table name. Default is ``gaia_sources``.
    """

    database: str
    table: str = "gaia_sources"


class Filters(Enum):
    """Allowed Gaia and 2MASS flux filter_band strings.

    Examples
    --------
    >>> from cabaret.queries import Filters
    >>> Filters.G
    <Filters.G: 'phot_g_mean_mag'>
    >>> Filters.from_string('RP')
    <Filters.RP: 'phot_rp_mean_mag'>
    >>> Filters.is_tmass('J')
    True
    >>> Filters.options()
    ('G', 'BP', 'RP', 'J', 'H', 'KS')
    """

    G = "phot_g_mean_mag"
    """Gaia G band magnitude [Gaia Vega system]"""
    BP = "phot_bp_mean_mag"
    """Gaia BP band magnitude [Gaia Vega system]"""
    RP = "phot_rp_mean_mag"
    """Gaia RP band magnitude [Gaia Vega system]"""
    J = "j_m"
    """2MASS J-band magnitude [2MASS Vega system]"""
    H = "h_m"
    """2MASS H-band magnitude [2MASS Vega system]"""
    KS = "ks_m"
    """2MASS KS-band magnitude [2MASS Vega system]"""

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
            name = value.name
        elif isinstance(value, str):
            name = value.upper()
        else:
            raise ValueError(
                f"Value must be a Filters enum or string, got {type(value)}"
            )
        return name in ("J", "H", "KS")

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
    def all(cls) -> list["Filters"]:
        """Return all filter bands in definition order."""
        return list(cls)

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
        "gaia_id": "gaia.source_id",
        "ra": "gaia.ra",
        "dec": "gaia.dec",
        "pmra": "gaia.pmra",
        "pmdec": "gaia.pmdec",
        "g_mag": "phot_g_mean_mag",
        "bp_mag": "phot_bp_mean_mag",
        "rp_mag": "phot_rp_mean_mag",
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
        "gaia_id": 'g."Source"',
        "ra": "g.RA_ICRS",
        "dec": "g.DE_ICRS",
        "pmra": "g.pmRA",
        "pmdec": "g.pmDE",
        "g_mag": "g.Gmag",
        "bp_mag": "g.BPmag",
        "rp_mag": "g.RPmag",
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

# Sentinel used for backends that cannot provide a Gaia source id (SQLite).
_NO_GAIA_ID = -1


class GaiaQuery:
    """Class to query Gaia DR3 data and retrieve sources.

    The class provides methods to query a configurable TAP service and return
    either the raw Astropy Table, a flux-normalised multi-band Table, or a
    Sources instance with RA-DEC coordinates and fluxes.

    The query backend is selected via ``GaiaQuery.DEFAULT_TAP_SOURCE`` (class
    level) or the per-call ``tap_source`` argument.  Accepts a
    ``GaiaTAPSource`` for an online TAP service, a ``GaiaSQLiteSource`` for an
    offline SQLite database, a TAP endpoint URL, a ``sqlite:///...`` URI, or a
    direct path to a ``.db``/``.sqlite``/``.sqlite3`` file.  The default is
    ``GaiaTAPSource.VIZIER`` (CDS VizieR), which hosts a copy of Gaia DR3.

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

    DEFAULT_TAP_SOURCE: GaiaTAPSource | GaiaSQLiteSource | str = GaiaTAPSource.VIZIER
    """Default TAP service used when ``tap_source=None`` is passed to query methods."""

    @staticmethod
    def query(
        center: tuple[float, float] | SkyCoord | None = None,
        radius: float | Angle | None = None,
        bounds: tuple[float, float, float, float] | None = None,
        filter_bands: Filters | str | Sequence[Filters | str] = Filters.G,
        limit: int = 100000,
        timeout: float | None = None,
        tap_source: GaiaTAPSource | GaiaSQLiteSource | str | None = None,
        allow_nulls: bool = False,
        filter_band: Filters | str | None = None,
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
        bounds : tuple or None
            Field-of-view bounds specified as ``(ra_min, ra_max, dec_min, dec_max)``
            in degrees. RA values are interpreted modulo 360 (degrees) and may
            be given so that ``ra_min <= ra_max`` for a contiguous RA interval,
            or ``ra_min > ra_max`` to indicate an interval that wraps across
            RA=0 (for example ``(350.0, 10.0, -5.0, 5.0)``). Declination bounds
            must satisfy ``-90 <= dec_min <= dec_max <= 90``.
        filter_bands : Filters, str, or sequence thereof, optional
            One or more filter bands to include as magnitude columns. Accepts a single
            ``Filters`` member or its name as a string, or a list of either.
            Pass ``"all"`` to request every available band (``Filters.all()``).
            Default is ``Filters.G``. When multiple bands are requested, the
            ``ORDER BY`` is determined by the first band (brightest-first,
            ASC for all magnitude columns).
            If any 2MASS band is included the 2MASS cross-match join is added.
        limit : int, optional
            The maximum number of sources to retrieve from the Gaia archive.
            By default, it is set to 100000.
        timeout : float, optional
            The maximum time to wait for the Gaia query to complete, in seconds.
            If None, there is no timeout. By default, it is set to None.
        tap_source : GaiaTAPSource, GaiaSQLiteSource, str, or None, optional
            Query backend to use. Accepts a ``GaiaTAPSource`` for an online TAP
            service, a ``GaiaSQLiteSource`` for an offline SQLite database, a TAP
            endpoint URL as a string, a ``sqlite:///...`` URI, or a direct path to
            a SQLite database file ending in ``.db``, ``.sqlite``, or
            ``.sqlite3``. If None, the default TAP backend is used.
        allow_nulls : bool, optional
            If False (default), only rows where all requested band columns are
            non-NULL are returned (``IS NOT NULL`` filter per band). Set to True
            to allow rows with missing magnitude values through.
        filter_band : Filters or str, optional
            Deprecated. Use ``filter_bands`` instead.

        Returns
        -------
        astropy.table.Table
            The raw Astropy Table returned by the TAP service, with columns
            normalised to ``gaia_id``, ``ra``, ``dec``, ``pmra``, ``pmdec``, and
            one magnitude column per requested band named after
            ``filter_band.value`` (e.g. ``phot_g_mean_mag``, ``h_m``). The
            SQLite backend has no Gaia identifier available and fills
            ``gaia_id`` with ``-1``.

        Examples
        --------
        >>> from cabaret.queries import GaiaQuery
        >>> from astropy.coordinates import SkyCoord
        >>> center = SkyCoord(ra=10.68458, dec=41.26917, unit='deg')
        >>> table = GaiaQuery.query(center, radius=0.1, limit=10, timeout=30)
        """
        if filter_band is not None:
            warnings.warn(
                "filter_band is deprecated; use filter_bands instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            filter_bands = filter_band
        requested_all = isinstance(filter_bands, str) and filter_bands.upper() == "ALL"
        bands = GaiaQuery._normalize_bands(filter_bands)
        center_coords, radius_deg, bounds_tuple = GaiaQuery._normalize_region(
            center=center,
            radius=radius,
            bounds=bounds,
        )
        backend = GaiaQuery._resolve_query_source(tap_source)

        if isinstance(backend, GaiaTAPSource):
            table = GaiaQuery._query_tap(
                tap_source=backend,
                bands=bands,
                center=center_coords,
                radius_deg=radius_deg,
                bounds=bounds_tuple,
                limit=limit,
                timeout=timeout,
                allow_nulls=allow_nulls,
            )
        else:
            table = GaiaQuery._query_sqlite(
                sqlite_source=backend,
                requested_bands=bands,
                requested_all=requested_all,
                center=center_coords,
                radius_deg=radius_deg,
                bounds=bounds_tuple,
                limit=limit,
                timeout=timeout,
                allow_nulls=allow_nulls,
            )
        return table

    @staticmethod
    def get_flux_table(
        center: tuple[float, float] | SkyCoord | None = None,
        radius: float | Angle | None = None,
        bounds: tuple[float, float, float, float] | None = None,
        filter_bands: Filters | str | Sequence[Filters | str] = Filters.G,
        dateobs: datetime | None = None,
        limit: int = 100000,
        timeout: float | None = None,
        tap_source: GaiaTAPSource | GaiaSQLiteSource | str | None = None,
        allow_nulls: bool = False,
        keep_mag: bool = False,
    ) -> Table:
        """Query and return a Table with all columns expressed as physical fluxes.

        Identical to :meth:`query` but additionally:

        * applies proper-motion correction when ``dateobs`` is given, and
        * converts supported magnitude columns to photons/s/m² using
          :meth:`_mag_to_photons` and renames them (e.g. ``"j_m"`` →
          ``"j_flux"`` and ``"phot_g_mean_mag"`` → ``"g_flux"``).

        Parameters
        ----------
        center : tuple or astropy.coordinates.SkyCoord
            The sky coordinates of the center of the FOV.
        radius : float or astropy.units.Quantity
            The radius of the FOV in degrees.
        bounds : tuple or None
            Field-of-view bounds specified as ``(ra_min, ra_max, dec_min, dec_max)``
            in degrees. RA values are interpreted modulo 360 (degrees) and may
            be given so that ``ra_min <= ra_max`` for a contiguous RA interval,
            or ``ra_min > ra_max`` to indicate an interval that wraps across
            RA=0 (for example ``(350.0, 10.0, -5.0, 5.0)``). Declination bounds
            must satisfy ``-90 <= dec_min <= dec_max <= 90``.
        filter_bands : Filters, str, or sequence thereof, optional
            One or more filter bands. Default is ``Filters.G``.
        dateobs : datetime.datetime or None, optional
            Observation date for proper-motion correction. Default is None.
        limit : int, optional
            Maximum number of sources to retrieve. Default is 100000.
        timeout : float or None, optional
            Query timeout in seconds. Default is None (no timeout).
        tap_source : GaiaTAPSource, GaiaSQLiteSource, str, or None, optional
            Query backend to use. Accepts a ``GaiaTAPSource``, ``GaiaSQLiteSource``,
            TAP URL, ``sqlite:///...`` URI, or ``.db``/``.sqlite``/``.sqlite3`` path.
            Default is ``GaiaQuery.DEFAULT_TAP_SOURCE``.
        allow_nulls : bool, optional
            Forwarded to :meth:`query`. If True, rows with NULL magnitude values
            are included in the result. Default is False.
        keep_mag : bool, optional
            If True, retain the original magnitude column (e.g.
            ``"phot_g_mean_mag"``) alongside the converted flux column
            (e.g. ``"g_flux"``). Default is False.

        Returns
        -------
        astropy.table.Table
            Table with columns ``ra``, ``dec``, ``pmra``, ``pmdec``, and one
            flux column per requested band, all in photons/s/m². Gaia band
            columns are renamed from ``"phot_<b>_mean_mag"`` to
            ``"<b>_flux"`` (e.g. ``"g_flux"``); 2MASS band columns are renamed
            from ``"<band>_m"`` to ``"<band>_flux"`` (e.g. ``"h_flux"``).
            When ``keep_mag=True``, the original magnitude columns are also
            present.

        Examples
        --------
        >>> from cabaret.queries import GaiaQuery, Filters
        >>> from astropy.coordinates import SkyCoord
        >>> center = SkyCoord(ra=10.68458, dec=41.26917, unit='deg')
        >>> table = GaiaQuery.get_flux_table(
        ...     center, radius=0.1, filter_bands=[Filters.G, Filters.H, Filters.KS],
        ...     limit=10, timeout=30,
        ... )  # doctest: +SKIP
        """
        bands = GaiaQuery._normalize_bands(filter_bands)

        table = GaiaQuery.query(
            center=center,
            radius=radius,
            bounds=bounds,
            filter_bands=bands,
            limit=limit,
            timeout=timeout,
            tap_source=tap_source,
            allow_nulls=allow_nulls,
        )

        if dateobs is not None:
            table = GaiaQuery._apply_proper_motion(table, dateobs)

        seen: set[Filters] = set()
        for band in bands:
            if band in seen:
                continue
            seen.add(band)
            col = band.value
            if Filters.is_tmass(band.name):
                new_name = str(col).removesuffix("_m") + "_flux"
            else:
                new_name = band.name.lower() + "_flux"  # e.g. "g_flux", "bp_flux"
            table[new_name] = GaiaQuery._mag_to_photons(
                np.ma.filled(table[col].value, np.nan),  # type: ignore
                band,
            )
            if not keep_mag:
                table.remove_column(col)

        return table

    @staticmethod
    def get_sources(
        center: tuple[float, float] | SkyCoord | None = None,
        radius: float | Angle | None = None,
        bounds: tuple[float, float, float, float] | None = None,
        filter_band: Filters | str = Filters.G,
        dateobs: datetime | None = None,
        limit: int = 100000,
        timeout: float | None = None,
        tap_source: GaiaTAPSource | GaiaSQLiteSource | str | None = None,
    ) -> Sources:
        """
        Query a Gaia DR3 TAP service to retrieve the RA-DEC coordinates of stars
        within a given radius of a center position, along with their fluxes.

        Parameters
        ----------
        center : tuple or astropy.coordinates.SkyCoord
            The sky coordinates of the center of the FOV.
            If a tuple is given, it should contain the RA and DEC in degrees.
        radius : float or astropy.units.Quantity or None
            The radius of the FOV in degrees. If a Quantity is given, it must be
            convertible to degrees.
        bounds : tuple or None
            Field-of-view bounds specified as ``(ra_min, ra_max, dec_min, dec_max)``
            in degrees. RA values are interpreted modulo 360 (degrees) and may
            be given so that ``ra_min <= ra_max`` for a contiguous RA interval,
            or ``ra_min > ra_max`` to indicate an interval that wraps across
            RA=0 (for example ``(350.0, 10.0, -5.0, 5.0)``). Declination bounds
            must satisfy ``-90 <= dec_min <= dec_max <= 90``.
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
        tap_source : GaiaTAPSource, GaiaSQLiteSource, str, or None, optional
            Query backend to use. Accepts a ``GaiaTAPSource`` for an online TAP
            service, a ``GaiaSQLiteSource`` for an offline SQLite database, a TAP
            endpoint URL as a string, a ``sqlite:///...`` URI, or a direct path to
            a SQLite database file ending in ``.db``, ``.sqlite``, or
            ``.sqlite3``. If None, the default TAP backend is used.

        Returns
        -------
        Sources
            A Sources instance containing the coordinates, fluxes and Gaia source
            ids of the retrieved sources.


        Notes
        -----
        Fluxes are always returned in photons/s/m² via :meth:`_mag_to_photons`.

        ``gaia_ids`` holds the Gaia DR3 ``source_id`` when a TAP backend is used.
        The SQLite backend carries no such identifier, so every id is ``-1``.

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
        ...     center, radius=0.1, timeout=30, limit=10
        ... )  # doctest: +SKIP

        """
        filter_band = Filters.ensure_enum(filter_band)

        table = GaiaQuery.query(
            center=center,
            radius=radius,
            bounds=bounds,
            limit=limit,
            timeout=timeout,
            filter_bands=filter_band,
            tap_source=tap_source,
        )

        if dateobs is not None:
            table = GaiaQuery._apply_proper_motion(table, dateobs)

        fluxes = GaiaQuery._mag_to_photons(
            np.ma.filled(table[filter_band.value].value, np.nan),  # type: ignore
            filter_band,
        )
        table.remove_rows(np.isnan(fluxes))
        fluxes = fluxes[~np.isnan(fluxes)]

        return Sources.from_arrays(
            ra=table["ra"].value.data,  # type: ignore
            dec=table["dec"].value.data,  # type: ignore
            fluxes=fluxes,
            gaia_ids=np.asarray(table["gaia_id"], dtype=np.int64),
        )

    @staticmethod
    def _resolve_query_source(
        tap_source: GaiaTAPSource | GaiaSQLiteSource | str | None,
    ) -> GaiaTAPSource | GaiaSQLiteSource:
        """Resolve a query source into either TAP or SQLite backend config."""
        if tap_source is None:
            return GaiaQuery.DEFAULT_TAP_SOURCE
        if isinstance(tap_source, GaiaSQLiteSource | GaiaTAPSource):
            return tap_source
        if not isinstance(tap_source, str):
            raise ValueError(
                "tap_source must be a GaiaTAPSource, GaiaSQLiteSource, or str, "
                f"got {type(tap_source)}"
            )

        source_text = tap_source.strip()
        if source_text.lower().startswith("sqlite:///"):
            db_path = source_text[len("sqlite:///") :]
            if not db_path:
                raise ValueError("SQLite URI must include a database path.")
            return GaiaSQLiteSource(database=db_path)
        if source_text.lower().endswith((".db", ".sqlite", ".sqlite3")):
            return GaiaSQLiteSource(database=source_text)
        return GaiaTAPSource.ensure_enum(source_text)

    @staticmethod
    def _normalize_region(
        center: tuple[float, float] | SkyCoord | None,
        radius: float | Angle | None,
        bounds: tuple[float, float, float, float] | None,
    ) -> tuple[
        tuple[float, float] | None,
        float | None,
        tuple[float, float, float, float] | None,
    ]:
        """Normalize circle/bounds geometry inputs.

        Exactly one of ``radius`` or ``bounds`` must be provided.

        The ``bounds`` tuple should be ``(ra_min, ra_max, dec_min, dec_max)``
        in degrees. RA values are normalised modulo 360; if ``ra_min <= ra_max``
        the RA interval is contiguous, otherwise ``ra_min > ra_max`` indicates
        an interval that wraps across RA=0. Declination bounds must satisfy
        ``-90 <= dec_min <= dec_max <= 90``. The function returns a tuple of
        ``(center_coords, radius_deg, bounds_tuple)`` where only the
        applicable geometry is populated.
        """
        if (radius is None) == (bounds is None):
            raise ValueError(
                "Exactly one geometry must be specified: either radius or bounds."
            )

        center_coords: tuple[float, float] | None = None
        radius_deg: float | None = None
        bounds_tuple: tuple[float, float, float, float] | None = None

        if radius is not None:
            if center is None:
                raise ValueError("center must be provided for radius-based queries.")
            radius_value = (
                radius.value if isinstance(radius, Quantity) else float(radius)
            )
            if radius_value <= 0:
                raise ValueError("radius must be > 0.")
            radius_deg = radius_value
            if isinstance(center, SkyCoord):
                center_coords = (
                    float(center.ra.deg % 360.0),  # type: ignore
                    float(center.dec.deg),  # type: ignore
                )
            else:
                center_coords = (float(center[0] % 360.0), float(center[1]))
        else:
            assert bounds is not None
            if len(bounds) != 4:
                raise ValueError(
                    "bounds must be a tuple of (ra_min, ra_max, dec_min, dec_max)."
                )
            ra_min, ra_max, dec_min, dec_max = bounds
            dec_min = float(dec_min)
            dec_max = float(dec_max)
            if dec_min > dec_max:
                raise ValueError("dec_min must be <= dec_max.")
            if dec_min < -90.0 or dec_max > 90.0:
                raise ValueError("declination bounds must stay within [-90, 90].")
            bounds_tuple = (
                float(ra_min % 360.0),
                float(ra_max % 360.0),
                dec_min,
                dec_max,
            )

        return center_coords, radius_deg, bounds_tuple

    @staticmethod
    def _query_tap(
        tap_source: GaiaTAPSource,
        bands: list[Filters],
        center: tuple[float, float] | None,
        radius_deg: float | None,
        bounds: tuple[float, float, float, float] | None,
        limit: int,
        timeout: float | None,
        allow_nulls: bool,
    ) -> Table:
        """Execute a Gaia query against a TAP backend."""
        cfg = _TAP_CONFIG[tap_source]

        select_cols = [
            f"{cfg['gaia_id']} AS gaia_id",
            f"{cfg['ra']} AS ra",
            f"{cfg['dec']} AS dec",
            f"{cfg['pmra']} AS pmra",
            f"{cfg['pmdec']} AS pmdec",
        ]
        where: list[str] = []
        joins: list[str] = []
        need_tmass_join = False

        for band in bands:
            if Filters.is_tmass(band.name):
                col_expr = cfg[_TMASS_COL_KEY[band.name]]
                need_tmass_join = True
            else:
                col_expr = cfg[band.name.lower() + "_mag"]
            select_cols.append(f"{col_expr} AS {band.value}")
            if not allow_nulls:
                where.append(f"{col_expr} IS NOT NULL")

        if need_tmass_join:
            joins.extend(cfg["tmass_joins"])

        first = bands[0]
        order_by = f"{first.value} ASC"

        if radius_deg is not None:
            assert center is not None
            where.append(
                f"1=CONTAINS("
                f"POINT('ICRS', {cfg['ra']}, {cfg['dec']}), "
                f"CIRCLE('ICRS', {center[0]}, {center[1]}, {radius_deg}))"
            )
        else:
            assert bounds is not None
            where.append(
                GaiaQuery._build_tap_bounds_where(
                    cfg["ra"],
                    cfg["dec"],
                    bounds,
                )
            )

        select_clause = ", ".join(select_cols)
        joins_clause = "\n".join(joins)
        where_clause = " AND ".join(where)

        adql = f"""
        SELECT TOP {limit} {select_clause}
        FROM {cfg["from"]}
        {joins_clause}
        WHERE {where_clause}
        ORDER BY {order_by}
        """

        return GaiaQuery._launch_job_with_timeout(
            adql, tap_source=tap_source, timeout=timeout
        )

    @staticmethod
    def _query_sqlite(
        sqlite_source: GaiaSQLiteSource,
        requested_bands: list[Filters],
        requested_all: bool,
        center: tuple[float, float] | None,
        radius_deg: float | None,
        bounds: tuple[float, float, float, float] | None,
        limit: int,
        timeout: float | None,
        allow_nulls: bool,
    ) -> Table:
        """Execute a Gaia-like query against a local SQLite catalog."""

        connection = sqlite3.connect(sqlite_source.database, check_same_thread=False)
        connection.row_factory = sqlite3.Row

        def _run_query() -> Table:
            with closing(connection):
                cursor = connection.cursor()
                selected_tables, schema_table = GaiaQuery._select_sqlite_tables(
                    cursor=cursor,
                    sqlite_source=sqlite_source,
                    center=center,
                    radius_deg=radius_deg,
                    bounds=bounds,
                )
                schema_table_name = GaiaQuery._quote_sql_identifier(schema_table)

                available_columns = {
                    str(row[1])
                    for row in cursor.execute(f"PRAGMA table_info({schema_table_name})")
                }
                required_position_columns = {"ra", "dec"}
                missing_required = sorted(required_position_columns - available_columns)
                if missing_required:
                    raise ValueError(
                        "SQLite source is missing required columns: "
                        f"{', '.join(missing_required)}"
                    )

                if requested_all:
                    selected_bands = [
                        band
                        for band in requested_bands
                        if band.value in available_columns
                    ]
                    if not selected_bands:
                        raise ValueError(
                            "SQLite source does not contain any supported "
                            "magnitude columns."
                        )
                else:
                    missing_requested = [
                        band.value
                        for band in requested_bands
                        if band.value not in available_columns
                    ]
                    if missing_requested:
                        raise ValueError(
                            "SQLite source is missing requested band columns: "
                            f"{', '.join(sorted(missing_requested))}"
                        )
                    selected_bands = requested_bands

                # Local catalogs are not guaranteed to carry Gaia source ids,
                # so they are reported as the -1 sentinel.
                select_cols = [
                    f"{_NO_GAIA_ID} AS gaia_id",
                    '"ra" AS ra',
                    '"dec" AS dec',
                    (
                        '"pmra" AS pmra'
                        if "pmra" in available_columns
                        else "NULL AS pmra"
                    ),
                    (
                        '"pmdec" AS pmdec'
                        if "pmdec" in available_columns
                        else "NULL AS pmdec"
                    ),
                ]
                for band in selected_bands:
                    quoted_col = GaiaQuery._quote_sql_identifier(band.value)
                    select_cols.append(f"{quoted_col} AS {band.value}")

                where_sql: list[str] = []
                params: list[float] = []
                if not allow_nulls:
                    for band in selected_bands:
                        quoted_col = GaiaQuery._quote_sql_identifier(band.value)
                        where_sql.append(f"{quoted_col} IS NOT NULL")

                if bounds is not None:
                    bounds_sql, bounds_params = GaiaQuery._build_sqlite_bounds_where(
                        bounds
                    )
                    where_sql.append(bounds_sql)
                    params.extend(bounds_params)
                else:
                    assert center is not None
                    assert radius_deg is not None
                    circle_sql, circle_params = (
                        GaiaQuery._build_sqlite_circle_prefilter_where(
                            center=center,
                            radius_deg=radius_deg,
                        )
                    )
                    where_sql.append(circle_sql)
                    params.extend(circle_params)

                first_band = selected_bands[0]

                output_names = ["gaia_id", "ra", "dec", "pmra", "pmdec"] + [
                    band.value for band in selected_bands
                ]
                if not selected_tables:
                    return Table(rows=[], names=output_names)

                output_rows: list[tuple] = []
                for table in selected_tables:
                    table_name = GaiaQuery._quote_sql_identifier(table)
                    query_parts = [
                        f"SELECT {', '.join(select_cols)}",
                        f"FROM {table_name}",
                    ]
                    if where_sql:
                        query_parts.append("WHERE " + " AND ".join(where_sql))
                    sql = "\n".join(query_parts)

                    rows = list(cursor.execute(sql, tuple(params)))
                    output_rows.extend(
                        tuple(
                            float("nan")
                            if row[name] is None or row[name] == ""
                            else row[name]
                            for name in output_names
                        )
                        for row in rows
                    )

                table = Table(rows=output_rows, names=output_names)

                if radius_deg is not None:
                    assert center is not None
                    inside_circle = GaiaQuery._on_sky_circle_mask(
                        ra=np.asarray(table["ra"], dtype=float),
                        dec=np.asarray(table["dec"], dtype=float),
                        center_ra=center[0],
                        center_dec=center[1],
                        radius_deg=radius_deg,
                    )
                    table = table[inside_circle]

                if len(table) > 0:
                    sort_values = np.asarray(
                        np.ma.filled(table[first_band.value], np.inf),  # type: ignore[index]
                        dtype=float,
                    )
                    table = table[np.argsort(sort_values, kind="stable")]
                return Table(table[:limit])

        return GaiaQuery._run_callable_with_timeout(
            func=_run_query,
            timeout=timeout,
            timeout_message=(
                "SQLite Gaia query timed out. "
                "You may want to increase the timeout or reduce the query size."
            ),
            on_timeout=connection.interrupt,
        )

    @staticmethod
    def _build_tap_bounds_where(
        ra_col: str,
        dec_col: str,
        bounds: tuple[float, float, float, float],
    ) -> str:
        """Build an ADQL bounds predicate with RA wrap handling."""
        ra_min, ra_max, dec_min, dec_max = bounds
        dec_clause = f"({dec_col} >= {dec_min} AND {dec_col} <= {dec_max})"
        if ra_min <= ra_max:
            ra_clause = f"({ra_col} >= {ra_min} AND {ra_col} <= {ra_max})"
        else:
            ra_clause = f"({ra_col} >= {ra_min} OR {ra_col} <= {ra_max})"
        return f"{dec_clause} AND {ra_clause}"

    @staticmethod
    def _select_sqlite_tables(
        cursor: sqlite3.Cursor,
        sqlite_source: GaiaSQLiteSource,
        center: tuple[float, float] | None,
        radius_deg: float | None,
        bounds: tuple[float, float, float, float] | None,
    ) -> tuple[list[str], str]:
        """Resolve SQLite tables for a query, including dec-ring shards."""
        all_tables = GaiaQuery._list_sqlite_tables(cursor)

        if sqlite_source.table in all_tables:
            return [sqlite_source.table], sqlite_source.table

        ring_tables: list[tuple[str, float, float]] = []
        for table in all_tables:
            dec_range = GaiaQuery._parse_dec_ring_table_name(table)
            if dec_range is not None:
                ring_tables.append((table, dec_range[0], dec_range[1]))

        if not ring_tables:
            raise ValueError(
                f"SQLite table '{sqlite_source.table}' not found and no sharded "
                "declination-ring tables were detected."
            )

        if bounds is not None:
            dec_min = bounds[2]
            dec_max = bounds[3]
        else:
            assert center is not None
            assert radius_deg is not None
            dec_min = max(-90.0, center[1] - radius_deg)
            dec_max = min(90.0, center[1] + radius_deg)

        ring_dec_min = {table: tmin for table, tmin, _ in ring_tables}

        selected = sorted(
            [
                table
                for table, ring_min, ring_max in ring_tables
                if ring_max >= dec_min and ring_min <= dec_max
            ],
            key=lambda name: ring_dec_min[name],
        )

        schema_table = min(ring_tables, key=lambda item: item[1])[0]
        return selected, schema_table

    @staticmethod
    def _list_sqlite_tables(cursor: sqlite3.Cursor) -> list[str]:
        """List user tables in a SQLite database."""
        return [
            str(row[0])
            for row in cursor.execute(
                "SELECT name FROM sqlite_master "
                "WHERE type='table' AND name NOT LIKE 'sqlite_%'"
            )
        ]

    @staticmethod
    def _parse_dec_ring_table_name(table_name: str) -> tuple[float, float] | None:
        """Parse sharded dec-ring table names such as '-25_-24'."""
        match = re.match(r"^(-?\d+(?:\.\d+)?)_(-?\d+(?:\.\d+)?)$", table_name)
        if match is None:
            return None
        dec_min = float(match.group(1))
        dec_max = float(match.group(2))
        if dec_min > dec_max:
            return dec_max, dec_min
        return dec_min, dec_max

    @staticmethod
    def _build_sqlite_bounds_where(
        bounds: tuple[float, float, float, float],
    ) -> tuple[str, tuple[float, ...]]:
        """Build SQL bounds predicate with parameters and RA wrap handling."""
        ra_min, ra_max, dec_min, dec_max = bounds
        if ra_min <= ra_max:
            return (
                '("dec" >= ? AND "dec" <= ?) AND ("ra" >= ? AND "ra" <= ?)',
                (dec_min, dec_max, ra_min, ra_max),
            )
        return (
            '("dec" >= ? AND "dec" <= ?) AND (("ra" >= ?) OR ("ra" <= ?))',
            (dec_min, dec_max, ra_min, ra_max),
        )

    @staticmethod
    def _build_sqlite_circle_prefilter_where(
        center: tuple[float, float],
        radius_deg: float,
    ) -> tuple[str, tuple[float, ...]]:
        """Build a rectangular SQL prefilter around a circle query."""
        center_ra, center_dec = center
        dec_min = max(-90.0, center_dec - radius_deg)
        dec_max = min(90.0, center_dec + radius_deg)

        cos_dec = max(abs(math.cos(math.radians(center_dec))), 1e-6)
        ra_half_width = radius_deg / cos_dec

        # When the circle is wide enough to span all RA values (e.g. near a pole),
        # skip the RA filter — only the dec bounds apply.
        if ra_half_width >= 180.0 or dec_min <= -90.0 or dec_max >= 90.0:
            return '"dec" >= ? AND "dec" <= ?', (dec_min, dec_max)

        ra_min = (center_ra - ra_half_width) % 360.0
        ra_max = (center_ra + ra_half_width) % 360.0

        return GaiaQuery._build_sqlite_bounds_where((ra_min, ra_max, dec_min, dec_max))

    @staticmethod
    def _on_sky_circle_mask(
        ra: np.ndarray,
        dec: np.ndarray,
        center_ra: float,
        center_dec: float,
        radius_deg: float,
    ) -> np.ndarray:
        """Return mask for points within an angular radius on the sphere."""
        ra_rad = np.radians(ra)
        dec_rad = np.radians(dec)
        center_ra_rad = math.radians(center_ra)
        center_dec_rad = math.radians(center_dec)

        cos_sep = np.sin(dec_rad) * math.sin(center_dec_rad) + np.cos(
            dec_rad
        ) * math.cos(center_dec_rad) * np.cos(ra_rad - center_ra_rad)
        cos_sep = np.clip(cos_sep, -1.0, 1.0)
        sep_deg = np.degrees(np.arccos(cos_sep))
        return sep_deg <= radius_deg

    @staticmethod
    def _quote_sql_identifier(identifier: str) -> str:
        """Quote a SQLite identifier safely."""
        if not identifier:
            raise ValueError("SQL identifier cannot be empty.")
        return '"' + identifier.replace('"', '""') + '"'

    @staticmethod
    def _run_callable_with_timeout(
        func: Callable[[], _T],
        timeout: float | None,
        timeout_message: str,
        on_timeout: Callable[[], None] | None = None,
    ) -> _T:
        """Run a callable with optional timeout using a daemon thread."""
        if timeout is None:
            return func()

        result: list[_T] = []
        exc: list = []

        def _run_safe():
            try:
                result.append(func())
            except Exception as e:
                exc.append(e)

        t = threading.Thread(target=_run_safe, daemon=True)
        t.start()
        t.join(timeout=timeout)
        if t.is_alive():
            if on_timeout is not None:
                on_timeout()
            raise TimeoutError(timeout_message)
        if exc:
            raise exc[0]
        return result[0]

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

        def _run() -> Table:
            tap = TapPlus(url=tap_source.value)
            job = tap.launch_job(query, **kwargs)
            return Table(job.get_results())  # type: ignore[arg-type]

        return GaiaQuery._run_callable_with_timeout(
            func=_run,
            timeout=timeout,
            timeout_message=(
                "Gaia query timed out. "
                "You may want to increase the timeout or reduce the query size. "
                f"Query was: {query}"
            ),
        )

    # Vega zero-points for all supported bands.
    # See _derive_band_properties for how these were obtained.
    _BAND_PROPS = {
        "G": {"dlam_lam": 0.7117, "flux_m0_Jy": 4031.5},
        "BP": {"dlam_lam": 0.5131, "flux_m0_Jy": 3683.21},
        "RP": {"dlam_lam": 0.3741, "flux_m0_Jy": 5040.41},
        "J": {"dlam_lam": 0.1312, "flux_m0_Jy": 1594.0},
        "H": {"dlam_lam": 0.151, "flux_m0_Jy": 1024.0},
        "KS": {"dlam_lam": 0.1214, "flux_m0_Jy": 666.8},
    }

    @staticmethod
    def _mag_to_photons(mags: np.ndarray, filter_band: Filters) -> np.ndarray:
        """Convert Vega magnitudes to photon flux in photons/s/m².

        Applies to all supported bands (Gaia G/BP/RP and 2MASS J/H/KS).
        Formula: Δλ/λ × F₀ × 1.51×10⁷ × 10^(−0.4 × mag),
        where F₀ is the Vega zero-point flux density in Jy.

        Parameters
        ----------
        mags : np.ndarray
            Vega magnitudes.
        filter_band : Filters
            Any supported filter band.

        Returns
        -------
        np.ndarray
            Flux in photons/s/m².
        """
        try:
            props = GaiaQuery._BAND_PROPS[filter_band.name]
        except KeyError:
            raise ValueError(
                f"_mag_to_photons: unsupported filter {filter_band}. "
                f"Supported: {list(GaiaQuery._BAND_PROPS)}."
            )
        Jy = 1.51e7  # [photons s^-1 m^-2 (Δλ/λ)^-1]
        return props["dlam_lam"] * props["flux_m0_Jy"] * Jy * 10 ** (-0.4 * mags)

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
        # Zero-fill missing proper motions so sources without pm data keep their
        # catalogue position rather than being silently NaN-corrupted.
        # np.where returns a plain ndarray (no units/mask), so use column
        # assignment (=) rather than in-place (+=) to avoid Astropy dropping
        # column metadata when mixing Column and ndarray operands.
        pmra = np.where(np.isnan(table["pmra"]), 0.0, table["pmra"])
        pmdec = np.where(np.isnan(table["pmdec"]), 0.0, table["pmdec"])
        table["ra"] = table["ra"] + years * pmra / 1000 / 3600  # type: ignore
        table["dec"] = table["dec"] + years * pmdec / 1000 / 3600  # type: ignore

        return table

    @staticmethod
    def _normalize_bands(
        filter_bands: Filters | str | Sequence[Filters | str],
    ) -> list["Filters"]:
        """Normalize filter_bands argument to a deduplicated list[Filters].

        Passing the string ``"all"`` (case-insensitive) expands to every
        available filter, equivalent to ``Filters.all()``.
        """
        if isinstance(filter_bands, str) and filter_bands.upper() == "ALL":
            return Filters.all()
        if isinstance(filter_bands, Filters | str):
            filter_bands = [filter_bands]
        if not isinstance(filter_bands, Sequence):
            raise ValueError(
                f"filter_bands must be a Filters, str, or sequence thereof, "
                f"got {type(filter_bands)}"
            )
        if not filter_bands:
            raise ValueError("At least one filter_band must be specified.")
        return list(dict.fromkeys(Filters.ensure_enum(b) for b in filter_bands))

    @staticmethod
    def _derive_band_properties():
        """
        Derives physical band properties for Gaia DR3 and 2MASS.

        This method serves as a reference for how the Vega zero-point fluxes and Δλ/λ
        values were obtained for the supported bands.

        Sources
        -------

        2MASS: Derived from Isophotal Fluxes
        Cohen et al. 2003 AJ 126 1090 (Table 2) | Bibcode: 2003AJ....126.1090C
        https://doi.org/10.1086/376474

        Gaia: Derived from Magnitude Zero Points
        Riello et al. 2021 A&A 649 A3 (Table 3) | Bibcode: 2021A&A...649A...3R
        https://doi.org/10.1051/0004-6361/202039587

        Example
        -------
        >>> from cabaret.queries import GaiaQuery
        >>> properties = GaiaQuery._derive_band_properties()
        >>> properties == GaiaQuery._BAND_PROPS
        True
        """
        AB_ZERO_POINT_JY = 3631.0

        # GAIA DR3 (Riello et al. 2021 A&A 649 A3, Table 3)
        # Format: {Band: [ZP_VEG, ZP_AB, FWHM_nm, LAM_0_nm]}
        gaia_table = {
            "G": [25.6874, 25.8010, 454.82, 639.07],
            "BP": [25.3385, 25.3540, 265.90, 518.26],
            "RP": [24.7479, 25.1040, 292.75, 782.51],
        }

        # 2MASS (Martin Cohen et al. 2003 AJ 126 1090, Table 2)
        # Format: {Band: [Flux_Jy, Bandwidth_um, Lam_Iso_um]}
        # Note: KS refers to the 2MASS "Short" K-band, which has a narrower
        # bandwidth and shorter pivot wavelength than standard Johnson K.
        tmass_table = {
            "J": [1594.0, 0.162, 1.235],
            "H": [1024.0, 0.251, 1.662],
            "KS": [666.8, 0.262, 2.159],
        }

        results = {}

        for band, (zp_v, zp_a, fwhm, l0) in gaia_table.items():
            # dlam_lam = FWHM / Pivot Wavelength
            dlam_lam = fwhm / l0
            # Flux density of Vega = 3631 * 10^(0.4 * (ZP_AB - ZP_VEG))
            flux_jy = AB_ZERO_POINT_JY * (10 ** (0.4 * (zp_a - zp_v)))

            results[band] = {
                "dlam_lam": round(dlam_lam, 4),
                "flux_m0_Jy": round(flux_jy, 2),
            }

        for band, (flux_jy, width, l_iso) in tmass_table.items():
            # dlam_lam = Width / Isophotal Wavelength
            dlam_lam = width / l_iso

            results[band] = {
                "dlam_lam": round(dlam_lam, 4),
                "flux_m0_Jy": round(flux_jy, 2),
            }

        return results
