import sqlite3

import numpy as np
import pytest
from astropy.coordinates import SkyCoord

from cabaret.queries import Filters, GaiaQuery, GaiaSQLiteSource, GaiaTAPSource
from cabaret.sources import Sources

from .utils import has_tap_source

skip_no_vizier = pytest.mark.skipif(
    not has_tap_source(GaiaTAPSource.VIZIER), reason="VizieR TAP unavailable"
)
skip_no_gaia = pytest.mark.skipif(
    not has_tap_source(GaiaTAPSource.GAIA), reason="Gaia Archive TAP unavailable"
)


@skip_no_vizier
def test_get_sources_default():
    center = SkyCoord(ra=10.68458, dec=41.26917, unit="deg")
    sources = GaiaQuery.get_sources(center, radius=0.05, limit=10, timeout=30)
    assert len(sources) <= 10
    assert sources is not None


@skip_no_vizier
def test_get_sources_vizier():
    center = SkyCoord(ra=10.68458, dec=41.26917, unit="deg")
    sources = GaiaQuery.get_sources(
        center, radius=0.05, limit=10, timeout=30, tap_source=GaiaTAPSource.VIZIER
    )
    assert len(sources) <= 10
    assert sources is not None


@skip_no_gaia
def test_get_sources_gaia():
    center = SkyCoord(ra=10.68458, dec=41.26917, unit="deg")
    sources = GaiaQuery.get_sources(
        center, radius=0.05, limit=10, timeout=30, tap_source=GaiaTAPSource.GAIA
    )
    assert len(sources) <= 10
    assert sources is not None


@pytest.mark.parametrize("tap_source", [GaiaTAPSource.VIZIER, GaiaTAPSource.GAIA])
def test_get_sources_gaia_ids(tap_source):
    """Both TAP backends resolve the same Gaia DR3 source ids."""
    if not has_tap_source(tap_source):
        pytest.skip(f"{tap_source.name} TAP unavailable")
    center = SkyCoord(ra=10.68458, dec=41.269, unit="deg")
    sources = GaiaQuery.get_sources(
        center, radius=0.05, limit=3, timeout=60, tap_source=tap_source
    )
    assert sources.gaia_ids.dtype == np.int64
    assert len(sources.gaia_ids) == len(sources)
    assert np.all(sources.gaia_ids > 0)
    # brightest-first ordering, so the same ids come back from either service
    assert sources.gaia_ids[0] == 381261910408440576


@skip_no_vizier
def test_get_sources_timeout():
    center = SkyCoord(ra=10.68458, dec=41.26917, unit="deg")
    with pytest.raises(TimeoutError):
        GaiaQuery.get_sources(center, radius=0.05, limit=10, timeout=0.0001)


@skip_no_vizier
def test_query_single_band():
    """query() still accepts a single band as a string."""
    center = SkyCoord(ra=10.68458, dec=41.26917, unit="deg")
    table = GaiaQuery.query(center, radius=0.05, filter_bands="G", limit=5, timeout=30)
    assert "phot_g_mean_mag" in table.colnames
    assert len(table) <= 5


@skip_no_vizier
def test_query_multi_band_tmass():
    """query() with two 2MASS bands returns both columns, no spurious Gaia flux."""
    center = SkyCoord(ra=10.68458, dec=41.26917, unit="deg")
    table = GaiaQuery.query(
        center,
        radius=0.05,
        filter_bands=[Filters.H, Filters.KS],
        limit=10,
        timeout=30,
    )
    assert "h_m" in table.colnames
    assert "ks_m" in table.colnames
    assert "phot_g_mean_mag" not in table.colnames
    assert len(table) <= 10
    # Ordered by H ASC (smallest magnitude = brightest first)
    assert table["h_m"][0] <= table["h_m"][-1]


@skip_no_vizier
def test_get_flux_table_multi_band():
    """get_flux_table() converts 2MASS mags to positive photon fluxes."""
    center = SkyCoord(ra=10.68458, dec=41.26917, unit="deg")
    table = GaiaQuery.get_flux_table(
        center,
        radius=0.05,
        filter_bands=[Filters.H, Filters.KS],
        limit=10,
        timeout=30,
    )
    assert "h_flux" in table.colnames
    assert "ks_flux" in table.colnames
    # Values should be positive photon fluxes, not raw magnitudes (~10–15)
    assert np.all(table["h_flux"] > 0)
    assert np.all(table["ks_flux"] > 0)
    assert np.median(table["h_flux"]) > 1000


@skip_no_vizier
def test_get_flux_table_gaia_bands():
    """get_flux_table() converts Gaia mags to positive photon fluxes."""
    center = SkyCoord(ra=10.68458, dec=41.26917, unit="deg")
    table = GaiaQuery.get_flux_table(
        center,
        radius=0.05,
        filter_bands=[Filters.G, Filters.BP, Filters.RP],
        limit=10,
        timeout=30,
    )
    assert "g_flux" in table.colnames
    assert "bp_flux" in table.colnames
    assert "rp_flux" in table.colnames
    assert "phot_g_mean_mag" not in table.colnames
    assert np.all(table["g_flux"] > 0)
    assert np.median(table["g_flux"]) > 1000


@skip_no_vizier
def test_get_flux_table_keep_mag():
    """keep_mag=True retains magnitude columns alongside flux columns."""
    center = SkyCoord(ra=10.68458, dec=41.26917, unit="deg")
    table = GaiaQuery.get_flux_table(
        center,
        radius=0.05,
        filter_bands=[Filters.G, Filters.H],
        limit=10,
        timeout=30,
        keep_mag=True,
    )
    assert "g_flux" in table.colnames
    assert "h_flux" in table.colnames
    assert "phot_g_mean_mag" in table.colnames
    assert "h_m" in table.colnames


@skip_no_vizier
def test_get_flux_table_allow_nulls():
    """allow_nulls=True lets rows with missing band values through."""
    center = SkyCoord(ra=10.68458, dec=41.26917, unit="deg")
    # Query Gaia-only: without allow_nulls the result has no NaNs
    strict = GaiaQuery.get_flux_table(
        center,
        radius=0.05,
        filter_bands="G",
        limit=50,
        timeout=30,
        allow_nulls=False,
    )
    lenient = GaiaQuery.get_flux_table(
        center,
        radius=0.05,
        filter_bands="G",
        limit=50,
        timeout=30,
        allow_nulls=True,
    )
    # Allowing nulls should return at least as many rows
    assert len(lenient) >= len(strict)


@skip_no_vizier
def test_get_flux_table_filter_bands_all():
    """filter_bands='all' returns flux columns for every supported band."""
    center = SkyCoord(ra=10.68458, dec=41.26917, unit="deg")
    table = GaiaQuery.get_flux_table(
        center,
        radius=0.05,
        filter_bands="all",
        limit=10,
        timeout=30,
    )
    for col in ("g_flux", "bp_flux", "rp_flux", "j_flux", "h_flux", "ks_flux"):
        assert col in table.colnames


def test_drop_nan_fluxes_removes_nans():
    """drop_nan_fluxes() removes sources with NaN flux values."""
    sources = Sources.from_arrays(
        ra=np.array([1.0, 2.0, 3.0]),
        dec=np.array([0.0, 0.0, 0.0]),
        fluxes=np.array([100.0, np.nan, 200.0]),
    )
    clean = sources.drop_nan_fluxes()
    assert len(clean) == 2
    assert not np.any(np.isnan(clean.fluxes))


def test_drop_nan_fluxes_no_nans_returns_self():
    """drop_nan_fluxes() returns self when no NaNs are present."""
    sources = Sources.from_arrays(
        ra=np.array([1.0, 2.0]),
        dec=np.array([0.0, 0.0]),
        fluxes=np.array([100.0, 200.0]),
    )
    assert sources.drop_nan_fluxes() is sources


def _make_sqlite_catalog(path: str):
    conn = sqlite3.connect(path)
    try:
        conn.execute(
            """
            CREATE TABLE gaia_sources (
                ra REAL,
                dec REAL,
                pmra REAL,
                pmdec REAL,
                phot_g_mean_mag REAL,
                h_m REAL
            )
            """
        )
        rows = [
            (359.9, 0.0, 1.0, 2.0, 10.0, 9.0),
            (0.2, 0.1, 0.0, 0.0, 11.0, 10.0),
            (15.0, 2.0, None, None, 13.0, 12.0),
            (45.0, 20.0, -1.0, 0.5, 12.0, None),
        ]
        conn.executemany(
            "INSERT INTO gaia_sources (ra, dec, pmra, pmdec, phot_g_mean_mag, h_m) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            rows,
        )
        conn.commit()
    finally:
        conn.close()


def _make_sharded_sqlite_catalog(path: str):
    conn = sqlite3.connect(path)
    try:
        conn.execute(
            'CREATE TABLE "-1_0" ('
            "ra REAL, dec REAL, pmra REAL, pmdec REAL, phot_g_mean_mag REAL"
            ")"
        )
        conn.execute(
            'CREATE TABLE "0_1" ('
            "ra REAL, dec REAL, pmra REAL, pmdec REAL, phot_g_mean_mag REAL"
            ")"
        )
        conn.executemany(
            'INSERT INTO "-1_0" (ra, dec, pmra, pmdec, phot_g_mean_mag) '
            "VALUES (?, ?, ?, ?, ?)",
            [
                (20.0, -0.6, 0.0, 0.0, 13.0),
                (10.0, -0.2, 0.0, 0.0, 12.0),
            ],
        )
        conn.executemany(
            'INSERT INTO "0_1" (ra, dec, pmra, pmdec, phot_g_mean_mag) '
            "VALUES (?, ?, ?, ?, ?)",
            [
                (40.0, 0.2, 0.0, 0.0, 11.0),
                (30.0, 0.6, 0.0, 0.0, 10.0),
            ],
        )
        conn.commit()
    finally:
        conn.close()


def test_query_sqlite_by_radius(tmp_path):
    db_path = tmp_path / "gaia_offline.db"
    _make_sqlite_catalog(str(db_path))

    table = GaiaQuery.query(
        center=(0.0, 0.0),
        radius=1.0,
        filter_bands="G",
        limit=10,
        tap_source=GaiaSQLiteSource(database=str(db_path)),
    )

    assert len(table) == 2
    assert "phot_g_mean_mag" in table.colnames
    assert np.all(table["phot_g_mean_mag"] >= 10.0)


def test_query_sqlite_by_bounds_wraparound(tmp_path):
    db_path = tmp_path / "gaia_offline.db"
    _make_sqlite_catalog(str(db_path))

    table = GaiaQuery.query(
        bounds=(359.5, 0.5, -1.0, 1.0),
        filter_bands="G",
        limit=10,
        tap_source=f"sqlite:///{db_path}",
    )

    assert len(table) == 2
    assert np.all((table["ra"] > 359.5) | (table["ra"] < 0.5))


def test_query_sqlite_filter_bands_all_partial_schema(tmp_path):
    db_path = tmp_path / "gaia_offline.db"
    _make_sqlite_catalog(str(db_path))

    table = GaiaQuery.query(
        center=(15.0, 2.0),
        radius=2.0,
        filter_bands="all",
        limit=10,
        tap_source=GaiaSQLiteSource(database=str(db_path)),
    )

    assert "phot_g_mean_mag" in table.colnames
    assert "h_m" in table.colnames
    assert "phot_bp_mean_mag" not in table.colnames


def test_query_sqlite_missing_explicit_band_raises(tmp_path):
    db_path = tmp_path / "gaia_offline.db"
    _make_sqlite_catalog(str(db_path))

    with pytest.raises(ValueError, match="missing requested band columns"):
        GaiaQuery.query(
            center=(15.0, 2.0),
            radius=2.0,
            filter_bands=[Filters.BP],
            limit=10,
            tap_source=GaiaSQLiteSource(database=str(db_path)),
        )


def test_get_sources_sqlite_bounds(tmp_path):
    db_path = tmp_path / "gaia_offline.db"
    _make_sqlite_catalog(str(db_path))

    sources = GaiaQuery.get_sources(
        bounds=(359.5, 0.5, -1.0, 1.0),
        filter_band=Filters.G,
        limit=10,
        tap_source=GaiaSQLiteSource(database=str(db_path)),
    )

    assert len(sources) == 2
    assert np.all(sources.fluxes > 0)
    # SQLite catalogs carry no Gaia source id
    assert np.all(sources.gaia_ids == -1)


def test_query_sqlite_sharded_auto_detect_bounds(tmp_path):
    db_path = tmp_path / "gaia_sharded.db"
    _make_sharded_sqlite_catalog(str(db_path))

    table = GaiaQuery.query(
        bounds=(0.0, 60.0, -0.3, 0.3),
        filter_bands=Filters.G,
        limit=10,
        tap_source=GaiaSQLiteSource(database=str(db_path)),
    )

    assert len(table) == 2
    assert np.all((-0.3 <= table["dec"]) & (table["dec"] <= 0.3))


def _make_sky_catalog(path: str, rows: list[tuple]):
    """rows: (ra, dec, phot_g_mean_mag)"""
    conn = sqlite3.connect(path)
    try:
        conn.execute(
            "CREATE TABLE gaia_sources "
            "(ra REAL, dec REAL, pmra REAL, pmdec REAL, phot_g_mean_mag REAL)"
        )
        conn.executemany(
            "INSERT INTO gaia_sources VALUES (?, ?, NULL, NULL, ?)",
            rows,
        )
        conn.commit()
    finally:
        conn.close()


def _query_ra_set(db_path: str, center: tuple, radius: float) -> set[float]:
    table = GaiaQuery.query(
        center=center,
        radius=radius,
        filter_bands="G",
        limit=100,
        tap_source=GaiaSQLiteSource(database=db_path),
    )
    return set(np.round(table["ra"], 1).tolist())


def test_sqlite_radius_ra_wrap_straddles_zero(tmp_path):
    """Stars straddling RA=0 must both appear when the circle crosses the wrap."""
    db = str(tmp_path / "cat.db")
    _make_sky_catalog(
        db,
        [
            (359.5, 0.0, 10.0),  # 0.5° west of RA=0 — inside
            (0.5, 0.0, 11.0),  # 0.5° east of RA=0 — inside
            (5.0, 0.0, 12.0),  # 5° east — outside
            (355.0, 0.0, 13.0),  # 5° west — outside
        ],
    )
    assert _query_ra_set(db, center=(0.0, 0.0), radius=2.0) == {359.5, 0.5}


def test_sqlite_radius_ra_wrap_center_near_360(tmp_path):
    """Circle centered at RA=359 must capture a star that wraps past RA=0."""
    db = str(tmp_path / "cat.db")
    _make_sky_catalog(
        db,
        [
            (1.0, 0.0, 10.0),  # 2° east via wrap — inside
            (357.0, 0.0, 11.0),  # 2° west — inside
            (5.0, 0.0, 12.0),  # 6° east via wrap — outside
        ],
    )
    assert _query_ra_set(db, center=(359.0, 0.0), radius=3.0) == {1.0, 357.0}


def test_sqlite_radius_near_north_pole_all_ra(tmp_path):
    """Stars near the pole at opposite RA must be found."""
    db = str(tmp_path / "cat.db")
    # (RA=0, Dec=89) → (RA=180, Dec=89.5): great-circle separation ≈ 1.5° < 2°
    _make_sky_catalog(
        db,
        [
            (0.0, 89.5, 10.0),  # same RA as center — inside
            (180.0, 89.5, 11.0),  # opposite RA, ~1.5° away — inside
            (90.0, 87.5, 12.0),  # ~2.5° from center — outside
        ],
    )
    assert _query_ra_set(db, center=(0.0, 89.0), radius=2.0) == {0.0, 180.0}


def test_sqlite_radius_near_south_pole_all_ra(tmp_path):
    """Same as north-pole test, mirrored to Dec=-89."""
    db = str(tmp_path / "cat.db")
    _make_sky_catalog(
        db,
        [
            (0.0, -89.5, 10.0),
            (180.0, -89.5, 11.0),
            (90.0, -87.5, 12.0),
        ],
    )
    assert _query_ra_set(db, center=(0.0, -89.0), radius=2.0) == {0.0, 180.0}


def test_query_sqlite_sharded_auto_detect_radius(tmp_path):
    db_path = tmp_path / "gaia_sharded.db"
    _make_sharded_sqlite_catalog(str(db_path))

    table = GaiaQuery.query(
        center=(30.0, 0.0),
        radius=30.0,
        filter_bands="G",
        limit=10,
        tap_source=f"sqlite:///{db_path}",
    )

    assert len(table) == 4
    assert np.all(table["dec"] >= -1.0)
    assert np.all(table["dec"] <= 1.0)
