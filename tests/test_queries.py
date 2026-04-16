import numpy as np
import pytest
from astropy.coordinates import SkyCoord

from cabaret.queries import Filters, GaiaQuery, GaiaTAPSource

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


@skip_no_vizier
def test_get_sources_timeout():
    center = SkyCoord(ra=10.68458, dec=41.26917, unit="deg")
    with pytest.raises(TimeoutError):
        GaiaQuery.get_sources(center, radius=0.05, limit=10, timeout=0.0001)


@skip_no_vizier
def test_query_single_band_backward_compat():
    """query() still accepts a single band as a string."""
    center = SkyCoord(ra=10.68458, dec=41.26917, unit="deg")
    table = GaiaQuery.query(center, radius=0.05, filter_bands="G", limit=5, timeout=30)
    assert "phot_g_mean_flux" in table.colnames
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
    assert "phot_g_mean_flux" not in table.colnames
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
