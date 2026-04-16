import pytest
from astropy.coordinates import SkyCoord

from cabaret.queries import GaiaQuery, GaiaTAPSource

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
