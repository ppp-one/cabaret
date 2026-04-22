import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

from cabaret.sources import Sources

COORDS = np.array([[10.64, 41.26], [10.68, 41.22]])
FLUXES = np.array([169435.6, 52203.9])


@pytest.fixture
def example_sources():
    return Sources(SkyCoord(COORDS, unit="deg"), FLUXES)


def test_sources_init_and_len(example_sources):
    sources = example_sources
    assert len(sources) == 2
    assert np.allclose(sources.fluxes, FLUXES)
    assert np.allclose(sources.coords.ra.deg, COORDS[:, 0])
    assert np.allclose(sources.coords.dec.deg, COORDS[:, 1])


def test_sources_ra_dec_properties(example_sources):
    sources = example_sources
    assert np.allclose(sources.ra.deg, COORDS[:, 0])
    assert np.allclose(sources.dec.deg, COORDS[:, 1])


def test_sources_to_pixel(example_sources):
    sources = example_sources
    wcs = WCS(naxis=2)
    wcs.wcs.crval = [10.65, 41.25]
    wcs.wcs.crpix = [100, 100]
    wcs.wcs.cdelt = [-0.0002777778, 0.0002777778]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    pixels = sources.to_pixel(wcs)
    assert pixels.shape == (2, 2)


def test_sources_invalid_init():
    # coords is not of type SkyCoord
    with pytest.raises(ValueError):
        Sources(COORDS, FLUXES)

    # Length mismatch
    with pytest.raises(ValueError):
        Sources(SkyCoord(COORDS, unit="deg"), FLUXES[:-1])


def test_sources_from_array():
    sources = Sources.from_arrays(COORDS[:, 0], COORDS[:, 1], FLUXES)
    assert np.allclose(sources.ra.deg, COORDS[:, 0])
    assert np.allclose(sources.dec.deg, COORDS[:, 1])
    assert np.allclose(sources.fluxes, FLUXES)


def test_center():
    coords = SkyCoord(ra=[179.999, 180.001], dec=[10.68, 10.689], unit="deg")
    sources = Sources(coords, fluxes=np.array([169_435.6, 92_203.9]))

    ra_center, dec_center = sources.center
    assert np.isclose(ra_center, 180.0)
    assert np.isclose(dec_center, 10.6845)


def test_center_with_wrap():
    coords = SkyCoord(ra=[359.999, 0.001], dec=[10.68, 10.689], unit="deg")
    sources = Sources(coords, fluxes=np.array([169_435.6, 92_203.9]))

    ra_center, _ = sources.center
    assert np.isclose(ra_center, 0.0)


# --- rates ---


def test_rates_default_to_zeros(example_sources):
    assert example_sources.ra_rates.shape == (2,)
    assert example_sources.dec_rates.shape == (2,)
    assert np.all(example_sources.ra_rates == 0.0)
    assert np.all(example_sources.dec_rates == 0.0)


def test_rates_stored_correctly():
    ra_rates = np.array([1.0, -0.5])
    dec_rates = np.array([0.2, 0.0])
    sources = Sources(
        SkyCoord(COORDS, unit="deg"), FLUXES, ra_rates=ra_rates, dec_rates=dec_rates
    )
    assert np.allclose(sources.ra_rates, ra_rates)
    assert np.allclose(sources.dec_rates, dec_rates)


def test_rates_wrong_shape_raises():
    with pytest.raises(ValueError):
        Sources(SkyCoord(COORDS, unit="deg"), FLUXES, ra_rates=np.array([1.0]))
    with pytest.raises(ValueError):
        Sources(SkyCoord(COORDS, unit="deg"), FLUXES, dec_rates=np.array([1.0]))


def test_getitem_slices_rates():
    ra_rates = np.array([1.0, 2.0])
    dec_rates = np.array([0.1, 0.2])
    sources = Sources(
        SkyCoord(COORDS, unit="deg"), FLUXES, ra_rates=ra_rates, dec_rates=dec_rates
    )
    sliced = sources[0]
    assert sliced.ra_rates.shape == (1,)
    assert np.isclose(sliced.ra_rates[0], 1.0)
    assert np.isclose(sliced.dec_rates[0], 0.1)

    masked = sources[np.array([False, True])]
    assert np.isclose(masked.ra_rates[0], 2.0)


def test_from_arrays_with_rates():
    ra_rates = np.array([0.5, -0.3])
    dec_rates = np.array([0.0, 0.1])
    sources = Sources.from_arrays(
        COORDS[:, 0], COORDS[:, 1], FLUXES, ra_rates=ra_rates, dec_rates=dec_rates
    )
    assert np.allclose(sources.ra_rates, ra_rates)
    assert np.allclose(sources.dec_rates, dec_rates)


def test_drop_nan_fluxes_preserves_rates():
    fluxes = np.array([1.0, float("nan"), 3.0])
    ra_rates = np.array([1.0, 2.0, 3.0])
    coords = SkyCoord(ra=[10.0, 10.1, 10.2], dec=[41.0, 41.1, 41.2], unit="deg")
    sources = Sources(coords, fluxes, ra_rates=ra_rates)
    cleaned = sources.drop_nan_fluxes()
    assert len(cleaned) == 2
    assert np.allclose(cleaned.ra_rates, [1.0, 3.0])


def test_concat_merges_rates():
    ra_rates_a = np.array([1.0, 2.0])
    sources_a = Sources(SkyCoord(COORDS, unit="deg"), FLUXES, ra_rates=ra_rates_a)
    sources_b = Sources.get_test_sources()  # zero rates
    merged = Sources.concat(sources_a, sources_b)
    assert len(merged) == len(sources_a) + len(sources_b)
    assert np.allclose(merged.ra_rates[:2], ra_rates_a)
    assert np.all(merged.ra_rates[2:] == 0.0)


def test_add_operator_delegates_to_concat(example_sources):
    combined = example_sources + example_sources
    assert len(combined) == 2 * len(example_sources)
    assert np.allclose(combined.fluxes[:2], example_sources.fluxes)


def test_rates_quantity_arcsec_per_hour():
    rates = np.array([3600.0, -3600.0]) * u.arcsec / u.hour
    sources = Sources(SkyCoord(COORDS, unit="deg"), FLUXES, ra_rates=rates)
    assert np.allclose(sources.ra_rates, [1.0, -1.0])


def test_rates_quantity_arcsec_per_second():
    rates = np.array([1.0, 0.5]) * u.arcsec / u.s
    sources = Sources(SkyCoord(COORDS, unit="deg"), FLUXES, ra_rates=rates)
    assert np.allclose(sources.ra_rates, [1.0, 0.5])


def test_rates_quantity_incompatible_unit_raises():
    with pytest.raises(u.UnitConversionError):
        Sources(
            SkyCoord(COORDS, unit="deg"),
            FLUXES,
            ra_rates=np.array([1.0, 2.0]) * u.meter,
        )
