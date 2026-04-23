import numpy as np
import pytest

from cabaret.image import _sample_trail, generate_star_image

RNG = np.random.default_rng(42)


# --- _sample_trail ---


def test_sample_trail_shape():
    result = _sample_trail(50.0, 60.0, 10.0, 5.0, 20, 0.0, RNG)
    assert result.shape == (2, 20)


def test_sample_trail_linear_no_jitter():
    result = _sample_trail(0.0, 0.0, 10.0, 0.0, 5, 0.0, RNG)
    assert np.allclose(result[0], [0.0, 2.5, 5.0, 7.5, 10.0])
    assert np.allclose(result[1], 0.0)


def test_sample_trail_end_position():
    x0, y0, dx, dy = 3.0, 7.0, 12.0, -4.0
    result = _sample_trail(x0, y0, dx, dy, 10, 0.0, RNG)
    assert np.isclose(result[0, 0], x0)
    assert np.isclose(result[1, 0], y0)
    assert np.isclose(result[0, -1], x0 + dx)
    assert np.isclose(result[1, -1], y0 + dy)


def test_sample_trail_jitter_changes_output():
    rng1 = np.random.default_rng(0)
    rng2 = np.random.default_rng(0)
    no_jitter = _sample_trail(50.0, 50.0, 10.0, 0.0, 20, 0.0, rng1)
    with_jitter = _sample_trail(50.0, 50.0, 10.0, 0.0, 20, 1.0, rng2)
    assert not np.allclose(no_jitter, with_jitter)


# --- generate_star_image ---


FRAME = (200, 200)


def test_zero_drift_regression():
    """Zero-drift path produces same result as calling without drift_pixels."""
    pos = np.array([[100.0], [100.0]])
    fluxes = [1e6]
    rng = np.random.default_rng(7)
    kwargs = dict(
        fluxes=fluxes,
        FWHM=3.0,
        frame_size=FRAME,
        rng=rng,
        telescope_aperture=1.0,
        site_elevation=0.0,
        exp_time=1.0,
    )
    img_no_drift = generate_star_image(pos=pos, **kwargs)
    rng2 = np.random.default_rng(7)
    img_zero_drift = generate_star_image(
        pos=pos, drift_pixels=np.zeros((2, 1)), **{**kwargs, "rng": rng2}
    )
    assert np.allclose(img_no_drift, img_zero_drift)


@pytest.mark.parametrize("dx,dy", [(0, 0), (10, 0), (0, 20), (10, 20)])
def test_flux_conservation(dx, dy):
    """Total flux in the image should approximately equal the input flux."""
    flux = 1e7
    pos = np.array([[100.0], [100.0]])
    drift = np.array([[float(dx)], [float(dy)]])
    rng = np.random.default_rng(1)
    img = generate_star_image(
        pos=pos,
        fluxes=[flux],
        FWHM=3.0,
        frame_size=FRAME,
        rng=rng,
        telescope_aperture=1.0,
        site_elevation=0.0,
        exp_time=1.0,
        drift_pixels=drift,
        n_trail_samples=50,
    )
    # Allow 2% tolerance (Poisson noise on large flux)
    assert abs(img.sum() - flux) / flux < 0.02


def test_drifting_star_spans_expected_range():
    """A star drifting 30px in x should illuminate pixels spread over ~30px."""
    pos = np.array([[85.0], [100.0]])
    drift = np.array([[30.0], [0.0]])
    rng = np.random.default_rng(2)
    img = generate_star_image(
        pos=pos,
        fluxes=[1e7],
        FWHM=2.0,
        frame_size=FRAME,
        rng=rng,
        telescope_aperture=1.0,
        site_elevation=0.0,
        exp_time=1.0,
        drift_pixels=drift,
        n_trail_samples=50,
    )
    # Find columns with significant flux
    col_sums = img.sum(axis=0)
    bright_cols = np.where(col_sums > col_sums.max() * 0.01)[0]
    span = bright_cols.max() - bright_cols.min()
    assert span >= 25, f"Expected trail span >= 25px, got {span}"


def test_off_frame_trail_start_renders_end():
    """Star starts off-frame left but drifts on-frame — should be rendered."""
    pos = np.array([[-5.0], [100.0]])
    drift = np.array([[30.0], [0.0]])
    rng = np.random.default_rng(3)
    img = generate_star_image(
        pos=pos,
        fluxes=[1e6],
        FWHM=3.0,
        frame_size=FRAME,
        rng=rng,
        telescope_aperture=1.0,
        site_elevation=0.0,
        exp_time=1.0,
        drift_pixels=drift,
        n_trail_samples=50,
    )
    assert img.sum() > 0, "Expected some flux on-frame from partially visible trail"


def test_crossing_trail_both_endpoints_off_frame():
    """Trail crossing the frame with both endpoints outside must produce flux.

    Regression for the endpoint-only on_camera_mask that dropped such trails.
    Start x=-20, end x=220 on a 200-wide frame: both off-frame, trail crosses fully.
    """
    pos = np.array([[-20.0], [100.0]])
    drift = np.array([[240.0], [0.0]])
    rng = np.random.default_rng(4)
    img = generate_star_image(
        pos=pos,
        fluxes=[1e6],
        FWHM=3.0,
        frame_size=FRAME,
        rng=rng,
        telescope_aperture=1.0,
        site_elevation=0.0,
        exp_time=1.0,
        drift_pixels=drift,
        n_trail_samples=50,
    )
    assert img.sum() > 0, (
        "Trail crossing the frame must render even when both endpoints are off-frame"
    )


def test_adaptive_sampling_short_trail_conserves_flux():
    """Short trail (length < FWHM) reduces to 1 sample but must conserve flux."""
    flux = 1e7
    pos = np.array([[100.0], [100.0]])
    drift = np.array([[1.0], [0.0]])  # 1 px drift with FWHM=3 → ceil(1/3*2)=1 sample
    rng = np.random.default_rng(5)
    img = generate_star_image(
        pos=pos,
        fluxes=[flux],
        FWHM=3.0,
        frame_size=FRAME,
        rng=rng,
        telescope_aperture=1.0,
        site_elevation=0.0,
        exp_time=1.0,
        drift_pixels=drift,
    )
    assert abs(img.sum() - flux) / flux < 0.02
