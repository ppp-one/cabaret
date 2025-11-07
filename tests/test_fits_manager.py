from datetime import datetime

import numpy as np
import pytest

from cabaret.camera import Camera
from cabaret.fits_manager import FITSManager
from cabaret.observatory import Observatory


@pytest.fixture
def observatory():
    return Observatory(
        name="TestObs",
        camera=Camera(width=10, height=10),
    )


def test_get_header_from_observatory(observatory):
    header = FITSManager.get_header_from_observatory(observatory)
    assert header["OBSNAME"] == "TestObs"
    assert "INSTRUME" in header
    assert "FOCUSER" in header
    assert "TELESCOP" in header
    assert "SITE" in header
    assert "CABARET" in header


def test_add_image_info_to_header():
    from astropy.io import fits

    header = fits.Header()
    dateobs = datetime(2024, 1, 1, 12, 0, 0)
    FITSManager.add_image_info_to_header(
        header, exp_time=10, ra=123.4, dec=56.7, dateobs=dateobs
    )
    assert header["EXPTIME"] == 10.0
    assert header["RA"] == 123.4
    assert header["DEC"] == 56.7
    assert header["DATE-OBS"].startswith("2024-01-01T12:00:00")


def test_to_hdu_list_and_save(tmp_path, observatory):
    image = np.ones((10, 10), dtype=np.uint16)
    ra, dec = 10.0, 20.0
    hdul = FITSManager.to_hdu_list(
        image=image,
        observatory=observatory,
        exp_time=30,
        ra=ra,
        dec=dec,
        dateobs=datetime(2024, 1, 1, 12, 0, 0),
        user_header={"TESTKEY": ("testval", "test comment")},
    )
    assert hdul[0].header["EXPTIME"] == 30.0
    assert hdul[0].header["RA"] == 10.0
    assert hdul[0].header["DEC"] == 20.0
    assert hdul[0].header["TESTKEY"] == "testval"

    # Test saving to file
    file_path = tmp_path / "test_image.fits"
    FITSManager.save(
        image=image,
        file_path=str(file_path),
        observatory=observatory,
        exp_time=30,
        ra=ra,
        dec=dec,
        dateobs=datetime(2024, 1, 1, 12, 0, 0),
        user_header={"TESTKEY": ("testval", "test comment")},
        overwrite=True,
    )
    from astropy.io import fits as afits

    with afits.open(file_path) as hdul2:
        assert hdul2[0].header["EXPTIME"] == 30.0
        assert hdul2[0].header["RA"] == 10.0
        assert hdul2[0].header["DEC"] == 20.0
        assert hdul2[0].header["TESTKEY"] == "testval"
        np.testing.assert_array_equal(hdul2[0].data, image)
