from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


class FITSManager:
    """A class to manage FITS file operations, including creating headers with
    observatory metadata.

    Examples
    --------
    >>> from cabaret.fits_manager import FITSManager
    >>> from cabaret import Observatory, Sources
    >>> observatory = Observatory(name="My Observatory")
    >>> sources = Sources.get_test_sources()
    >>> ra, dec = sources.ra.deg.mean(), sources.dec.deg.mean()
    >>> image = observatory.generate_image(
    ...     ra=ra, dec=dec, exp_time=30, sources=sources, seed=0
    ... )
    >>> hdu_list = FITSManager.get_fits_from_array(
    ...     image=image, observatory=observatory, exp_time=30, ra=ra, dec=dec
    ... )
    >>> hdu_list[0].header['DEC']

    """

    @staticmethod
    def save(
        observatory,
        file_path: str | Path,
        hdu_list: fits.HDUList | None = None,
        image: np.ndarray | None = None,
        user_header: dict[str, Any] | fits.Header | None = None,
        exp_time: float | None = None,
        ra: float | None = None,
        dec: float | None = None,
        dateobs: datetime | None = None,
        wcs: WCS | None = None,
        overwrite: bool = True,
    ) -> fits.HDUList:
        """
        Save a numpy array image to a FITS file using observatory metadata.

        Parameters
        ----------
        observatory : Observatory
            The observatory instance for header metadata.
        file_path : str
            The path to the FITS file to write.
        hdu_list : fits.HDUList, optional
            An existing HDUList to save. If provided, image and header parameters
            are ignored.
        image : numpy.ndarray
            The image data to save.
        user_header : dict or fits.Header, optional
            Additional header keywords to add.
        exp_time : float, optional
            Exposure time in seconds. Passed to the header.
        ra : float, optional
            Right ascension of the image center in degrees. Passed to the header.
        dec : float, optional
            Declination of the image center in degrees. Passed to the header.
        dateobs : datetime, optional
            The observation date and time. Passed to the header.
        wcs : WCS, optional
            World Coordinate System information for the image. Passed to the header.
        overwrite : bool, optional
            Whether to overwrite existing file (default: True).
        """
        if hdu_list is None:
            hdu_list = FITSManager.to_hdu_list(
                image=image,
                observatory=observatory,
                exp_time=exp_time,
                ra=ra,
                dec=dec,
                dateobs=dateobs,
                wcs=wcs,
                user_header=user_header,
            )
        elif not isinstance(hdu_list, fits.HDUList):
            raise TypeError("hdu_list must be an instance of fits.HDUList.")

        hdu_list.writeto(file_path, overwrite=overwrite)

        return hdu_list

    @staticmethod
    def to_hdu_list(
        observatory,
        image,
        exp_time: float | None = None,
        ra: float | None = None,
        dec: float | None = None,
        dateobs: datetime | None = None,
        wcs: WCS | None = None,
        user_header: dict[str, Any] | fits.Header | None = None,
    ) -> fits.HDUList:
        """
        Create a FITS HDUList from a numpy array, optionally with a header.

        Parameters
        ----------
        observatory : Observatory
            The observatory instance for header metadata.
        image : numpy.ndarray
            The image data to convert to FITS.
        exp_time : float, optional
            Exposure time in seconds.
        ra : float, optional
            Right ascension of the image center in degrees.
        dec : float, optional
            Declination of the image center in degrees.
        dateobs : datetime, optional
            The observation date and time.
        wcs : WCS, optional
            World Coordinate System information for the image.
        user_header : dict or fits.Header, optional
            Additional header keywords to add.

        Returns
        -------
        fits.HDUList
            The created FITS HDUList.
        """

        header = FITSManager.get_header_from_observatory(
            observatory,
            user_header=user_header,
        )
        FITSManager.add_image_info_to_header(
            header, exp_time=exp_time, ra=ra, dec=dec, dateobs=dateobs, wcs=wcs
        )
        hdu = fits.PrimaryHDU(data=image, header=header)
        hdul = fits.HDUList([hdu])
        return hdul

    @staticmethod
    def add_image_info_to_header(
        header: fits.Header,
        exp_time: float | None = None,
        ra: float | None = None,
        dec: float | None = None,
        dateobs: datetime | None = None,
        wcs: WCS | None = None,
    ):
        """Add image-specific info to a FITS header."""
        if exp_time is not None:
            header["EXPTIME"] = (float(exp_time), "Exposure time in seconds")
        if ra is not None:
            header["RA"] = (float(ra), "Right ascension of image center [deg]")
        if dec is not None:
            header["DEC"] = (float(dec), "Declination of image center [deg]")
        if dateobs is None:
            dateobs = datetime.now()
        if wcs is not None:
            wcs_header = wcs.to_header()
            for card in wcs_header.cards:
                header[card.keyword] = (card.value, card.comment)
        header["DATE-OBS"] = (dateobs.isoformat(), "UTC datetime start of exposure")

    @staticmethod
    def get_header_from_observatory(
        observatory,
        user_header: dict[str, Any] | fits.Header | None = None,
    ) -> "fits.Header":
        """
        Create a FITS header populated with metadata from an Observatory instance.

        Parameters
        ----------
        observatory : Observatory
            The observatory instance to extract metadata from.
        extra : dict, optional
            Additional header keywords to add.

        Returns
        -------
        fits.Header
            The populated FITS header.
        """
        header = fits.Header()
        header["OBSNAME"] = (observatory.name, "Observatory name")

        # Camera-related fields
        camera = observatory.camera
        header["INSTRUME"] = (getattr(camera, "name", "Camera"), "Instrument used")
        header["XBINNING"] = (
            getattr(camera, "bin_x", 1),
            "Binning level along the X-axis",
        )
        header["YBINNING"] = (
            getattr(camera, "bin_y", 1),
            "Binning level along the Y-axis",
        )
        header["XPIXSZ"] = (
            getattr(camera, "pitch", None),
            "Pixel Width in microns (after binning)",
        )
        header["YPIXSZ"] = (
            getattr(camera, "pitch", None),
            "Pixel Height in microns (after binning)",
        )
        header["CAM-DNAM"] = (
            getattr(camera, "name", "cabaret_camera"),
            "Short name of Camera driver",
        )

        header["FOCUSER"] = (
            getattr(observatory.focuser, "name", "cabaret_focuser"),
            "Focuser name",
        )
        header["FOCUSPOS"] = (observatory.focuser.position, "Focuser position")
        header["TELESCOP"] = (
            getattr(observatory.telescope, "name", "cabaret_telescope"),
            "Telescope name",
        )
        telescope = observatory.telescope
        header["FOCALLEN"] = (
            telescope.focal_length,
            "[m] Focal length of telescope",
        )
        header["APTDIA"] = (
            telescope.diameter,
            "[m] Aperture diameter of telescope",
        )
        header["APTAREA"] = (
            telescope.collecting_area,
            "[m^2] Aperture area of telescope",
        )
        header["SITE"] = (
            getattr(observatory.site, "name", "cabaret_site"),
            "Site name",
        )
        # Add more fields as needed from the observatory/camera/telescope/site

        if user_header:
            # Accept both dict and fits.Header
            if isinstance(user_header, dict):
                for k, v in user_header.items():
                    header[k] = v
            elif isinstance(user_header, fits.Header):
                for card in user_header.cards:
                    header[card.keyword] = (card.value, card.comment)
            else:
                raise TypeError("user_header must be a dict or fits.Header")

        try:
            from importlib.metadata import version

            header["CABARET"] = (version("cabaret"), "Version of Cabaret")
        except Exception:
            header["CABARET"] = ("unknown", "Version of Cabaret")
        return header
