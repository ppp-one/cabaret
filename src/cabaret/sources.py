from dataclasses import dataclass

import numpy as np
from astropy.coordinates import Longitude, SkyCoord
from astropy.wcs import WCS


@dataclass
class Sources:
    """A collection of sources with their sky coordinates and fluxes.

    Examples
    --------

    Create a Sources instance from arrays:

    >>> from cabaret.queries import Sources
    >>> import numpy as np
    >>> from astropy.coordinates import SkyCoord
    >>> coords = np.array([[10.64, 41.26], [10.68, 41.22]])
    >>> fluxes = np.array([169435.6, 52203.9])
    >>> sources = Sources(SkyCoord(coords, unit='deg'), fluxes)
    >>> sources
    Sources(coords=<SkyCoord (ICRS): (ra, dec) in deg
        [(10.64, 41.26), (10.68, 41.22)]>, fluxes=array([169435.6,  52203.9]))

    """

    coords: SkyCoord
    """SkyCoords instance with the RA and DEC coordinates of the sources."""
    fluxes: np.ndarray
    """An array of shape (n,) containing the fluxes of the sources."""

    def __post_init__(self):
        if not isinstance(self.coords, SkyCoord):
            raise ValueError("coords must be an instance of SkyCoord.")
        if not isinstance(self.fluxes, np.ndarray):
            try:
                self.fluxes = np.array(self.fluxes)
            except Exception:
                raise ValueError("fluxes must be an instance of np.ndarray.")
        if self.coords.size != self.fluxes.size:
            raise ValueError("coords and fluxes must have the same length.")

    @property
    def ra(self) -> Longitude:
        """Right Ascension coordinates of the sources."""
        return self.coords.ra  # type: ignore

    @property
    def dec(self) -> Longitude:
        """Declination coordinates of the sources."""
        return self.coords.dec  # type: ignore

    def to_pixel(self, wcs: WCS) -> np.ndarray:
        """Convert the RA-DEC coordinates to pixel coordinates using the given WCS.

        Parameters
        ----------
        wcs : astropy.wcs.WCS
            The WCS object used for the conversion.

        Returns
        -------
        np.ndarray
            An array of shape (n, 2) containing the pixel coordinates of the sources.
        """
        return np.array(self.coords.to_pixel(wcs))

    def __len__(self) -> int:
        return len(self.fluxes)

    @classmethod
    def from_arrays(
        cls,
        ra: np.ndarray | list,
        dec: np.ndarray | list,
        fluxes: np.ndarray | list,
        units: str = "deg",
    ) -> "Sources":
        """Create a Sources instance from separate RA and DEC arrays.

        Parameters
        ----------
        ra : np.ndarray
            An array of shape (n,) containing the RA coordinates of the sources in deg.
        dec : np.ndarray
            An array of shape (n,) containing the DEC coordinates of the sources in deg.
        fluxes : np.ndarray
            An array of shape (n,) containing the fluxes of the sources.
        **kwargs
            Additional keyword arguments passed to the Sources constructor.

        Returns
        -------
        Sources
            A Sources instance.
        """
        if not isinstance(ra, np.ndarray):
            try:
                ra = np.array(ra)
            except Exception:
                raise ValueError("ra must be an instance of np.ndarray.")
        if not isinstance(dec, np.ndarray):
            try:
                dec = np.array(dec)
            except Exception:
                raise ValueError("dec must be an instance of np.ndarray.")
        if ra.shape != dec.shape:
            raise ValueError("ra and dec must have the same shape.")

        parameters = {
            "coords": SkyCoord(ra=ra, dec=dec, unit=units),
            "fluxes": fluxes,
        }
        return cls(**(parameters))

    @classmethod
    def get_test_source(cls) -> "Sources":
        """Return a simple test Sources instance."""
        coords = SkyCoord(ra=[10.64], dec=[10.68], unit="deg")
        fluxes = np.array([169_435.6])
        return cls(coords, fluxes)
