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

    >>> from cabaret.sources import Sources
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
    """An array of shape (n,) containing the fluxes of the sources in units of
    photons/s/m²."""

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

    @property
    def center(self) -> tuple[float, float]:
        """Midpoint RA and DEC of the sources in degrees."""
        ra_min, ra_max = self.ra.deg.min(), self.ra.deg.max()  # type: ignore
        ra_range = ra_max - ra_min

        if ra_range > 180:
            ra_shifted = (self.ra.deg + 180) % 360  # type: ignore
            ra_min, ra_max = ra_shifted.min(), ra_shifted.max()  # type: ignore
            ra_center = (ra_min + ra_max) / 2 - 180
            ra_center %= 360
        else:
            ra_center = (ra_min + ra_max) / 2

        dec_center = (self.dec.deg.min() + self.dec.deg.max()) / 2  # type: ignore
        return ra_center, dec_center

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

    def __getitem__(self, key) -> "Sources":
        new_fluxes = np.asarray(self.fluxes)[key]
        if np.isscalar(new_fluxes):
            new_fluxes = np.array([new_fluxes])

        new_coords = self.coords[key]

        return Sources(new_coords, new_fluxes)  # type: ignore

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
            These should be of units photons/s/m².
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

    def drop_nan_fluxes(self) -> "Sources":
        """Return a new Sources instance with any sources that have NaN fluxes removed.

        If no sources have NaN fluxes, returns self.

        Returns
        -------
        Sources
            A Sources instance with NaN flux entries removed.
        """
        mask = np.isnan(self.fluxes)
        if mask.any():
            return self[~mask]
        return self

    @classmethod
    def get_test_sources(cls) -> "Sources":
        """Return a simple test Sources instance."""
        coords = SkyCoord(
            ra=[12.29611593, 12.29929654, 12.33757534, 12.34247842, 12.29354464],
            dec=[30.45675318, 30.44855405, 30.42613357, 30.48059276, 30.47310728],
            unit="deg",
        )
        fluxes = np.array([307220.0, 64271.0, 61002.0, 43466.0, 9239.0])
        return cls(coords, fluxes=fluxes)
