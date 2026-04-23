from dataclasses import dataclass, field

import numpy as np
from astropy import units as u
from astropy.coordinates import Longitude, SkyCoord
from astropy.wcs import WCS


def _normalize_rates(
    rates: np.ndarray | u.Quantity | None, n: int, name: str
) -> np.ndarray:
    """Convert rates to a plain arcsec/s numpy array of shape (n,)."""
    if rates is None:
        return np.zeros(n, dtype=np.float64)
    if isinstance(rates, u.Quantity):
        rates = rates.to(u.arcsec / u.s).value
    rates = np.asarray(rates, dtype=np.float64)
    if rates.shape != (n,):
        raise ValueError(f"{name} must have shape ({n},).")
    return rates


@dataclass(init=False)
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
    ra_rates: np.ndarray = field(init=False, repr=False)
    """On-sky RA motion rates dα·cos(δ)/dt, shape (n,) in arcsec/s. Accepts an
    astropy Quantity with angular velocity units (e.g. ``u.arcsec/u.hour``) or a
    plain array assumed to be in arcsec/s. Stored internally as arcsec/s. Positive
    values move east. This matches the JPL Horizons ``RA_rate`` convention."""
    dec_rates: np.ndarray = field(init=False, repr=False)
    """Dec motion rates, shape (n,) in arcsec/s. Accepts an astropy Quantity with
    angular velocity units or a plain array assumed to be in arcsec/s. Stored
    internally as arcsec/s. Positive values move north."""

    def __init__(
        self,
        coords: SkyCoord,
        fluxes: np.ndarray,
        ra_rates: np.ndarray | u.Quantity | None = None,
        dec_rates: np.ndarray | u.Quantity | None = None,
    ) -> None:
        if not isinstance(coords, SkyCoord):
            raise ValueError("coords must be an instance of SkyCoord.")
        if not isinstance(fluxes, np.ndarray):
            try:
                fluxes = np.array(fluxes)
            except Exception:
                raise ValueError("fluxes must be an instance of np.ndarray.")
        if coords.size != fluxes.size:
            raise ValueError("coords and fluxes must have the same length.")

        self.coords = coords
        self.fluxes = fluxes
        n = coords.size
        self.ra_rates = _normalize_rates(ra_rates, n, "ra_rates")
        self.dec_rates = _normalize_rates(dec_rates, n, "dec_rates")

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
        new_ra_rates = np.asarray(self.ra_rates)[key]
        new_dec_rates = np.asarray(self.dec_rates)[key]
        if np.isscalar(new_fluxes):
            new_fluxes = np.array([new_fluxes])
            new_ra_rates = np.array([new_ra_rates])
            new_dec_rates = np.array([new_dec_rates])

        new_coords = self.coords[key]

        return Sources(new_coords, new_fluxes, new_ra_rates, new_dec_rates)

    def __add__(self, other: "Sources") -> "Sources":
        return Sources.concat(self, other)

    @classmethod
    def concat(cls, *sources_list: "Sources") -> "Sources":
        """Concatenate multiple Sources objects into one.

        Parameters
        ----------
        *sources_list : Sources
            One or more Sources objects to concatenate.

        Returns
        -------
        Sources
            A new Sources instance with all sources combined.
        """
        if not sources_list:
            raise ValueError("Must provide at least one Sources object to concatenate.")
        return cls(
            coords=np.concatenate([s.coords for s in sources_list]),  # type: ignore[arg-type]
            fluxes=np.concatenate([s.fluxes for s in sources_list]),
            ra_rates=np.concatenate([s.ra_rates for s in sources_list]),
            dec_rates=np.concatenate([s.dec_rates for s in sources_list]),
        )

    @classmethod
    def from_arrays(
        cls,
        ra: np.ndarray | list,
        dec: np.ndarray | list,
        fluxes: np.ndarray | list,
        units: str = "deg",
        ra_rates: np.ndarray | list | None = None,
        dec_rates: np.ndarray | list | None = None,
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
        units : str, optional
            Units for ra and dec, by default "deg".
        ra_rates : np.ndarray or list, optional
            RA motion rates in arcsec/s, shape (n,). Defaults to zeros.
        dec_rates : np.ndarray or list, optional
            Dec motion rates in arcsec/s, shape (n,). Defaults to zeros.

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
        if not isinstance(fluxes, np.ndarray):
            try:
                fluxes = np.array(fluxes)
            except Exception:
                raise ValueError("fluxes must be an instance of np.ndarray.")
        if not (ra.shape == fluxes.shape):
            raise ValueError("ra, dec, and fluxes must all have the same shape.")

        return cls(
            coords=SkyCoord(ra=ra, dec=dec, unit=units),
            fluxes=fluxes,
            ra_rates=ra_rates,
            dec_rates=dec_rates,
        )

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
