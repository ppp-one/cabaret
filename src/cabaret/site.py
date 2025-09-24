from dataclasses import dataclass


@dataclass
class Site:
    """
    Observatory site configuration.

    Attributes
    ----------
    sky_background : float
        Sky background in e-/m^2/arcsec^2/s.
    seeing : float
        Atmospheric seeing in arcseconds.
    latitude : float | None
        Site latitude in degrees.
    longitude : float | None
        Site longitude in degrees.
    """

    sky_background: float = 150
    seeing: float = 1.3
    latitude: float | None = None
    longitude: float | None = None
