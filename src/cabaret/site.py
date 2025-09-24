from dataclasses import dataclass


@dataclass
class Site:
    """
    Observatory site configuration.
    """

    sky_background: float = 150
    """Sky background in e-/m^2/arcsec^2/s."""

    seeing: float = 1.3
    """Atmospheric seeing in arcseconds."""

    latitude: float | None = None
    """Site latitude in degrees."""

    longitude: float | None = None
    """Site longitude in degrees."""
    sky_background: float = 150
    seeing: float = 1.3
    latitude: float | None = None
    longitude: float | None = None
