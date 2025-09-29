# cabaret
![](example.jpg)
*cabaret* is a Python package to simulate astronomical images using the [Gaia catalog](https://en.wikipedia.org/wiki/Gaia_catalogues) of stars.

Documentation can be found at [cabaret.readthedocs.io](https://cabaret.readthedocs.io/).

## Installation

You can install *cabaret* in a Python (`>=3.11`) environment with

```bash
pip install cabaret
```

or from a local clone

```bash
git clone https://github.com/ppp-one/cabaret
pip install -e cabaret
```

You can test the package has been properly installed with

```bash
python -c "import cabaret"
```

## Quick Start

### Basic image 

To generate an image from RA/DEC coordinates, run:
```python
import cabaret

image = cabaret.Observatory().generate_image(
    ra=12.33230,  # right ascension in degrees
    dec=30.4343,  # declination in degrees
    exp_time=10,  # exposure time in seconds
)
```

To display the image:

```python
from cabaret.plot import plot_image

plot_image(image, contrast=0.25)
```

### Custom Observatory Configuration

Create a fully customized observatory by defining each component:

```python
from datetime import UTC, datetime
import cabaret

# Define observatory components
site = cabaret.Site(
    sky_background=150,  # e-/m^2/arcsec^2/s
    seeing=1.5,  # arcseconds
    elevation=2500,  # meters
)

telescope = cabaret.Telescope(
    focal_length=8,  # meters
    diameter=1.0,  # meters
)

camera = cabaret.Camera(
    name="Example Camera",
    width=2048,  # pixels
    height=2048,  # pixels
    pitch=10,  # microns
    gain=1,  # electrons per ADU
    read_noise=6.2,  # electrons
)

# Create the observatory
observatory = cabaret.Observatory(
    name="My Observatory",
    site=site,
    telescope=telescope,
    camera=camera,
)

# Generate a FITS image with metadata
hdu = observatory.generate_fits_image(
    ra=323.362583,
    dec=-0.82325,
    exp_time=0.5,
    dateobs=datetime.now(UTC),
    filter_band=cabaret.Filters.G,
    seed=42,
)

# Save as FITS file
hdu.writeto("simulated_image.fits", overwrite=True)
```

### Custom Source Lists

You can manipulate the source catalog before generating images:

```python
from astropy.coordinates import SkyCoord
import numpy as np

# Query Gaia catalog
center = SkyCoord(ra=323.362583, dec=-0.82325, unit="deg")
table = cabaret.GaiaQuery.query(
    center=center,
    radius=camera.get_fov_radius() * 1.5,
    filter_band=cabaret.Filters.G,
)

# Filter out bright sources
fluxes = table[cabaret.Filters.G.value].value.data
mask = fluxes < 1e4
filtered_table = table[mask]

# Create custom source list
sources = cabaret.Sources.from_arrays(
    ra=filtered_table["ra"].value.data,
    dec=filtered_table["dec"].value.data,
    fluxes=filtered_table[cabaret.Filters.G.value].value.data,
)

# Generate image with custom sources
image = observatory.generate_image(
    ra=center.ra.deg,
    dec=center.dec.deg,
    exp_time=0.5,
    sources=sources,
)
```

Explore the full documentation, API reference, and advanced usage examples at [cabaret.readthedocs.io](https://cabaret.readthedocs.io/).
