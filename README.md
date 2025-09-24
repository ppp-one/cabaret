# cabaret
![](example.jpg)
*cabaret* is a Python package to simulate astronomical images using the [GAIA catalog](https://en.wikipedia.org/wiki/Gaia_catalogues) of stars.

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

## Example

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

To display the image (`matplotlib` required here):

```python
import matplotlib.pyplot as plt
from cabaret.plot import plot_image

plot_image(image)
plt.show()
```

## Documentation

Explore the full documentation, API reference, and advanced usage examples at [ppp-one.github.io/cabaret/](https://ppp-one.github.io/cabaret/)
