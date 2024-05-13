# cabaret

*cabaret* is a Python package to simulate astronomical images using the [GAIA catalog](https://en.wikipedia.org/wiki/Gaia_catalogues) of stars.
## Installation

You can install *cabaret* with

```bash
pip install cabaret
```

or from a local clone

```bash
git clone https://github.com/ppp-one/cabaret
pip install -e cabaret
```

## Example

### Basic image 

To generate an image from RA/DEC coordinates and a field of vue specified in degrees:

```python
import cabaret

field_of_view = (0.1, 0.1) # in degrees 
ra, dec = 12.3323, 30.4343 # in degrees

image = cabaret.image(ra, dec, field_of_view)
```

and to display the image (`matplotlib` required here):

```python
import matplotlib.pyplot as plt

plt.imshow(image)
```

### Using the camera characteristics

To  adjust the physical characteristics of the camera, you can define and pass a `Camera` object

```python
import cabaret
from cabaret import Camera

camera = Camera(read_noise=10, gain=1, exposure_time=1)

field_of_view = (0.1, 0.1) # in degrees
ra, dec = 12.3323, 30.4343 # in degrees

image = cabaret.image(ra, dec, field_of_view, camera=camera)
``` 