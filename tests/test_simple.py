import pytest
from cabaret import generate_image
from cabaret import Camera


def test_simple():
    camera = Camera(width=100, height=100)
    ra, dec = 323.36152, -0.82325
    image = generate_image(ra, dec, 1, camera=camera)
