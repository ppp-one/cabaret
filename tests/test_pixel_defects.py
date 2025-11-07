import numpy as np
import pytest

from cabaret.camera import (
    Camera,
    ColumnPixelDefect,
    ConstantPixelDefect,
    RandomNoisePixelDefect,
)


@pytest.fixture
def camera():
    return Camera(height=128, width=128)


@pytest.fixture
def image(camera):
    return np.zeros((camera.height, camera.width), dtype=np.float32)


def test_constant_pixel_defect_initialization():
    defect = ConstantPixelDefect(name="constant_defect", rate=0.1, value=1024)
    assert defect.name == "constant_defect"
    assert defect.rate == 0.1
    assert defect.value == 1024


def test_pixel_defect_number_of_defect_pixels():
    camera = Camera(width=1024, height=1024)
    defect = ConstantPixelDefect(rate=0.1)
    assert defect.number_of_defect_pixels(camera) == np.round(0.1 * 1024 * 1024)


def test_pixel_defect_select_random_pixels(camera):
    defect = ConstantPixelDefect(name="test_defect", rate=0.1, seed=42)
    defect.set_pixels(defect._select_random_pixels(camera), camera)

    expected_number_of_pixels = defect.number_of_defect_pixels(camera)
    assert defect.pixels.shape[0] == expected_number_of_pixels


def test_constant_pixel_defect_introduce_pixel_defect_with_set_pixels(image):
    camera = Camera(width=5, height=5)
    defect = ConstantPixelDefect(value=128, seed=42)

    # Set specific pixels
    specific_pixels = np.array([[1, 1], [2, 2], [3, 3]])
    defect.set_pixels(specific_pixels, camera)

    defect.introduce_pixel_defect(image, camera)

    # Check that the specific pixels have been set to the defect value
    for pixel in specific_pixels:
        assert image[pixel[0], pixel[1]] == 128


def test_pixel_defect_overwrite_pixel_values(camera, image):
    defect = ConstantPixelDefect(name="test_defect", rate=0.1, seed=42)
    defect.set_pixels(defect._select_random_pixels(camera), camera)

    # Create a copy of the image to modify
    modified_image = image.copy()
    defect._overwrite_pixel_values(modified_image, defect.pixels, 255)

    # Check that the specified pixels have been modified
    for pixel in defect.pixels:
        assert modified_image[pixel[0], pixel[1]] == 255


def test_pixel_defect_check_pixel_bounds_valid(camera):
    defect = ConstantPixelDefect(name="test_defect", rate=0.1, seed=42)
    valid_pixels = np.array([[0, 0], [camera.height - 1, camera.width - 1]])
    defect._check_pixel_bounds(
        valid_pixels, camera.height, camera.width, defect.name
    )  # Should not raise


def test_pixel_defect_check_pixel_bounds_invalid(camera):
    defect = ConstantPixelDefect(name="test_defect", rate=0.1, seed=42)
    invalid_pixels = np.array([[camera.height, camera.width]])  # Out of bounds
    with pytest.raises(ValueError):
        defect._check_pixel_bounds(
            invalid_pixels, camera.height, camera.width, defect.name
        )


def test_constant_pixel_defect_introduce(camera, image):
    defect = ConstantPixelDefect(name="constant_defect", rate=0.1, value=1024)
    defect.set_pixels(defect._select_random_pixels(camera), camera)
    modified_image = defect.introduce_pixel_defect(image.copy(), camera)

    # Check that the defect pixels have the correct value
    for pixel in defect.pixels:
        assert modified_image[pixel[0], pixel[1]] == 1024


def test_constant_pixel_defect_no_defect_pixels(camera, image):
    defect = ConstantPixelDefect(name="constant_defect", rate=0.0, value=255)
    defect.set_pixels(defect._select_random_pixels(camera), camera)
    modified_image = defect.introduce_pixel_defect(image.copy(), camera)

    # Check that the image remains unchanged
    assert np.array_equal(modified_image, image)


def test_constant_pixel_defect_invalid_pixel_selection(camera):
    defect = ConstantPixelDefect(name="constant_defect", rate=0.1, value=255)
    with pytest.raises(ValueError):
        defect.set_pixels(np.array([[camera.height, camera.width]]), camera)


def test_column_pixel_defect_initialization():
    defect = ColumnPixelDefect(name="test_defect", rate=0.1, value=1024, dim=0)
    assert defect.name == "test_defect"
    assert defect.rate == 0.1
    assert defect.value == 1024
    assert defect.dim == 0


def test_column_pixel_defect_introduce(camera, image):
    defect = ColumnPixelDefect(name="test_defect", rate=0.1, value=1024, dim=0)
    # defect.set_pixels(defect._select_random_pixel(camera), camera)
    modified_image = defect.introduce_pixel_defect(image.copy(), camera)

    # Check that the defect pixels have the correct value
    for pixel in defect.pixels:
        assert modified_image[pixel[0], pixel[1]] == 1024


def test_column_pixel_defect_no_defect_pixels(camera, image):
    defect = ColumnPixelDefect(name="test_defect", rate=0.0, value=255, dim=0)
    modified_image = defect.introduce_pixel_defect(image.copy(), camera)

    assert np.array_equal(modified_image, image)


def test_column_pixel_defect_zero_defect_pixels(camera, image):
    defect = ColumnPixelDefect(name="test_defect", rate=0.0, value=255, dim=0)
    modified_image = defect.introduce_pixel_defect(image.copy(), camera)

    assert np.array_equal(modified_image, image)


def test_random_noise_pixel_defect_initialization():
    defect = RandomNoisePixelDefect(
        name="test_noise", rate=0.1, noise_level=10.0, distribution="normal"
    )
    assert defect.name == "test_noise"
    assert defect.rate == 0.1
    assert defect.noise_level == 10.0
    assert defect.distribution == "normal"


def test_random_noise_pixel_defect_introduce(camera, image):
    defect = RandomNoisePixelDefect(
        name="test_noise", rate=0.1, noise_level=10.0, distribution="normal"
    )
    modified_image = defect.introduce_pixel_defect(image.copy(), camera)

    assert np.sum(np.abs(modified_image - image)) > 0


def test_random_noise_pixel_defect_distribution(camera, image):
    defect = RandomNoisePixelDefect(
        name="test_noise", rate=0.1, noise_level=10.0, distribution="normal"
    )
    defect.set_pixels(defect._select_random_pixels(camera), camera)
    modified_image = defect.introduce_pixel_defect(image.copy(), camera)

    assert np.sum(np.abs(modified_image - image)) > 0


def test_random_noise_pixel_defect_clipping(camera, image):
    defect = RandomNoisePixelDefect(
        name="test_noise", rate=0.1, noise_level=1000.0, distribution="normal"
    )
    defect.set_pixels(defect._select_random_pixels(camera), camera)
    modified_image = defect.introduce_pixel_defect(image.copy(), camera)

    assert np.all(modified_image <= camera.max_adu)
