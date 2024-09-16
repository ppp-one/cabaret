def test_readme_1():
    import cabaret

    ra, dec = 12.3323, 30.4343  # in degrees
    exposure_time = 10  # in seconds

    image = cabaret.generate_image(ra, dec, exposure_time)

    assert image.shape == (
        1024,
        1024,
    ), f"Expected image shape (1024, 1024), but got {image.shape}"


def test_readme_2():
    import cabaret
    from cabaret import Camera

    camera = Camera(read_noise=10, gain=1, width=100, height=100)

    ra, dec = 12.3323, 30.4343  # in degrees
    exposure_time = 10  # in seconds

    image = cabaret.generate_image(ra, dec, exposure_time, camera=camera)

    assert image.shape == (
        100,
        100,
    ), f"Expected image shape (1000, 1000), but got {image.shape}"
