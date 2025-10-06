---
title: 'cabaret: A Python package for simulating astronomical images'
tags:
  - Python
  - astronomy
  - image simulation
  - Gaia catalog
  - stellar fields
  - observatory instrumentation
authors:
  - name: Peter Pedersen
    orcid: 0000-0002-5220-609X
    affiliation: 1
    equal-contrib: true
  - name: David Degen
    orcid: 0009-0008-1068-481X
    affiliation: 1
    equal-contrib: true
  - name: Lionel Garcia
    orcid: 0000-0002-4296-2246
    affiliation: 2
    equal-contrib: true

affiliations:
 - name: ETH Zürich, Department of Physics, ETH Zurich, Wolfgang-Pauli-Strasse 2, 8093 Zurich, Switzerland
   index: 1
 - name: Institution Name, Address, Country
   index: 2

date: 1 October 2025
bibliography: paper.bib
---

# Summary

Astronomical research increasingly relies on realistic simulations to interpret observations, test hypotheses, and develop new analysis techniques. `cabaret` is a Python package designed to simulate astronomical images using the [Gaia catalog](https://www.cosmos.esa.int/web/gaia) of stars, providing researchers and educators with a fast, flexible tool for generating synthetic stellar field images. The package seamlessly integrates real astronomical data from Gaia with customizable observatory configurations, enabling users to simulate images that accurately reflect observational conditions including telescope optics, camera properties, atmospheric seeing, and detector effects. `cabaret` is particularly well-suited for validating data reduction pipelines, training machine learning models, developing observatory control software, and educational applications where realistic astronomical data is needed.

# Statement of need

The development of modern astronomical instrumentation and analysis pipelines requires extensive testing with realistic data. While real observational data is invaluable, it has limitations: observations are expensive, time-consuming, and often lack the ground truth necessary for validation. Synthetic images provide a controlled environment where all parameters are known, making them essential for testing algorithms, training machine learning models, and validating software systems.

Existing image simulation tools often fall into two categories: highly specialized packages designed for specific surveys or instruments, and general-purpose simulators that require extensive configuration and domain expertise [@ufig]. `cabaret` fills a gap by providing an accessible, easy-to-use package that generates realistic stellar field images with minimal setup while maintaining the flexibility to customize observatory parameters for specific use cases.

The package is designed around several key principles:

- **Accessibility**: With just a few lines of Python code, users can generate their first synthetic image using Gaia catalog data.
- **Realism**: Images incorporate stellar positions and flux values from Gaia, combined with Moffat point spread functions (PSF) and detector effects.
- **Flexibility**: All observatory components (telescope, camera, focuser, site) are configurable through a simple API, allowing users to match specific instruments or explore parameter spaces.
- **Reproducibility**: All simulations are deterministic when a random seed is provided, ensuring reproducible results for scientific workflows.

`cabaret` has already proven valuable in the development of `alpaca-simulators` [@alpaca], a comprehensive simulator for ASCOM Alpaca devices. By providing realistic device simulation and image generation, `alpaca-simulators` enables thorough testing of robotic telescope control software and observatory automation systems without requiring access to physical hardware.

The package is also particularly relevant for projects like SPECULOOS (Search for habitable Planets EClipsing ULtra-cOOl Stars) [@speculoos], which monitors ultracool dwarf stars to detect transiting exoplanets. For such surveys, the ability to generate synthetic images that closely match real observations is crucial for validating photometric pipelines, characterizing systematic uncertainties, and optimizing observing strategies.

# Key Features

## Gaia Catalog Integration

`cabaret` provides access to the Gaia catalog through the `astroquery` [@astroquery] package, automatically querying and retrieving stellar positions, proper motions, fluxes, for any field of view. Users simply specify the sky coordinates and field size, and `cabaret` handles the catalog retrieval and flux conversion. Alternatively, users can provide their own source catalogs for full control over the simulated stellar population.

## Configurable Observatory Model

The package implements a modular observatory model with four main components:

- **Telescope**: Configurable aperture, focal length, and optical properties
- **Camera**: Customizable detector dimensions, pixel scale, gain, dark current, readout noise, and exposure time
- **Focuser**: Position and offset parameters for focus effects
- **Site**: Atmospheric seeing and sky background conditions

All components can be instantiated with sensible defaults or customized to match real instruments. For example, simulating images from a specific telescope requires only specifying its aperture and focal length:

```python
import cabaret

telescope = cabaret.Telescope(aperture=1000, focal_length=8000)
observatory = cabaret.Observatory(telescope=telescope)
image = observatory.generate_image(ra=12.33, dec=30.43, exp_time=10)
```

## Point Spread Function Modeling

Stars are rendered using a Moffat profile, a physically-motivated functional form commonly used to model atmospheric seeing. The full width at half maximum (FWHM) can be specified directly or computed from seeing conditions and airmass. The package supports:

- Configurable Moffat β parameter for PSF shape
- Adaptive rendering radius based on FWHM
- Proper flux normalization
- Poisson noise in photon generation

## Detector Effects and Noise

`cabaret` includes a comprehensive set of detector effects to create realistic images:

- **Photon shot noise**: Poisson-distributed noise from stellar photons
- **Readout noise**: Gaussian noise from detector readout
- **Dark current**: Temperature-dependent dark current accumulation
- **Bias level**: Configurable detector bias
- **Gain**: Conversion from electrons to analog-to-digital units (ADU)
- **Saturation**: Proper handling of well depth and ADU limits
- **Pixel defects**: Support for cold pixels, hot pixels, quantum efficiency variations, and readout smear

These effects can be enabled or disabled individually, allowing users to isolate specific contributions or create idealized images for testing.

## World Coordinate System

All generated images include proper World Coordinate System (WCS) headers, enabling seamless integration with astronomical analysis tools. The WCS implementation correctly handles:

- Tangent plane projection
- Pixel scale from camera and telescope properties
- Coordinate system transformations
- FITS header generation

# Implementation

`cabaret` is implemented in pure Python with minimal dependencies (`numpy` [@numpy], `astropy` [@astropy], `astroquery` [@astroquery]), making it easy to install and integrate into existing workflows. The package follows modern Python best practices:

- Type hints throughout the codebase
- Dataclass-based configuration for clean, validated interfaces
- Comprehensive docstrings with usage examples
- Automated testing with `pytest` [@pytest]
- Continuous integration and documentation hosting

The image generation algorithm uses an efficient windowed rendering approach, where stars are rendered only within a small region around their position (typically 5× the FWHM). This dramatically reduces computation time compared to global convolution methods while maintaining accuracy.

The package is designed to be extensible, with clear interfaces for adding new detector effects, PSF models, or catalog sources.

# Example Usage and Validation

\autoref{fig:comparison} shows a comparison between a real SPECULOOS observation and a simulated image generated with `cabaret` using matching observatory parameters. The simulated image accurately reproduces the stellar field, PSF characteristics, background noise, and overall appearance of the real observation, demonstrating the package's ability to create realistic synthetic data for pipeline validation and testing.

![Comparison between a real SPECULOOS observation (left) and a simulated image generated with `cabaret` (right). Both images show the same field of view with matching exposure time and seeing conditions. The simulated image successfully reproduces the stellar positions, brightness distribution, and detector characteristics of the real observation.\label{fig:comparison}](figures/speculoos_comparison.png)

The basic workflow for generating a synthetic image is straightforward:

```python
import cabaret

# Create observatory with default configuration
observatory = cabaret.Observatory()

# Generate image from Gaia catalog
image = observatory.generate_image(
    ra=12.33230,      # right ascension in degrees
    dec=30.4343,      # declination in degrees
    exp_time=10,      # exposure time in seconds
)
```

For more advanced use cases, all observatory components can be customized:

```python
from cabaret import Observatory, Camera, Telescope, Site

# Configure a specific instrument
camera = Camera(
    width=2048,
    height=2048,
    pitch=13.5,         # microns
    gain=1.5,           # e-/ADU
    read_noise=10,      # electrons
    dark_current=0.1,   # e-/pixel/s
)

telescope = Telescope(
    aperture=1000,      # mm
    focal_length=8000,  # mm
)

site = Site(
    latitude=28.3,      # degrees
    longitude=-16.5,    # degrees
    elevation=2400,     # meters
)

observatory = Observatory(
    camera=camera,
    telescope=telescope,
    site=site,
)
```

# Applications

`cabaret` is designed for a wide range of applications:

- **Pipeline validation**: Test photometry, astrometry, and image processing algorithms with known ground truth
- **Machine learning**: Generate training datasets for neural networks doing source detection, classification, or deblending
- **Observatory control**: Simulate telescope and camera behavior for software development and testing
- **Education**: Teach students about observational astronomy, detector physics, and data analysis
- **Survey planning**: Explore parameter spaces to optimize observing strategies
- **Systematic uncertainty characterization**: Quantify the impact of various instrumental and atmospheric effects

The package has been successfully used in the development of `alpaca-simulators`, which provides a complete simulation environment for ASCOM Alpaca devices, enabling developers to test observatory control software without physical hardware.

# Acknowledgements

We acknowledge the contributions of several key libraries to the functionality of cabaret, specifically `astropy` [@astropy],`astroquery` [@astroquery], and `numpy` [@numpy]. Additionally, we utilized `Matplotlib` [@matplotlib] for the plots in this paper and the package's documentation. We also utilized `prose` [@prose], `photutils` [@photutils], and `pandas` [@reback2020pandas] in the comparison with real observations. Furthermore, testing was conducted using `pytest` [@pytest] to ensure the reliability of our code. This work made use of the Gaia catalog, ESA's space mission for stellar astrometry.

# References
