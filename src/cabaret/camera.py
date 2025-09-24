from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Literal

import numpy as np
import numpy.random
from astropy.wcs import WCS


@dataclass
class Camera:
    """
    Camera configuration and properties.

    Attributes
    ----------
    name : str
        Name of the camera.
    width : int
        Number of pixels in the x-direction (width).
    height : int
        Number of pixels in the y-direction (height).
    bin_x : int
        Binning factor in x.
    bin_y : int
        Binning factor in y.
    pitch : float
        Pixel pitch in microns.
    plate_scale : float or None
        Arcseconds per pixel (if None, calculated from pitch and telescope).
    max_adu : int
        Maximum ADU value.
    well_depth : int
        Full well capacity in electrons.
    bias : int
        Bias level in ADU.
    gain : float
        Gain in electrons per ADU.
    read_noise : float
        Read noise in electrons.
    dark_current : float
        Dark current in electrons per second.
    average_quantum_efficiency : float
        Average quantum efficiency (fraction).
    pixel_defects : dict
        Dictionary of pixel defect configurations.

    Examples
    --------
    Create a camera with multiple pixel defects:

    >>> from cabaret import Camera, Observatory, Sources
    >>> pixel_defects = {
    ...     "hot": {"type": "constant", "value": 10_000, "rate": 0.1, "seed": 0},
    ...     "cold": {"type": "constant", "value": 0, "rate": 0.01, "seed": 1},
    ...     "noise": {"type": "noise", "rate": 0.02, "noise_level": 100, "seed": 2},
    ...     "qe": {"type": "quantum_efficiency_map", "seed": 3}
    ... }
    >>> camera = Camera(width=100, height=100, pixel_defects=pixel_defects)
    >>> observatory = Observatory(camera=camera)
    >>> sources = Sources.get_test_source()
    >>> ra, dec = sources.ra.mean().deg, sources.dec.mean().deg,

    >>> image = observatory.generate_image(
    ...     exp_time=10, ra=ra, dec=dec, seed=42, sources=sources
    ... )
    >>> clean_image = observatory.generate_image(
    ...     exp_time=10, ra=ra, dec=dec,
    ...     apply_pixel_defects=False, seed=42, sources=sources
    ... )

    >>> import matplotlib.pyplot as plt
    >>> fig, axes = plt.subplots(1, 2, figsize=(7, 5), sharex=True, sharey=True)
    >>> im0 = axes[0].imshow(clean_image, cmap='gray')
    >>> cbar0 = fig.colorbar(im0, ax=axes[0], orientation="horizontal", pad=0.2)
    >>> cbar0.set_label('ADU')
    >>> im1 = axes[1].imshow(image, cmap='gray')
    >>> cbar1 = fig.colorbar(im1, ax=axes[1], orientation="horizontal", pad=0.2)
    >>> cbar1.set_label('ADU')
    >>> axes[0].set_title('Image without defects')
    >>> axes[1].set_title('Image with multiple defects')
    >>> axes[0].set_ylabel('Pixels')
    >>> for ax in axes:
    ...     ax.set_xlabel('Pixels')
    >>> plt.subplots_adjust(wspace=0.1)
    >>> plt.show()
    """

    name: str = "gaia-camera-simulated"
    width: int = 1024  # pixels
    height: int = 1024  # pixels
    bin_x: int = 1  # binning factor in x
    bin_y: int = 1  # binning factor in y
    pitch: float = 13.5  # pixel pitch, microns
    plate_scale: float | None = (
        None  # arcsec/pixel (calculated from pitch+telescope if None)
    )
    max_adu: int = 2**16 - 1  # maximum ADU value
    well_depth: int = 2**16 - 1  # electrons
    bias: int = 300  # ADU
    gain: float = 1.0  # e-/ADU
    read_noise: float = 6.2  # e-
    dark_current: float = 0.2  # e-/s
    average_quantum_efficiency: float = 0.8  # fraction
    pixel_defects: dict = field(default_factory=dict)

    @property
    def shape(self) -> tuple[int, int]:
        """Tuple of (height, width) of the camera in pixels."""
        return (self.height, self.width)

    @property
    def size(self) -> int:
        """Total number of pixels in the camera."""
        return self.height * self.width

    def __post_init__(self):
        if self.pixel_defects:
            self.pixel_defects = {
                key: (
                    self._create_pixel_defect(key, **value)
                    if isinstance(value, dict)
                    else value
                )
                for key, value in self.pixel_defects.items()
            }

    def _create_pixel_defect(
        self,
        name: str,
        type: Literal["constant", "telegraphic", "noise", "fixed_pattern"] = "constant",
        **kwargs,
    ) -> "PixelDefect":
        defect_classes = {
            "constant": ConstantPixelDefect,
            "telegraphic": TelegraphicPixelDefect,
            "noise": RandomNoisePixelDefect,
            "quantum_efficiency_map": QuantumEfficiencyMapPixelDefect,
            "cti": ChargeTransferInefficiencyPixelDefect,
            "readout_smear": ReadoutSmearPixelDefect,
        }
        defect_type = type

        if defect_type not in defect_classes:
            raise ValueError(f"Unknown pixel defect type for {name}")

        return defect_classes[defect_type](name=name, **kwargs)

    @classmethod
    def create_ideal_camera(cls, **kwargs) -> "Camera":
        """Create a defect-free camera."""
        parameters = {
            "read_noise": 0,
            "dark_current": 0,
            "average_quantum_efficiency": 1.0,
            "bias": 0,
            "gain": 1.0,
        }
        return cls(**(parameters | kwargs))

    def get_wcs(self, center):
        """Get a WCS object for the camera centered at the given sky coordinate."""
        if self.plate_scale is None:
            raise ValueError("plate_scale must be set to compute WCS.")

        wcs = WCS(naxis=2)
        wcs.wcs.cdelt = [-self.plate_scale / 3600, -self.plate_scale / 3600]
        wcs.wcs.cunit = ["deg", "deg"]
        wcs.wcs.crpix = [int(self.width / 2), int(self.height / 2)]
        wcs.wcs.crval = [center.ra.deg, center.dec.deg]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        return wcs

    def _create_blank_image(self) -> np.ndarray:
        """Create a blank image filled with zeros for testing."""
        return np.zeros((self.height, self.width), dtype=np.uint16)


@dataclass
class PixelDefect(ABC):
    name: str = "defect"
    rate: float = 0.0
    seed: int = 0
    _rng: numpy.random.Generator | None = None
    _pixels: np.ndarray | None = None

    @property
    def pixels(self) -> np.ndarray:
        if self._pixels is None:
            raise ValueError(f"{self.name} pixels are not defined.")
        return self._pixels

    @property
    def rng(self) -> numpy.random.Generator:
        if self._rng is None:
            self._rng = numpy.random.default_rng(self.seed)

        return self._rng

    @abstractmethod
    def introduce_pixel_defect(self, image: np.ndarray, camera: Camera):
        """Introduce the defect into the image.

        Parameters
        ----------
        image : np.ndarray
            The image to which the defect will be introduced.
        camera : Camera
            The camera to which the defect applies.
        """
        raise NotImplementedError

    def set_pixels(self, pixels: np.ndarray, camera: Camera):
        """Set the pixels for the defect.

        Parameters
        ----------
        pixels : np.ndarray
            The pixel coordinates of the defect.
        camera : Camera
            The camera to which the defect applies.
        """
        self._check_pixel_bounds(pixels, camera.height, camera.width, self.name)
        self._pixels = pixels

    def number_of_defect_pixels(self, camera: Camera) -> int:
        return int(round(self.rate * camera.width * camera.height))

    def _select_random_pixels(self, camera: Camera) -> np.ndarray:
        return self.rng.integers(
            [camera.height, camera.width],
            size=(self.number_of_defect_pixels(camera), 2),
        )

    @staticmethod
    def _overwrite_pixel_values(
        image: np.ndarray, pixels: np.ndarray, pixel_values: int | np.ndarray
    ):
        if isinstance(pixel_values, np.ndarray):
            if not pixel_values.size == pixels.shape[0]:
                raise ValueError("Pixel values must match the number of pixels.")

        image[pixels[:, 0], pixels[:, 1]] = pixel_values

    @staticmethod
    def _check_pixel_bounds(pixels, height: int, width: int, name: str):
        if pixels is None:
            raise ValueError(f"{name} pixels are not defined.")
        if isinstance(pixels, list):
            pixels = np.array(pixels)
        if not pixels.ndim == 2 or not pixels.shape[1] == 2:
            raise ValueError(
                f"{name} pixels must be a numpy array of (x, y) tuples."
                f"Got shape {pixels.shape} instead."
            )
        if np.any(pixels[:, 0] >= height) or np.any(pixels[:, 1] >= width):
            raise ValueError(f"{name} pixels are outside the frame.")


@dataclass
class ConstantPixelDefect(PixelDefect):
    """A pixel defect that sets selected pixels to a constant value.

    Examples
    --------
    >>> from cabaret import Observatory, Sources
    >>> pixel_defects = {
        "hot": {"type": "constant", "value": 5000, "rate": 0.01, "seed": 0},
        "cold": {"type": "constant", "value": 0, "rate": 0.01, "seed": 1}
    }
    >>> observatory = Observatory(
    ...     camera={"width": 80, "height": 80, "pixel_defects": pixel_defects}
    ... )
    >>> sources = Sources.get_test_source()
    >>> ra, dec = sources.ra.mean().deg, sources.dec.mean().deg,
    >>> image = observatory.generate_image(
    ...     exp_time=3, ra=ra, dec=dec, sources=sources, seed=42
    ... )
    >>> clean_image = observatory.generate_image(
    ...     exp_time=3, ra=ra, dec=dec, sources=sources,
    ...     apply_pixel_defects=False, seed=42
    ... )

    >>> import matplotlib.pyplot as plt
    >>> fig, axes = plt.subplots(1, 2, figsize=(7, 5), sharex=True, sharey=True)
    >>> im0 = axes[0].imshow(clean_image, cmap='gray')
    >>> cbar0 = fig.colorbar(im0, ax=axes[0], orientation="horizontal", pad=0.2)
    >>> cbar0.set_label('ADU')
    >>> im1 = axes[1].imshow(image, cmap='gray')
    >>> cbar1 = fig.colorbar(im1, ax=axes[1], orientation="horizontal", pad=0.2)
    >>> cbar1.set_label('ADU')
    >>> axes[0].set_title('Image without defects')
    >>> axes[1].set_title('Image with constant defects')
    >>> axes[0].set_ylabel('Pixels')
    >>> for ax in axes:
    ...     ax.set_xlabel('Pixels')
    >>> plt.subplots_adjust(wspace=0.1)
    >>> plt.show()


    """

    value: int = 0
    seed: int = 0

    def introduce_pixel_defect(self, image, camera, **kwargs):
        if self._pixels is None:
            self.set_pixels(self._select_random_pixels(camera), camera)

        self._overwrite_pixel_values(image, self.pixels, self.value)

        return image


@dataclass
class TelegraphicPixelDefect(PixelDefect):
    """A pixel defect that simulates telegraphic noise by selecting entire rows or
    columns of pixels to be defective.

    Examples
    --------
    >>> from cabaret import Observatory, Sources
    >>> pixel_defects = {
    ...     "telegraphic": {"type": "telegraphic", "value": 5000, "rate": 0.5, "dim": 0}
    ... }
    >>> observatory = Observatory(
    ...     camera={"width": 80, "height": 80, "pixel_defects": pixel_defects}
    ... )
    >>> sources = Sources.get_test_source()
    >>> ra, dec = sources.ra.mean().deg, sources.dec.mean().deg,
    >>> image = observatory.generate_image(
    ...     exp_time=3, ra=ra, dec=dec, sources=sources, seed=42
    ... )
    >>> clean_image = observatory.generate_image(
    ...     exp_time=3, ra=ra, dec=dec, sources=sources,
    ...     apply_pixel_defects=False, seed=42
    ... )

    >>> import matplotlib.pyplot as plt
    >>> fig, axes = plt.subplots(1, 2, figsize=(7, 5), sharex=True, sharey=True)
    >>> im0 = axes[0].imshow(clean_image, cmap='gray')
    >>> cbar0 = fig.colorbar(im0, ax=axes[0], orientation="horizontal", pad=0.2)
    >>> cbar0.set_label('ADU')
    >>> im1 = axes[1].imshow(image, cmap='gray')
    >>> cbar1 = fig.colorbar(im1, ax=axes[1], orientation="horizontal", pad=0.2)
    >>> cbar1.set_label('ADU')
    >>> axes[0].set_title('Image without defects')
    >>> axes[1].set_title('Image with telegraphic defects')
    >>> axes[0].set_ylabel('Pixels')
    >>> for ax in axes:
    ...     ax.set_xlabel('Pixels')
    >>> plt.subplots_adjust(wspace=0.1)
    >>> plt.show()

    """

    value: int = 0
    seed: int = 0
    dim: int = 0
    _lines: np.ndarray | None = None

    def introduce_pixel_defect(self, image, camera, **kwargs):
        if self._lines is None:
            self.set_lines(self._select_random_lines(camera, self.dim), camera)

        self._overwrite_pixel_values(image, self.pixels, self.value)

        return image

    @property
    def lines(self) -> np.ndarray:
        if self._lines is None:
            raise ValueError(f"{self.name} lines are not defined.")
        return self._lines

    def set_lines(self, lines: np.ndarray | list, camera: Camera):
        self._lines = np.array(lines)

        # Set pixels based on the selected lines
        if self.dim == 0:
            if not np.all(self._lines < camera.height):
                raise ValueError("Selected lines are outside the frame.")
            X, Y = np.meshgrid(np.arange(camera.width), self._lines)
        else:
            if not np.all(self._lines < camera.width):
                raise ValueError("Selected lines are outside the frame.")
            X, Y = np.meshgrid(self._lines, np.arange(camera.width))

        self.set_pixels(np.column_stack((X.ravel(), Y.ravel())), camera)

    def _select_random_lines(self, camera: Camera, dim: int = 0) -> np.ndarray:
        line_length = camera.width if dim == 0 else camera.height
        number_of_lines = self.number_of_defect_pixels(camera) // line_length
        selected_lines = self.rng.integers(line_length, size=(number_of_lines))
        return selected_lines


@dataclass
class RandomNoisePixelDefect(PixelDefect):
    """A pixel defect that introduces random noise to selected pixels.

    This is different from read noise, which is applied to all pixels and represents
    the electronic noise of the sensor. RandomNoisePixelDefect simulates pixels that
    are abnormally noisy (e.g., "noisy" or "unstable" pixels), and the noise is
    added only to those selected pixels.

    Examples
    --------
    >>> from cabaret import Observatory, Sources
    >>> pixel_defects = {"noise": {"type": "noise", "rate": 0.01, "noise_level": 5e3}}
    >>> observatory = Observatory(
    ...     camera={"width": 80, "height": 80, "pixel_defects": pixel_defects}
    ... )
    >>> sources = Sources.get_test_source()

    >>> ra, dec = sources.ra.mean().deg, sources.dec.mean().deg,
    >>> image = observatory.generate_image(
    ...     exp_time=3, ra=ra, dec=dec, sources=sources
    ... )
    >>> clean_image = observatory.generate_image(
    ...     exp_time=3, ra=ra, dec=dec, sources=sources, apply_pixel_defects=False
    ... )

    >>> import matplotlib.pyplot as plt
    >>> fig, axes = plt.subplots(1, 2, figsize=(7, 5), sharex=True, sharey=True)
    >>> im0 = axes[0].imshow(clean_image, cmap='gray')
    >>> cbar0 = fig.colorbar(im0, ax=axes[0], orientation="horizontal", pad=0.2)
    >>> cbar0.set_label('ADU')
    >>> im1 = axes[1].imshow(image, cmap='gray')
    >>> cbar1 = fig.colorbar(im1, ax=axes[1], orientation="horizontal", pad=0.2)
    >>> cbar1.set_label('ADU')
    >>> axes[0].set_title('Image without defects')
    >>> axes[1].set_title('Image with random noise defects')
    >>> axes[0].set_ylabel('Pixels')
    >>> for ax in axes:
    ...     ax.set_xlabel('Pixels')
    >>> plt.subplots_adjust(wspace=0.1)
    >>> plt.show()

    """

    noise_level: float = 10.0  # Standard deviation for Gaussian noise
    distribution: Literal["normal", "poisson"] = "normal"

    def generate_noise(self, size: int) -> np.ndarray:
        if self.distribution == "poisson":
            return self.noise_level * self.rng.poisson(size=size)
        elif self.distribution == "normal":
            return self.noise_level * self.rng.normal(size=size)
        else:
            raise ValueError("Unknown noise distribution.")

    def introduce_pixel_defect(self, image, camera, seed: int | None = None):
        if seed is not None:
            self._rng = numpy.random.default_rng(seed)

        if self._pixels is None:
            self.set_pixels(self._select_random_pixels(camera), camera)

        noise = self.generate_noise(self.pixels.shape[0])

        # Add noise only to the selected defect pixels
        image[self.pixels[:, 0], self.pixels[:, 1]] += noise
        image = np.clip(image, 0, camera.max_adu)

        return image


@dataclass
class QuantumEfficiencyMapPixelDefect(PixelDefect):
    """
    A pixel defect that simulates variations in quantum efficiency across the sensor.

    Examples
    --------
    >>> from cabaret import Observatory, Sources
    >>> pixel_defects = {
    ...     "qe": {"type": "quantum_efficiency_map", "quantum_efficiency_std": 0.25}
    ... }
    >>> observatory = Observatory(
    ...     camera={"width": 40, "height": 40, "pixel_defects": pixel_defects}
    ... )
    >>> sources = Sources.get_test_source()

    >>> ra, dec = sources.ra.mean().deg, sources.dec.mean().deg,
    >>> image = observatory.generate_image(
    ...     exp_time=3, ra=ra, dec=dec, sources=sources
    ... )
    >>> clean_image = observatory.generate_image(
    ...     exp_time=3, ra=ra, dec=dec, sources=sources, apply_pixel_defects=False
    ... )


    >>> import matplotlib.pyplot as plt
    >>> fig, axes = plt.subplots(1, 3, figsize=(10, 5), sharex=True, sharey=True)
    >>> im0 = axes[0].imshow(clean_image, cmap='gray')
    >>> cbar0 = fig.colorbar(im0, ax=axes[0], orientation="horizontal", pad=0.2)
    >>> cbar0.set_label('ADU')
    >>> im1 = axes[1].imshow(image, cmap='gray')
    >>> cbar1 = fig.colorbar(im1, ax=axes[1], orientation="horizontal", pad=0.2)
    >>> cbar1.set_label('ADU')
    >>> im2 = axes[2].imshow(
    ...     observatory.camera.pixel_defects["qe"].quantum_efficiency_map, cmap='gray'
    ... )
    >>> cbar2 = fig.colorbar(im2, ax=axes[2], orientation="horizontal", pad=0.2)
    >>> cbar2.set_label('Quantum Efficiency')
    >>> axes[0].set_title('Image without defects')
    >>> axes[1].set_title('Image with QE defects')
    >>> axes[2].set_title('Quantum Efficiency Map')
    >>> axes[0].set_ylabel('Pixels')
    >>> for ax in axes:
    ...     ax.set_xlabel('Pixels')
    >>> plt.subplots_adjust(wspace=0.3, hspace=0.3)
    >>> plt.show()

    """

    quantum_efficiency_std: float = 0.1
    quantum_efficiency_map: np.ndarray | None = None
    seed: int = 0

    def __post_init__(self):
        if self.quantum_efficiency_map is not None:
            self.quantum_efficiency_map = np.array(self.quantum_efficiency_map)

    def generate_quantum_efficiency_map(self, camera: Camera):
        self.quantum_efficiency_map = np.clip(
            self.rng.normal(
                loc=camera.average_quantum_efficiency,
                scale=self.quantum_efficiency_std,
                size=(camera.height, camera.width),
            ),
            0,
            1,
        )

    def introduce_pixel_defect(
        self, image: np.ndarray, camera: Camera, seed: int | None = None
    ):
        if seed is not None:
            self._rng = numpy.random.default_rng(seed)
        if self.quantum_efficiency_map is None:
            self.generate_quantum_efficiency_map(camera)

        image = image * self.quantum_efficiency_map / camera.average_quantum_efficiency
        image = np.clip(image, 0, camera.max_adu)
        return image


@dataclass
class ChargeTransferInefficiencyPixelDefect(PixelDefect):
    """
    Simulates charge transfer inefficiency (CTI) by trailing charge along columns or
    rows.

    Each affected pixel leaks a fraction of its value to downstream pixels along the
    readout direction.

    Parameters
    ----------
    cti_fraction : float
        Fraction of charge transferred per pixel.
    dimension : int
        Direction of readout: 0 for columns (default), 1 for rows.
    seed : int
        Random seed for reproducibility.

    Examples
    --------
    >>> from cabaret import Observatory, Sources
    >>> pixel_defects = {
    ...     "cti": {"type": "cti", "cti_fraction": 0.05, "rate": 0.5, "dimension": 0}
    ... }
    >>> observatory = Observatory(
    ...     camera={"width": 40, "height": 40, "pixel_defects": pixel_defects}
    ... )
    >>> sources = Sources.get_test_source()
    >>> ra, dec = sources.ra.mean().deg, sources.dec.mean().deg,
    >>> image = observatory.generate_image(
    ...     exp_time=3, ra=ra, dec=dec, sources=sources
    ... )
    >>> clean_image = observatory.generate_image(
    ...     exp_time=3, ra=ra, dec=dec, sources=sources, apply_pixel_defects=False
    ... )

    >>> import matplotlib.pyplot as plt
    >>> fig, axes = plt.subplots(1, 2, figsize=(7, 5), sharex=True, sharey=True)
    >>> im0 = axes[0].imshow(clean_image, cmap='gray')
    >>> cbar0 = fig.colorbar(im0, ax=axes[0], orientation="horizontal", pad=0.2)
    >>> cbar0.set_label('ADU')
    >>> im1 = axes[1].imshow(image, cmap='gray')
    >>> cbar1 = fig.colorbar(im1, ax=axes[1], orientation="horizontal", pad=0.2)
    >>> cbar1.set_label('ADU')
    >>> axes[0].set_title('Image without defects')
    >>> axes[1].set_title('Image with CTI defects')
    >>> axes[0].set_ylabel('Pixels')
    >>> for ax in axes:
    ...     ax.set_xlabel('Pixels')
    >>> plt.subplots_adjust(wspace=0.1)
    >>> plt.show()
    """

    cti_fraction: float = 0.01  # Fraction of charge transferred per pixel
    dimension: int = 0  # 0: columns (default), 1: rows
    seed: int = 0

    def introduce_pixel_defect(self, image, camera, **kwargs):
        if self._pixels is None:
            self.set_pixels(self._select_random_pixels(camera), camera)
        # For each affected pixel, trail charge along the selected dimension
        if self.dimension == 0:
            # Trail down columns (y increases)
            for y, x in self.pixels:
                for y2 in range(y + 1, camera.height):
                    image[y2, x] += image[y, x] * self.cti_fraction
        elif self.dimension == 1:
            # Trail along rows (x increases)
            for y, x in self.pixels:
                for x2 in range(x + 1, camera.width):
                    image[y, x2] += image[y, x] * self.cti_fraction
        else:
            raise ValueError("dimension must be 0 (columns) or 1 (rows)")
        return image


@dataclass
class ReadoutSmearPixelDefect(PixelDefect):
    """
    Simulates readout smear (frame transfer smear) caused by continued exposure during
    readout.

    Parameters
    ----------
    smear_fraction : float
        Fraction of the exposure time spent during readout (smear strength).
    dim : int
        Direction of readout: 0 for vertical (default), 1 for horizontal.
    readout_time : float
        Time taken to read out the entire frame (seconds).
    seed : int
        Random seed for reproducibility.

    Examples
    --------
    >>> from cabaret import Observatory, Sources
    >>> pixel_defects = {
    ...     "smear": {"type": "readout_smear", "readout_time": 1}
    ... }
    >>> observatory = Observatory(
    ...     camera={"width": 40, "height": 40, "pixel_defects": pixel_defects}
    ... )
    >>> sources = Sources.get_test_source()
    >>> ra, dec = sources.ra.mean().deg, sources.dec.mean().deg,
    >>> image = observatory.generate_image(
    ...     exp_time=1, ra=ra, dec=dec, sources=sources
    ... )
    >>> clean_image = observatory.generate_image(
    ...     exp_time=1, ra=ra, dec=dec, sources=sources, apply_pixel_defects=False
    ... )

    >>> import matplotlib.pyplot as plt
    >>> fig, axes = plt.subplots(1, 2, figsize=(7, 5), sharex=True, sharey=True)
    >>> im0 = axes[0].imshow(clean_image, cmap='gray')
    >>> cbar0 = fig.colorbar(im0, ax=axes[0], orientation="horizontal", pad=0.2)
    >>> cbar0.set_label('ADU')
    >>> im1 = axes[1].imshow(image, cmap='gray')
    >>> cbar1 = fig.colorbar(im1, ax=axes[1], orientation="horizontal", pad=0.2)
    >>> cbar1.set_label('ADU')
    >>> axes[0].set_title('Image without defects')
    >>> axes[1].set_title('Image with readout smear')
    >>> axes[0].set_ylabel('Pixels')
    >>> for ax in axes:
    ...     ax.set_xlabel('Pixels')
    >>> plt.subplots_adjust(wspace=0.1)
    >>> plt.show()
    """

    readout_time: float = 0.1
    dim: int = 0
    seed: int = 0

    def __post__init__(self):
        if self.dim not in (0, 1):
            raise ValueError("dim must be 0 (vertical) or 1 (horizontal)")

    def introduce_pixel_defect(self, image, camera, exp_time, **kwargs):
        # Fraction of the exposure time that each pixel is exposed
        # during the readout of a single row/column, assuming uniform readout speed
        smear_fraction_per_pixel_readout = (
            self.readout_time / exp_time / image.shape[self.dim]
        )

        # A fraction of each pixels values to those read out after it
        smear = smear_fraction_per_pixel_readout * np.cumsum(image, axis=self.dim)

        image = np.clip(image + smear, 0, camera.max_adu)
        return image
