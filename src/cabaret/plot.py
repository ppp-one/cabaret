import importlib.util

if not importlib.util.find_spec("matplotlib"):
    raise ImportError("Please install matplotlib to use cabaret.utils.")


import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval


def plot_image(
    image,
    title: str | None = None,
    cmap="gray",
    ax=None,
    add_colorbar=True,
    colorbar_kwargs={},
):
    """Plot a 2D image with zscale normalization."""
    if ax is None:
        _, ax = plt.subplots()
    interval = ZScaleInterval()
    vmin, vmax = interval.get_limits(image)
    img = ax.imshow(image, vmin=vmin, vmax=vmax, cmap=cmap)
    if title is not None:
        ax.set_title(title)
    ax.set_xlabel("Pixels")
    ax.set_ylabel("Pixels")

    if add_colorbar:
        cbar = plt.colorbar(
            img, ax=ax, **{"fraction": 0.046, "pad": 0.04} | colorbar_kwargs
        )
        cbar.set_label("ADU")

    return img
