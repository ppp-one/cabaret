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
    transparent=True,
    contrast=0.25,
):
    """Plot a 2D image with zscale normalization.

    Parameters
    ----------
    image : np.ndarray
        The image to plot.
    title : str, optional
        Title for the plot.
    cmap : str, optional
        Colormap to use.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on.
    add_colorbar : bool, optional
        Whether to add a colorbar.
    colorbar_kwargs : dict, optional
        Additional kwargs for colorbar.
    transparent : bool, optional
        If True, set figure background to transparent.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))
    else:
        fig = plt.gcf()
    interval = ZScaleInterval(contrast=contrast)
    vmin, vmax = interval.get_limits(image)
    img = ax.imshow(image, vmin=vmin, vmax=vmax, cmap=cmap, rasterized=True)
    if title is not None:
        ax.set_title(title)
    ax.set_xlabel("Pixels")
    ax.set_ylabel("Pixels")

    if add_colorbar:
        cbar = plt.colorbar(
            img, ax=ax, **{"fraction": 0.046, "pad": 0.04} | colorbar_kwargs
        )
        cbar.set_label("ADU")

    if transparent:
        fig.patch.set_alpha(0)
        ax.patch.set_alpha(0)

    return img
