"""
Image Display Utilities

Functions for displaying astronomical images with appropriate scaling,
colormaps, and styling.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
from typing import Optional, Tuple


def display_image(image: np.ndarray,
                 zmin: Optional[float] = None,
                 zmax: Optional[float] = None,
                 title: str = "",
                 cmap: str = 'gray',
                 figsize: Tuple[float, float] = (14, 8),
                 colorbar: bool = True,
                 origin: str = 'lower') -> None:
    """
    Display an astronomical image with proper scaling and styling.

    Parameters
    ----------
    image : np.ndarray
        2D image array to display
    zmin : float, optional
        Minimum value for display scaling. If None, uses image min
    zmax : float, optional
        Maximum value for display scaling. If None, uses image max
    title : str, optional
        Plot title (default: "")
    cmap : str, optional
        Matplotlib colormap name (default: 'gray')
    figsize : tuple of float, optional
        Figure size (width, height) in inches (default: (14, 8))
    colorbar : bool, optional
        Whether to show colorbar (default: True)
    origin : str, optional
        Image origin ('lower' or 'upper') (default: 'lower')

    Examples
    --------
    >>> from samos.utils import display
    >>> display.display_image(image, zmin=500, zmax=3000, title='Science Frame')

    >>> # Auto-scale based on percentiles
    >>> vmin, vmax = display.get_percentile_limits(image, 1, 99)
    >>> display.display_image(image, zmin=vmin, zmax=vmax)
    """
    # Use astropy styling
    plt.style.use(astropy_mpl_style)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Display image
    im = ax.imshow(image, origin=origin, cmap=cmap, vmin=zmin, vmax=zmax)

    # Add title
    if title:
        ax.set_title(title)

    # Add colorbar scaled to figure height
    if colorbar:
        im_ratio = image.shape[0] / image.shape[1]
        plt.colorbar(im, ax=ax, fraction=0.046 * im_ratio, pad=0.04)

    plt.tight_layout()
    plt.show()


def get_percentile_limits(image: np.ndarray,
                         lower_percentile: float = 1.0,
                         upper_percentile: float = 99.0) -> Tuple[float, float]:
    """
    Calculate display limits based on image percentiles.

    This is useful for auto-scaling images with outliers or varying
    dynamic range.

    Parameters
    ----------
    image : np.ndarray
        2D image array
    lower_percentile : float, optional
        Lower percentile for minimum (default: 1.0)
    upper_percentile : float, optional
        Upper percentile for maximum (default: 99.0)

    Returns
    -------
    vmin : float
        Lower display limit
    vmax : float
        Upper display limit

    Examples
    --------
    >>> vmin, vmax = get_percentile_limits(image, 2, 98)
    >>> display_image(image, zmin=vmin, zmax=vmax)
    """
    vmin = np.percentile(image, lower_percentile)
    vmax = np.percentile(image, upper_percentile)

    return vmin, vmax


def get_sigma_limits(image: np.ndarray,
                    n_sigma: float = 3.0) -> Tuple[float, float]:
    """
    Calculate display limits based on sigma clipping.

    Parameters
    ----------
    image : np.ndarray
        2D image array
    n_sigma : float, optional
        Number of standard deviations from median (default: 3.0)

    Returns
    -------
    vmin : float
        Lower display limit (median - n_sigma * std)
    vmax : float
        Upper display limit (median + n_sigma * std)

    Examples
    --------
    >>> vmin, vmax = get_sigma_limits(image, n_sigma=5.0)
    >>> display_image(image, zmin=vmin, zmax=vmax)
    """
    median = np.median(image)
    std = np.std(image)

    vmin = median - n_sigma * std
    vmax = median + n_sigma * std

    return vmin, vmax


def display_comparison(images: list,
                      titles: list,
                      zmin: Optional[float] = None,
                      zmax: Optional[float] = None,
                      cmap: str = 'gray',
                      figsize: Tuple[float, float] = (16, 6)) -> None:
    """
    Display multiple images side-by-side for comparison.

    Parameters
    ----------
    images : list of np.ndarray
        List of 2D image arrays to display
    titles : list of str
        List of titles for each image
    zmin : float, optional
        Minimum value for display scaling (same for all images)
    zmax : float, optional
        Maximum value for display scaling (same for all images)
    cmap : str, optional
        Matplotlib colormap name (default: 'gray')
    figsize : tuple of float, optional
        Figure size (default: (16, 6))

    Examples
    --------
    >>> display_comparison([raw, cleaned], ['Raw', 'CR Cleaned'])
    >>> display_comparison([bias, flat, science],
    ...                   ['Bias', 'Flat', 'Science'],
    ...                   zmin=0, zmax=5000)
    """
    n_images = len(images)

    if len(titles) != n_images:
        raise ValueError("Number of titles must match number of images")

    plt.style.use(astropy_mpl_style)
    fig, axes = plt.subplots(1, n_images, figsize=figsize)

    if n_images == 1:
        axes = [axes]

    for ax, image, title in zip(axes, images, titles):
        im = ax.imshow(image, origin='lower', cmap=cmap, vmin=zmin, vmax=zmax)
        ax.set_title(title)

        im_ratio = image.shape[0] / image.shape[1]
        plt.colorbar(im, ax=ax, fraction=0.046 * im_ratio, pad=0.04)

    plt.tight_layout()
    plt.show()


def save_image(image: np.ndarray,
              filepath: str,
              zmin: Optional[float] = None,
              zmax: Optional[float] = None,
              title: str = "",
              cmap: str = 'gray',
              figsize: Tuple[float, float] = (14, 8),
              dpi: int = 300) -> None:
    """
    Save an image to a file.

    Parameters
    ----------
    image : np.ndarray
        2D image array
    filepath : str
        Output file path (e.g., 'output.png', 'output.pdf')
    zmin : float, optional
        Minimum display value
    zmax : float, optional
        Maximum display value
    title : str, optional
        Plot title
    cmap : str, optional
        Colormap (default: 'gray')
    figsize : tuple of float, optional
        Figure size in inches (default: (14, 8))
    dpi : int, optional
        Resolution in dots per inch (default: 300)

    Examples
    --------
    >>> save_image(image, 'science_frame.png', zmin=500, zmax=3000)
    """
    plt.style.use(astropy_mpl_style)
    fig, ax = plt.subplots(figsize=figsize)

    im = ax.imshow(image, origin='lower', cmap=cmap, vmin=zmin, vmax=zmax)

    if title:
        ax.set_title(title)

    im_ratio = image.shape[0] / image.shape[1]
    plt.colorbar(im, ax=ax, fraction=0.046 * im_ratio, pad=0.04)

    plt.tight_layout()
    plt.savefig(filepath, dpi=dpi, bbox_inches='tight')
    plt.close(fig)

    print(f"Saved image to: {filepath}")
