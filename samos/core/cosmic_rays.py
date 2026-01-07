"""
Cosmic Ray Detection and Removal

This module provides functions for detecting and removing cosmic ray hits
from CCD images using the L.A.Cosmic algorithm.
"""

import numpy as np
from lacosmic.core import lacosmic
from typing import Tuple


def remove_cosmic_rays(image: np.ndarray,
                      contrast: float = 5.0,
                      cr_threshold: float = 20.0,
                      neighbor_threshold: float = 5.0,
                      readnoise: float = 3.0,
                      effective_gain: float = 2.0,
                      verbose: bool = False) -> np.ndarray:
    """
    Remove cosmic rays from an image using the L.A.Cosmic algorithm.

    This function wraps the lacosmic package with sensible defaults for
    SAMOS data. The L.A.Cosmic algorithm identifies cosmic rays by their
    sharp, isolated peaks that differ from astronomical sources.

    Parameters
    ----------
    image : np.ndarray
        2D image array to clean
    contrast : float, optional
        Contrast threshold between CR and underlying object (default: 5.0)
        Higher values = more conservative (fewer false positives)
    cr_threshold : float, optional
        Detection threshold in standard deviations (default: 20.0)
        Higher values = fewer detections
    neighbor_threshold : float, optional
        Threshold for neighboring pixels (default: 5.0)
    readnoise : float, optional
        CCD readout noise in electrons (default: 3.0)
    effective_gain : float, optional
        Effective gain in e-/ADU (default: 2.0)
    verbose : bool, optional
        Print detection statistics (default: False)

    Returns
    -------
    np.ndarray
        Cleaned image with cosmic rays replaced by interpolated values

    Examples
    --------
    >>> from samos.core import cosmic_rays
    >>> cleaned = cosmic_rays.remove_cosmic_rays(image)

    >>> # More aggressive cleaning
    >>> cleaned = cosmic_rays.remove_cosmic_rays(image, cr_threshold=10.0)

    Notes
    -----
    The L.A.Cosmic algorithm (van Dokkum 2001, PASP 113, 1420) detects
    cosmic rays using Laplacian edge detection. Detected pixels are
    replaced by the median of surrounding unaffected pixels.

    References
    ----------
    van Dokkum, P. 2001, PASP, 113, 1420
    """
    cleaned_data, mask = lacosmic(
        image,
        contrast=contrast,
        cr_threshold=cr_threshold,
        neighbor_threshold=neighbor_threshold,
        readnoise=readnoise,
        effective_gain=effective_gain
    )

    if verbose:
        n_pixels = np.sum(mask)
        fraction = n_pixels / mask.size * 100
        print(f"Cosmic rays detected: {n_pixels} pixels ({fraction:.2f}%)")

    return cleaned_data


def remove_cosmic_rays_with_mask(image: np.ndarray,
                                 contrast: float = 5.0,
                                 cr_threshold: float = 20.0,
                                 neighbor_threshold: float = 5.0,
                                 readnoise: float = 3.0,
                                 effective_gain: float = 2.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Remove cosmic rays and return both cleaned image and detection mask.

    Parameters
    ----------
    image : np.ndarray
        2D image array to clean
    contrast : float, optional
        Contrast threshold (default: 5.0)
    cr_threshold : float, optional
        Detection threshold in sigma (default: 20.0)
    neighbor_threshold : float, optional
        Threshold for neighboring pixels (default: 5.0)
    readnoise : float, optional
        CCD readout noise in electrons (default: 3.0)
    effective_gain : float, optional
        Effective gain in e-/ADU (default: 2.0)

    Returns
    -------
    cleaned_data : np.ndarray
        Cleaned image
    mask : np.ndarray
        Boolean mask where True indicates cosmic ray pixels

    Examples
    --------
    >>> cleaned, mask = cosmic_rays.remove_cosmic_rays_with_mask(image)
    >>> print(f"Found {mask.sum()} cosmic ray pixels")
    """
    cleaned_data, mask = lacosmic(
        image,
        contrast=contrast,
        cr_threshold=cr_threshold,
        neighbor_threshold=neighbor_threshold,
        readnoise=readnoise,
        effective_gain=effective_gain
    )

    return cleaned_data, mask


def get_recommended_parameters(instrument: str = 'samos') -> dict:
    """
    Get recommended cosmic ray removal parameters for specific instruments.

    Parameters
    ----------
    instrument : str, optional
        Instrument name (default: 'samos')
        Options: 'samos', 'conservative', 'aggressive'

    Returns
    -------
    dict
        Dictionary of recommended parameters

    Examples
    --------
    >>> params = cosmic_rays.get_recommended_parameters('samos')
    >>> cleaned = cosmic_rays.remove_cosmic_rays(image, **params)
    """
    presets = {
        'samos': {
            'contrast': 5.0,
            'cr_threshold': 20.0,
            'neighbor_threshold': 5.0,
            'readnoise': 3.0,
            'effective_gain': 2.0
        },
        'conservative': {
            'contrast': 7.0,
            'cr_threshold': 25.0,
            'neighbor_threshold': 7.0,
            'readnoise': 3.0,
            'effective_gain': 2.0
        },
        'aggressive': {
            'contrast': 3.0,
            'cr_threshold': 10.0,
            'neighbor_threshold': 3.0,
            'readnoise': 3.0,
            'effective_gain': 2.0
        }
    }

    if instrument not in presets:
        raise ValueError(f"Unknown instrument: {instrument}. "
                        f"Available: {list(presets.keys())}")

    return presets[instrument]
