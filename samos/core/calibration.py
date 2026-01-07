"""
Calibration Frame Processing

Functions for creating and applying calibration frames (bias, dark, flat).
"""

import numpy as np
from typing import List, Optional
from astropy.io import fits


def create_master_bias(bias_frames: List[np.ndarray],
                      method: str = 'median') -> np.ndarray:
    """
    Create master bias frame from multiple bias exposures.

    Parameters
    ----------
    bias_frames : list of np.ndarray
        List of bias frame arrays
    method : str, optional
        Combination method: 'median' or 'mean' (default: 'median')

    Returns
    -------
    np.ndarray
        Master bias frame

    Examples
    --------
    >>> bias_frames = [read_fits(f) for f in bias_files]
    >>> master_bias = create_master_bias(bias_frames)
    """
    stack = np.array(bias_frames)

    if method == 'median':
        master = np.median(stack, axis=0)
    elif method == 'mean':
        master = np.mean(stack, axis=0)
    else:
        raise ValueError(f"Unknown method: {method}")

    return master


def create_master_flat(flat_frames: List[np.ndarray],
                      master_bias: Optional[np.ndarray] = None,
                      method: str = 'median',
                      normalize: bool = True) -> np.ndarray:
    """
    Create master flat field frame.

    Parameters
    ----------
    flat_frames : list of np.ndarray
        List of flat field arrays
    master_bias : np.ndarray, optional
        Master bias to subtract before combining
    method : str, optional
        Combination method (default: 'median')
    normalize : bool, optional
        Normalize flat to median=1 (default: True)

    Returns
    -------
    np.ndarray
        Master flat field frame

    Examples
    --------
    >>> flat_frames = [read_fits(f) for f in flat_files]
    >>> master_flat = create_master_flat(flat_frames, master_bias)
    """
    # Bias subtract if provided
    if master_bias is not None:
        flats_corrected = [flat - master_bias for flat in flat_frames]
    else:
        flats_corrected = flat_frames

    # Combine
    stack = np.array(flats_corrected)

    if method == 'median':
        master = np.median(stack, axis=0)
    elif method == 'mean':
        master = np.mean(stack, axis=0)
    else:
        raise ValueError(f"Unknown method: {method}")

    # Normalize
    if normalize:
        median_value = np.median(master)
        master = master / median_value

    return master


def apply_bias_correction(image: np.ndarray,
                         master_bias: np.ndarray) -> np.ndarray:
    """
    Apply bias correction to an image.

    Parameters
    ----------
    image : np.ndarray
        Raw image
    master_bias : np.ndarray
        Master bias frame

    Returns
    -------
    np.ndarray
        Bias-corrected image

    Examples
    --------
    >>> corrected = apply_bias_correction(science_frame, master_bias)
    """
    return image - master_bias


def apply_flat_correction(image: np.ndarray,
                         master_flat: np.ndarray,
                         min_value: float = 0.1) -> np.ndarray:
    """
    Apply flat field correction to an image.

    Parameters
    ----------
    image : np.ndarray
        Bias-corrected image
    master_flat : np.ndarray
        Normalized master flat frame
    min_value : float, optional
        Minimum flat value to avoid division by near-zero (default: 0.1)

    Returns
    -------
    np.ndarray
        Flat-corrected image

    Examples
    --------
    >>> corrected = apply_flat_correction(bias_corrected, master_flat)
    """
    # Avoid division by very small values
    flat_safe = np.where(master_flat < min_value, 1.0, master_flat)

    return image / flat_safe


def reduce_image(image: np.ndarray,
                master_bias: Optional[np.ndarray] = None,
                master_flat: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Apply full calibration pipeline to an image.

    Parameters
    ----------
    image : np.ndarray
        Raw image
    master_bias : np.ndarray, optional
        Master bias frame
    master_flat : np.ndarray, optional
        Master flat frame

    Returns
    -------
    np.ndarray
        Calibrated image

    Examples
    --------
    >>> calibrated = reduce_image(raw_image, master_bias, master_flat)
    """
    result = image.copy()

    if master_bias is not None:
        result = apply_bias_correction(result, master_bias)

    if master_flat is not None:
        result = apply_flat_correction(result, master_flat)

    return result
