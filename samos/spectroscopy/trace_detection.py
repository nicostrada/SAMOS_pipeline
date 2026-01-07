"""
Spectroscopic Trace Detection

Functions for identifying spectral traces (slits) in multi-object spectroscopy data.
"""

import numpy as np
from scipy.signal import find_peaks
import copy
from typing import List, Tuple, Optional


def find_slit_edges(image: np.ndarray,
                   x_range: Tuple[int, int] = (1500, 2000),
                   peak_threshold: float = 4e5,
                   distance: int = 3) -> Tuple[np.ndarray, np.ndarray]:
    """
    Find slit edges by detecting peaks in the collapsed spatial profile.

    This function collapses the image along the dispersion direction within
    a specified range, then identifies slit edges as peaks in the derivative
    of the spatial profile.

    Parameters
    ----------
    image : np.ndarray
        2D image (typically arc lamp exposure)
    x_range : tuple of int, optional
        (x_min, x_max) range for collapsing in dispersion direction
        (default: (1500, 2000))
    peak_threshold : float, optional
        Minimum peak height for edge detection (default: 4e5)
    distance : int, optional
        Minimum separation between peaks in pixels (default: 3)

    Returns
    -------
    slit_up : np.ndarray
        Y-coordinates of upper slit edges
    slit_down : np.ndarray
        Y-coordinates of lower slit edges

    Examples
    --------
    >>> from samos.spectroscopy import trace_detection
    >>> slit_up, slit_down = trace_detection.find_slit_edges(arc_image)
    >>> print(f"Found {len(slit_up)} slits")
    """
    # Collapse image along dispersion direction
    cut = image[:, x_range[0]:x_range[1]].sum(axis=1)

    # Compute derivative to find edges
    cut1 = np.roll(cut, +1)
    diff = cut - cut1

    # Separate positive and negative edges
    diff_neg = copy.deepcopy(diff)
    diff_neg[diff > 0] = 0

    diff_pos = copy.deepcopy(diff)
    diff_pos[diff < 0] = 0

    # Find peaks
    v_slit_up, _ = find_peaks(diff_pos, distance=distance, height=peak_threshold)
    v_slit_down, _ = find_peaks(-diff_neg, distance=distance, height=peak_threshold)

    return v_slit_up, v_slit_down


def resolve_overlapping_slits(slit_up: np.ndarray,
                              slit_down: np.ndarray,
                              verbose: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """
    Remove overlapping or conflicting slit definitions.

    Identifies and removes slits where the edges overlap or are too close
    to each other, which can occur when slits are closely packed.

    Parameters
    ----------
    slit_up : np.ndarray
        Y-coordinates of upper slit edges
    slit_down : np.ndarray
        Y-coordinates of lower slit edges
    verbose : bool, optional
        Print information about removed slits (default: True)

    Returns
    -------
    slit_up : np.ndarray
        Cleaned upper slit edges
    slit_down : np.ndarray
        Cleaned lower slit edges

    Examples
    --------
    >>> slit_up, slit_down = resolve_overlapping_slits(slit_up, slit_down)
    """
    to_be_removed = []

    for i in range(1, len(slit_up)):
        try:
            # Check for overlap with next slit
            if i + 1 < len(slit_up) and slit_down[i] > slit_up[i + 1]:
                if verbose:
                    print(f"Slit {i}: {slit_up[i]}-{slit_down[i]} conflicts with next slit, removing")
                to_be_removed.append(i)

            # Check for overlap with previous slit
            elif slit_down[i - 1] > slit_up[i]:
                if verbose:
                    print(f"Slit {i}: {slit_up[i]}-{slit_down[i]} conflicts with previous slit, removing")
                to_be_removed.append(i)

        except IndexError:
            continue

    # Remove conflicting entries
    slit_up = np.delete(slit_up, to_be_removed, axis=0)
    slit_down = np.delete(slit_down, to_be_removed, axis=0)

    if verbose:
        print(f"\nRemoved {len(to_be_removed)} overlapping slits")
        print(f"Final count: {len(slit_up)} slits")

    return slit_up, slit_down


def define_slit_regions(slit_up: np.ndarray,
                       slit_down: np.ndarray,
                       margin: int = 12) -> List[List[int]]:
    """
    Define slit extraction regions including margins for sky subtraction.

    For each slit, defines four y-coordinates:
    - s0: Lower edge of extraction region (including sky margin)
    - s1: Upper edge of spectrum
    - s2: Lower edge of spectrum
    - s3: Upper edge of extraction region (including sky margin)

    Parameters
    ----------
    slit_up : np.ndarray
        Y-coordinates of upper slit edges
    slit_down : np.ndarray
        Y-coordinates of lower slit edges
    margin : int, optional
        Margin in pixels for sky regions above/below spectrum (default: 12)

    Returns
    -------
    list of list of int
        List of slit regions, each element is [s0, s1, s2, s3]

    Examples
    --------
    >>> slits = define_slit_regions(slit_up, slit_down, margin=12)
    >>> print(f"Slit 0: {slits[0]}")  # [s0, s1, s2, s3]
    >>> # Extract slit 0 from image:
    >>> slit_cutout = image[slits[0][0]:slits[0][3], :]
    """
    slits = []

    for i in range(len(slit_up)):
        # Determine lower margin
        if i == 0:
            s0 = slit_up[i] - margin
        else:
            s0 = max(slit_down[i - 1], slit_up[i] - margin)

        # Determine upper margin
        if i == len(slit_down) - 1:
            s3 = slit_down[i] + margin
        else:
            s3 = min(slit_up[i + 1], slit_down[i] + margin)

        slits.append([s0, slit_up[i], slit_down[i], s3])

    return slits


def extract_slit_cutouts(image: np.ndarray,
                        slits: List[List[int]],
                        extension_names: Optional[List[str]] = None) -> List[np.ndarray]:
    """
    Extract 2D cutouts for each slit from a full frame image.

    Parameters
    ----------
    image : np.ndarray
        Full frame 2D image
    slits : list of list of int
        Slit definitions from define_slit_regions()
    extension_names : list of str, optional
        Names for each slit (for logging)

    Returns
    -------
    list of np.ndarray
        List of 2D cutouts, one per slit

    Examples
    --------
    >>> cutouts = extract_slit_cutouts(science_frame, slits)
    >>> print(f"Slit 0 shape: {cutouts[0].shape}")
    """
    cutouts = []

    for i, slit in enumerate(slits):
        cutout = image[slit[0]:slit[3], :]
        cutouts.append(cutout)

        if extension_names:
            print(f"{extension_names[i]}: shape {cutout.shape}")
        else:
            print(f"Slit {i}: shape {cutout.shape}, y-range [{slit[0]}:{slit[3]}]")

    return cutouts


def detect_traces_from_arc(arc_image: np.ndarray,
                          arc_dark: Optional[np.ndarray] = None,
                          x_range: Tuple[int, int] = (1500, 2000),
                          peak_threshold: float = 4e5,
                          margin: int = 12,
                          verbose: bool = True) -> List[List[int]]:
    """
    Complete pipeline to detect spectral traces from an arc lamp image.

    Parameters
    ----------
    arc_image : np.ndarray
        Arc lamp exposure
    arc_dark : np.ndarray, optional
        Arc lamp dark frame to subtract (e.g., DMD off exposure)
    x_range : tuple of int, optional
        Region for edge detection (default: (1500, 2000))
    peak_threshold : float, optional
        Minimum peak height (default: 4e5)
    margin : int, optional
        Sky margin in pixels (default: 12)
    verbose : bool, optional
        Print diagnostic information (default: True)

    Returns
    -------
    list of list of int
        Slit definitions [s0, s1, s2, s3] for each detected trace

    Examples
    --------
    >>> # Simple case: arc only
    >>> slits = detect_traces_from_arc(arc_image)

    >>> # With dark subtraction
    >>> slits = detect_traces_from_arc(arc_image, arc_dark=arc_dark_image)
    """
    # Subtract dark if provided
    if arc_dark is not None:
        arc_clean = arc_image - arc_dark
        if verbose:
            print("Arc dark subtracted")
    else:
        arc_clean = arc_image

    # Find edges
    if verbose:
        print(f"\nDetecting edges in x-range {x_range}...")

    slit_up, slit_down = find_slit_edges(arc_clean, x_range, peak_threshold)

    if verbose:
        print(f"Found {len(slit_up)} upper edges and {len(slit_down)} lower edges")

    # Resolve overlaps
    slit_up, slit_down = resolve_overlapping_slits(slit_up, slit_down, verbose=verbose)

    # Define regions
    slits = define_slit_regions(slit_up, slit_down, margin=margin)

    if verbose:
        print(f"\nFinal slit definitions:")
        for i, slit in enumerate(slits):
            dy = slit[3] - slit[0]
            spectrum_dy = slit[2] - slit[1]
            print(f"Slit {i:2d}: y=[{slit[0]:4d}:{slit[3]:4d}] "
                  f"(total dy={dy:2d}, spectrum dy={spectrum_dy:2d})")

    return slits
