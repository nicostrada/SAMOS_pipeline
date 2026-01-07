"""
Spectral Extraction and Trace Fitting

This module provides functions for extracting 1D spectra from 2D spectral images,
including trace fitting, rectification, and optimal extraction.
"""

import numpy as np
from scipy import ndimage, signal, interpolate
from typing import Tuple, Optional, Dict
import warnings


def find_trace_edges(image: np.ndarray,
                     x_range: Optional[Tuple[int, int]] = None,
                     threshold_fraction: float = 0.2) -> Tuple[np.ndarray, np.ndarray]:
    """
    Find slit edges in a 2D spectrum using edge detection.

    Parameters
    ----------
    image : np.ndarray
        2D spectral image (usually flat field or arc lamp)
    x_range : tuple, optional
        (x_start, x_end) region to analyze for edge detection
    threshold_fraction : float, optional
        Fraction of max gradient to consider as edge (default: 0.2)

    Returns
    -------
    x_pixels : np.ndarray
        X pixel positions where edges were detected
    y_centers : np.ndarray
        Y pixel positions of trace center

    Examples
    --------
    >>> x_pix, y_ctr = find_trace_edges(flat_field, x_range=(1500, 2000))
    """
    if x_range is None:
        x_range = (0, image.shape[1])

    x_start, x_end = x_range
    dx, dy = image.shape[1], image.shape[0]

    x_pixels, y_centers = [], []

    for ix in range(x_start, x_end):
        profile = image[:, ix]

        # Find edges using gradient
        gradient = np.gradient(profile)
        threshold = (np.max(np.abs(gradient))) * threshold_fraction

        edges = np.where(np.abs(gradient) > threshold)[0]

        if len(edges) >= 2:
            # Use first and last edge
            e0, e1 = edges[0], edges[-1]
            x_pixels.append(ix)
            y_centers.append((e0 + e1) / 2.0)

    return np.array(x_pixels), np.array(y_centers)


def fit_trace(x_pixels: np.ndarray,
             y_pixels: np.ndarray,
             method: str = 'polynomial',
             degree: int = 3,
             sigma_clip: float = 3.0) -> np.poly1d:
    """
    Fit a polynomial or spline to trace curvature.

    Parameters
    ----------
    x_pixels : np.ndarray
        X pixel positions
    y_pixels : np.ndarray
        Y pixel positions (trace center)
    method : str, optional
        Fitting method: 'polynomial' or 'spline' (default: 'polynomial')
    degree : int, optional
        Polynomial degree or spline order (default: 3)
    sigma_clip : float, optional
        Sigma clipping threshold for outlier rejection (default: 3.0)

    Returns
    -------
    trace_poly : np.poly1d
        Polynomial function describing the trace

    Examples
    --------
    >>> x, y = find_trace_edges(flat, x_range=(1500, 2000))
    >>> trace = fit_trace(x, y, degree=3)
    >>> y_fitted = trace(x)
    """
    # Sigma clipping to remove outliers
    mask = np.ones(len(y_pixels), dtype=bool)

    if sigma_clip > 0:
        for _ in range(3):  # Iterate sigma clipping
            coeffs = np.polyfit(x_pixels[mask], y_pixels[mask], degree)
            poly = np.poly1d(coeffs)
            residuals = y_pixels[mask] - poly(x_pixels[mask])
            sigma = np.std(residuals)
            mask[mask] = np.abs(residuals) < sigma_clip * sigma

    # Final fit
    if method == 'polynomial':
        coeffs = np.polyfit(x_pixels[mask], y_pixels[mask], degree)
        trace_poly = np.poly1d(coeffs)
    elif method == 'spline':
        # Use UnivariateSpline
        spline = interpolate.UnivariateSpline(
            x_pixels[mask], y_pixels[mask], k=min(degree, 5), s=1.0
        )
        # Convert to polynomial approximation for consistency
        x_fit = np.linspace(x_pixels[0], x_pixels[-1], 1000)
        y_fit = spline(x_fit)
        coeffs = np.polyfit(x_fit, y_fit, degree)
        trace_poly = np.poly1d(coeffs)
    else:
        raise ValueError(f"Unknown method: {method}. Use 'polynomial' or 'spline'")

    return trace_poly


def create_spatial_mask(image_shape: Tuple[int, int],
                       trace_poly: np.poly1d,
                       slit_edges: Tuple[int, int, int, int],
                       margin: int = 12) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create masks separating slit (illuminated) from sky (background) regions.

    Parameters
    ----------
    image_shape : tuple
        (ny, nx) shape of the image
    trace_poly : np.poly1d
        Polynomial describing trace center
    slit_edges : tuple
        (y0, y1_trace, y2_trace, y3) slit boundaries
    margin : int, optional
        Sky margin in pixels (default: 12)

    Returns
    -------
    mask_slit : np.ndarray
        Binary mask where True = slit (illuminated region)
    mask_sky : np.ndarray
        Binary mask where True = sky (background region)

    Examples
    --------
    >>> mask_slit, mask_sky = create_spatial_mask(
    ...     flat.shape, trace_poly, slit_edges=[895, 907, 925, 937]
    ... )
    """
    ny, nx = image_shape
    y0, y1, y2, y3 = slit_edges

    mask_slit = np.zeros(image_shape, dtype=bool)
    mask_sky = np.zeros(image_shape, dtype=bool)

    # Define slit width relative to trace
    below = y1 - y0 - margin
    above = y3 - y2 - margin

    # Fill masks column by column following the trace
    for ix in range(nx):
        trace_y = int(trace_poly(ix))

        # Slit mask (illuminated region)
        y_start = max(0, trace_y - below)
        y_end = min(ny, trace_y + above)
        mask_slit[y_start:y_end, ix] = True

        # Sky mask (background regions)
        mask_sky[0:max(0, trace_y - below), ix] = True
        mask_sky[min(ny, trace_y + above):, ix] = True

    return mask_slit, mask_sky


def rectify_spectrum(data_2d: np.ndarray,
                     trace_poly: np.poly1d,
                     output_height: Optional[int] = None,
                     center_y: int = 20) -> np.ndarray:
    """
    Rectify a curved 2D spectrum to straighten the trace.

    Parameters
    ----------
    data_2d : np.ndarray
        2D spectral image with curved trace
    trace_poly : np.poly1d
        Polynomial describing trace curvature
    output_height : int, optional
        Height of rectified spectrum (default: same as input)
    center_y : int, optional
        Y pixel to center the straightened trace (default: 20)

    Returns
    -------
    rectified : np.ndarray
        Rectified 2D spectrum with straight trace

    Examples
    --------
    >>> rect_spec = rectify_spectrum(data_2d, trace_poly, center_y=20)
    """
    ny, nx = data_2d.shape

    if output_height is None:
        output_height = ny

    rectified = np.zeros((output_height, nx))

    # Shift each column to straighten the trace
    for ix in range(nx):
        trace_y = int(trace_poly(ix))
        shift = -trace_y + center_y

        # Use np.roll for shifting
        rectified[:, ix] = np.roll(data_2d[:, ix], shift)

    return rectified


def subtract_background(data_2d: np.ndarray,
                       mask_sky: np.ndarray,
                       method: str = 'median',
                       poly_order: int = 1) -> np.ndarray:
    """
    Subtract sky/background from 2D spectrum.

    Parameters
    ----------
    data_2d : np.ndarray
        2D spectral image
    mask_sky : np.ndarray
        Boolean mask where True indicates sky regions
    method : str, optional
        Background estimation method: 'median', 'mean', or 'poly' (default: 'median')
    poly_order : int, optional
        Polynomial order for 'poly' method (default: 1)

    Returns
    -------
    data_sub : np.ndarray
        Background-subtracted 2D spectrum

    Examples
    --------
    >>> data_clean = subtract_background(data_2d, mask_sky, method='median')
    """
    ny, nx = data_2d.shape
    data_sub = data_2d.copy()

    for ix in range(nx):
        sky_pixels = data_2d[mask_sky[:, ix], ix]

        if len(sky_pixels) == 0:
            continue

        if method == 'median':
            background = np.median(sky_pixels)
        elif method == 'mean':
            background = np.mean(sky_pixels)
        elif method == 'poly':
            # Fit polynomial to sky regions
            y_sky = np.where(mask_sky[:, ix])[0]
            if len(y_sky) > poly_order + 1:
                coeffs = np.polyfit(y_sky, sky_pixels, poly_order)
                poly = np.poly1d(coeffs)
                background_profile = poly(np.arange(ny))
                data_sub[:, ix] -= background_profile
                continue
            else:
                background = np.median(sky_pixels)
        else:
            raise ValueError(f"Unknown method: {method}")

        data_sub[:, ix] -= background

    return data_sub


def extract_1d_boxcar(data_2d: np.ndarray,
                      mask_slit: Optional[np.ndarray] = None,
                      width: Optional[int] = None,
                      center_y: Optional[int] = None) -> np.ndarray:
    """
    Extract 1D spectrum using boxcar (simple summation).

    Parameters
    ----------
    data_2d : np.ndarray
        2D spectral image (should be rectified and background-subtracted)
    mask_slit : np.ndarray, optional
        Boolean mask defining extraction region
    width : int, optional
        Extraction width in pixels (used if mask_slit not provided)
    center_y : int, optional
        Center of extraction (used if mask_slit not provided)

    Returns
    -------
    spectrum_1d : np.ndarray
        Extracted 1D spectrum

    Examples
    --------
    >>> spec_1d = extract_1d_boxcar(rectified, mask_slit=mask)
    >>> # Or specify width and center
    >>> spec_1d = extract_1d_boxcar(rectified, width=10, center_y=20)
    """
    ny, nx = data_2d.shape

    if mask_slit is not None:
        # Use provided mask
        spectrum_1d = np.sum(data_2d * mask_slit, axis=0)
    elif width is not None and center_y is not None:
        # Extract based on width and center
        y_start = max(0, center_y - width // 2)
        y_end = min(ny, center_y + width // 2)
        spectrum_1d = np.sum(data_2d[y_start:y_end, :], axis=0)
    else:
        raise ValueError("Must provide either mask_slit or (width, center_y)")

    return spectrum_1d


def extract_1d_optimal(data_2d: np.ndarray,
                      variance: Optional[np.ndarray] = None,
                      mask_slit: Optional[np.ndarray] = None,
                      spatial_profile: Optional[np.ndarray] = None,
                      readnoise: float = 3.0,
                      gain: float = 2.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract 1D spectrum using optimal (variance-weighted) extraction.

    Based on Horne 1986, PASP 98, 609.

    Parameters
    ----------
    data_2d : np.ndarray
        2D spectral image (rectified, background-subtracted)
    variance : np.ndarray, optional
        Variance image (if None, estimated from data)
    mask_slit : np.ndarray, optional
        Boolean mask defining extraction region
    spatial_profile : np.ndarray, optional
        Normalized spatial profile (if None, estimated from data)
    readnoise : float, optional
        CCD readout noise in electrons (default: 3.0)
    gain : float, optional
        CCD gain in e-/ADU (default: 2.0)

    Returns
    -------
    spectrum_1d : np.ndarray
        Optimally extracted 1D spectrum
    variance_1d : np.ndarray
        Variance of extracted spectrum

    Examples
    --------
    >>> spec_1d, var_1d = extract_1d_optimal(rectified, mask_slit=mask)
    >>> snr = spec_1d / np.sqrt(var_1d)

    Notes
    -----
    Optimal extraction weights each pixel by its inverse variance and
    the spatial profile, maximizing S/N for the extracted spectrum.

    References
    ----------
    Horne, K. 1986, PASP, 98, 609
    """
    ny, nx = data_2d.shape

    # Estimate variance if not provided
    if variance is None:
        # Variance = (data * gain + readnoise^2) / gain^2
        variance = np.maximum(data_2d * gain + readnoise**2, 1.0) / gain**2

    # Estimate spatial profile if not provided
    if spatial_profile is None:
        # Median spatial profile from bright regions
        spectrum_1d_temp = np.sum(data_2d, axis=0)
        bright_cols = spectrum_1d_temp > np.percentile(spectrum_1d_temp, 75)

        if mask_slit is not None:
            profile_data = data_2d[:, bright_cols] * mask_slit[:, bright_cols]
        else:
            profile_data = data_2d[:, bright_cols]

        spatial_profile = np.median(profile_data, axis=1)
        spatial_profile /= np.sum(spatial_profile)  # Normalize

        # Smooth the profile
        spatial_profile = ndimage.gaussian_filter1d(spatial_profile, sigma=1.0)

    # Apply mask to profile
    if mask_slit is not None:
        profile_2d = spatial_profile[:, np.newaxis] * mask_slit
    else:
        profile_2d = spatial_profile[:, np.newaxis] * np.ones((ny, nx))

    # Optimal extraction weights
    # w = P / V where P is profile, V is variance
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        weights = profile_2d / variance
        weights = np.nan_to_num(weights, 0.0)

    # Extracted spectrum
    numerator = np.sum(weights * data_2d, axis=0)
    denominator = np.sum(weights * profile_2d, axis=0)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        spectrum_1d = np.divide(numerator, denominator,
                               out=np.zeros_like(numerator),
                               where=denominator!=0)

    # Variance of extracted spectrum
    variance_1d = np.sum(profile_2d, axis=0) / np.sum(weights * profile_2d, axis=0)
    variance_1d = np.nan_to_num(variance_1d, np.inf)

    return spectrum_1d, variance_1d


def full_extraction_pipeline(spec_file: str,
                             config: Optional[Dict] = None) -> Dict[str, np.ndarray]:
    """
    Complete extraction pipeline for a single slit.

    Parameters
    ----------
    spec_file : str
        Path to spec_*.fits file from initial reduction
    config : dict, optional
        Configuration parameters (see examples)

    Returns
    -------
    results : dict
        Dictionary containing:
        - 'spec1d_boxcar': Boxcar extracted spectrum
        - 'spec1d_optimal': Optimally extracted spectrum
        - 'variance_1d': Variance of optimal extraction
        - 'spec2d_rect': Rectified 2D spectrum
        - 'trace_poly': Trace polynomial
        - 'wavelengths': Wavelength array (if calibrated)

    Examples
    --------
    >>> config = {
    ...     'trace_fit': {'method': 'polynomial', 'degree': 3},
    ...     'extraction': {'width': 10},
    ...     'background': {'method': 'median'}
    ... }
    >>> results = full_extraction_pipeline('spec_010.fits', config)
    """
    from astropy.io import fits

    # Default configuration
    if config is None:
        config = {
            'trace_fit': {'method': 'polynomial', 'degree': 3},
            'extraction': {'width': 10},
            'background': {'method': 'median'}
        }

    # Load data
    with fits.open(spec_file) as hdul:
        data = hdul['DATA'].data
        flat = hdul['FLAT'].data
        lines = hdul['LINES'].data  # Arc lamp
        header = hdul[0].header

        slit_edges = (header['Y0'], header['Y1'], header['Y2'], header['Y3'])

    # Step 1: Fit trace from flat
    x_pix, y_ctr = find_trace_edges(flat, x_range=(1500, 2000))
    trace_poly = fit_trace(x_pix, y_ctr, **config['trace_fit'])

    # Step 2: Create masks
    mask_slit, mask_sky = create_spatial_mask(data.shape, trace_poly, slit_edges)

    # Step 3: Rectify
    data_rect = rectify_spectrum(data, trace_poly)

    # Step 4: Background subtraction
    data_clean = subtract_background(data_rect, mask_sky, **config['background'])

    # Step 5: Extract 1D
    spec1d_boxcar = extract_1d_boxcar(data_clean, mask_slit=mask_slit)
    spec1d_optimal, variance_1d = extract_1d_optimal(data_clean, mask_slit=mask_slit)

    results = {
        'spec1d_boxcar': spec1d_boxcar,
        'spec1d_optimal': spec1d_optimal,
        'variance_1d': variance_1d,
        'spec2d_rect': data_clean,
        'trace_poly': trace_poly,
    }

    return results
