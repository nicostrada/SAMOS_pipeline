"""
Input/Output Utilities

Functions for reading and writing FITS files, managing file paths,
and handling data organization.
"""

import os
from pathlib import Path
from typing import List, Optional, Union
from astropy.io import fits
import numpy as np


def read_fits(filepath: str, ext: int = 0) -> np.ndarray:
    """
    Read data from a FITS file.

    Parameters
    ----------
    filepath : str
        Path to FITS file
    ext : int, optional
        Extension number to read (default: 0)

    Returns
    -------
    np.ndarray
        Image data array

    Examples
    --------
    >>> data = read_fits('science_001.fits')
    >>> data = read_fits('multi_ext.fits', ext=1)
    """
    with fits.open(filepath) as hdul:
        data = hdul[ext].data

    return data


def write_fits(data: np.ndarray,
              filepath: str,
              header: Optional[fits.Header] = None,
              overwrite: bool = True) -> None:
    """
    Write data to a FITS file.

    Parameters
    ----------
    data : np.ndarray
        Image data to write
    filepath : str
        Output file path
    header : astropy.io.fits.Header, optional
        FITS header to include
    overwrite : bool, optional
        Overwrite existing file (default: True)

    Examples
    --------
    >>> write_fits(reduced_data, 'reduced_001.fits')
    >>> write_fits(spectrum, 'spec1d.fits', header=hdr)
    """
    hdu = fits.PrimaryHDU(data=data, header=header)
    hdu.writeto(filepath, overwrite=overwrite)
    print(f"Wrote: {filepath}")


def find_files(directory: str,
              pattern: str = "*.fits",
              recursive: bool = False) -> List[str]:
    """
    Find files matching a pattern in a directory.

    Parameters
    ----------
    directory : str
        Directory to search
    pattern : str, optional
        File pattern (supports wildcards) (default: "*.fits")
    recursive : bool, optional
        Search recursively in subdirectories (default: False)

    Returns
    -------
    list of str
        Sorted list of matching file paths

    Examples
    --------
    >>> files = find_files('/data/raw', 'bias*.fits')
    >>> files = find_files('/data', 'sci_*.fits', recursive=True)
    """
    path = Path(directory)

    if recursive:
        files = list(path.rglob(pattern))
    else:
        files = list(path.glob(pattern))

    # Return sorted absolute paths as strings
    return sorted([str(f.resolve()) for f in files])


def organize_files_by_type(directory: str,
                           file_types: dict = None) -> dict:
    """
    Organize FITS files by observation type from headers.

    Parameters
    ----------
    directory : str
        Directory containing FITS files
    file_types : dict, optional
        Dictionary mapping OBSTYPE header values to categories
        Default: {'bias': 'BIAS', 'flat': 'FLAT', 'arc': 'ARC',
                 'science': 'OBJECT'}

    Returns
    -------
    dict
        Dictionary with keys as categories and values as file lists

    Examples
    --------
    >>> files = organize_files_by_type('/data/raw')
    >>> print(files['bias'])
    ['/data/raw/bias_001.fits', '/data/raw/bias_002.fits']
    """
    if file_types is None:
        file_types = {
            'bias': 'BIAS',
            'flat': 'FLAT',
            'arc': 'ARC',
            'science': 'OBJECT'
        }

    # Initialize result dictionary
    organized = {key: [] for key in file_types.keys()}
    organized['unknown'] = []

    # Find all FITS files
    fits_files = find_files(directory, "*.fits")

    for filepath in fits_files:
        try:
            with fits.open(filepath) as hdul:
                header = hdul[0].header
                obstype = header.get('OBSTYPE', '').upper()

                # Find matching category
                matched = False
                for category, obstype_value in file_types.items():
                    if obstype == obstype_value:
                        organized[category].append(filepath)
                        matched = True
                        break

                if not matched:
                    organized['unknown'].append(filepath)

        except Exception as e:
            print(f"Warning: Could not read {filepath}: {e}")
            organized['unknown'].append(filepath)

    return organized


def create_output_directory(base_path: str,
                           subdirs: Optional[List[str]] = None) -> Path:
    """
    Create output directory structure.

    Parameters
    ----------
    base_path : str
        Base output directory path
    subdirs : list of str, optional
        Subdirectories to create within base_path

    Returns
    -------
    pathlib.Path
        Path object for base directory

    Examples
    --------
    >>> outdir = create_output_directory('/data/reduced',
    ...                                  ['calibrated', 'extracted', 'final'])
    """
    base = Path(base_path)
    base.mkdir(parents=True, exist_ok=True)

    if subdirs:
        for subdir in subdirs:
            (base / subdir).mkdir(exist_ok=True)

    return base


def get_header_value(filepath: str,
                    keyword: str,
                    ext: int = 0,
                    default: Optional[any] = None) -> any:
    """
    Get a header keyword value from a FITS file.

    Parameters
    ----------
    filepath : str
        Path to FITS file
    keyword : str
        Header keyword to retrieve
    ext : int, optional
        Extension number (default: 0)
    default : any, optional
        Default value if keyword not found

    Returns
    -------
    any
        Header keyword value or default

    Examples
    --------
    >>> exptime = get_header_value('science.fits', 'EXPTIME')
    >>> target = get_header_value('science.fits', 'OBJECT', default='Unknown')
    """
    try:
        with fits.open(filepath) as hdul:
            value = hdul[ext].header.get(keyword, default)
        return value
    except Exception as e:
        print(f"Warning: Could not read {keyword} from {filepath}: {e}")
        return default


def combine_images(filepaths: List[str],
                  method: str = 'median',
                  output_file: Optional[str] = None) -> np.ndarray:
    """
    Combine multiple FITS images using median or mean.

    Parameters
    ----------
    filepaths : list of str
        List of FITS file paths to combine
    method : str, optional
        Combination method: 'median' or 'mean' (default: 'median')
    output_file : str, optional
        If provided, save combined image to this file

    Returns
    -------
    np.ndarray
        Combined image array

    Examples
    --------
    >>> bias_files = find_files('/data/raw', 'bias*.fits')
    >>> master_bias = combine_images(bias_files, method='median',
    ...                              output_file='master_bias.fits')
    """
    if not filepaths:
        raise ValueError("No files provided for combination")

    # Read all images
    images = []
    for filepath in filepaths:
        data = read_fits(filepath)
        images.append(data)

    # Stack into 3D array
    image_stack = np.array(images)

    # Combine
    if method == 'median':
        combined = np.median(image_stack, axis=0)
    elif method == 'mean':
        combined = np.mean(image_stack, axis=0)
    else:
        raise ValueError(f"Unknown method: {method}. Use 'median' or 'mean'")

    # Save if requested
    if output_file:
        write_fits(combined, output_file)

    return combined
