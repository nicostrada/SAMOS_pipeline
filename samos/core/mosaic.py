"""
SAMOS Mosaic Creation Module

This module handles the reading and assembly of multi-CCD mosaic images
from SAMOS multi-extension FITS files.
"""

import numpy as np
from astropy.io import fits
from typing import Tuple


def read_sami_mosaic(filepath: str, fix_bad_columns: bool = True) -> fits.PrimaryHDU:
    """
    Read and assemble SAMOS 4-CCD mosaic from multi-extension FITS file.

    SAMOS uses 4 CCDs arranged in a 2x2 pattern. This function reads the
    individual CCD data, trims overscan regions, flips arrays to correct
    orientation, and assembles them into a single mosaic image.

    Parameters
    ----------
    filepath : str
        Path to the SAMOS multi-extension FITS file
    fix_bad_columns : bool, optional
        If True, interpolate over known bad columns (default: True)

    Returns
    -------
    astropy.io.fits.PrimaryHDU
        HDU containing the assembled mosaic data with original header

    Notes
    -----
    CCD Layout (as read from file):
        Extension 1: Top left CCD (flipped)
        Extension 2: Top right CCD (flipped)
        Extension 3: Bottom left CCD (flipped)
        Extension 4: Bottom right CCD (flipped)

    Trimming:
        - Extensions 1, 3: columns 54:2101 (2047 pixels)
        - Extensions 2, 4: columns 64:2111 (2047 pixels)

    Examples
    --------
    >>> from samos.core import mosaic
    >>> hdu = mosaic.read_sami_mosaic('science_001.fits')
    >>> mosaic_data = hdu.data
    >>> print(mosaic_data.shape)
    (4096, 4094)
    """
    with fits.open(filepath) as hdul:
        hdr = hdul[0].header

        # Read and trim each CCD
        # Top left CCD (extension 1)
        data1 = hdul[1].data[:, 54:2101]

        # Top right CCD (extension 2)
        data2 = hdul[2].data[:, 64:2111]

        # Bottom left CCD (extension 3)
        data3 = hdul[3].data[:, 54:2101]

        # Bottom right CCD (extension 4)
        data4 = hdul[4].data[:, 64:2111]

    # Create mosaic array
    mosaic_data = assemble_mosaic(data1, data2, data3, data4)

    # Fix known bad columns
    if fix_bad_columns:
        mosaic_data = fix_bad_column(mosaic_data, column=2046)

    # Create HDU with mosaic
    hdu = fits.PrimaryHDU(data=mosaic_data, header=hdr)

    return hdu


def assemble_mosaic(data1: np.ndarray, data2: np.ndarray,
                    data3: np.ndarray, data4: np.ndarray) -> np.ndarray:
    """
    Assemble 4 CCD arrays into a single mosaic.

    Parameters
    ----------
    data1 : np.ndarray
        Top left CCD data
    data2 : np.ndarray
        Top right CCD data
    data3 : np.ndarray
        Bottom left CCD data
    data4 : np.ndarray
        Bottom right CCD data

    Returns
    -------
    np.ndarray
        Assembled mosaic array with shape (2*height, width1+width2)
    """
    height = data1.shape[0]
    width1 = data1.shape[1]
    width2 = data2.shape[1]

    # Create blank mosaic
    mosaic = np.zeros((height * 2, width1 + width2))

    # Place CCDs with proper orientation (flip along axis 0)
    mosaic[height:, 0:width1] = np.flip(data1, axis=0)          # Top left
    mosaic[height:, width1:] = np.flip(data2, axis=0)           # Top right
    mosaic[0:height, 0:width1] = np.flip(data3, axis=0)         # Bottom left
    mosaic[0:height, width1:] = np.flip(data4, axis=0)          # Bottom right

    # Final flip to get correct orientation
    mosaic = np.flip(mosaic, axis=0)

    return mosaic


def fix_bad_column(data: np.ndarray, column: int) -> np.ndarray:
    """
    Interpolate over a bad column using adjacent columns.

    Parameters
    ----------
    data : np.ndarray
        2D image array
    column : int
        Column index to fix

    Returns
    -------
    np.ndarray
        Image with interpolated column
    """
    if column < 1 or column >= data.shape[1] - 1:
        raise ValueError(f"Column {column} cannot be interpolated (edge of array)")

    data[:, column] = (data[:, column - 1] + data[:, column + 1]) / 2

    return data


def get_mosaic_dimensions() -> Tuple[int, int]:
    """
    Return the expected dimensions of a SAMOS mosaic.

    Returns
    -------
    tuple of int
        (height, width) of assembled mosaic
    """
    # Each CCD: 2048 rows, ~2047 columns after trimming
    # Total: 4096 rows, 4094 columns
    return (4096, 4094)
