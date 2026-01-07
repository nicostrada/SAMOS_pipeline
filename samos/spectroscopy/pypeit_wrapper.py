"""
PypeIt Wrapper for SAMOS

This module provides a simplified interface to PypeIt for wavelength
calibration of SAMOS spectroscopic data.
"""

import numpy as np
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from astropy.io import fits
import warnings


def create_samos_spectrograph_file():
    """
    Create a custom PypeIt spectrograph configuration for SAMOS.

    This defines SAMOS as a custom instrument for PypeIt.

    Returns
    -------
    str
        Python code defining SAMOS spectrograph class

    Notes
    -----
    SAMOS specifications:
    - Detector: 4k x 4k CCD mosaic
    - Pixel scale: ~0.18 arcsec/pixel
    - Dispersion: ~1.5 Å/pixel (red arm)
    - Wavelength range: 5500-9000 Å (red arm)
    """
    samos_spec = """
from pypeit import spectrographs
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container

class SAMOSRedSpectrograph(spectrograph.Spectrograph):
    '''
    Child class for SAMOS Red spectrograph.
    '''
    ndet = 1
    name = 'samos_red'
    telescope = spectrograph.TelescopePar(
        name='Soar',
        longitude=-70.73,
        latitude=-30.24,
        elevation=2738.0,
        eff_aperture=4.1
    )
    camera = 'red'
    url = 'http://www.ctio.noirlab.edu/soar/content/soar-adaptive-module-samos'
    header_name = 'SAMOS'
    pypeline = 'MultiSlit'

    def get_detector_par(self, det, hdu=None):
        '''
        Return metadata for SAMOS detector.
        '''
        detector_dict = dict(
            binning='1,1',
            det=1,
            dataext=0,
            specaxis=0,
            specflip=False,
            spatflip=False,
            platescale=0.18,  # arcsec/pixel
            darkcurr=0.0,
            saturation=65535.0,
            nonlinear=0.99,
            mincounts=-1e10,
            numamplifiers=1,
            gain=np.asarray([2.0]),
            ronoise=np.asarray([3.0]),
            datasec=np.asarray(['[:,:]']),
            oscansec=None
        )
        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        '''
        Return default parameters for SAMOS.
        '''
        par = super().default_pypeit_par()

        # Wavelength calibration
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['lamps'] = ['ArI', 'NeI', 'HgI']
        par['calibrations']['wavelengths']['rms_threshold'] = 0.5
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        par['calibrations']['wavelengths']['fwhm'] = 4.0

        # Arc line detection
        par['calibrations']['wavelengths']['n_first'] = 2
        par['calibrations']['wavelengths']['n_final'] = 4
        par['calibrations']['wavelengths']['sigrej_first'] = 2.0
        par['calibrations']['wavelengths']['sigrej_final'] = 3.0

        # Spectral range
        par['calibrations']['wavelengths']['spec_min_max'] = [5500.0, 9000.0]

        return par

    def config_specific_par(self, scifile, inp_par=None):
        '''
        Modify parameters for specific configuration.
        '''
        par = super().config_specific_par(scifile, inp_par=inp_par)
        return par

    def init_meta(self):
        '''
        Define metadata for SAMOS files.
        '''
        self.meta = {}
        self.meta['ra'] = dict(ext=0, card='RA')
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['decker'] = dict(ext=0, card='SLITNAME')
        self.meta['binning'] = dict(ext=0, card=None, default='1,1')
        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        self.meta['dispname'] = dict(ext=0, card='GRATING')

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        '''
        Check for frames of a given type.
        '''
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype == 'science':
            return good_exp & (fitstbl['exptime'] > 60)
        if ftype == 'bias':
            return good_exp & (fitstbl['exptime'] < 1)
        if ftype == 'pixelflat' or ftype == 'trace':
            return good_exp & (fitstbl['exptime'] > 1) & (fitstbl['exptime'] < 30)
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['exptime'] > 1) & (fitstbl['exptime'] < 30)
        return np.zeros(len(fitstbl), dtype=bool)
"""
    return samos_spec


def configure_pypeit_for_samos(output_dir: str = '.') -> Path:
    """
    Configure PypeIt for SAMOS data reduction.

    Parameters
    ----------
    output_dir : str
        Directory to write configuration files

    Returns
    -------
    config_file : Path
        Path to created configuration file

    Examples
    --------
    >>> config = configure_pypeit_for_samos(output_dir='./pypeit')
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Write spectrograph file
    spec_file = output_dir / 'samos_spectrograph.py'
    spec_file.write_text(create_samos_spectrograph_file())

    print(f"✓ Created SAMOS spectrograph configuration: {spec_file}")

    return spec_file


def prepare_arc_for_pypeit(arc_fits: str,
                           output_file: str,
                           slit_index: int = 10) -> str:
    """
    Prepare SAMOS arc lamp data for PypeIt wavelength calibration.

    Parameters
    ----------
    arc_fits : str
        Path to spec_*.fits file containing arc lamp in LINES extension
    output_file : str
        Output FITS file for PypeIt
    slit_index : int
        Slit number to extract (default: 10, good reference slit)

    Returns
    -------
    output_file : str
        Path to created FITS file

    Examples
    --------
    >>> arc_file = prepare_arc_for_pypeit('spec_010.fits', 'pypeit_arc.fits')
    """
    # Load arc data
    with fits.open(arc_fits) as hdul:
        arc_data = hdul['LINES'].data
        header = hdul[0].header.copy()

    # Create simple FITS for PypeIt
    hdu = fits.PrimaryHDU(arc_data, header=header)
    hdu.header['OBJECT'] = 'ARC'
    hdu.header['IMAGETYP'] = 'ARC'

    hdu.writeto(output_file, overwrite=True)

    print(f"✓ Created PypeIt arc file: {output_file}")
    return output_file


def run_wavelength_calibration_simple(arc_1d: np.ndarray,
                                      line_list: str = 'HgArNe',
                                      wv_range: Tuple[float, float] = (5500, 9000),
                                      ) -> Dict:
    """
    Simple wavelength calibration using scipy (fallback if PypeIt too complex).

    This provides a simpler alternative using manual line matching.

    Parameters
    ----------
    arc_1d : np.ndarray
        1D arc spectrum
    line_list : str
        Lamp type: 'HgArNe', 'Ne', 'Ar', etc.
    wv_range : tuple
        Expected wavelength range (angstroms)

    Returns
    -------
    solution : dict
        Wavelength solution with keys:
        - 'coeffs': Polynomial coefficients
        - 'rms': RMS residual
        - 'matched_lines': Matched pixel/wavelength pairs

    Examples
    --------
    >>> solution = run_wavelength_calibration_simple(arc_1d)
    >>> wavelengths = np.polyval(solution['coeffs'], pixels)
    """
    from scipy.signal import find_peaks

    # Reference wavelengths for HgArNe
    reference_lines = {
        'HgArNe': np.array([
            5460.7500, 5769.6100, 5790.6700, 5852.4879, 5944.8342,
            6074.3377, 6096.1631, 6143.5939, 6163.5939, 6217.2812,
            6266.4950, 6334.4278, 6382.9917, 6402.2480, 6506.5281,
            6532.8822, 6598.9529, 6678.2762, 6717.0430, 6929.4673,
            6965.4310, 7032.4131, 7173.9381, 7245.1666, 7635.1060,
            7724.2070, 7948.1760, 8115.3110, 8264.5220, 8377.6080,
            8424.6480
        ]),
        'Ne': np.array([
            5852.4879, 5944.8342, 6074.3377, 6096.1631, 6143.5939,
            6163.5939, 6217.2812, 6266.4950, 6334.4278, 6382.9917,
            6402.2480, 6506.5281, 6532.8822, 6598.9529, 6678.2762,
            6717.0430, 6929.4673, 6965.4310, 7032.4131, 7173.9381,
            7245.1666
        ]),
        'Ar': np.array([
            7635.1060, 7724.2070, 7948.1760, 8115.3110, 8264.5220,
            8377.6080, 8424.6480
        ])
    }

    ref_wavelengths = reference_lines.get(line_list, reference_lines['HgArNe'])

    # Filter to expected range
    ref_wavelengths = ref_wavelengths[
        (ref_wavelengths >= wv_range[0]) & (ref_wavelengths <= wv_range[1])
    ]

    # Detect peaks in arc
    peaks, properties = find_peaks(arc_1d, height=100, prominence=50)

    if len(peaks) < 10:
        raise ValueError(f"Only {len(peaks)} peaks found. Need at least 10 for calibration.")

    # Simple initial guess: linear mapping
    pixel_range = (0, len(arc_1d))
    initial_dispersion = (wv_range[1] - wv_range[0]) / len(arc_1d)

    # Match peaks to reference lines (simple nearest neighbor)
    matched_pixels = []
    matched_wavelengths = []

    for ref_wv in ref_wavelengths:
        # Expected pixel position
        expected_pixel = (ref_wv - wv_range[0]) / initial_dispersion

        # Find nearest peak
        distances = np.abs(peaks - expected_pixel)
        if np.min(distances) < 50:  # Within 50 pixels
            nearest_peak = peaks[np.argmin(distances)]
            matched_pixels.append(nearest_peak)
            matched_wavelengths.append(ref_wv)

    if len(matched_pixels) < 10:
        raise ValueError(f"Only matched {len(matched_pixels)} lines. Need at least 10.")

    matched_pixels = np.array(matched_pixels)
    matched_wavelengths = np.array(matched_wavelengths)

    # Fit polynomial
    degree = 3
    coeffs = np.polyfit(matched_pixels, matched_wavelengths, degree)

    # Calculate RMS
    fitted_wavelengths = np.polyval(coeffs, matched_pixels)
    residuals = fitted_wavelengths - matched_wavelengths
    rms = np.sqrt(np.mean(residuals**2))

    solution = {
        'coeffs': coeffs,
        'rms': rms,
        'matched_pixels': matched_pixels,
        'matched_wavelengths': matched_wavelengths,
        'n_lines': len(matched_pixels),
        'residuals': residuals
    }

    print(f"✓ Wavelength calibration: {len(matched_pixels)} lines, RMS = {rms:.3f} Å")

    return solution


def apply_wavelength_solution(pixel_array: np.ndarray,
                              solution: Dict) -> np.ndarray:
    """
    Apply wavelength solution to pixel array.

    Parameters
    ----------
    pixel_array : np.ndarray
        Array of pixel positions
    solution : dict
        Wavelength solution from run_wavelength_calibration_simple

    Returns
    -------
    wavelengths : np.ndarray
        Wavelength array

    Examples
    --------
    >>> pixels = np.arange(4094)
    >>> wavelengths = apply_wavelength_solution(pixels, solution)
    """
    coeffs = solution['coeffs']
    wavelengths = np.polyval(coeffs, pixel_array)
    return wavelengths


def create_wavelength_fits(spec_file: str,
                          solution: Dict,
                          output_1d: str,
                          output_2d: str,
                          extraction_results: Dict):
    """
    Create wavelength-calibrated FITS files.

    Parameters
    ----------
    spec_file : str
        Original spec_*.fits file
    solution : dict
        Wavelength solution
    output_1d : str
        Output file for 1D spectrum
    output_2d : str
        Output file for 2D spectrum
    extraction_results : dict
        Results from extraction.full_extraction_pipeline()

    Returns
    -------
    tuple
        (output_1d, output_2d) paths

    Examples
    --------
    >>> create_wavelength_fits('spec_010.fits', solution,
    ...                        'spec1d_010.fits', 'spec2d_010.fits',
    ...                        extraction_results)
    """
    from astropy import units as u
    from specutils import Spectrum1D

    # Get wavelength array
    nx = len(extraction_results['spec1d_optimal'])
    pixels = np.arange(nx)
    wavelengths = apply_wavelength_solution(pixels, solution)

    # Trim to valid wavelength range
    valid = (wavelengths >= 5800) & (wavelengths <= 9500)

    # Create 1D spectrum
    spec1d = Spectrum1D(
        spectral_axis=wavelengths[valid] * u.AA,
        flux=extraction_results['spec1d_optimal'][valid] * u.adu
    )

    # Create 2D spectrum
    spec2d = Spectrum1D(
        spectral_axis=wavelengths[valid] * u.AA,
        flux=extraction_results['spec2d_rect'][:, valid] * u.adu
    )

    # Write FITS
    spec1d.write(output_1d, overwrite=True)
    spec2d.write(output_2d, overwrite=True)

    print(f"✓ Created wavelength-calibrated spectra:")
    print(f"  1D: {output_1d}")
    print(f"  2D: {output_2d}")

    return output_1d, output_2d


def batch_wavelength_calibration(spec_files: List[str],
                                 reference_slit: int = 10,
                                 output_dir: str = '.') -> Dict:
    """
    Apply wavelength calibration to all slits.

    Parameters
    ----------
    spec_files : list
        List of spec_*.fits files
    reference_slit : int
        Slit to use for wavelength solution (default: 10)
    output_dir : str
        Output directory

    Returns
    -------
    results : dict
        Dictionary with calibration results for each slit

    Examples
    --------
    >>> from pathlib import Path
    >>> spec_files = sorted(Path('.').glob('spec_*.fits'))
    >>> results = batch_wavelength_calibration(spec_files)
    """
    from . import extraction

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get wavelength solution from reference slit
    ref_file = f'spec_{reference_slit:03d}.fits'

    print(f"Calibrating reference slit {reference_slit}...")

    # Extract and get 1D arc
    ref_results = extraction.full_extraction_pipeline(ref_file)

    # Load arc
    with fits.open(ref_file) as hdul:
        arc_rect = extraction.rectify_spectrum(
            hdul['LINES'].data,
            ref_results['trace_poly']
        )
    arc_1d = np.sum(arc_rect, axis=0)

    # Get wavelength solution
    solution = run_wavelength_calibration_simple(arc_1d)

    print(f"✓ Reference wavelength solution: RMS = {solution['rms']:.3f} Å")

    # Apply to all slits
    results = {}

    for i, spec_file in enumerate(spec_files):
        slit_num = int(Path(spec_file).stem.split('_')[1])

        print(f"Processing slit {slit_num} ({i+1}/{len(spec_files)})...")

        # Extract
        extract_results = extraction.full_extraction_pipeline(str(spec_file))

        # Create wavelength-calibrated FITS
        out_1d = output_dir / f'spec1d_{slit_num:03d}.fits'
        out_2d = output_dir / f'spec2d_{slit_num:03d}.fits'

        create_wavelength_fits(
            str(spec_file), solution, str(out_1d), str(out_2d), extract_results
        )

        results[slit_num] = {
            'spec1d': str(out_1d),
            'spec2d': str(out_2d),
            'extraction': extract_results
        }

    print(f"\n✓ Batch calibration complete: {len(results)} slits processed")

    return results
