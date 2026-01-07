#!/usr/bin/env python3
"""
SAMOS Spectroscopy Pipeline - Automated Reduction Script

This script runs the complete SAMOS spectroscopic data reduction pipeline
in automated mode using parameters from a configuration file.

Usage:
    python run_spectroscopy_pipeline.py --config config.yaml [--output /path/to/output]

The pipeline performs:
    1. Read and combine calibration frames (bias, flat, arc)
    2. Detect spectral traces from arc lamp
    3. Extract individual slit cutouts
    4. Apply calibrations (bias, flat, CR removal)
    5. Save reduced 2D spectra and intermediate products
"""

import argparse
import yaml
import sys
import os
from pathlib import Path
import numpy as np
from astropy.io import fits

# Add parent directory to path to import samos package
sys.path.insert(0, str(Path(__file__).parent.parent))

from samos.core import mosaic, cosmic_rays, calibration
from samos.spectroscopy import trace_detection
from samos.utils import io, display


def load_config(config_file: str) -> dict:
    """Load configuration from YAML file."""
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    return config


def setup_output_directory(config: dict) -> Path:
    """Create output directory structure."""
    output_base = Path(config['output']['base_directory'])

    subdirs = [
        'calibrations',
        'traces',
        'extracted_slits',
        'final_products',
        'qa_plots'
    ]

    for subdir in subdirs:
        (output_base / subdir).mkdir(parents=True, exist_ok=True)

    print(f"\nOutput directory: {output_base}")
    return output_base


def process_calibration_frames(config: dict, output_dir: Path) -> dict:
    """
    Process bias and flat calibration frames.

    Returns dictionary with master calibration frames.
    """
    print("\n" + "="*60)
    print("STEP 1: Processing Calibration Frames")
    print("="*60)

    calib_config = config['calibrations']
    data_dir = Path(config['data']['raw_directory'])

    calibrations = {}

    # Process bias frames
    if calib_config.get('use_bias', True):
        print("\nCreating master bias...")
        bias_files = io.find_files(str(data_dir), calib_config['bias_pattern'])
        print(f"Found {len(bias_files)} bias frames")

        bias_mosaics = []
        for bias_file in bias_files:
            hdu = mosaic.read_sami_mosaic(bias_file)
            bias_mosaics.append(hdu.data)

        master_bias = calibration.create_master_bias(bias_mosaics, method='median')
        calibrations['bias'] = master_bias

        # Save master bias
        io.write_fits(master_bias, str(output_dir / 'calibrations' / 'master_bias.fits'))
        print("✓ Master bias created and saved")
    else:
        calibrations['bias'] = None

    # Process flat frames
    if calib_config.get('use_flat', True):
        print("\nCreating master flat...")
        flat_files = io.find_files(str(data_dir), calib_config['flat_pattern'])
        print(f"Found {len(flat_files)} flat frames")

        flat_mosaics = []
        for flat_file in flat_files:
            hdu = mosaic.read_sami_mosaic(flat_file)
            flat_mosaics.append(hdu.data)

        master_flat = calibration.create_master_flat(
            flat_mosaics,
            master_bias=calibrations['bias'],
            method='median',
            normalize=True
        )
        calibrations['flat'] = master_flat

        # Save master flat
        io.write_fits(master_flat, str(output_dir / 'calibrations' / 'master_flat.fits'))
        print("✓ Master flat created and saved")
    else:
        calibrations['flat'] = None

    return calibrations


def detect_traces(config: dict, calibrations: dict, output_dir: Path) -> list:
    """
    Detect spectral traces using arc lamp exposure.

    Returns list of slit definitions.
    """
    print("\n" + "="*60)
    print("STEP 2: Detecting Spectral Traces")
    print("="*60)

    arc_config = config['arc_lamp']
    data_dir = Path(config['data']['raw_directory'])

    # Read arc lamp frame
    arc_file = data_dir / arc_config['arc_file']
    print(f"\nReading arc lamp: {arc_file.name}")
    arc_hdu = mosaic.read_sami_mosaic(str(arc_file))
    arc_data = arc_hdu.data

    # Read arc dark if specified
    arc_dark = None
    if arc_config.get('arc_dark_file'):
        arc_dark_file = data_dir / arc_config['arc_dark_file']
        print(f"Reading arc dark: {arc_dark_file.name}")
        arc_dark_hdu = mosaic.read_sami_mosaic(str(arc_dark_file))
        arc_dark = arc_dark_hdu.data

    # Remove cosmic rays from arc
    print("Removing cosmic rays from arc...")
    arc_clean = cosmic_rays.remove_cosmic_rays(arc_data, verbose=True)

    # Subtract dark if present
    if arc_dark is not None:
        arc_dark_clean = cosmic_rays.remove_cosmic_rays(arc_dark, verbose=False)
        arc_final = arc_clean - arc_dark_clean
        print("Arc dark subtracted")
    else:
        arc_final = arc_clean

    # Save cleaned arc
    io.write_fits(arc_final, str(output_dir / 'calibrations' / 'arc_clean.fits'))

    # Detect traces
    trace_config = config['trace_detection']
    slits = trace_detection.detect_traces_from_arc(
        arc_final,
        x_range=tuple(trace_config['x_range']),
        peak_threshold=trace_config['peak_threshold'],
        margin=trace_config['margin'],
        verbose=True
    )

    print(f"\n✓ Detected {len(slits)} spectral traces")

    # Save slit definitions
    slit_info = {
        'n_slits': len(slits),
        'slits': slits
    }
    with open(output_dir / 'traces' / 'slit_definitions.yaml', 'w') as f:
        yaml.dump(slit_info, f)

    return slits


def reduce_science_frames(config: dict, calibrations: dict, slits: list,
                         output_dir: Path) -> None:
    """
    Reduce science exposures and extract individual slits.
    """
    print("\n" + "="*60)
    print("STEP 3: Reducing Science Frames")
    print("="*60)

    data_dir = Path(config['data']['raw_directory'])
    science_config = config['science']

    # Find science files
    science_files = io.find_files(str(data_dir), science_config['science_pattern'])
    print(f"\nFound {len(science_files)} science frames")

    for science_file in science_files:
        print(f"\nProcessing: {Path(science_file).name}")

        # Read and create mosaic
        science_hdu = mosaic.read_sami_mosaic(science_file)
        science_data = science_hdu.data

        # Remove cosmic rays
        science_clean = cosmic_rays.remove_cosmic_rays(science_data, verbose=False)
        print("  ✓ Cosmic rays removed")

        # Apply calibrations
        science_reduced = calibration.reduce_image(
            science_clean,
            master_bias=calibrations.get('bias'),
            master_flat=calibrations.get('flat')
        )
        print("  ✓ Calibrations applied")

        # Extract and save individual slits
        extract_and_save_slits(science_reduced, calibrations, slits,
                              science_file, output_dir)


def extract_and_save_slits(science_data: np.ndarray, calibrations: dict,
                          slits: list, science_file: str, output_dir: Path) -> None:
    """
    Extract individual slits and save as multi-extension FITS files.
    """
    science_name = Path(science_file).stem

    for i_slit, slit in enumerate(slits):
        # Extract slit cutout
        slit_data = science_data[slit[0]:slit[3], :]

        # Extract flat cutout
        if calibrations.get('flat') is not None:
            flat_cutout = calibrations['flat'][slit[0]:slit[3], :]
        else:
            flat_cutout = np.ones_like(slit_data)

        # Create multi-extension FITS
        hdr = fits.Header()
        hdr['SLIT'] = i_slit
        hdr['Y0'] = (slit[0], 'y[0] pixel of the 2D cut')
        hdr['Y1'] = (slit[1], 'y[0] pixel of the trace')
        hdr['Y2'] = (slit[2], 'y[1] pixel of the trace')
        hdr['Y3'] = (slit[3], 'y[1] pixel of the 2D cut')
        hdr['SOURCE'] = Path(science_file).name

        primary_hdu = fits.PrimaryHDU(header=hdr)
        data_hdu = fits.ImageHDU(data=slit_data, name='DATA')
        flat_hdu = fits.ImageHDU(data=flat_cutout, name='FLAT')

        hdulist = fits.HDUList([primary_hdu, data_hdu, flat_hdu])

        # Save
        output_file = output_dir / 'extracted_slits' / f'{science_name}_slit_{i_slit:03d}.fits'
        hdulist.writeto(output_file, overwrite=True)

    print(f"  ✓ Extracted {len(slits)} slits")


def main():
    parser = argparse.ArgumentParser(
        description='Run SAMOS spectroscopy reduction pipeline'
    )
    parser.add_argument('--config', required=True,
                       help='Path to configuration YAML file')
    parser.add_argument('--output', default=None,
                       help='Output directory (overrides config file)')

    args = parser.parse_args()

    # Load configuration
    print("\n" + "="*60)
    print("SAMOS Spectroscopy Pipeline")
    print("="*60)
    print(f"\nLoading configuration: {args.config}")
    config = load_config(args.config)

    # Override output directory if specified
    if args.output:
        config['output']['base_directory'] = args.output

    # Setup output
    output_dir = setup_output_directory(config)

    # Run pipeline
    try:
        # Step 1: Process calibrations
        calibrations = process_calibration_frames(config, output_dir)

        # Step 2: Detect traces
        slits = detect_traces(config, calibrations, output_dir)

        # Step 3: Reduce science frames
        reduce_science_frames(config, calibrations, slits, output_dir)

        print("\n" + "="*60)
        print("Pipeline completed successfully!")
        print("="*60)
        print(f"\nResults saved to: {output_dir}")

    except Exception as e:
        print(f"\n ERROR: Pipeline failed with exception:")
        print(f"{type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == '__main__':
    sys.exit(main())
