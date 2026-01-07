# SAMOS Data Reduction Pipeline

A modular Python pipeline for reducing multi-object spectroscopy and imaging data from the SAMOS instrument.

## Features

- **Dual-mode operation**: Interactive Jupyter notebooks for parameter tuning + automated scripts for production
- **Multi-instrument support**: Separate modules for spectroscopy and imaging (imaging in development)
- **Modular architecture**: Easy to extend and integrate with other tools (e.g., PypeIt)
- **Professional structure**: Standard Python package with proper documentation and testing

## Installation

### Quick Start

```bash
# Clone or navigate to the Pipeline directory
cd /path/to/SAMOS/Pipeline

# Create a dedicated conda environment (recommended)
conda create -n samos python=3.13
conda activate samos

# Install the pipeline in development mode
pip install -e .
```

### Full Installation with Dependencies

```bash
# Install from requirements file
pip install -r requirements.txt

# Or install with optional dependencies
pip install -e .[dev]      # Development tools
pip install -e .[pypeit]   # PypeIt integration
```

## Quick Start Guide

### Interactive Mode (Recommended for First-Time Users)

1. **Explore your data** using Jupyter notebooks:
   ```bash
   cd notebooks/spectroscopy
   jupyter lab
   ```

2. **Run notebooks in order**:
   - `01_initial_inspection.ipynb` - Visualize raw data
   - `02_calibration_frames.ipynb` - Create master bias/flat
   - `03_trace_identification.ipynb` - Detect spectral traces
   - `04_trace_extraction.ipynb` - Extract 2D spectra
   - `05_wavelength_calibration.ipynb` - Wavelength solution
   - `06_apply_calibration.ipynb` - Apply to all slits
   - `07_visualization.ipynb` - Inspect final products

3. **Save your parameters** to a configuration file for automated runs

### Automated Mode (Production Use)

1. **Copy and edit configuration template**:
   ```bash
   cp configs/spectroscopy_default.yaml my_observation.yaml
   # Edit my_observation.yaml with your file paths and parameters
   ```

2. **Run the automated pipeline**:
   ```bash
   python scripts/run_spectroscopy_pipeline.py \
       --config my_observation.yaml \
       --output /path/to/output
   ```

3. **Results** will be saved in:
   - `calibrations/` - Master bias, flat, arc frames
   - `extracted_slits/` - Individual 2D slit cutouts
   - `final_products/` - 1D spectra (future)
   - `qa_plots/` - Quality assessment plots

## Project Structure

```
Pipeline/
├── samos/                      # Main Python package
│   ├── core/                   # Core processing (mosaic, CR removal, calibration)
│   ├── spectroscopy/           # Spectroscopy tools
│   ├── imaging/                # Imaging tools (future)
│   ├── utils/                  # Utilities (I/O, display, config)
│   └── pipeline/               # Automated pipeline orchestration
│
├── notebooks/                  # Interactive Jupyter notebooks
│   ├── spectroscopy/          # Spectroscopy workflow (01-07)
│   └── imaging/               # Imaging workflow (future)
│
├── scripts/                    # Command-line executables
│   └── run_spectroscopy_pipeline.py
│
├── configs/                    # Configuration templates
│   └── spectroscopy_default.yaml
│
├── calibration_data/           # Reference calibration files
│   └── line_lists/            # Wavelength calibration line lists
│
├── docs/                       # Documentation
│   ├── user_guide/
│   └── developer_guide/
│
└── examples/                   # Example datasets and workflows
```

## Usage Examples

### Example 1: Process a Single Observation

```python
from samos.core import mosaic, cosmic_rays, calibration
from samos.spectroscopy import trace_detection
from samos.utils import display

# Read and display a raw frame
hdu = mosaic.read_sami_mosaic('science_001.fits')
display.display_image(hdu.data, zmin=500, zmax=3000)

# Remove cosmic rays
cleaned = cosmic_rays.remove_cosmic_rays(hdu.data)

# Apply calibrations
reduced = calibration.reduce_image(cleaned, master_bias, master_flat)
```

### Example 2: Detect Traces from Arc Lamp

```python
from samos.spectroscopy import trace_detection

# Detect slits automatically
slits = trace_detection.detect_traces_from_arc(
    arc_image,
    arc_dark=arc_dark_image,
    peak_threshold=4e5,
    margin=12
)

print(f"Found {len(slits)} spectral traces")
```

### Example 3: Batch Processing

```python
from samos.utils import io

# Find all science frames
science_files = io.find_files('/data/raw', 'target*.fits')

# Process each one
for science_file in science_files:
    # ... reduction steps ...
    pass
```

## Data Organization

Keep your data organized outside the Pipeline directory:

```
~/SAMOS_Data/
├── raw/                        # Raw observations
│   └── 2025-01-15_ABELL3120/
│       ├── bias/
│       ├── flat/
│       ├── arc/
│       └── science/
│
├── reduced/                    # Pipeline outputs
│   └── 2025-01-15_ABELL3120/
│       ├── calibrations/
│       ├── extracted_slits/
│       └── final_products/
│
└── configs/                    # Your configurations
    └── abell3120_config.yaml
```

## Configuration

The pipeline uses YAML configuration files. Key sections:

```yaml
data:
  raw_directory: "/path/to/raw/data"
  target_name: "ABELL3120"

calibrations:
  bias_pattern: "bias*.fits"
  flat_pattern: "flat*.fits"

arc_lamp:
  arc_file: "calibration.032.fits"
  arc_dark_file: "calibration.035.fits"

trace_detection:
  x_range: [1500, 2000]
  peak_threshold: 400000.0
  margin: 12
```

See [configs/spectroscopy_default.yaml](configs/spectroscopy_default.yaml) for all options.

## Future Integration: PypeIt

The pipeline is designed to integrate with [PypeIt](https://pypeit.readthedocs.io/) for advanced spectroscopic reduction:

```python
# Future: Use PypeIt for wavelength calibration
from samos.spectroscopy import pypeit_interface

wavelength_solution = pypeit_interface.calibrate_wavelength(
    arc_spectrum,
    use_pypeit=True  # Enable PypeIt
)
```

To install PypeIt support:
```bash
pip install -e .[pypeit]
```

## Development

### Running Tests

```bash
# Install development dependencies
pip install -e .[dev]

# Run tests
pytest tests/

# With coverage
pytest --cov=samos tests/
```

### Contributing

1. Create a feature branch
2. Make your changes
3. Add tests
4. Submit a pull request

### Code Style

We use:
- `black` for code formatting
- `flake8` for linting
- NumPy-style docstrings

```bash
# Format code
black samos/

# Check style
flake8 samos/
```

## Documentation

Full documentation is available in the `docs/` directory:

- [User Guide](docs/user_guide/) - How to use the pipeline
- [Developer Guide](docs/developer_guide/) - Contributing and extending
- [API Reference](docs/api/) - Function and class documentation

To build documentation:
```bash
cd docs
make html
```

## Support

- **Issues**: Report bugs at [GitHub Issues](https://github.com/samos/pipeline/issues)
- **Questions**: Contact the SAMOS team
- **Documentation**: See `docs/` directory

## License

MIT License (update as needed)

## Citation

If you use this pipeline in your research, please cite:

```
@software{samos_pipeline,
  author = {SAMOS Team},
  title = {SAMOS Data Reduction Pipeline},
  year = {2025},
  version = {2.0.0},
  url = {https://github.com/samos/pipeline}
}
```

## Changelog

### Version 2.0.0 (2025-01-15)
- Complete restructuring into modular Python package
- Added automated pipeline scripts
- Separated interactive notebooks from source code
- Added configuration file support
- Prepared for PypeIt integration

### Version 1.0.0 (2024)
- Initial notebook-based pipeline
