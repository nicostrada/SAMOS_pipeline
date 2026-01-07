# SAMOS Pipeline Installation Guide

## Prerequisites

- Python 3.10 or newer
- conda or mamba (recommended) or pip
- ~5 GB disk space for dependencies

## Recommended Installation (Conda)

### Step 1: Create Clean Environment

```bash
# Create new conda environment with Python 3.13
conda create -n samos python=3.13
conda activate samos
```

### Step 2: Navigate to Pipeline Directory

```bash
cd /Users/nestrada/Documents/SAMOS/Pipeline
```

### Step 3: Install Pipeline

```bash
# Install in development mode (recommended)
pip install -e .
```

This will:
- Install all required dependencies from `requirements.txt`
- Make the `samos` package importable from anywhere
- Allow you to edit the code without reinstalling

### Step 4: Verify Installation

```bash
# Test imports
python -c "import samos; from samos.core import mosaic; print('Success!')"

# Check version
python -c "import samos; print(samos.__version__)"
```

## Alternative: Using Existing SAMOS Environment

If you already have the `samos` conda environment from PipelineV2:

```bash
# Activate existing environment
conda activate samos

# Navigate to Pipeline directory
cd /Users/nestrada/Documents/SAMOS/Pipeline

# Install pipeline
pip install -e .
```

## Optional: Install Development Tools

For developers who want to contribute:

```bash
pip install -e .[dev]
```

This adds:
- pytest (testing)
- sphinx (documentation)
- black (code formatting)
- flake8 (linting)

## Optional: Install PypeIt Integration

For advanced wavelength calibration using PypeIt:

```bash
pip install -e .[pypeit]
```

**Note**: PypeIt has many dependencies and may take time to install.

## Jupyter Notebook Setup

### Register Kernel

To use the pipeline in Jupyter notebooks:

```bash
# Install Jupyter (if not already installed)
pip install jupyter jupyterlab

# Register kernel
python -m ipykernel install --user --name samos --display-name "Python 3.13 (SAMOS Pipeline)"
```

### Launch Jupyter

```bash
cd /Users/nestrada/Documents/SAMOS/Pipeline/notebooks/spectroscopy
jupyter lab
```

Select the "Python 3.13 (SAMOS Pipeline)" kernel when running notebooks.

## Troubleshooting

### Issue: ModuleNotFoundError: No module named 'cv2'

OpenCV may need system dependencies on some platforms:

```bash
# macOS
brew install opencv

# Then reinstall opencv-python
pip uninstall opencv-python
pip install opencv-python
```

### Issue: Import errors in notebooks

Make sure you're using the correct kernel:
- In Jupyter: Kernel → Change Kernel → "Python 3.13 (SAMOS Pipeline)"
- In VSCode: Select kernel in top-right corner

### Issue: Permission denied when installing

Don't use `sudo pip install`. Instead:
- Use conda environments (recommended)
- Or use `--user` flag: `pip install --user -e .`

### Issue: Conflicting package versions

```bash
# Create fresh environment
conda deactivate
conda remove -n samos --all
conda create -n samos python=3.13
conda activate samos
cd /Users/nestrada/Documents/SAMOS/Pipeline
pip install -e .
```

## Verification

After installation, run this verification script:

```python
# verification.py
import sys
print(f"Python: {sys.version}")

import numpy as np
print(f"NumPy: {np.__version__}")

import astropy
print(f"Astropy: {astropy.__version__}")

import samos
print(f"SAMOS Pipeline: {samos.__version__}")

from samos.core import mosaic, cosmic_rays, calibration
from samos.spectroscopy import trace_detection
from samos.utils import io, display

print("\n✓ All imports successful!")
print(f"\nInstallation location: {samos.__file__}")
```

Save and run:
```bash
python verification.py
```

## Next Steps

After successful installation:

1. **Interactive mode**: Navigate to `notebooks/spectroscopy/` and open `01_initial_inspection.ipynb`
2. **Automated mode**: Copy and edit `configs/spectroscopy_default.yaml`
3. **Read the docs**: Check `README.md` for usage examples

## Updating the Pipeline

To get latest changes:

```bash
cd /Users/nestrada/Documents/SAMOS/Pipeline
git pull  # If using git
# No need to reinstall - changes are live due to editable install (-e flag)
```

## Uninstalling

```bash
# Remove package
pip uninstall samos-pipeline

# Remove environment (if desired)
conda deactivate
conda remove -n samos --all
```
