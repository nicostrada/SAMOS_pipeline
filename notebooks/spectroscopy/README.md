# SAMOS Spectroscopy Notebooks

**Current Workflow - Production Ready**

---

## Quick Start

### Setup

```bash
# Activate environment
conda activate samos

# Launch Jupyter
cd /Users/nestrada/Documents/SAMOS/Pipeline/notebooks/spectroscopy
jupyter lab
```

### Select Kernel
Choose: **"Python 3.13 (SAMOS)"** in Jupyter

---

## Current Workflow

### Step 1: Initial Reduction ✅

**Notebook:** [01_initial_inspection.ipynb](01_initial_inspection.ipynb)

**What it does:**
- Reads raw SAMOS multi-CCD data
- Creates master calibration frames (bias, flat, arc)
- Removes cosmic rays
- Automatically detects spectral traces
- Extracts individual slit cutouts

**Inputs:**
- Raw science frames
- Bias frames
- Flat frames
- Arc lamp frames

**Outputs:**
- `BIAS.fits`, `FLAT.fits`, `ARC.fits`
- `spec_000.fits` to `spec_NNN.fits` (one per slit)

**Time:** ~2-5 minutes

**Status:** ✅ Production ready

---

### Step 2: Visual QA & Inspection ✅

**Notebook:** [02_visual_qa.ipynb](02_visual_qa.ipynb)

**What it does:**
- Interactive visualization of all slits
- Quality assessment metrics (S/N, arc brightness)
- Identify reference slit for wavelength calibration
- Flag bad slits
- Generate QA report (PDF)
- Export configuration for next step

**Outputs:**
- `qa_config.yaml` - Configuration for Step 3
- `qa_report.pdf` - QA summary (optional)

**Time:** ~5-10 minutes

**Status:** ✅ Ready to use

---

### Step 3: Spectroscopic Reduction ✅

**Notebook:** [03_spectroscopy_pypeit.ipynb](03_spectroscopy_pypeit.ipynb)

**What it does:**

**Part A: Trace Fitting & Extraction**
- Fit trace curvature from flat field
- Create spatial masks (slit vs sky)
- Rectify 2D spectra
- Background/sky subtraction
- Extract 1D spectra (boxcar + optimal)

**Part B: Wavelength Calibration**
- Prepare arc lamp spectrum
- Identify HgArNe emission lines
- Fit wavelength solution (polynomial)
- Assess calibration quality

**Part C: Batch Processing**
- Apply to all slits
- Create wavelength-calibrated FITS files
- Generate summary statistics

**Outputs:**
- `calibrated/spec1d_*.fits` - Wavelength-calibrated 1D spectra
- `calibrated/spec2d_*.fits` - Wavelength-calibrated 2D spectra
- `wavelength_solution_*.txt` - Solution coefficients

**Time:** ~10-20 minutes

**Status:** ✅ Production ready

---

### Step 4: Visualization ✅

**Notebook:** [04_visualization.ipynb](04_visualization.ipynb)

**What it does:**
- Interactive visualization using jdaviz
- Specviz for 1D spectra
- Mosviz for multi-object spectroscopy

**Inputs:**
- `calibrated/spec1d_*.fits` files (from Step 3)
- `calibrated/spec2d_*.fits` files (from Step 3)

**Status:** ✅ Ready to use

---

## Using the Modules Directly

While notebooks 02 and 03 are in development, you can use the modules directly:

### Example: Complete Spectroscopic Reduction

```python
from pathlib import Path
from samos.spectroscopy import extraction, pypeit_wrapper

# Step 1: Extract all slits
spec_files = sorted(Path('.').glob('spec_*.fits'))

# Step 2: Batch wavelength calibration and extraction
results = pypeit_wrapper.batch_wavelength_calibration(
    spec_files,
    reference_slit=10,
    output_dir='./calibrated'
)

print(f"✓ Processed {len(results)} slits")
print(f"✓ Output: spec1d_*.fits and spec2d_*.fits")
```

### Example: Custom Extraction

```python
from samos.spectroscopy import extraction

# Custom extraction for single slit
config = {
    'trace_fit': {'method': 'polynomial', 'degree': 3},
    'extraction': {'width': 10},
    'background': {'method': 'median'}
}

results = extraction.full_extraction_pipeline('spec_010.fits', config)

# Access results
spec_1d_boxcar = results['spec1d_boxcar']
spec_1d_optimal = results['spec1d_optimal']
spec_2d_rect = results['spec2d_rect']
```

---

## Available Modules

### Core Modules

```python
from samos.core import mosaic, cosmic_rays, calibration
from samos.spectroscopy import trace_detection, extraction, pypeit_wrapper
from samos.utils import display, io
```

### Key Functions

**Extraction:**
- `extraction.fit_trace()` - Fit trace curvature
- `extraction.rectify_spectrum()` - Straighten spectra
- `extraction.subtract_background()` - Sky subtraction
- `extraction.extract_1d_boxcar()` - Simple extraction
- `extraction.extract_1d_optimal()` - Optimal extraction
- `extraction.full_extraction_pipeline()` - Complete pipeline

**Wavelength Calibration:**
- `pypeit_wrapper.run_wavelength_calibration_simple()` - Wavelength solution
- `pypeit_wrapper.batch_wavelength_calibration()` - Batch processing

See module docstrings for full API documentation:
```python
help(extraction.fit_trace)
```

---

## Data Organization

### Recommended Structure

```
~/SAMOS_Data/
├── _Run2/Data/SAMI/          # Raw data
│   ├── 20241016/             # Bias frames
│   └── 20241017/             # Science, flats, arcs
│
└── SAMOS_REDUCED/            # Pipeline outputs
    └── ABELL3120/
        └── SAMI_manual_Mask_T00_Low_Red/
            ├── BIAS.fits
            ├── FLAT.fits
            ├── ARC.fits
            ├── spec_*.fits           # From notebook 01
            ├── spec1d_*.fits         # From notebook 03
            └── spec2d_*.fits         # From notebook 03
```

---

## Archive

Previous notebooks and documentation are archived in:
- `archive/old_notebooks/` - Previous notebook versions
- `archive/old_documentation/` - Transition documentation

See [archive/ARCHIVE_SUMMARY.md](archive/ARCHIVE_SUMMARY.md) for details.

**Note:** The archived OLD notebooks contain complete, tested algorithms that still work with manual import updates.

---

## Troubleshooting

### Import Errors

```bash
cd /Users/nestrada/Documents/SAMOS/Pipeline
conda activate samos
pip install -e .
```

### Kernel Not Found

```bash
conda activate samos
python -m ipykernel install --user --name samos --display-name "Python 3.13 (SAMOS)"
```

### Missing Dependencies

```bash
conda activate samos
pip install -r requirements.txt
```

---

## What's Working Now

✅ **Notebook 01** - Complete initial reduction
✅ **Notebook 02** - Visual QA and inspection
✅ **Notebook 03** - Complete spectroscopic reduction with wavelength calibration
✅ **Notebook 04** - Interactive visualization
✅ **Modules** - All spectroscopy modules functional
✅ **Python API** - Direct module usage for automation

**Complete workflow ready for production use!**

---

## Getting Help

### Documentation
- Module docstrings: `help(function_name)`
- Archive documentation: [archive/](archive/)
- Main pipeline docs: `../../README.md`

### Examples
- Check notebook 01 for complete reduction example
- See module docstrings for API examples
- Review archived notebooks for algorithm details

---

## Summary

**Current Status:**
- Initial reduction: ✅ Working
- Visualization: ✅ Working
- Spectroscopy modules: ✅ Ready
- Consolidated notebooks: 🔄 In development

**Recommended Path:**
1. Run notebook 01 for initial reduction
2. Use modules directly for spectroscopic reduction (or wait for notebook 03)
3. Use notebook 04 for visualization

---

Last updated: December 15, 2025
