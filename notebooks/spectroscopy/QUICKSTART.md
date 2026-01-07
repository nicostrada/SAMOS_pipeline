# SAMOS Spectroscopy - Quick Reference

**One-page guide for common tasks**

---

## 🚀 Quick Start (5 minutes)

```bash
# 1. Activate environment
conda activate samos

# 2. Launch Jupyter
cd /Users/nestrada/Documents/SAMOS/Pipeline/notebooks/spectroscopy
jupyter lab

# 3. Run notebook 01_initial_inspection.ipynb
#    (Follow cells in order)

# 4. Results: spec_*.fits files created
```

---

## 📓 Notebook Workflow

| # | Notebook | Status | Time | Purpose |
|---|----------|--------|------|---------|
| 01 | initial_inspection | ✅ Ready | 5 min | Reduce raw data → spec_*.fits |
| 02 | visual_qa | 🔄 Coming | 10 min | QA and inspection |
| 03 | spectroscopy_pypeit | 🔄 Coming | 15 min | Extract + wavelength cal |
| 04 | visualization | ✅ Ready | Variable | Interactive jdaviz |

**Total time (when complete): ~30-45 minutes**

---

## 🐍 Python API (Use Now)

### Complete Reduction in Python

```python
from pathlib import Path
from samos.spectroscopy import extraction, pypeit_wrapper

# Batch process all slits
spec_files = sorted(Path('.').glob('spec_*.fits'))
results = pypeit_wrapper.batch_wavelength_calibration(
    spec_files,
    reference_slit=10,
    output_dir='./calibrated'
)

# Output: spec1d_*.fits and spec2d_*.fits
```

### Custom Extraction

```python
from samos.spectroscopy import extraction

results = extraction.full_extraction_pipeline('spec_010.fits')

# Access results
spec_1d = results['spec1d_optimal']  # Best S/N
spec_2d = results['spec2d_rect']     # Rectified 2D
trace = results['trace_poly']         # Trace polynomial
```

---

## 📦 Key Modules

```python
# Import everything you need
from samos.core import mosaic, cosmic_rays, calibration
from samos.spectroscopy import trace_detection, extraction, pypeit_wrapper
from samos.utils import display, io
```

### Most Used Functions

| Function | Purpose |
|----------|---------|
| `mosaic.read_sami_mosaic()` | Read multi-CCD FITS |
| `cosmic_rays.remove_cosmic_rays()` | Clean cosmic rays |
| `trace_detection.detect_traces_from_arc()` | Find slits |
| `extraction.full_extraction_pipeline()` | Extract 1D spectrum |
| `pypeit_wrapper.batch_wavelength_calibration()` | Wavelength cal |

---

## 📂 Data Paths

### Configure in Notebook 01

```python
# Raw data location
data_directory = "/Users/nestrada/Documents/SAMOS/_Run2/Data/SAMI"

# Output location
analysis_top_directory = "/Users/nestrada/Documents/SAMOS/SAMOS_REDUCED"

# Target information
target_name = "ABELL3120"
target_mode = "SAMI_manual_Mask_T00_Low_Red"
```

### Output Structure

```
SAMOS_REDUCED/
└── ABELL3120/
    └── SAMI_manual_Mask_T00_Low_Red/
        ├── BIAS.fits              # Master bias
        ├── FLAT.fits              # Master flat
        ├── ARC.fits               # Arc lamp
        ├── spec_000.fits          # Slit cutouts
        ├── spec1d_000.fits        # 1D spectra
        └── spec2d_000.fits        # 2D spectra
```

---

## 🔧 Common Tasks

### View a 2D Spectrum

```python
from astropy.io import fits
from samos.utils import display

hdul = fits.open('spec_010.fits')
display.display_image(hdul['DATA'].data, title='Slit 10')
```

### Extract Custom Width

```python
from samos.spectroscopy import extraction

results = extraction.full_extraction_pipeline(
    'spec_010.fits',
    config={'extraction': {'width': 15}}  # Wider extraction
)
```

### Check Wavelength Solution

```python
from samos.spectroscopy import pypeit_wrapper
import numpy as np

# Get solution
arc_1d = np.load('arc_spectrum.npy')  # Your arc
solution = pypeit_wrapper.run_wavelength_calibration_simple(arc_1d)

print(f"RMS: {solution['rms']:.3f} Å")
print(f"Lines matched: {solution['n_lines']}")
```

---

## ❗ Troubleshooting

### Problem: Import errors

```bash
cd /Users/nestrada/Documents/SAMOS/Pipeline
pip install -e .
```

### Problem: No cosmic ray module

**Fixed!** Module uses correct import now.

### Problem: Kernel not found

```bash
conda activate samos
python -m ipykernel install --user --name samos
```

### Problem: Notebooks look different

Check you're using current notebooks, not archived ones:
- ✅ Use: `01_initial_inspection.ipynb`
- ❌ Avoid: `archive/old_notebooks/01_initial_inspection_OLD.ipynb`

---

## 📚 Learn More

| Resource | Location |
|----------|----------|
| Full README | [README.md](README.md) |
| Module docs | `help(function_name)` in Python |
| Archive | [archive/ARCHIVE_SUMMARY.md](archive/ARCHIVE_SUMMARY.md) |
| Original algorithms | `archive/old_notebooks/*_OLD.ipynb` |

---

## 🎯 Typical Session

```bash
# Start
conda activate samos
jupyter lab

# Run notebooks
# 01_initial_inspection.ipynb
# → Creates spec_*.fits files

# Python reduction
python
>>> from samos.spectroscopy import pypeit_wrapper
>>> from pathlib import Path
>>> files = sorted(Path('.').glob('spec_*.fits'))
>>> results = pypeit_wrapper.batch_wavelength_calibration(files)
>>> exit()

# Visualize
# 04_visualization.ipynb
```

---

## ⚡ Tips

1. **Always run notebook 01 first** - Creates all needed files
2. **Use reference slit 10** - Usually good quality
3. **Check archive for algorithms** - OLD notebooks have complete code
4. **Module > Notebook for automation** - Faster, reproducible
5. **jdaviz for inspection** - Notebook 04, interactive tools

---

Last updated: December 15, 2025
