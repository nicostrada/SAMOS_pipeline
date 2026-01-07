# SAMOS Spectroscopy Notebooks

Interactive workflow for SAMOS spectroscopic data reduction.

## Quick Start

1. **Activate environment**:
   ```bash
   conda activate samos
   ```

2. **Launch Jupyter**:
   ```bash
   cd /Users/nestrada/Documents/SAMOS/Pipeline/notebooks/spectroscopy
   jupyter lab
   ```

3. **Select kernel**: "Python 3.13 (SAMOS)"

4. **Run notebooks in order** (see below)

---

## Notebook Workflow

### ✅ 01_initial_inspection.ipynb **[UPDATED - READY TO USE]**

**What it does:**
- Reads raw SAMOS multi-CCD data
- Creates master calibration frames (bias, flat, arc)
- Detects spectral traces automatically
- Extracts individual slit cutouts
- Saves multi-extension FITS files for each slit

**Inputs needed:**
- Raw science frames (e.g., `target.027.fits`)
- Bias frames (e.g., `bias.055-074.fits`)
- Flat frames (e.g., `target.028.fits`, `calibration.029.fits`)
- Arc lamp frames (e.g., `calibration.032.fits`)
- Arc dark frame (optional, e.g., `calibration.035.fits`)

**Outputs:**
- `BIAS.fits` - Master bias frame
- `FLAT.fits` - Master flat frame
- `ARC.fits` - Arc lamp frame
- `spec_000.fits` to `spec_NNN.fits` - Individual slit cutouts

**Time to run:** ~2-5 minutes (depending on number of frames)

**Status:** ✅ Fully updated to use new `samos` package

---

### 📝 02_calibration_frames.ipynb **[PLACEHOLDER CREATED]**

**Future enhancement** for detailed calibration frame quality control.

For now, calibration is fully covered in notebook 01.

**What's there**: Placeholder with examples and API docs

---

### 📝 03_trace_identification.ipynb **[PLACEHOLDER CREATED]**

**Future enhancement** for interactive trace identification and fine-tuning.

For now, automated trace detection in notebook 01 works well.

**What's there**: Placeholder with advanced examples and manual override options

---

### ⚠️ 04_trace_extraction.ipynb **[UPDATED - CONCEPTUAL FRAMEWORK]**

**What it does:**
- Fits polynomial to trace curvature
- Performs flat fielding and background subtraction
- Rectifies curved spectra to straight traces
- Extracts 1D spectra

**Status:** ⚠️ Updated with simplified workflow and concepts

**Note:** For full algorithms, see `04_trace_extraction_OLD.ipynb` or use PypeIt

**What's available:** Conceptual examples showing the key steps with working code snippets

---

### ⚠️ 05_wavelength_calibration.ipynb **[UPDATED - CONCEPTUAL FRAMEWORK]**

**What it does:**
- Wavelength calibration for reference slit
- Identifies HgArNe emission lines
- Fits wavelength solution
- Saves polynomial coefficients

**Status:** ⚠️ Updated with simplified workflow and concepts

**Note:** For full algorithms, see `05_wavelength_calibration_OLD.ipynb` or use PypeIt

**What's available:** Arc line detection, wavelength fitting examples, quality assessment

---

### ⚠️ 06_apply_calibration.ipynb **[UPDATED - CONCEPTUAL FRAMEWORK]**

**What it does:**
- Applies wavelength calibration to all slits
- Creates wavelength-calibrated 1D and 2D spectra
- Saves as `spec1d_*.fits` and `spec2d_*.fits`

**Status:** ⚠️ Updated with simplified workflow and concepts

**Note:** For full algorithms, see `06_apply_calibration_OLD.ipynb` or use PypeIt

**What's available:** Batch processing example, Spectrum1D output, verification plots

---

### ✅ 07_visualization.ipynb **[UPDATED - READY TO USE]**

**What it does:**
- Interactive visualization using jdaviz
- Uses Specviz for viewing individual 1D spectra
- Uses Mosviz for multi-object spectroscopy (2D + images)

**Inputs needed:**
- `spec1d_*.fits` files (from notebooks 04-06, when updated)
- `spec2d_*.fits` files (from notebooks 04-06, when updated)
- `tile_*.fits` files (optional cutout images)

**Status:** ✅ Fully updated, ready to use

**Note**: Will work best after notebooks 04-06 are updated. For now, you can visualize the raw 2D spectra from notebook 01.

---

## Updated Notebook Structure

The updated notebook 01 uses the **new `samos` package**:

```python
# NEW IMPORTS
from samos.core import mosaic, cosmic_rays, calibration
from samos.spectroscopy import trace_detection
from samos.utils import display, io

# Example: Read SAMOS mosaic
hdu = mosaic.read_sami_mosaic('science.fits')

# Example: Remove cosmic rays
cleaned = cosmic_rays.remove_cosmic_rays(hdu.data, verbose=True)

# Example: Detect traces
slits = trace_detection.detect_traces_from_arc(
    arc_image,
    peak_threshold=4e5,
    margin=12,
    verbose=True
)
```

### Old vs New

**Old way** (Class_SAMOS):
```python
from Class_SAMOS import SAMOS
SAMOS = SAMOS(data_path)
hdu = SAMOS.read_SAMI_mosaic(file)
cleaned = SAMOS.CR_correct(image)
SAMOS.display_image(image, zmin, zmax)
```

**New way** (samos package):
```python
from samos.core import mosaic, cosmic_rays
from samos.utils import display

hdu = mosaic.read_sami_mosaic(file)
cleaned = cosmic_rays.remove_cosmic_rays(image)
display.display_image(image, zmin=zmin, zmax=zmax)
```

---

## Data Organization

**Recommended structure:**

```
~/SAMOS_Data/
├── _Run2/Data/SAMI/          # Your raw data (unchanged)
│   ├── 20241016/
│   │   └── bias.*.fits
│   ├── 20241017/
│   │   ├── target.*.fits
│   │   └── calibration.*.fits
│
└── SAMOS_REDUCED/            # Pipeline outputs
    └── ABELL3120/
        └── SAMI_manual_Mask_T00_Low_Red/
            ├── BIAS.fits
            ├── FLAT.fits
            ├── ARC.fits
            └── spec_*.fits
```

**Configure in notebook:**
```python
# Set these paths in notebook 01, cell 2
data_directory = "/Users/nestrada/Documents/SAMOS/_Run2/Data/SAMI"
analysis_top_directory = "/Users/nestrada/Documents/SAMOS/SAMOS_REDUCED"
target_name = "ABELL3120"
```

---

## Troubleshooting

### ModuleNotFoundError: No module named 'samos'

**Solution:**
```bash
cd /Users/nestrada/Documents/SAMOS/Pipeline
conda activate samos
pip install -e .
```

### Kernel not found

**Solution:**
```bash
conda activate samos
python -m ipykernel install --user --name samos --display-name "Python 3.13 (SAMOS)"
```

Then in Jupyter: Kernel → Change Kernel → "Python 3.13 (SAMOS)"

### Import errors in notebooks

Make sure you're using the **updated** notebook (`01_initial_inspection.ipynb`), not the OLD version.

If you accidentally opened the OLD notebook, it's saved as `01_initial_inspection_OLD.ipynb`.

---

## What's Next?

### Immediate (Ready Now):

1. **Run notebook 01** - Complete initial data reduction
2. **Inspect outputs** - Check `spec_*.fits` files
3. **Try notebook 07** - Visualize results with jdaviz

### Coming Soon (Needs Work):

4. **Notebooks 04-06** - Will be updated to use new package
5. **Wavelength calibration** - Dedicated module in progress
6. **1D extraction** - Optimal extraction algorithm

### Future Enhancements:

7. **PypeIt integration** - Advanced wavelength calibration
8. **Automated QA** - Quality control plots
9. **Batch processing** - Process multiple targets at once

---

## Need Help?

- **Documentation**: See main `README.md` in Pipeline root
- **Examples**: Check `QUICKSTART.md` for common tasks
- **API docs**: Use `help(function_name)` in Python

```python
# Get help on any function
from samos.core import mosaic
help(mosaic.read_sami_mosaic)
```

---

Happy reducing! 🔭
