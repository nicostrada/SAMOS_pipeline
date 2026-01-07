# Status of Notebooks 04-06

## Summary

Notebooks 04, 05, and 06 contain valuable spectroscopic reduction code but **need to be updated** to use the new `samos` package structure.

Currently, they use the old `Class_SAMOS` which is no longer available.

---

## Current Status

| Notebook | Function | Status | Priority |
|----------|----------|---------|----------|
| 04_trace_extraction.ipynb | Trace fitting, rectification, 1D extraction | ⚠️ **Needs update** | Medium |
| 05_wavelength_calibration.ipynb | Arc line ID, wavelength solution | ⚠️ **Needs update** | **HIGH** |
| 06_apply_calibration.ipynb | Apply wavelength cal to all slits | ⚠️ **Needs update** | High |

---

## What Each Notebook Does

### 04_trace_extraction.ipynb

**Purpose**: Extract and rectify 2D spectra

**Key steps**:
1. Fit polynomial to trace curvature
2. Perform flat fielding
3. Background subtraction
4. Rectify curved spectra to straight traces
5. Extract 1D spectra (boxcar or optimal)

**Needed for**: Creating `spec2d_*.fits` and `spec1d_*.fits` files

**Update needed**: Replace `Class_SAMOS` calls with `samos` package functions

---

### 05_wavelength_calibration.ipynb

**Purpose**: Wavelength calibration from arc lamp

**Key steps**:
1. Load arc lamp spectrum for reference slit
2. Identify emission lines (HgArNe)
3. Match to line list
4. Fit wavelength solution (polynomial)
5. Evaluate quality of fit
6. Save wavelength solution

**Needed for**: Assigning wavelengths to pixel positions

**Update needed**:
- Replace `Class_SAMOS` calls
- Create `samos.spectroscopy.wavelength_cal` module
- Possibly integrate PypeIt for robust wavelength calibration

**Priority**: **HIGH** - This is the next critical step after trace detection

---

### 06_apply_calibration.ipynb

**Purpose**: Apply wavelength calibration to all slits

**Key steps**:
1. Load wavelength solution from slit #10 (reference)
2. Apply to all other slits
3. Create wavelength-calibrated 1D spectra
4. Save as standard FITS format with WCS headers

**Needed for**: Final science-ready spectra

**Update needed**: Replace `Class_SAMOS` calls

---

## Workaround: What You Can Do Now

### Option 1: Wait for Updates

These notebooks will be updated soon. Check back or request priority update.

### Option 2: Use Current Notebooks with Modifications

The notebooks still contain valuable algorithms. You can:

1. **Read the old notebooks** to understand the workflow
2. **Extract the core algorithms** (polynomial fitting, line matching, etc.)
3. **Rewrite** using new `samos` package functions

### Option 3: Use External Tools

For wavelength calibration, consider:
- **PypeIt**: Professional-grade wavelength calibration
  - https://pypeit.readthedocs.io/
  - Install: `pip install pypeit`
  - More robust than custom solutions

- **IRAF/PyRAF**: Traditional but reliable
  - identify, reidentify, dispcor

- **specutils**: Python package for spectral analysis
  - Already installed in your environment
  - Has some wavelength calibration utilities

---

## Update Plan

### Phase 1: Wavelength Calibration (High Priority)

**Goal**: Update notebook 05 to use new structure

**Tasks**:
1. Create `samos/spectroscopy/wavelength_cal.py` module
2. Extract line identification algorithms from old notebook
3. Implement wavelength solution fitting
4. Create interactive plots for quality control
5. Update notebook 05 to use new module

**Timeline**: Can be done in 1-2 hours of focused work

### Phase 2: Trace Extraction (Medium Priority)

**Goal**: Update notebook 04

**Tasks**:
1. Create `samos/spectroscopy/extraction.py` module
2. Implement trace fitting (polynomial or spline)
3. Implement background subtraction
4. Implement rectification
5. Implement 1D extraction (boxcar and optimal)
6. Update notebook 04

**Timeline**: 2-3 hours

### Phase 3: Apply Calibration (Lower Priority)

**Goal**: Update notebook 06

**Tasks**:
1. Simple application of wavelength solution
2. WCS header creation
3. Standard FITS format output

**Timeline**: 1 hour

---

## Request Updates

If you need these notebooks updated urgently:

1. **Contact the development team**
2. **Specify which notebook** is most critical for your work
3. **Share your use case** to help prioritize

Or:

1. **Try updating them yourself** using the patterns from notebook 01
2. **Share your updated version** with the team
3. **Contribute back** to help other users

---

## Example: How to Update

### Old Code (Class_SAMOS):
```python
from Class_SAMOS import SAMOS
SAMOS = SAMOS(data_path)

# Read file
hdu = SAMOS.read_SAMI_mosaic(file)

# Remove cosmic rays
cleaned = SAMOS.CR_correct(image)

# Display
SAMOS.display_image(image, zmin, zmax)
```

### New Code (samos package):
```python
from samos.core import mosaic, cosmic_rays
from samos.utils import display

# Read file
hdu = mosaic.read_sami_mosaic(file)

# Remove cosmic rays
cleaned = cosmic_rays.remove_cosmic_rays(image)

# Display
display.display_image(image, zmin=zmin, zmax=zmax)
```

---

## Current Workflow Recommendation

**Until notebooks 04-06 are updated**, here's the recommended workflow:

### ✅ **You CAN do** (Working now):

1. **Run Notebook 01** - Complete initial reduction
   - Creates master calibrations
   - Detects traces
   - Extracts slit cutouts as `spec_*.fits`

2. **Run Notebook 07** - Visualization
   - View 2D spectra from notebook 01
   - (Won't have 1D or wavelength calibration yet)

### ⏳ **Coming Soon** (Needs updates):

3. **Notebook 04** - Extract 1D spectra from 2D cutouts
4. **Notebook 05** - Wavelength calibration
5. **Notebook 06** - Apply calibration to all slits

---

## Bottom Line

**For now**: Notebook 01 gets you 80% of the way to reduced data.

**Next critical step**: Wavelength calibration (notebook 05)

**Timeline**: Updates can be done, just need prioritization based on user needs.

---

Last updated: December 15, 2025
