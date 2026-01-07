# SAMOS Spectroscopy Workflow with PypeIt Integration

**Design Document**
**Date:** December 15, 2025

---

## Overview

This document outlines the new modular workflow integrating PypeIt for production-quality spectroscopic reduction while maintaining visual inspection capabilities.

---

## Design Principles

1. **Modular**: Each step can be run independently
2. **Visual First**: Interactive inspection and quality control
3. **Automated Ready**: All steps can be automated in pipeline
4. **Production Quality**: Uses PypeIt for wavelength calibration
5. **Flexible**: Can skip steps or use custom parameters

---

## Proposed Workflow Structure

### Phase 1: Initial Reduction (WORKING)
**Notebook:** `01_initial_inspection.ipynb`

**What it does:**
- Multi-CCD mosaic assembly
- Bias/flat calibration
- Cosmic ray removal
- Automatic trace detection
- Slit extraction

**Output:** `spec_*.fits` (multi-extension FITS per slit)

**Status:** ✅ Production ready

---

### Phase 2: Visual Inspection & QA (NEW)
**Notebook:** `02_visual_qa.ipynb` (consolidate 02+03)

**What it does:**
- Interactive visualization of all slits
- Quality assessment plots
- Identify reference slit for wavelength cal
- Flag bad slits
- Adjust trace parameters if needed
- Export QA report

**Output:**
- QA report (PDF/HTML)
- `qa_config.yaml` (flagged slits, reference slit)

**Status:** 🆕 To be created

---

### Phase 3: Spectroscopic Reduction with PypeIt (CONSOLIDATED)
**Notebook:** `03_spectroscopy_pypeit.ipynb` (consolidate 04+05+06)

**What it does:**

**3a. Trace Fitting & Rectification**
- Load slit cutouts from notebook 01
- Fit trace curvature (polynomial/spline)
- Create spatial mask (slit vs sky)
- Rectify 2D spectra
- Background/sky subtraction
- Flat fielding
- Extract 1D spectra (boxcar + optimal)

**3b. Wavelength Calibration (PypeIt)**
- Configure PypeIt for SAMOS
- Load arc lamp data
- Automatic line identification (HgArNe)
- Robust wavelength solution fitting
- Quality assessment
- Apply to all slits

**3c. Final Calibration**
- Apply wavelength solution to all slits
- Create wavelength-calibrated 1D/2D spectra
- Write standard FITS with WCS headers
- Generate calibration QA plots

**Output:**
- `spec1d_*.fits` - Wavelength-calibrated 1D spectra
- `spec2d_*.fits` - Wavelength-calibrated 2D spectra
- `wavelength_solution.fits` - PypeIt solution
- `reduction_qa.pdf` - Quality assessment

**Status:** 🆕 To be created

---

### Phase 4: Visualization (WORKING)
**Notebook:** `04_visualization.ipynb` (rename from 07)

**What it does:**
- Interactive jdaviz (Specviz/Mosviz)
- Spectral analysis tools
- Export for publication

**Status:** ✅ Working, rename only

---

## File Organization

### Notebooks (New Structure)

```
notebooks/spectroscopy/
├── 01_initial_inspection.ipynb     ✅ Ready (existing)
├── 02_visual_qa.ipynb              🆕 Create (consolidate 02+03)
├── 03_spectroscopy_pypeit.ipynb    🆕 Create (consolidate 04+05+06)
├── 04_visualization.ipynb          ✅ Rename (from 07)
│
├── archive/                        📁 Move old notebooks
│   ├── 02_calibration_frames.ipynb
│   ├── 03_trace_identification.ipynb
│   ├── 04_trace_extraction.ipynb
│   ├── 05_wavelength_calibration.ipynb
│   ├── 06_apply_calibration.ipynb
│   └── 07_visualization.ipynb
│
└── reference/                      📁 Keep OLD notebooks
    ├── 04_trace_extraction_OLD.ipynb
    ├── 05_wavelength_calibration_OLD.ipynb
    └── 06_apply_calibration_OLD.ipynb
```

---

## Module Development

### New Modules Required

**1. `samos/spectroscopy/extraction.py`**

```python
def fit_trace(flat, method='polynomial', degree=3):
    """Fit trace curvature."""

def rectify_spectrum(data_2d, trace_poly):
    """Straighten curved spectrum."""

def create_spatial_mask(flat, trace_poly, margin=12):
    """Separate slit from sky regions."""

def subtract_background(data, mask_sky):
    """Sky/background subtraction."""

def extract_1d_boxcar(data_2d, trace_poly, width=10):
    """Boxcar extraction."""

def extract_1d_optimal(data_2d, trace_poly, variance=None):
    """Optimal extraction with variance weighting."""
```

**2. `samos/spectroscopy/pypeit_wrapper.py`**

```python
def configure_pypeit_for_samos():
    """Create PypeIt configuration for SAMOS."""

def prepare_pypeit_files(spec_files, output_dir):
    """Convert SAMOS FITS to PypeIt format."""

def run_wavelength_calibration(arc_file, instrument='samos_red'):
    """Run PypeIt wavelength calibration."""

def apply_wavelength_solution(spec_files, solution):
    """Apply PypeIt solution to all slits."""

def extract_pypeit_qa(pypeit_dir):
    """Extract QA plots and metrics."""
```

**3. `samos/quality/qa.py`**

```python
def create_slit_qa_plots(spec_files, output_pdf):
    """Create comprehensive QA report."""

def assess_trace_quality(spec_file):
    """Assess trace detection quality."""

def assess_wavelength_solution(solution, residuals):
    """Assess wavelength calibration quality."""

def flag_bad_slits(spec_files, criteria):
    """Identify problematic slits."""
```

---

## Automated Pipeline Integration

### Pipeline Configuration

**File:** `config/spectroscopy_config.yaml`

```yaml
# SAMOS Spectroscopy Pipeline Configuration

reduction:
  # Initial reduction (notebook 01)
  cosmic_ray_removal: true
  trace_detection:
    method: 'automatic'
    peak_threshold: 4e5
    margin: 12

  # Spectroscopic reduction (notebook 03)
  trace_fitting:
    method: 'polynomial'
    degree: 3

  extraction:
    method: 'optimal'  # or 'boxcar'
    width: 10

  # Wavelength calibration
  wavelength:
    use_pypeit: true
    reference_slit: 10  # or 'auto'
    line_list: 'HgArNe'
    solution_order: 3

quality:
  # Quality assessment
  generate_qa: true
  qa_format: 'pdf'
  flag_bad_slits: true
  min_snr: 5.0

output:
  # Output format
  write_1d: true
  write_2d: true
  wcs_headers: true
  format: 'standard_fits'
```

### Automated Script

**File:** `scripts/reduce_spectroscopy.py`

```python
#!/usr/bin/env python
"""
Automated SAMOS spectroscopy reduction pipeline.

Usage:
    python reduce_spectroscopy.py --target ABELL3120 --mode SAMI_manual_Mask_T00_Low_Red
"""

from samos.pipeline import spectroscopy_pipeline
from samos.utils import config_loader

def main(target, mode, config_file='config/spectroscopy_config.yaml'):
    # Load configuration
    config = config_loader.load_config(config_file)

    # Run pipeline
    pipeline = spectroscopy_pipeline.SpectroscopyPipeline(config)

    # Step 1: Initial reduction
    pipeline.run_initial_reduction(target, mode)

    # Step 2: Generate QA (optional visual inspection)
    pipeline.generate_qa_report()

    # Step 3: Spectroscopic reduction with PypeIt
    pipeline.run_spectroscopic_reduction()

    # Step 4: Create summary
    pipeline.create_summary()

    print("✓ Spectroscopy pipeline complete!")
```

---

## Implementation Plan

### Phase 1: Setup (1-2 hours)
- ✅ Install PypeIt
- ✅ Test PypeIt with SAMOS data
- 🆕 Create PypeIt configuration for SAMOS

### Phase 2: Modules (3-4 hours)
- 🆕 Create `samos/spectroscopy/extraction.py`
- 🆕 Create `samos/spectroscopy/pypeit_wrapper.py`
- 🆕 Create `samos/quality/qa.py`

### Phase 3: Notebooks (2-3 hours)
- 🆕 Create `02_visual_qa.ipynb`
- 🆕 Create `03_spectroscopy_pypeit.ipynb`
- ✅ Rename `07→04_visualization.ipynb`
- 📁 Archive old notebooks 02-06

### Phase 4: Pipeline Integration (2-3 hours)
- 🆕 Create `config/spectroscopy_config.yaml`
- 🆕 Create `scripts/reduce_spectroscopy.py`
- 🆕 Update `samos/pipeline/spectroscopy_pipeline.py`

### Phase 5: Testing & Documentation (1-2 hours)
- 🧪 Test full workflow with real data
- 📝 Update README
- 📝 Create quickstart guide

**Total Estimated Time:** 9-14 hours

---

## Benefits of New Workflow

### For Interactive Use:
1. **Fewer notebooks** (4 instead of 7)
2. **Clearer workflow** (inspect → reduce → visualize)
3. **Better QA** (dedicated QA notebook)
4. **Production quality** (PypeIt integration)

### For Automated Use:
1. **Single command** reduction
2. **Configurable** via YAML
3. **Reproducible** results
4. **Automatic QA** generation

### For Development:
1. **Modular** code in separate modules
2. **Reusable** functions
3. **Testable** components
4. **Maintainable** structure

---

## Migration Path

Users can choose:

**Option A: Use New Workflow (Recommended)**
- Run notebooks 01 → 02 → 03 → 04
- Use PypeIt for wavelength calibration
- Modern, production-ready

**Option B: Use Old Workflow (Reference)**
- Run notebooks 01 → 04_OLD → 05_OLD → 06_OLD → 07
- Original algorithms
- Still works

**Option C: Mix and Match**
- Use notebook 01 for initial reduction
- Use PypeIt directly (command line)
- Use notebook 04 for visualization

---

## Next Steps

1. **Verify PypeIt installation**
2. **Test PypeIt with SAMOS arc lamp**
3. **Create extraction module**
4. **Create pypeit_wrapper module**
5. **Create consolidated notebook 03**
6. **Create QA notebook 02**
7. **Test full workflow**
8. **Update documentation**

---

Last updated: December 15, 2025
