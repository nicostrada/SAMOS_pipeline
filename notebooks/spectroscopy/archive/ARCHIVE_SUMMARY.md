# Archive Summary - SAMOS Spectroscopy Notebooks

**Archive Date:** December 15, 2025

---

## What's In This Archive

This archive contains the previous iteration of SAMOS spectroscopy notebooks and associated documentation from the transition period when we migrated from the old `Class_SAMOS` structure to the new modular `samos` package.

---

## Archived Notebooks

### `old_notebooks/` Directory

**Working Notebooks (OLD versions):**
- `01_initial_inspection_OLD.ipynb` - Original initial reduction (before modular rewrite)
- `04_trace_extraction_OLD.ipynb` - Complete trace extraction algorithms (production-ready)
- `05_wavelength_calibration_OLD.ipynb` - Complete wavelength calibration (production-ready)
- `06_apply_calibration_OLD.ipynb` - Complete calibration application (production-ready)

**Placeholder/Conceptual Notebooks:**
- `02_calibration_frames.ipynb` - Placeholder for detailed calibration QC
- `03_trace_identification.ipynb` - Placeholder for interactive trace tuning
- `04_trace_extraction.ipynb` - Conceptual framework (simplified version)
- `05_wavelength_calibration.ipynb` - Conceptual framework (simplified version)
- `06_apply_calibration.ipynb` - Conceptual framework (simplified version)

### Why These Were Archived

**OLD notebooks (01, 04, 05, 06):**
- Contain complete, working algorithms
- Use old `Class_SAMOS` structure (no longer maintained)
- Preserved as reference for algorithm implementation
- **Can still be used** if you update the imports manually

**Placeholder notebooks (02, 03):**
- Created during transition
- Provided examples and API documentation
- Functionality now covered in main notebooks or modules

**Conceptual notebooks (04, 05, 06 - non-OLD):**
- Created as simplified examples during transition
- Showed workflow concepts
- Replaced by production-ready modules in `samos/spectroscopy/`

---

## Archived Documentation

### `old_documentation/` Directory

**Transition Documentation:**
- `README.md` - Previous README (before workflow reorganization)
- `NOTEBOOK_UPDATE_SUMMARY.md` - Summary of initial notebook updates
- `NOTEBOOKS_04-06_STATUS.md` - Status report on notebooks 04-06
- `UPDATE_COMPLETE.md` - Completion report for notebook updates
- `VERIFICATION_REPORT.md` - Verification of notebook consistency

**PypeIt Integration Documentation:**
- `PYPEIT_WORKFLOW_DESIGN.md` - Design document for PypeIt integration
- `PYPEIT_INTEGRATION_STATUS.md` - Status of PypeIt implementation

**Tools:**
- `verify_notebooks.py` - Script to verify notebook consistency

### Why This Documentation Was Archived

- Represents transition phase
- Historical record of migration process
- Useful for understanding design decisions
- No longer needed for current workflow

---

## Current Workflow (After Archive)

### Active Notebooks

```
notebooks/spectroscopy/
├── 01_initial_inspection.ipynb    - Initial reduction (READY)
├── 04_visualization.ipynb         - jdaviz visualization (READY)
└── archive/                       - This archive
```

### Next Steps (In Development)

```
├── 02_visual_qa.ipynb            - Quality assessment (TO CREATE)
└── 03_spectroscopy_pypeit.ipynb  - Spectroscopic reduction (TO CREATE)
```

### Production Modules

All core functionality now in modules:

```python
from samos.core import mosaic, cosmic_rays, calibration
from samos.spectroscopy import trace_detection, extraction, pypeit_wrapper
from samos.utils import display, io
```

---

## How To Use Archived Content

### Option 1: Reference OLD Algorithms

If you want to use the complete, tested algorithms:

1. **Open an OLD notebook** (e.g., `04_trace_extraction_OLD.ipynb`)
2. **Extract the algorithm code** (it's well-documented)
3. **Update imports** to use new `samos` package:
   ```python
   # Old
   from Class_SAMOS import SAMOS
   SAMOS = SAMOS(data_path)

   # New
   from samos.core import mosaic, cosmic_rays
   from samos.utils import display, io
   ```
4. **Run with your data**

### Option 2: Use as Learning Resource

The OLD notebooks contain detailed explanations of:
- Trace fitting algorithms
- Wavelength calibration procedures
- Line matching techniques
- Optimal extraction methods

These are useful for understanding the theory even if you use the new modules.

### Option 3: Compare Implementations

Compare OLD notebook algorithms with new module implementations:
- OLD: `04_trace_extraction_OLD.ipynb`
- NEW: `samos/spectroscopy/extraction.py`

This helps understand how the algorithms were modernized.

---

## Migration History

### Phase 1: Initial Migration (Dec 15, 2025)
- Rewrote notebook 01 using new `samos` package
- Created placeholder notebooks 02, 03
- Created conceptual notebooks 04, 05, 06
- Updated notebook 07 (visualization)

### Phase 2: Module Development (Dec 15, 2025)
- Created `samos.spectroscopy.extraction` module
- Created `samos.spectroscopy.pypeit_wrapper` module
- Installed and integrated PypeIt

### Phase 3: Cleanup (Dec 15, 2025)
- Archived transition notebooks
- Consolidated documentation
- Streamlined to 2 working notebooks + 2 planned

---

## Key Algorithms Preserved

The OLD notebooks contain production-ready implementations of:

**Trace Extraction (04_OLD):**
- Edge detection using derivative analysis
- Polynomial fitting with sigma clipping
- Spatial mask creation
- Spectrum rectification
- Background subtraction
- Boxcar and optimal extraction

**Wavelength Calibration (05_OLD):**
- Arc line detection
- HgArNe line identification
- Cross-correlation for initial alignment
- Polynomial wavelength solution fitting
- Quality assessment (residuals, RMS)

**Calibration Application (06_OLD):**
- Batch processing of all slits
- Wavelength solution propagation
- FITS file creation with WCS headers
- Text file output (wavelength, flux)

All these algorithms are now implemented in the new modules but the OLD notebooks provide:
- Original tested code
- Detailed comments and explanations
- Step-by-step visualization
- Real-world examples

---

## Important Notes

### These Notebooks Still Work

The OLD notebooks can still be used if you:
1. Have the old `Class_SAMOS.py` file, OR
2. Update the imports to use the new `samos` package

### Don't Delete This Archive

Keep this archive because:
- OLD notebooks contain complete, tested algorithms
- Transition documentation shows design evolution
- Useful reference for understanding spectroscopic reduction
- May be needed to compare results

### Recommended Approach

For new work:
1. Use `01_initial_inspection.ipynb` (updated, working)
2. Use new modules (`samos.spectroscopy.extraction`, etc.)
3. Wait for `03_spectroscopy_pypeit.ipynb` (consolidated workflow)

For reference:
1. Check OLD notebooks for algorithm details
2. Read transition documentation to understand decisions
3. Compare with new module implementations

---

## File Inventory

### Notebooks Archived: 9 files
- 4 OLD versions (complete algorithms)
- 5 transition versions (placeholders/concepts)

### Documentation Archived: 8 files
- 5 transition documents
- 2 PypeIt integration documents
- 1 verification tool

### Total Archive Size: ~11 MB
- Mostly notebook output (plots, images)
- Preserves complete working examples

---

## Contact

For questions about archived content:
- See main `README.md` in parent directory
- Check module documentation in `samos/spectroscopy/`
- Review current workflow documentation

---

Last updated: December 15, 2025
