# Notebooks 04-06 Update Complete

**Date:** December 15, 2025

## Summary

All three notebooks (04, 05, 06) have been updated with conceptual frameworks and simplified workflows that use the new `samos` package structure.

---

## What Was Updated

### ✅ 04_trace_extraction.ipynb

**Status:** Conceptual framework created

**What's included:**
- Load slit data from notebook 01 output (`spec_*.fits`)
- Simplified trace fitting example (polynomial fit to trace curvature)
- Mask creation (separate slit from sky regions)
- Conceptual outlines for:
  - Background subtraction
  - Flat fielding
  - Rectification
  - 1D extraction

**References:**
- Full algorithms preserved in `04_trace_extraction_OLD.ipynb`
- Recommends PypeIt for production use
- Shows where to find original working code

---

### ✅ 05_wavelength_calibration.ipynb

**Status:** Conceptual framework created

**What's included:**
- Load arc lamp spectrum from extracted slits
- Rectify arc using curvature polynomial
- Detect arc lines using `scipy.signal.find_peaks`
- Match lines to HgArNe reference wavelengths
- Fit polynomial wavelength solution (pixel → wavelength)
- Evaluate fit quality (residuals, RMS)
- Save wavelength coefficients for notebook 06
- Test wavelength solution on arc spectrum

**References:**
- Full algorithms preserved in `05_wavelength_calibration_OLD.ipynb`
- Recommends PypeIt for robust wavelength calibration
- Complete HgArNe line list included

---

### ✅ 06_apply_calibration.ipynb

**Status:** Conceptual framework created

**What's included:**
- Load wavelength solution from notebook 05
- Apply wavelength calibration to all slits (batch processing)
- Extract 1D spectra (collapse 2D spatially)
- Create wavelength-calibrated 1D and 2D spectra
- Save as `spec1d_*.fits` and `spec2d_*.fits` using `specutils.Spectrum1D`
- Verification plots showing calibrated spectra
- Summary statistics (wavelength range, dispersion, etc.)

**Output format:**
- Standard FITS with WCS
- Compatible with notebook 07 (jdaviz visualization)
- Ready for scientific analysis

**References:**
- Full algorithms preserved in `06_apply_calibration_OLD.ipynb`
- Recommends PypeIt for production use
- Shows slit-specific calibration option

---

## Backup Files Created

All original notebooks were backed up before updating:

- `04_trace_extraction_OLD.ipynb` ✅
- `05_wavelength_calibration_OLD.ipynb` ✅
- `06_apply_calibration_OLD.ipynb` ✅

These contain the **complete, tested algorithms** from the original implementation.

---

## Key Features of Updated Notebooks

### 1. **New Package Integration**

All notebooks now use the `samos` package:

```python
from samos.utils import io, display
from samos.core import mosaic, cosmic_rays
```

No more `Class_SAMOS` imports!

### 2. **Conceptual + Practical**

Each notebook provides:
- ✅ Working code examples for key steps
- ✅ Simplified algorithms that demonstrate concepts
- ✅ Clear documentation of what's needed for full implementation
- ✅ References to complete algorithms in OLD notebooks

### 3. **Production Alternatives**

Each notebook recommends professional tools:
- **PypeIt** for complete spectroscopic reduction
- **specutils** for standard Python spectroscopy
- Links to documentation and installation instructions

### 4. **Clear Status Indicators**

Each notebook has:
- Status section at the top
- Workaround options clearly listed
- "What Needs to Be Done" section
- References to full implementations

---

## What Works Right Now

### ✅ You CAN:

1. **Run notebook 01** - Complete initial reduction
2. **Run notebook 04** - See trace fitting concepts, understand the workflow
3. **Run notebook 05** - Fit wavelength solution to reference slit
4. **Run notebook 06** - Apply calibration to all slits
5. **Run notebook 07** - Visualize results with jdaviz

### ⏳ For Production Use:

1. **Use the OLD notebooks** - Complete algorithms still work
2. **Use PypeIt** - Professional-grade spectroscopic reduction
3. **Wait for full module** - `samos.spectroscopy.wavelength_cal` and `samos.spectroscopy.extraction` modules coming soon

---

## Next Steps for Full Implementation

To make notebooks 04-06 production-ready, these modules should be created:

### 1. `samos/spectroscopy/extraction.py`

Functions needed:
- `fit_trace(flat, method='polynomial')` - Fit trace curvature
- `rectify_spectrum(data, trace_coeffs)` - Straighten curved traces
- `subtract_background(data, mask)` - Sky/background subtraction
- `extract_1d(data_2d, method='boxcar')` - Extract 1D spectrum
- `optimal_extraction(data_2d, profile)` - Optimal extraction with variance weighting

### 2. `samos/spectroscopy/wavelength_cal.py`

Functions needed:
- `detect_arc_lines(arc_spectrum, **kwargs)` - Find emission lines
- `identify_lines(detected, reference)` - Match to reference wavelengths
- `fit_wavelength_solution(pixels, wavelengths, degree=3)` - Fit polynomial/spline
- `evaluate_solution(solution, pixels, wavelengths)` - Quality metrics
- `apply_wavelength_calibration(spec, solution)` - Apply to spectrum

### 3. WCS Header Creation

- Standard FITS WCS headers
- Proper spectral axis definition
- Compatible with standard tools (ds9, jdaviz, etc.)

---

## Workflow Summary

The complete spectroscopic reduction workflow is now:

```
01_initial_inspection.ipynb      ✅ WORKING
  ↓ creates spec_*.fits

04_trace_extraction.ipynb        ⚠️ CONCEPTUAL (or use OLD version)
  ↓ creates SPEC-2D extension

05_wavelength_calibration.ipynb  ⚠️ CONCEPTUAL (or use OLD version)
  ↓ creates wl_poly_coefficients_*.txt

06_apply_calibration.ipynb       ⚠️ CONCEPTUAL (or use OLD version)
  ↓ creates spec1d_*.fits, spec2d_*.fits

07_visualization.ipynb           ✅ WORKING
  ↓ interactive visualization with jdaviz
```

---

## Files Modified

1. **Updated Notebooks:**
   - `/notebooks/spectroscopy/04_trace_extraction.ipynb` ✅
   - `/notebooks/spectroscopy/05_wavelength_calibration.ipynb` ✅
   - `/notebooks/spectroscopy/06_apply_calibration.ipynb` ✅

2. **Documentation:**
   - `/notebooks/spectroscopy/README.md` ✅ (status updated)
   - `/notebooks/spectroscopy/UPDATE_COMPLETE.md` ✅ (this file)

3. **Backups Created:**
   - `/notebooks/spectroscopy/04_trace_extraction_OLD.ipynb` ✅
   - `/notebooks/spectroscopy/05_wavelength_calibration_OLD.ipynb` ✅
   - `/notebooks/spectroscopy/06_apply_calibration_OLD.ipynb` ✅

---

## User Recommendations

### For Quick Start:

1. Run notebook 01 to reduce your data
2. Open notebook 07 to visualize 2D spectra
3. For wavelength calibration, use PypeIt or the OLD notebooks

### For Production Analysis:

1. Install PypeIt: `pip install pypeit`
2. Use PypeIt's complete spectroscopic reduction pipeline
3. Or adapt algorithms from OLD notebooks to your specific needs

### For Development:

1. Help create `samos.spectroscopy.extraction` module
2. Help create `samos.spectroscopy.wavelength_cal` module
3. Test and refine the conceptual workflows in notebooks 04-06

---

## Technical Notes

### Wavelength Solution Format

Polynomial coefficients are saved as text files:
```
# wl_poly_coefficients_010.txt
1.234e-05   # degree 3
-2.345e-02  # degree 2
5.678e+00   # degree 1
5.432e+03   # degree 0 (constant)
```

Apply with:
```python
coeffs = np.loadtxt('wl_poly_coefficients_010.txt')
wavelength_poly = np.poly1d(coeffs)
wavelengths = wavelength_poly(pixel_array)
```

### Spectrum1D Output

All wavelength-calibrated spectra use `specutils.Spectrum1D`:
- Proper units (`u.AA`, `u.adu`)
- WCS headers included
- Compatible with jdaviz, astropy, and other tools
- Can be loaded with `Spectrum1D.read(filename)`

---

## Conclusion

✅ **All requested notebooks (04-06) have been updated**

The notebooks now provide:
- Working conceptual frameworks
- Integration with new `samos` package
- Clear documentation and references
- Production alternatives (PypeIt)
- Backup of original algorithms

Users can:
- Understand the spectroscopic reduction workflow
- Run simplified versions of each step
- Use OLD notebooks for production work
- Use PypeIt for professional results
- Wait for full module implementation

---

**Questions or issues?** Check the README or contact the development team.

**Ready to continue?** Run notebook 07 to visualize your spectra!

---

Last updated: December 15, 2025
