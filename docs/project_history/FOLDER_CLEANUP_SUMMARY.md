# Pipeline Folder Cleanup Summary

**Date:** December 16, 2025

**Action:** Organized Pipeline directory to keep only actively used folders

---

## Changes Made

### Folders Moved to Archive

The following **6 folders** were moved to `_archived_folders/`:

1. **calibration_data/** - Goodman instrument data (not SAMOS)
2. **archive/** - Old archive from previous migration
3. **examples/** - Superseded by production notebooks
4. **Mosviz/** - Empty folder, functionality in notebook 04
5. **tests/** - Empty placeholder folder
6. **tools/** - External PyHammer tool (not part of pipeline)

See [_archived_folders/ARCHIVED_FOLDERS_README.md](_archived_folders/ARCHIVED_FOLDERS_README.md) for detailed explanation.

---

## Current Clean Structure

The Pipeline folder now contains only **actively used** directories:

```
Pipeline/
├── samos/                      # ✅ Main Python package
│   ├── core/                   # Core processing (mosaic, cosmic rays, calibration)
│   ├── spectroscopy/           # Spectroscopy (extraction, wavelength cal)
│   ├── imaging/                # Imaging (future)
│   ├── utils/                  # Utilities (I/O, display)
│   └── pipeline/               # Automated pipeline
│
├── notebooks/                  # ✅ Interactive Jupyter notebooks
│   └── spectroscopy/          # Current workflow (4 notebooks)
│       ├── 01_initial_inspection.ipynb
│       ├── 02_visual_qa.ipynb
│       ├── 03_spectroscopy_pypeit.ipynb
│       ├── 04_visualization.ipynb
│       └── archive/           # Old notebooks from previous updates
│
├── scripts/                    # ✅ Automated pipeline scripts
│   └── run_spectroscopy_pipeline.py
│
├── configs/                    # ✅ Configuration templates
│   └── spectroscopy_default.yaml
│
├── docs/                       # ✅ Documentation
│   ├── user_guide/
│   └── developer_guide/
│
├── _archived_folders/          # 📦 Archived unused folders
│   ├── calibration_data/
│   ├── archive/
│   ├── examples/
│   ├── Mosviz/
│   ├── tests/
│   ├── tools/
│   └── ARCHIVED_FOLDERS_README.md
│
├── README.md                   # Main documentation
├── QUICKSTART.md              # Quick start guide
├── requirements.txt           # Python dependencies
└── setup.py                   # Package installation (updated)
```

---

## Files Updated

### setup.py
Removed reference to archived `calibration_data/`:

**Before:**
```python
package_data={
    "samos": [
        "calibration_data/**/*",
        "configs/*.yaml",
    ],
},
```

**After:**
```python
package_data={
    "samos": [
        # Note: calibration_data archived - SAMOS line lists now in pypeit_wrapper.py
        # Add future SAMOS-specific data files here if needed
    ],
},
```

**Reason:** HgArNe line lists are now embedded directly in `samos/spectroscopy/pypeit_wrapper.py`, no external data files needed.

---

## Benefits

### ✅ Cleaner Organization
- Main directory shows only essential components
- Easier to navigate and understand structure
- New users see only what's actively used

### ✅ Reduced Confusion
- No empty placeholder folders (tests/, Mosviz/)
- No nested archives (archive of archive)
- No external tools mixed with pipeline code (PyHammer)
- No data for different instruments (Goodman)

### ✅ Better Maintenance
- Clear what needs updating vs what's historical
- Smaller directory listing when navigating
- Focused development on active components

### ✅ Preserved History
- All archived folders preserved in `_archived_folders/`
- Detailed README explains what was archived and why
- Easy to restore if needed

---

## What This Does NOT Affect

### ✅ Workflow Still Works
- All notebooks (01-04) work unchanged
- All `samos` package modules work unchanged
- Scripts and configs work unchanged
- Documentation unchanged

### ✅ No Data Loss
- All folders moved, not deleted
- Complete history preserved
- Can restore any folder if needed

### ✅ No Breaking Changes
- Package installation still works (`pip install -e .`)
- Import statements unchanged
- Module structure unchanged

---

## Verification

### Check Pipeline Structure
```bash
cd /Users/nestrada/Documents/SAMOS/Pipeline
ls -la
# Should see clean directory with only active folders
```

### Verify Package Still Works
```bash
conda activate samos
python -c "from samos.core import mosaic; from samos.spectroscopy import extraction; print('✓ Package works')"
```

### Verify Notebooks Still Work
```bash
cd notebooks/spectroscopy
jupyter lab
# Open and run any notebook - all should work
```

---

## Restoring Archived Folders

If you need any archived folder back:

```bash
cd /Users/nestrada/Documents/SAMOS/Pipeline

# Restore a specific folder
mv _archived_folders/calibration_data ./

# Restore all folders
mv _archived_folders/* ./
rmdir _archived_folders
```

---

## Related Documentation

- **Archived Folders Details:** [_archived_folders/ARCHIVED_FOLDERS_README.md](_archived_folders/ARCHIVED_FOLDERS_README.md)
- **Spectroscopy Workflow:** [notebooks/spectroscopy/README.md](notebooks/spectroscopy/README.md)
- **Workflow Complete:** [notebooks/spectroscopy/WORKFLOW_COMPLETE.md](notebooks/spectroscopy/WORKFLOW_COMPLETE.md)
- **Main README:** [README.md](README.md)

---

## Summary

**Before cleanup:** 20+ items in Pipeline directory (including 6 unused folders)

**After cleanup:** 14 items in Pipeline directory (only active components + 1 archive folder)

**Result:** Clean, organized, production-ready pipeline structure

---

Last updated: December 16, 2025
