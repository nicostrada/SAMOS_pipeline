# Archived Folders - Pipeline Cleanup

**Archive Date:** December 16, 2025

**Reason:** Pipeline folder organization and cleanup to keep only actively used folders in the main directory.

---

## What Was Archived

The following folders and **historical documentation** were moved to `_archived_folders/` because they are **not actively used** in the current SAMOS spectroscopy workflow:

### Archived Documentation (2 files)

**Historical transition documents:**

1. **MIGRATION_SUMMARY.md**
   - Documents Dec 15 migration from Class_SAMOS to modular structure
   - Historical record of initial restructuring
   - Superseded by PROJECT_STATUS.md

2. **WORKFLOW_CLEANUP_COMPLETE.md**
   - Documents interim status when only 2 notebooks were ready
   - Now superseded by WORKFLOW_COMPLETE.md (4 notebooks ready)

**Why archived:**
- Historical documentation of migration process
- Useful reference for understanding past changes
- Superseded by current consolidated documentation (PROJECT_STATUS.md)

---

### Archived Folders (6 folders)

The following folders were moved to `_archived_folders/`:

### 1. calibration_data/
**Contents:**
- `line_lists/goodman_lines/` - Line list PDFs for Goodman spectrograph (different instrument)
- `sensitivity/` - Empty
- `standards/` - Empty

**Why archived:**
- Contains calibration data for Goodman instrument, not SAMOS
- Not referenced in current samos package code
- Not used by notebooks
- Line lists for SAMOS are embedded in `samos/spectroscopy/pypeit_wrapper.py` (HgArNe)

**Note:** If SAMOS-specific calibration data is needed in the future, create a new `calibration_data/` folder with SAMOS-specific files.

---

### 2. archive/
**Contents:**
- `deprecated/` - Old deprecated files
- `development/` - Development files
- `legacy_notebooks/` - Legacy notebooks

**Why archived:**
- This was an archive from a PREVIOUS migration
- Creates confusing nested archives
- Already have proper archive in `notebooks/spectroscopy/archive/`
- All relevant old notebooks properly archived in the spectroscopy workflow

**Note:** This is an "archive of an archive" - contains files from earlier development phases.

---

### 3. examples/
**Contents:**
- `spectroscopy_example/` - Single example folder with minimal content

**Why archived:**
- Not referenced in documentation
- Not used in current workflow
- Current notebooks (01-04) serve as complete working examples
- README.md already provides usage examples

**Note:** The production notebooks in `notebooks/spectroscopy/` are better examples.

---

### 4. Mosviz/
**Contents:**
- Only `.ipynb_checkpoints/` (Jupyter temp files)

**Why archived:**
- Essentially empty folder
- Mosviz functionality now integrated in notebook 04_visualization.ipynb
- No standalone Mosviz scripts needed

---

### 5. tests/
**Contents:**
- Empty folder

**Why archived:**
- No test files present
- When tests are added in the future, can recreate this folder
- Testing infrastructure not yet implemented

**Future:** When adding tests, recreate `tests/` folder with proper pytest structure.

---

### 6. tools/
**Contents:**
- `pyhammer/` - PyHammer spectral typing tool (separate project)

**Why archived:**
- PyHammer is an external tool, not part of SAMOS pipeline
- Not integrated with current workflow
- Not documented in pipeline README
- If needed, can be installed separately as standalone tool

**Note:** PyHammer is for stellar spectral classification, not reduction.

---

## Current Clean Pipeline Structure

After cleanup, the Pipeline folder contains only actively used directories:

```
Pipeline/
├── samos/                      # ✅ Main Python package (ESSENTIAL)
│   ├── core/                   # Core processing modules
│   ├── spectroscopy/           # Spectroscopy modules
│   ├── imaging/                # Imaging modules
│   ├── utils/                  # Utility functions
│   └── pipeline/               # Automated pipeline
│
├── notebooks/                  # ✅ Interactive Jupyter notebooks (ACTIVE)
│   └── spectroscopy/          # Current 4-notebook workflow
│       ├── 01_initial_inspection.ipynb
│       ├── 02_visual_qa.ipynb
│       ├── 03_spectroscopy_pypeit.ipynb
│       ├── 04_visualization.ipynb
│       └── archive/           # Properly organized archive
│
├── scripts/                    # ✅ Automated pipeline scripts (ACTIVE)
│   └── run_spectroscopy_pipeline.py
│
├── configs/                    # ✅ Configuration files (ACTIVE)
│   └── spectroscopy_default.yaml
│
├── docs/                       # ✅ Documentation (ACTIVE)
│   ├── user_guide/
│   └── developer_guide/
│
├── _archived_folders/          # 📦 This archive
│
├── README.md                   # Main documentation
├── QUICKSTART.md              # Quick reference
├── requirements.txt           # Dependencies
└── setup.py                   # Package installation
```

---

## What Changed in setup.py

The `setup.py` file referenced `calibration_data` in `package_data`:

```python
package_data={
    "samos": [
        "calibration_data/**/*",  # This path no longer exists
        "configs/*.yaml",
    ],
},
```

**Action needed:** This should be updated or removed since `calibration_data/` is archived.

**Alternative:** If SAMOS-specific calibration data is needed, create `samos/data/` inside the package:
```
samos/
└── data/
    └── line_lists/
        └── hgarne.dat  # SAMOS line lists
```

---

## Restoring Archived Folders

If you need any of these folders back:

```bash
cd /Users/nestrada/Documents/SAMOS/Pipeline

# Restore a specific folder
mv _archived_folders/calibration_data ./

# Restore all
mv _archived_folders/* ./
rmdir _archived_folders
```

---

## Benefits of This Cleanup

### ✅ Clearer Structure
- Main directory shows only active components
- Easier to understand what's essential vs optional
- New users won't be confused by unused folders

### ✅ Reduced Confusion
- No nested archives (archive of archive)
- No empty placeholder folders
- No external tools mixed with pipeline code

### ✅ Better Maintenance
- Clear separation of SAMOS code vs external tools
- Easier to identify what needs updating
- Less clutter when navigating the project

### ✅ Proper Documentation
- This README explains what was archived and why
- Easy to restore if needed
- Clear record of organizational decisions

---

## Summary

| Item | Size | Reason | Can Delete? |
|------|------|--------|-------------|
| **Documentation** | | | |
| MIGRATION_SUMMARY.md | <50 KB | Historical, superseded | Keep for history |
| WORKFLOW_CLEANUP_COMPLETE.md | <50 KB | Interim status, superseded | Keep for history |
| **Folders** | | | |
| calibration_data | ~500 KB | Goodman data, not SAMOS | Yes |
| archive | Small | Archive of previous archive | Yes |
| examples | Small | Superseded by notebooks | Yes |
| Mosviz | Minimal | Empty, functionality in nb 04 | Yes |
| tests | Empty | No tests yet | Yes |
| tools | ~10 MB | External tool (PyHammer) | Keep if using PyHammer |

**Recommendation:**
- **Documentation:** Keep for historical reference
- **Folders:** Can remain archived. If specific functionality is needed later (e.g., PyHammer), restore only that folder.

---

## Related Documentation

### Current Documentation
- **PROJECT_STATUS.md** - Consolidated current status: `../PROJECT_STATUS.md`
- **Main Pipeline README** - General info: `../README.md`
- **Spectroscopy Workflow** - Complete workflow: `../notebooks/spectroscopy/README.md`
- **Workflow Complete** - Detailed status: `../notebooks/spectroscopy/WORKFLOW_COMPLETE.md`

### Historical Documentation (Archived)
- **MIGRATION_SUMMARY.md** - Dec 15 migration: `./MIGRATION_SUMMARY.md`
- **WORKFLOW_CLEANUP_COMPLETE.md** - Interim status: `./WORKFLOW_CLEANUP_COMPLETE.md`
- **Spectroscopy Archive** - Old notebooks: `../notebooks/spectroscopy/archive/ARCHIVE_SUMMARY.md`

---

Last updated: December 16, 2025
