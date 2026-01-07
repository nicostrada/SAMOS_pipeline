# SAMOS Pipeline - Documentation Index

**Last Updated:** December 16, 2025

This document provides a complete index of all current documentation and where to find archived/historical documents.

---

## Quick Navigation

**Start here:**
- 📍 **[PROJECT_STATUS.md](PROJECT_STATUS.md)** - Complete current status (RECOMMENDED)
- 🚀 **[QUICKSTART.md](QUICKSTART.md)** - Get started in 5 minutes

**For specific needs:**
- 📖 **[README.md](README.md)** - General pipeline information
- 💾 **[INSTALL.md](INSTALL.md)** - Detailed installation instructions
- 📊 **[notebooks/spectroscopy/README.md](notebooks/spectroscopy/README.md)** - Complete workflow guide

---

## Current Documentation (Use These)

### Pipeline Root Documentation

Located in: `/Users/nestrada/Documents/SAMOS/Pipeline/`

| Document | Purpose | When to Use |
|----------|---------|-------------|
| **PROJECT_STATUS.md** | Complete current status | Want comprehensive overview |
| **README.md** | General pipeline info | First-time orientation |
| **INSTALL.md** | Installation guide | Setting up environment |
| **QUICKSTART.md** | 5-minute quick start | Want to start immediately |
| **FOLDER_CLEANUP_SUMMARY.md** | Recent folder cleanup | Understanding organization |
| **DOCUMENTATION_INDEX.md** | This file | Finding documentation |

### Spectroscopy Workflow Documentation

Located in: `notebooks/spectroscopy/`

| Document | Purpose | When to Use |
|----------|---------|-------------|
| **README.md** | Complete workflow guide | Running full workflow |
| **QUICKSTART.md** | Quick reference | Common tasks, troubleshooting |
| **WORKFLOW_COMPLETE.md** | Detailed completion status | Understanding what's ready |

### Archive Documentation

Located in: `_archived_folders/` and `notebooks/spectroscopy/archive/`

| Document | Purpose | When to Use |
|----------|---------|-------------|
| **_archived_folders/ARCHIVED_FOLDERS_README.md** | Explains archived folders | Understanding what was archived |
| **notebooks/spectroscopy/archive/ARCHIVE_SUMMARY.md** | Explains archived notebooks | Finding old algorithms |

---

## Documentation by Purpose

### 🎯 "I want to get started quickly"
→ **[QUICKSTART.md](QUICKSTART.md)** (5 minutes)

### 📚 "I want complete information"
→ **[PROJECT_STATUS.md](PROJECT_STATUS.md)** (comprehensive)

### 🔧 "I need to install the pipeline"
→ **[INSTALL.md](INSTALL.md)** (step-by-step)

### 📓 "I want to run the workflow"
→ **[notebooks/spectroscopy/README.md](notebooks/spectroscopy/README.md)** (complete guide)

### 🐛 "I have a problem"
→ **[QUICKSTART.md](QUICKSTART.md)** or **[notebooks/spectroscopy/QUICKSTART.md](notebooks/spectroscopy/QUICKSTART.md)** (troubleshooting sections)

### 🔍 "I want to understand the code"
→ Module docstrings: `help(module.function)` in Python

### 📜 "I want to see old algorithms"
→ **[notebooks/spectroscopy/archive/old_notebooks/](notebooks/spectroscopy/archive/old_notebooks/)** (OLD notebooks)

### 🏛️ "I want to understand the history"
→ **[_archived_folders/MIGRATION_SUMMARY.md](_archived_folders/MIGRATION_SUMMARY.md)** (historical)

---

## Documentation Structure

```
Pipeline/
│
├── 📄 PROJECT_STATUS.md              ← START HERE (comprehensive)
├── 📄 README.md                      General info
├── 📄 INSTALL.md                     Installation
├── 📄 QUICKSTART.md                  Quick start
├── 📄 FOLDER_CLEANUP_SUMMARY.md      Cleanup documentation
├── 📄 DOCUMENTATION_INDEX.md         ← This file
│
├── notebooks/spectroscopy/
│   ├── 📄 README.md                  Complete workflow
│   ├── 📄 QUICKSTART.md              Quick reference
│   ├── 📄 WORKFLOW_COMPLETE.md       Detailed status
│   │
│   └── archive/
│       ├── 📄 ARCHIVE_SUMMARY.md     Archive explanation
│       ├── 📁 old_notebooks/         9 archived notebooks
│       └── 📁 old_documentation/     8 archived docs
│
└── _archived_folders/
    ├── 📄 ARCHIVED_FOLDERS_README.md Archive explanation
    ├── 📄 MIGRATION_SUMMARY.md       Historical migration
    ├── 📄 WORKFLOW_CLEANUP_COMPLETE.md Historical status
    └── 📁 [6 archived folders]
```

---

## Current vs. Historical Documentation

### ✅ Current (Use These)

**These documents reflect the current state (Dec 16, 2025):**

1. **PROJECT_STATUS.md** - Complete current status
2. **README.md** - General pipeline information
3. **INSTALL.md** - Installation instructions
4. **QUICKSTART.md** (both locations) - Quick references
5. **FOLDER_CLEANUP_SUMMARY.md** - Recent cleanup
6. **notebooks/spectroscopy/README.md** - Workflow guide
7. **notebooks/spectroscopy/WORKFLOW_COMPLETE.md** - Detailed status

### 📦 Historical (Reference Only)

**These documents are archived for historical reference:**

1. **MIGRATION_SUMMARY.md** - Dec 15 initial migration
2. **WORKFLOW_CLEANUP_COMPLETE.md** - Interim status (2 notebooks)
3. **notebooks/spectroscopy/archive/old_documentation/** - Transition docs
   - Old README.md
   - NOTEBOOK_UPDATE_SUMMARY.md
   - PYPEIT_WORKFLOW_DESIGN.md
   - etc.

---

## What Each Document Contains

### PROJECT_STATUS.md (RECOMMENDED)

**Most comprehensive document. Contains:**
- Complete current structure
- All 4 notebooks described
- All modules documented
- Installation instructions
- Python API examples
- Recent changes (Dec 15-16)
- Dependencies
- Known issues
- Future development
- Quick start guide

**When to use:** Want complete overview of project

### README.md (Pipeline root)

**General pipeline information. Contains:**
- Features and capabilities
- Quick installation
- Basic usage examples
- Project structure
- Configuration basics
- Future PypeIt notes

**When to use:** First-time orientation

### INSTALL.md

**Detailed installation guide. Contains:**
- Prerequisites
- Step-by-step installation
- Conda environment setup
- Jupyter configuration
- Troubleshooting common issues
- Verification scripts

**When to use:** Setting up environment for first time

### QUICKSTART.md

**Get running in 5 minutes. Contains:**
- Minimal installation steps
- Quick workflow overview
- Common tasks with code
- Troubleshooting
- Pro tips

**When to use:** Want to start immediately

### notebooks/spectroscopy/README.md

**Complete workflow guide. Contains:**
- All 4 notebooks described in detail
- Module usage examples
- Data organization
- Troubleshooting
- What's working now

**When to use:** Running the spectroscopy workflow

### notebooks/spectroscopy/WORKFLOW_COMPLETE.md

**Detailed completion status. Contains:**
- Complete work summary
- All phases documented
- Module architecture
- Recent fixes
- Python API examples
- Archive information

**When to use:** Understanding what's been completed

### FOLDER_CLEANUP_SUMMARY.md

**Recent cleanup documentation. Contains:**
- What folders were archived
- Why they were archived
- New structure benefits
- How to restore if needed

**When to use:** Understanding folder organization

---

## Documentation Philosophy

### Consolidation Principle

We've consolidated documentation to avoid duplication:

**Before:** 12+ scattered documents with overlapping information

**After:** 7 current documents + organized archives

**Benefits:**
- ✅ Clear which docs are current
- ✅ No conflicting information
- ✅ Easy to find what you need
- ✅ Historical docs preserved for reference

### Documentation Levels

1. **Quick (5 min)** → QUICKSTART.md
2. **Complete (15 min)** → PROJECT_STATUS.md
3. **Specific (as needed)** → Individual docs (INSTALL, README, etc.)
4. **Historical (reference)** → Archived docs

---

## Maintenance Guidelines

### When to Update Each Document

**PROJECT_STATUS.md**
- Update when: Major features added, status changes
- Frequency: After significant milestones

**README.md**
- Update when: General project info changes
- Frequency: Infrequently (stable overview)

**INSTALL.md**
- Update when: Installation process changes
- Frequency: When dependencies or steps change

**QUICKSTART.md**
- Update when: Basic workflow changes
- Frequency: When simplifications found

**notebooks/spectroscopy/README.md**
- Update when: Workflow steps change
- Frequency: When notebooks modified

**WORKFLOW_COMPLETE.md**
- Update when: Major work completed
- Frequency: After completion of major phases

---

## Archive Policy

### What Gets Archived

Documents get archived when:
1. Superseded by newer documentation
2. Describe historical states/transitions
3. No longer reflect current workflow
4. Useful for reference but not active use

### What Never Gets Archived

These always stay current:
- PROJECT_STATUS.md
- INSTALL.md
- Current README files
- Current QUICKSTART files

### Where Archives Live

- **Folder archives:** `_archived_folders/`
- **Notebook archives:** `notebooks/spectroscopy/archive/`
- **Documentation archives:** In respective archive folders

---

## Finding Specific Information

### Installation
→ [INSTALL.md](INSTALL.md)

### Quick Start
→ [QUICKSTART.md](QUICKSTART.md)

### Complete Status
→ [PROJECT_STATUS.md](PROJECT_STATUS.md)

### Workflow Guide
→ [notebooks/spectroscopy/README.md](notebooks/spectroscopy/README.md)

### Module Documentation
→ `help(module.function)` in Python

### Old Algorithms
→ [notebooks/spectroscopy/archive/old_notebooks/](notebooks/spectroscopy/archive/old_notebooks/)

### Migration History
→ [_archived_folders/MIGRATION_SUMMARY.md](_archived_folders/MIGRATION_SUMMARY.md)

### Troubleshooting
→ QUICKSTART.md (both locations) → Troubleshooting section

---

## Summary

**Current active documentation:** 7 files
- PROJECT_STATUS.md (comprehensive)
- README.md (general)
- INSTALL.md (installation)
- QUICKSTART.md (quick start)
- FOLDER_CLEANUP_SUMMARY.md (cleanup)
- notebooks/spectroscopy/README.md (workflow)
- notebooks/spectroscopy/WORKFLOW_COMPLETE.md (detailed)

**Archived documentation:** 10+ files
- Historical migration documents
- Transition status reports
- Old notebook documentation
- Design documents

**Organization principle:** Current vs. archived, with clear navigation

**Recommended starting point:** **PROJECT_STATUS.md** for complete overview

---

Last updated: December 16, 2025
