#!/usr/bin/env python3
"""
Migration script to reorganize SAMOS Pipeline structure.

This script creates the new directory structure and moves files to their
appropriate locations while preserving the original files in an archive.

Usage:
    python migrate_structure.py [--dry-run] [--backup]

Options:
    --dry-run   Show what would be done without making changes
    --backup    Create backup of entire Pipeline folder before migration
"""

import os
import shutil
import argparse
from pathlib import Path
from datetime import datetime


class PipelineMigration:
    def __init__(self, base_path, dry_run=False):
        self.base_path = Path(base_path)
        self.dry_run = dry_run
        self.actions = []

    def log_action(self, action_type, source, dest=None):
        """Log migration actions."""
        if dest:
            self.actions.append(f"{action_type}: {source} -> {dest}")
        else:
            self.actions.append(f"{action_type}: {source}")
        print(f"  {action_type}: {source}" + (f" -> {dest}" if dest else ""))

    def create_directory(self, path):
        """Create directory if it doesn't exist."""
        full_path = self.base_path / path
        if not self.dry_run:
            full_path.mkdir(parents=True, exist_ok=True)
        self.log_action("CREATE DIR", path)

    def move_file(self, source, dest):
        """Move file from source to destination."""
        src_path = self.base_path / source
        dst_path = self.base_path / dest

        if not src_path.exists():
            print(f"  WARNING: Source not found: {source}")
            return

        if not self.dry_run:
            dst_path.parent.mkdir(parents=True, exist_ok=True)
            if dst_path.exists():
                print(f"  WARNING: Destination exists, skipping: {dest}")
                return
            shutil.move(str(src_path), str(dst_path))

        self.log_action("MOVE", source, dest)

    def copy_file(self, source, dest):
        """Copy file from source to destination."""
        src_path = self.base_path / source
        dst_path = self.base_path / dest

        if not src_path.exists():
            print(f"  WARNING: Source not found: {source}")
            return

        if not self.dry_run:
            dst_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(str(src_path), str(dst_path))

        self.log_action("COPY", source, dest)

    def move_directory(self, source, dest):
        """Move entire directory."""
        src_path = self.base_path / source
        dst_path = self.base_path / dest

        if not src_path.exists():
            print(f"  WARNING: Source not found: {source}")
            return

        if not self.dry_run:
            dst_path.parent.mkdir(parents=True, exist_ok=True)
            if dst_path.exists():
                print(f"  WARNING: Destination exists, skipping: {dest}")
                return
            shutil.move(str(src_path), str(dst_path))

        self.log_action("MOVE DIR", source, dest)

    def create_structure(self):
        """Create the new directory structure."""
        print("\n=== Creating Directory Structure ===")

        directories = [
            # Main package
            "samos",
            "samos/core",
            "samos/spectroscopy",
            "samos/imaging",
            "samos/utils",
            "samos/pipeline",

            # Notebooks
            "notebooks/spectroscopy",
            "notebooks/imaging",

            # Scripts
            "scripts",

            # Configs
            "configs",
            "configs/instrument_profiles",

            # Calibration data
            "calibration_data",
            "calibration_data/line_lists",
            "calibration_data/sensitivity",
            "calibration_data/standards",

            # Tools
            "tools",

            # Tests
            "tests",

            # Documentation
            "docs",
            "docs/user_guide",
            "docs/developer_guide",
            "docs/presentations",
            "docs/reference_plots",

            # Examples
            "examples",
            "examples/spectroscopy_example",

            # Archive
            "archive",
            "archive/development",
            "archive/legacy_notebooks",
            "archive/deprecated",
        ]

        for directory in directories:
            self.create_directory(directory)

    def migrate_files(self):
        """Migrate files to new locations."""
        print("\n=== Migrating Files ===")

        # Move main notebooks to notebooks/spectroscopy/
        print("\nMigrating main notebooks...")
        notebook_mapping = {
            "1.SAMOS_reduction_FullFrame_V1.ipynb": "notebooks/spectroscopy/01_initial_inspection.ipynb",
            "2.SAMOS_reduction_splittraces_V2.ipynb": "notebooks/spectroscopy/04_trace_extraction.ipynb",
            "3.SAMOS_reduction_HgArNe.ipynb": "notebooks/spectroscopy/05_wavelength_calibration.ipynb",
            "4.SAMOS_reduction_all_wlcal.ipynb": "notebooks/spectroscopy/06_apply_calibration.ipynb",
            "5.SAMOS_reduction_VIZTOOLS.ipynb": "notebooks/spectroscopy/07_visualization.ipynb",
        }

        for source, dest in notebook_mapping.items():
            self.move_file(source, dest)

        # Move Class_SAMOS.py (will be split later)
        print("\nMigrating Python source files...")
        if (self.base_path / "Class_SAMOS.py").exists():
            self.copy_file("Class_SAMOS.py", "samos/core/samos_legacy.py")
            self.move_file("Class_SAMOS.py", "archive/deprecated/Class_SAMOS.py")

        # Move special notebooks to examples
        print("\nMigrating example notebooks...")
        if (self.base_path / "SAMOS_SISI_V1.ipynb").exists():
            self.move_file("SAMOS_SISI_V1.ipynb", "examples/spectroscopy_example/samos_sisi_v1.ipynb")

        if (self.base_path / "Mosviz/MosvizExample.ipynb").exists():
            self.move_file("Mosviz/MosvizExample.ipynb", "examples/spectroscopy_example/mosviz_example.ipynb")

        # Move calibration data
        print("\nMigrating calibration data...")
        if (self.base_path / "Calibration_GoodmanLines").exists():
            self.move_directory("Calibration_GoodmanLines", "calibration_data/line_lists/goodman_lines")

        # Move tools
        print("\nMigrating external tools...")
        if (self.base_path / "PyHammer").exists():
            self.move_directory("PyHammer", "tools/pyhammer")

        # Move documentation
        print("\nMigrating documentation...")
        if (self.base_path / "Pipeline Presentation1.pdf").exists():
            self.move_file("Pipeline Presentation1.pdf", "docs/presentations/pipeline_presentation1.pdf")

        if (self.base_path / "Pipeline Presentation1.pptx").exists():
            self.move_file("Pipeline Presentation1.pptx", "docs/presentations/pipeline_presentation1.pptx")

        if (self.base_path / "Notes.md").exists():
            self.move_file("Notes.md", "docs/notes.md")

        # Move reference plots
        print("\nMigrating reference plots...")
        plot_files = ["plot_1200m5_light_theme_03.pdf", "plot_1200m6_light_theme_04.pdf"]
        for plot_file in plot_files:
            if (self.base_path / plot_file).exists():
                self.move_file(plot_file, f"docs/reference_plots/{plot_file}")

        # Archive old development materials
        print("\nArchiving development materials...")
        if (self.base_path / "earlytests").exists():
            self.move_directory("earlytests", "archive/development/earlytests")

        if (self.base_path / "current").exists():
            self.move_directory("current", "archive/development/current")

        if (self.base_path / "Nicolas_01").exists():
            self.move_directory("Nicolas_01", "archive/development/nicolas_01")

        # Clean up duplicate V2 notebook
        if (self.base_path / "4.SAMOS_reduction_all_wlcal_V2.ipynb").exists():
            self.move_file("4.SAMOS_reduction_all_wlcal_V2.ipynb",
                          "archive/legacy_notebooks/4.SAMOS_reduction_all_wlcal_V2.ipynb")

    def create_init_files(self):
        """Create __init__.py files for Python package."""
        print("\n=== Creating __init__.py Files ===")

        init_files = [
            "samos/__init__.py",
            "samos/core/__init__.py",
            "samos/spectroscopy/__init__.py",
            "samos/imaging/__init__.py",
            "samos/utils/__init__.py",
            "samos/pipeline/__init__.py",
        ]

        for init_file in init_files:
            init_path = self.base_path / init_file
            if not self.dry_run and not init_path.exists():
                init_path.parent.mkdir(parents=True, exist_ok=True)
                init_path.write_text('"""SAMOS Pipeline package."""\n')
            self.log_action("CREATE", init_file)

    def run(self):
        """Execute the migration."""
        print(f"\n{'='*60}")
        print(f"SAMOS Pipeline Migration")
        print(f"Base path: {self.base_path}")
        print(f"Mode: {'DRY RUN' if self.dry_run else 'LIVE'}")
        print(f"{'='*60}")

        self.create_structure()
        self.migrate_files()
        self.create_init_files()

        print(f"\n{'='*60}")
        print(f"Migration {'would be' if self.dry_run else 'is'} complete!")
        print(f"Total actions: {len(self.actions)}")
        print(f"{'='*60}")

        if self.dry_run:
            print("\nThis was a DRY RUN. No changes were made.")
            print("Run without --dry-run to execute the migration.")


def create_backup(base_path):
    """Create backup of Pipeline folder."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_name = f"Pipeline_backup_{timestamp}"
    backup_path = base_path.parent / backup_name

    print(f"\n=== Creating Backup ===")
    print(f"Backing up to: {backup_path}")

    shutil.copytree(base_path, backup_path, symlinks=True)
    print(f"Backup created successfully!")

    return backup_path


def main():
    parser = argparse.ArgumentParser(description="Migrate SAMOS Pipeline structure")
    parser.add_argument("--dry-run", action="store_true",
                       help="Show what would be done without making changes")
    parser.add_argument("--backup", action="store_true",
                       help="Create backup before migration")
    parser.add_argument("--path", default="/Users/nestrada/Documents/SAMOS/Pipeline",
                       help="Path to Pipeline directory")

    args = parser.parse_args()

    base_path = Path(args.path)

    if not base_path.exists():
        print(f"Error: Pipeline directory not found: {base_path}")
        return 1

    # Create backup if requested
    if args.backup and not args.dry_run:
        backup_path = create_backup(base_path)
        print(f"\nBackup created at: {backup_path}")

    # Run migration
    migration = PipelineMigration(base_path, dry_run=args.dry_run)
    migration.run()

    return 0


if __name__ == "__main__":
    exit(main())
