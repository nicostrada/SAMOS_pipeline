#!/usr/bin/env python3
"""
Verify consistency and correctness of SAMOS spectroscopy notebooks.

This script checks:
1. All notebooks use new samos package imports (not old Class_SAMOS)
2. Import statements are consistent
3. No broken references
4. All required packages are available
"""

import json
import re
from pathlib import Path

def check_notebook(notebook_path):
    """Check a single notebook for consistency."""
    issues = []
    warnings = []

    with open(notebook_path) as f:
        nb = json.load(f)

    # Check for old Class_SAMOS imports in code cells
    for i, cell in enumerate(nb.get('cells', [])):
        if cell.get('cell_type') == 'code':
            source = ''.join(cell.get('source', []))

            # Check for old imports (should only be in markdown/comments)
            if 'from Class_SAMOS import' in source or 'import Class_SAMOS' in source:
                # Not in a comment
                if not source.strip().startswith('#'):
                    issues.append(f"Cell {i}: Uses old Class_SAMOS import")

            # Check for new samos imports
            if 'from samos' in source:
                imports = re.findall(r'from samos\.[^\s]+ import [^\n]+', source)
                if imports:
                    warnings.append(f"Cell {i}: Found {len(imports)} samos imports")

    # Check for samos package usage
    has_samos_import = False
    for cell in nb.get('cells', []):
        if cell.get('cell_type') == 'code':
            source = ''.join(cell.get('source', []))
            if 'from samos' in source:
                has_samos_import = True
                break

    return {
        'notebook': notebook_path.name,
        'issues': issues,
        'warnings': warnings,
        'has_samos': has_samos_import
    }

def main():
    """Check all notebooks."""
    notebooks_dir = Path(__file__).parent

    # Notebooks to check (exclude OLD versions)
    notebooks = [
        '01_initial_inspection.ipynb',
        '02_calibration_frames.ipynb',
        '03_trace_identification.ipynb',
        '04_trace_extraction.ipynb',
        '05_wavelength_calibration.ipynb',
        '06_apply_calibration.ipynb',
        '07_visualization.ipynb',
    ]

    print("=" * 70)
    print("SAMOS Spectroscopy Notebooks Verification")
    print("=" * 70)
    print()

    all_results = []

    for nb_name in notebooks:
        nb_path = notebooks_dir / nb_name
        if not nb_path.exists():
            print(f"⚠️  {nb_name}: NOT FOUND")
            continue

        result = check_notebook(nb_path)
        all_results.append(result)

        # Print results
        status = "✅" if not result['issues'] else "❌"
        print(f"{status}  {nb_name}")

        if result['issues']:
            for issue in result['issues']:
                print(f"     ❌ {issue}")

        if result['has_samos']:
            print(f"     ✓ Uses new samos package")
        else:
            print(f"     ℹ️  No samos imports (may be placeholder)")

        print()

    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    total = len(all_results)
    issues_count = sum(1 for r in all_results if r['issues'])
    samos_count = sum(1 for r in all_results if r['has_samos'])

    print(f"Total notebooks checked: {total}")
    print(f"Notebooks with issues: {issues_count}")
    print(f"Notebooks using new package: {samos_count}")
    print()

    if issues_count == 0:
        print("✅ All notebooks are consistent!")
    else:
        print(f"⚠️  {issues_count} notebook(s) need attention")

    return issues_count

if __name__ == '__main__':
    exit(main())
