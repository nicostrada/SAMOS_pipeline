# SAMOS Pipeline Quick Start

Get up and running in 5 minutes!

## 1. Install the Pipeline

```bash
# Activate your conda environment
conda activate samos

# Navigate to Pipeline directory
cd /Users/nestrada/Documents/SAMOS/Pipeline

# Install in development mode
pip install -e .

# Verify installation
python -c "import samos; print('Success! Version:', samos.__version__)"
```

## 2. Choose Your Workflow

### Option A: Interactive Mode (First Time Users)

**Best for**: Learning the pipeline, exploring data, tuning parameters

```bash
# Launch Jupyter
cd notebooks/spectroscopy
jupyter lab
```

Work through notebooks in order:
1. `01_initial_inspection.ipynb` - Look at your raw data
2. `02_calibration_frames.ipynb` - Create master bias/flat (future)
3. `03_trace_identification.ipynb` - Find spectral traces (future)
4. Continue through the workflow...

### Option B: Automated Mode (Production)

**Best for**: Batch processing, reproducible reductions

```bash
# 1. Copy configuration template
cp configs/spectroscopy_default.yaml my_config.yaml

# 2. Edit my_config.yaml with your paths:
#    - data.raw_directory: where your FITS files are
#    - output.base_directory: where to save results
#    - File patterns (bias_pattern, flat_pattern, etc.)

# 3. Run pipeline
python scripts/run_spectroscopy_pipeline.py --config my_config.yaml
```

## 3. Quick Test with Example Data

Try the pipeline with your existing data:

```python
# test_pipeline.py
from samos.core import mosaic, cosmic_rays
from samos.utils import display

# Read a SAMOS frame
hdu = mosaic.read_sami_mosaic('/Users/nestrada/Documents/SAMOS/_Run2/Data/SAMI/20241017/target.027.fits')

# Display it
vmin, vmax = display.get_percentile_limits(hdu.data, 1, 99)
display.display_image(hdu.data, zmin=vmin, zmax=vmax, title='Raw Science Frame')

# Remove cosmic rays
cleaned = cosmic_rays.remove_cosmic_rays(hdu.data, verbose=True)
display.display_image(cleaned, zmin=vmin, zmax=vmax, title='CR Cleaned')

print("✓ Pipeline working!")
```

Run with:
```bash
python test_pipeline.py
```

## 4. Directory Structure at a Glance

```
Pipeline/
├── samos/              # Python package (import from anywhere)
├── notebooks/          # Interactive analysis
├── scripts/            # Automated pipelines
├── configs/            # Configuration templates
└── calibration_data/   # Reference line lists, etc.

Your Data (outside Pipeline/):
~/SAMOS_Data/
├── raw/                # Your FITS files
├── reduced/            # Pipeline outputs
└── configs/            # Your configurations
```

## 5. Common Tasks

### Create Master Bias

```python
from samos.utils import io
from samos.core import calibration, mosaic

# Find all bias files
bias_files = io.find_files('/path/to/data', 'bias*.fits')

# Read and create mosaics
bias_mosaics = [mosaic.read_sami_mosaic(f).data for f in bias_files]

# Combine
master_bias = calibration.create_master_bias(bias_mosaics, method='median')

# Save
io.write_fits(master_bias, 'master_bias.fits')
```

### Detect Traces

```python
from samos.spectroscopy import trace_detection

# Detect slits from arc lamp
slits = trace_detection.detect_traces_from_arc(
    arc_image,
    arc_dark=arc_dark_image,  # Optional
    peak_threshold=4e5,
    verbose=True
)

print(f"Found {len(slits)} traces")
```

### Process a Science Frame

```python
from samos.core import mosaic, cosmic_rays, calibration

# Read
hdu = mosaic.read_sami_mosaic('science.fits')

# Remove cosmic rays
cleaned = cosmic_rays.remove_cosmic_rays(hdu.data)

# Apply calibrations
reduced = calibration.reduce_image(cleaned, master_bias, master_flat)

# Save
io.write_fits(reduced, 'science_reduced.fits')
```

## 6. Getting Help

- **Documentation**: Check `README.md` for detailed info
- **Examples**: Look in `examples/` directory
- **Config help**: See `configs/spectroscopy_default.yaml` with comments
- **Function docs**: Use `help(function_name)` in Python

```python
# Get help on any function
from samos.core import mosaic
help(mosaic.read_sami_mosaic)
```

## 7. Next Steps

After the quick start:

1. **Read the full README**: [README.md](README.md)
2. **Explore notebooks**: Work through the spectroscopy workflow
3. **Configure for your data**: Edit `configs/spectroscopy_default.yaml`
4. **Run automated pipeline**: Process your full dataset

## Troubleshooting

**Can't import samos?**
```bash
# Make sure you installed it
pip install -e .

# Check it's installed
pip list | grep samos
```

**Notebooks can't find samos?**
- Change kernel to "Python 3.13 (SAMOS Pipeline)"
- Or install kernel: `python -m ipykernel install --user --name samos`

**Module not found (cv2, lacosmic, etc.)?**
```bash
# Reinstall requirements
pip install -r requirements.txt
```

## Pro Tips

1. **Always use conda environment**: `conda activate samos` before working
2. **Keep data outside Pipeline/**: Don't clutter the code directory
3. **Use configuration files**: Easier to reproduce and share
4. **Start interactive**: Learn with notebooks before going automated
5. **Check examples**: Real working code in `examples/` directory

Happy reducing! 🔭
