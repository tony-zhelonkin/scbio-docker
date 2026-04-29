# Python Virtual Environment Guide for DC Dictionary Project

## Overview
The `scdock-r-dev` Docker image uses a **layered virtual environment** approach to manage Python packages efficiently. Many usefull packages are already pre-installed. 

## Pre-installed Python Environments

### 1. Base Environment (`/opt/venvs/base`)
**Already includes** these core packages:
- `scanpy==1.10.4` - Single-cell analysis
- `scvi-tools==1.2.0` - Variational inference tools
- `cellrank==2.0.6` - Cell fate mapping
- `scvelo==0.3.2` - RNA velocity
- `muon==0.1.7` - Multi-modal data
- `harmony-pytorch` - Batch correction
- `scrublet==0.2.3` - Doublet detection
- `scib==1.1.5` - Benchmarking
- `scirpy==0.19.0` - TCR/BCR analysis
- Plus many more (numpy, pandas, matplotlib, seaborn, etc.)

**This is the default environment** - no installation needed!

### 2. Specialized Layered Environments
These are created on-demand when you first use them:

#### ATAC Environment
```bash
usepy atac  # Creates and activates ATAC environment
```
Includes: `snapatac2==2.7.1`, `episcanpy==0.4.0`

#### Spatial Environment
```bash
usepy squid  # Creates and activates spatial environment
```
Includes: `squidpy==1.6.5`, `spatialdata==0.4.0`, and spatial dependencies

#### Cell Communication Environment
```bash
usepy comms  # Creates and activates communication environment
```
Includes: Cell-cell communication and GRN tools

## How to Use

### Switching Between Environments
```bash
# Default (already active when container starts)
usepy base

# For ATAC-seq analysis
usepy atac

# For spatial transcriptomics
usepy squid

# For cell communication
usepy comms

# Check current environment
which python
python -V
```

### Installing Additional Packages

#### Option 1: Install in Current Environment (if you have permissions)
```bash
pip install --user package-name
```

#### Option 2: Create Your Own Virtual Environment
If you need packages not in the pre-installed environments:

```bash
# Create a new venv that inherits from base
python -m venv --system-site-packages ~/my-project-venv

# Activate it
source ~/my-project-venv/bin/activate

# Install your specific packages
pip install TEtranscripts==2.2.3
pip install any-other-package
```

#### Option 3: Use Conda (if needed)

I don\`t like conda but it is possible to install it in the same manner 
```bash
# Install miniconda in your home directory
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3
~/miniconda3/bin/conda init bash
source ~/.bashrc
```

## Package Conflicts Resolution

### If you see "Permission denied" errors:
This happens because `/opt/venvs/base` is owned by UID 1000, but VS Code remaps you to your actual UID. Solutions:
1. Use `pip install --user` to install in your home directory
2. Create your own venv (see Option 2 above)
3. Use specialized environments (`usepy atac/squid/comms`)

### If packages are already installed but wrong version:
The base environment has pinned versions. To override:
1. Create your own venv with `--system-site-packages`
2. Install the specific version you need

## Common Workflows

### Single-cell RNA-seq Analysis
```bash
# Already in base environment - just start working!
python
>>> import scanpy as sc
>>> import scvi
>>> adata = sc.read_h5ad("your_data.h5ad")
```

### ATAC-seq Analysis
```bash
usepy atac
python
>>> import snapatac2 as snap
>>> import episcanpy as epi
```

### Spatial Transcriptomics
```bash
usepy squid
python
>>> import squidpy as sq
>>> import spatialdata as sd
```

### Multi-modal Analysis
```bash
# Base environment has muon
python
>>> import muon as mu
>>> mdata = mu.read("multimodal.h5mu")
```

## Checking Installed Packages
```bash
# List all packages in current environment
pip list

# Check if a specific package is installed
python -c "import scanpy; print(scanpy.__version__)"

# See which venv is active
echo $VIRTUAL_ENV
```

## Troubleshooting

### "Module not found" error
1. Check you're in the right environment: `which python`
2. Check if package is installed: `pip list | grep package-name`
3. If not installed, see installation options above

### Permission errors during pip install
- Use `pip install --user`
- Or create your own venv (recommended for project-specific packages)

### Package version conflicts
- Create a project-specific venv to isolate dependencies
- Use `pip install --force-reinstall` if needed (in your own venv)

## Notes
- The Docker image is optimized for size - only essential packages are pre-installed
- Specialized packages are in layered venvs to avoid conflicts
- You can always create your own venv for full control
- The base environment covers 95% of single-cell analysis needs