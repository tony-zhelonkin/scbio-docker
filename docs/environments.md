# Environments

How the `scdock-r-dev` image organizes Python and R environments, and how to
move data between them. For where these fit in the overall image layering, see
[architecture.md](architecture.md).

## Python environments

Python is **3.11** (deadsnakes PPA, not the Ubuntu 22.04 default 3.10). The
image ships one fully-resolved **base venv** plus three **layered venvs** that
are created on demand.

### Base venv

`/opt/venvs/base` is the default environment and covers the bulk of single-cell
work: scanpy, anndata, scvi-tools, cellrank, scvelo, muon, MACS3, scIBD, plus
`numpy`/`pandas`/`matplotlib`/`seaborn`. `radian` (the R console) also lives
here. No activation is needed for base — it is active when the container starts.

### Switching environments

Two mechanisms are wired into the shell:

| Command | Effect |
|---------|--------|
| `usepy base\|squid\|atac\|comms` | Activate an environment in the current shell (creates it if missing) |
| `py-base`, `py-squid`, `py-atac`, `py-comms` | Run a one-off command in that environment without switching |

```bash
usepy atac                 # activate the ATAC env
py-squid python -c "import squidpy"   # one-off in the squid env
which python && python -V  # confirm the active interpreter
echo "$VIRTUAL_ENV"
```

### Layered venvs

The three specialized venvs inherit the base packages via
`--system-site-packages`, so each adds only its own extras on top of base rather
than duplicating the stack.

| Env | Focus | Key packages |
|-----|-------|--------------|
| `squid` | Spatial transcriptomics | squidpy, spatialdata |
| `atac` | scATAC-seq | snapatac2, episcanpy |
| `comms` | Cell communication / GRN | liana, cellphonedb, scglue, pyscenic |

They are built the first time you `usepy <name>` them. To (re)build one
explicitly, use the helper, which runs `python3.11 -m venv
--system-site-packages` and installs from the matching requirements file:

```bash
create_layered_venv.sh squid docker/requirements/squid.txt
create_layered_venv.sh atac  docker/requirements/atac.txt
create_layered_venv.sh comms docker/requirements/comms.txt
```

Requirements live in `docker/requirements/{base,squid,atac,comms}.txt`. Each
build also writes a frozen lockfile (`/opt/environments/<name>_frozen.txt`)
capturing inherited plus added packages.

### Your own venv

For project-specific packages that do not belong in a shared env, layer your
own from base:

```bash
python3.11 -m venv --system-site-packages ~/venvs/myproject
source ~/venvs/myproject/bin/activate
pip install some-package
```

### Runtime installs

Packages install into the **active** venv. Since `/opt/venvs/*` is not usually
persisted across container rebuilds, freeze anything you want to keep:

```bash
usepy base
pip install scikit-learn umap-learn
pip freeze > requirements-project.txt
```

Most packages install from wheels; the image keeps `build-essential` and dev
headers so anything needing compilation still builds. `sudo` (NOPASSWD) is
available for system libraries, e.g. `sudo apt-get install -y libhdf5-dev`.

### Jupyter kernels

Two kernels are registered out of the box:

| Kernel name | Language | Interpreter |
|-------------|----------|-------------|
| `python311-scagent` | Python | base venv (`ipykernel`) |
| `ir` | R | system R (IRkernel) |

`jupyter-scatter` is preinstalled in the base venv for interactive embeddings.

## R environments

R is **4.5.3** (built from source) on **Bioconductor 3.22**. The library layout
is **two-tier**:

| Tier | Path | Properties |
|------|------|-----------|
| System | `/usr/local/lib/R/library` | Read-only, ~80 core packages, renv-pinned, shared across containers |
| User | `~/R/x86_64-pc-linux-gnu-library/4.5` | Writable, takes precedence, per-user runtime installs |

`.libPaths()` lists the user library first, so runtime installs shadow the
system copies without touching them.

### Installing packages at runtime

No sudo is needed — packages land in the writable user library:

```r
install.packages("ggExtra")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
remotes::install_github("satijalab/seurat-data")

.libPaths()
# [1] "/home/devuser/R/x86_64-pc-linux-gnu-library/4.5"  # user (writable)
# [2] "/usr/local/lib/R/library"                          # system (core)
```

Heavy annotation packages (`BSgenome.*`, `EnsDb.*`, `org.*.eg.db`) are **not**
preinstalled — pull them on demand, one at a time if memory is tight.

### The "Installation paths not writeable" warning is normal

When installing with `BiocManager`, you will see:

```
Installation paths not writeable, unable to update packages
  path: /usr/local/lib/R/library
  packages:
    aplot, BiocGenerics, Matrix, Seurat, ...
```

This is **expected and harmless**. Your package installed fine into the user
library; BiocManager merely checked whether the read-only system packages need
updates and found it cannot write there (by design — the system tier is pinned
for reproducibility). Suppress it with:

```r
BiocManager::install("PACKAGE", update = FALSE)
```

### Missing system dependency during compilation

If a package fails to configure (e.g. `libxml2 not found`), install the dev
library with sudo and retry:

| R package | System dependency |
|-----------|-------------------|
| XML / xml2 | `libxml2-dev` |
| RCurl / httr | `libcurl4-openssl-dev` |
| sf | `libgdal-dev libproj-dev libgeos-dev` |
| rJava | `openjdk-11-jdk` |

```bash
sudo apt-get update && sudo apt-get install -y libxml2-dev
R -e 'install.packages("XML")'
```

### Reproducibility with renv

The system tier is pinned via `/opt/settings/renv.lock` (manifest at
`/opt/settings/R-packages-manifest.csv`). For per-project pinning, snapshot the
user-library additions:

```r
renv::init()       # once per project
renv::snapshot()   # after installing packages
renv::restore()    # reproduce elsewhere
```

Commit the resulting `renv.lock` to your project repo.

## R ↔ Python interoperability

For multi-modal work (e.g. scRNA + scATAC), move objects between Seurat and the
scanpy/muon ecosystem via on-disk `.h5ad`/`.h5mu`, or bridge in-process with
`reticulate`.

**Export from R (Seurat):**

```r
library(MuDataSeurat)
WriteH5AD(object = s, file = "rna.h5ad", assay = "RNA")  # single assay
WriteH5MU(s, "multiome.h5mu")                            # multiome
```

**Process in Python:**

```bash
usepy base
```

```python
import scanpy as sc, muon as mu
m = mu.read_h5mu("multiome.h5mu")
rna, atac = m.mod["RNA"], m.mod["ATAC"]
# run scVI/PeakVI/scGLUE; store embeddings in obsm["X_scvi"], etc.
mu.write_h5mu("multiome.h5mu", m)
```

**Import back to R:**

```r
library(MuDataSeurat)
s_multi <- ReadH5MU("multiome.h5mu")
names(s_multi@reductions)   # pca, umap, scvi, peakvi, ...
```

**In-process bridge:** `reticulate` is preinstalled and points at the base venv,
so you can call scanpy directly from R when a file round-trip is overkill. Other
converters in the image — `anndataR`, `sceasy`, `zellkonverter` — offer
alternative object translations.
