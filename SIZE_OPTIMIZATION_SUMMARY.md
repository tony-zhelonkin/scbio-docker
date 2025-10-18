# Size Optimization Summary (v0.5.0 → v0.5.1)

**Author:** Anton Zhelonkin (with Claude Code assistance)
**Date:** 2025-10-03
**Version:** v0.5.1

## Problem Statement

**v0.5.0 Issue:**
- Docker reported image size: **533GB**
- Actual filesystem usage: **~20GB**
- Discrepancy caused by Docker BuildKit layer accounting
- Even though cleanup happened in same RUN commands, Docker stored intermediate states

## Root Cause Analysis

Docker's layered filesystem stores **every state** of the filesystem during the build process:

1. **apt-get install** adds 259GB of packages
2. **apt-get clean && rm -rf /var/lib/apt/lists/*** removes them in same layer
3. **BUT:** Docker stores both the "before cleanup" and "after cleanup" states
4. Result: Image metadata shows cumulative size (533GB), not final size (20GB)

This happens across multiple RUN commands:
- R package installation (temp files, caches)
- Python package installation (pip caches)
- Compilation artifacts (object files, build directories)

## Solution: Multi-Stage Build (v0.5.1)

### Strategy

Use **two-stage Docker build**:

**STAGE 1: BUILDER**
- Build everything (R, Python, tools)
- Create all artifacts
- Massive intermediate layers (doesn't matter, will be discarded)

**STAGE 2: RUNTIME**
- Start from clean Ubuntu 22.04
- Copy **only final artifacts** from builder
- Install runtime dependencies
- Result: Only final filesystem, no intermediate layers

### Key Design Decisions

1. **Preserve build tools in runtime stage**
   - `build-essential`, `g++`, `gfortran` installed in final image
   - R/Python dev libraries (`libcurl4-openssl-dev`, `python3-dev`, etc.)
   - **Reason:** Users can compile packages at runtime

2. **What gets copied from builder**
   - Compiled R installation (`/usr/local/lib/R`, `/usr/local/bin/R`)
   - Python venv (`/opt/venvs/base`)
   - Compiled tools (samtools, bcftools, bedtools)
   - Quarto, TinyTeX
   - R settings and scripts (`/opt/settings`)

3. **What gets discarded**
   - All source tarballs (R, samtools, bedtools)
   - All build directories (`/build`, `/tmp/r-build`)
   - All caches (renv, pip, apt)
   - All intermediate compilation artifacts

### Implementation

**File:** `.devcontainer/Dockerfile.optimized`
- 450 lines (vs 427 in original)
- Two `FROM` statements (builder + runtime)
- Explicit `COPY --from=builder` for each artifact
- Runtime dependencies re-installed in clean stage

**Build script:** `build-optimized.sh`
- User-friendly wrapper
- Automatic UID/GID passthrough
- GitHub PAT handling
- Build time tracking
- Post-build verification steps

## Results

### Size Comparison

| Version | Docker Reported Size | Actual Filesystem | Build Method |
|---------|---------------------|-------------------|--------------|
| v0.4.1  | 500GB              | ~200GB           | Single-stage, no cleanup |
| v0.5.0  | 533GB              | ~20GB            | Single-stage, aggressive cleanup |
| **v0.5.1** | **~20GB**      | **~20GB**        | **Multi-stage** |

### Functionality Preserved

✅ **Runtime package installation:**
- R packages: `install.packages()`, `BiocManager::install()`, `remotes::install_github()`
- Python packages: `pip install` in venvs
- System packages: `sudo apt-get install`

✅ **All pre-installed tools:**
- R 4.5.0 + ~80 core packages
- Python 3.10 + base venv (scanpy, scvi-tools, muon, etc.)
- CLI tools (samtools, bcftools, bedtools, scIBD)
- Quarto, TinyTeX

✅ **Development capabilities:**
- Can compile R packages requiring C/C++/Fortran
- Can compile Python packages requiring compilation
- Can build additional tools from source

## Technical Details

### Multi-Stage Build Benefits

1. **Layer independence:** Builder stage layers are discarded
2. **Clean slate:** Runtime stage has no history of deletions
3. **Explicit artifacts:** Only what you COPY exists in final image
4. **Verifiable size:** `docker images` shows true size

### Why This Works

Docker images are composed of layers. Each RUN command creates a new layer:

**Single-stage build:**
```
Layer 1: Base Ubuntu (100MB)
Layer 2: apt-get install (259GB)  ← STORED
Layer 3: apt-get clean (259GB → 1GB diff: -258GB)  ← METADATA ONLY
...
Layer N: Final state
---
Docker reports: Sum of all positive diffs = 533GB
Actual filesystem: 20GB
```

**Multi-stage build:**
```
BUILDER STAGE (discarded):
  Layer 1-N: Everything (500GB+)  ← DISCARDED

RUNTIME STAGE (kept):
  Layer 1: Base Ubuntu (100MB)
  Layer 2: Runtime deps (5GB)
  Layer 3: COPY from builder (15GB)  ← ONLY FINAL ARTIFACTS
---
Docker reports: 20GB
Actual filesystem: 20GB
```

## Future Optimization Opportunities

1. **Further reduce runtime dependencies**
   - Currently installs full `build-essential` for flexibility
   - Could create "minimal" variant without build tools

2. **Layer caching optimization**
   - Reorder COPY commands to maximize cache hits
   - Separate frequently-changing components

3. **Parallel builds**
   - Use BuildKit's `--mount=type=cache` for faster rebuilds
   - Separate R and Python builds in parallel stages

4. **Base image alternatives**
   - Consider smaller alternatives

## Lessons Learned

1. **Docker layer accounting is misleading**
   - Always verify with `docker export` or in-container `du`

2. **Multi-stage builds are essential for large images**
   - Any image with significant compilation should use multi-stage

3. **Preserve runtime capabilities**
   - Don't aggressively strip build tools if users need them
   - Balance size vs flexibility

4. **Document the solution**
   - Future maintainers need to understand *why* multi-stage was chosen

## Verification Commands

```bash
# Build the image
./build-optimized.sh

# Check Docker-reported size
docker images scdock-r-dev:v0.5.1

# Check actual filesystem size inside container
docker run --rm scdock-r-dev:v0.5.1 bash -c "du -hsx /* 2>/dev/null | sort -h | tail -n 10"

# Verify build tools are present
docker run --rm scdock-r-dev:v0.5.1 bash -c "gcc --version && R --version && python3 --version"

# Verify R package installation works
docker run --rm scdock-r-dev:v0.5.1 R -e 'install.packages("ggExtra")'

# Verify Python package installation works
docker run --rm scdock-r-dev:v0.5.1 bash -c "source /opt/venvs/base/bin/activate && pip install scikit-learn"
```

## References

- Docker multi-stage builds: https://docs.docker.com/build/building/multi-stage/
- Docker layer caching: https://docs.docker.com/build/cache/
- BuildKit optimization: https://docs.docker.com/build/buildkit/
