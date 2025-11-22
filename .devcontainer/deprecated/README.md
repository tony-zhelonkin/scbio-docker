# Deprecated Files

**⚠️ WARNING: Files in this directory are DEPRECATED and will be removed in v0.6.0**

This directory contains files that have been moved to new locations as part of the repository restructuring in v0.5.2.

## What Changed

### Build Files
- **`build-optimized.sh`** → Use `scripts/build.sh` instead
- **`Dockerfile.optimized`** → Use `docker/base/Dockerfile` instead

### R Installation Scripts
- **`install_R_core.R`** → Moved to `docker/base/R/install_core.R`
- **`install_R_packages.R`** → Moved to `docker/base/R/install_core.R`
- **`install_R_archr.R`** → Deprecated (use official ArchR image instead)
- **`install_renv_project.R`** → Moved to `docker/base/R/install_renv_project.R`
- **`install_httpgd.R`** → Moved to `docker/base/R/install_httpgd.R`

### Helper Scripts
- **`create_layered_venv.sh`** → Moved to `scripts/create_layered_venv.sh`
- **`poststart_sanity.sh`** → Moved to `scripts/poststart_sanity.sh`

### Python Requirements
- **`.environments/`** → Moved to `docker/requirements/`
  - `.environments/base_requirements.txt` → `docker/requirements/base.txt`
  - `.environments/squid_requirements.txt` → `docker/requirements/squid.txt`
  - `.environments/atac_requirements.txt` → `docker/requirements/atac.txt`
  - `.environments/comms_requirements.txt` → `docker/requirements/comms.txt`

## Migration Guide

If you have scripts or documentation referencing these old paths, update them as follows:

### Build Commands
```bash
# OLD
./build-optimized.sh --tag my-image:v1

# NEW
scripts/build.sh --tag my-image:v1
```

### Dockerfile References
```bash
# OLD
docker build -f .devcontainer/Dockerfile.optimized .

# NEW
docker build -f docker/base/Dockerfile .
```

### Requirements Files
```dockerfile
# OLD
COPY .environments/base_requirements.txt /opt/requirements.txt

# NEW
COPY docker/requirements/base.txt /opt/requirements.txt
```

## Timeline

- **v0.5.2 (current)**: Files moved to `.devcontainer/deprecated/` with warnings
- **v0.5.3**: Final release with deprecated files present
- **v0.6.0**: This entire directory will be deleted

## Why the Change?

The restructuring improves:
1. **Clarity**: Clear separation between build assets (`docker/`, `scripts/`) and configuration (`.devcontainer/`)
2. **Discoverability**: Docker files in `docker/`, scripts in `scripts/`
3. **Standards**: Aligns with Docker repository best practices
4. **AI-friendliness**: Predictable file locations for tools

## Need Help?

- See [docs/migration.md](../../docs/migration.md) for detailed migration guide
- See [docs/repo-structure.md](../../docs/repo-structure.md) for new structure overview
- See [docs/build.md](../../docs/build.md) for updated build instructions
