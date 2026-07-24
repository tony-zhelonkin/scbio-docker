Migration Notes: dev-restructure

Summary
- Build assets moved from .devcontainer/ and .environments/ into docker/ and scripts/.
- New canonical build entry points:
  - Dockerfile: docker/base/Dockerfile
  - Build script: scripts/build.sh
  - Python requirements: docker/requirements/*.txt
  - R installers: docker/base/R/*.R

Compat
- Legacy files remain for one release to avoid breaking existing workflows.
- init-project.sh now copies sanity script from scripts/poststart_sanity.sh.

Action items for users
- Replace any direct references to .devcontainer/Dockerfile(.optimized) with docker/base/Dockerfile.
- Replace .environments/*.txt references with docker/requirements/*.txt.
- Use scripts/build.sh instead of build-optimized.sh.

