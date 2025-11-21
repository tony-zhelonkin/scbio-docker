# Handoff Notes (v0.5.2)

Short and practical. What changed, where things live, and how to build/test.

Summary
- Added libhdf5-dev and libgsl-dev so R packages can compile at runtime.
- Expanded the R stack: chromVAR, motifmatchr, TFBSTools, JASPAR2022, SingleR, celldex, EnsDb.Mmusculus.v79, etc.
- Fixed safe_install() to treat meta-packages correctly (tidyverse installs as expected).
- Version aligned to v0.5.2.

Where things live
- Dockerfile: docker/base/Dockerfile
- R installers: docker/base/R/*.R
- Python requirements: docker/requirements/*.txt
- Build script: scripts/build.sh
- Sanity script: scripts/poststart_sanity.sh

Build
```bash
# Generic build
scripts/build.sh --tag scdock-r-dev:v0.5.2

# With GitHub PAT (avoids rate limits)
scripts/build.sh --github-pat ghp_xxxxx --tag scdock-r-dev:v0.5.2
```

Quick test
```bash
# Sanity
docker run --rm scdock-r-dev:v0.5.2 bash -lc 'scripts/poststart_sanity.sh'

# Spot-check R packages
docker run --rm scdock-r-dev:v0.5.2 R -q -e "library(chromVAR); library(SingleR); library(TFBSTools); cat('OK\\n')"
```

Lockfiles (first build)
```bash
CID=$(docker create scdock-r-dev:v0.5.2)
docker cp $CID:/opt/settings/renv.lock ./renv.lock
docker cp $CID:/opt/settings/R-packages-manifest.csv ./R-packages-manifest.csv
docker rm $CID
```

Notes
- Multi-stage build preserved; runtime includes build tools for on-demand installs.
- Docker may report larger size than actual filesystem. Use in-container du to verify.
- If you hit GitHub rate limits during build, pass a PAT.

That’s it. Build, sanity, and go.

