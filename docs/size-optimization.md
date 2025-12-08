# Image Size: the short version

Problem
- Single-stage builds made Docker report huge sizes (layer accounting), even when actual filesystem use was ~20GB.

Fix
- Switch to a multi-stage build: compile in a builder stage, copy only final artifacts to a fresh runtime image.

What stays in runtime
- R 4.5 + core packages
- Python base venv
- CLI tools (samtools, bcftools, bedtools)
- TinyTeX
- Build tools needed for runtime installs (kept intentionally)

How to build now
```bash
scripts/build.sh --tag scdock-r-dev:v0.5.2
# or
docker build . -f docker/base/Dockerfile -t scdock-r-dev:v0.5.2
```

Verify
```bash
docker images scdock-r-dev:v0.5.2
docker run --rm scdock-r-dev:v0.5.2 bash -c "du -hsx /* 2>/dev/null | sort -h | tail -n 10"
```

Notes
- Expect Docker to sometimes report larger numbers than the in-container du. The latter reflects real disk use.

