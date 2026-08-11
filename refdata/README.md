# refdata — shared reference-data cache

The **mechanism** for fetching and versioning large public reference data.
The **bytes** live outside this repo, on a host path you choose.

That split is the whole point: one copy of each database on the machine,
mounted read-only into every container, refreshed by updating the cache — not
by rebuilding images or re-downloading per project.

> Nothing here is baked into `scdock-r-dev`. Reference data is orders of
> magnitude larger than the image and changes on a different clock.

## Placement, and when to extract it

This is a **satellite subsystem**: it runs on the upstream-data clock, not the
image clock, so it sits outside `VERSION` and `docs/changelog.md` and versions
itself here. It lives in this repo because the `:ro` mount seam
(`--refcache`, `REFCACHE_ROOT`) is genuinely devcontainer wiring and because
extraction stays cheap — there are no imports either way.

**Extract to its own repo when a third source lands in `sources/`.** Two is a
deployment detail; three is a subsystem with a growth curve, at which point the
credential matrix, the fetcher's own dependency stack and an independent release
cadence stop being hypothetical. Other triggers: a consumer that never renders a
devcontainer (an HPC or bare-metal job), or a collaborator who needs the
fetchers but not the image.

No host path belongs here. `/data2/users/shared/refcache` is one deployment,
passed in via `--refcache`, exactly like project data mounts.

## Layout

```
refdata/
├── refcache.sh          # generic driver: snapshot → verify → manifest → flip → prune
├── sources/
│   ├── cistarget.sh     # aertslab motif databases (feathers + motif2tf)
│   └── coresh.sh        # CoReSh GEO chunk compendium (Synapse syn66227307)
└── fetcher/Dockerfile   # one image that runs any source
```

Two different paths are involved, and they must not be confused:

| Variable | Side | Value | Set by |
|---|---|---|---|
| `REFCACHE_HOST` | host | e.g. `/data2/users/shared/refcache` | you, when fetching or rendering |
| `REFCACHE_ROOT` | container | always `/refcache` | `init-container.sh`, into the project `.env` |

Fetchers bind `$REFCACHE_HOST`; analysis code inside a container reads
`$REFCACHE_ROOT`. On the host, each source gets a directory of dated snapshots
plus a `current` symlink:

```
$REFCACHE_HOST/
├── cistarget/
│   ├── cistarget_20260810/
│   │   ├── region_based/  gene_based/  motif2tf/
│   │   └── MANIFEST.json
│   └── current -> cistarget_20260810
└── coresh/
    ├── syn66227307_20260721/
    └── current -> syn66227307_20260721
```

**Always resolve through `current/`.** The flip is a `rename(2)`, so a reader
never sees a missing or half-written snapshot, and a refresh does not disturb
a running analysis that already opened files under the old snapshot.

## Usage

Build the fetcher once:

```bash
docker build -t refdata-fetcher:latest -f refdata/fetcher/Dockerfile refdata
```

See what a refresh would cost before committing to it — `--dry-run` on
`cistarget` reports per-file sizes from upstream `Content-Length`:

```bash
docker run --rm -v "$REFCACHE_HOST:/cache" refdata-fetcher:latest cistarget --dry-run
```

Then fetch:

```bash
# cisTarget — set selection drives the size, see the table below
docker run --rm -v "$REFCACHE_HOST:/cache" \
  -e CISTARGET_SETS=gene_hg38,gene_mm10,motif2tf \
  refdata-fetcher:latest cistarget

# CoReSh — needs a Synapse PAT with view + download
docker run --rm -v "$REFCACHE_HOST:/cache" \
  -e SYNAPSE_AUTH_TOKEN="$SYNAPSE_AUTH_TOKEN" \
  refdata-fetcher:latest coresh
```

Common flags: `--dry-run`, `--keep N` (old snapshots to retain, default 1),
`--force` (resume into an existing snapshot dir instead of refusing),
`--verify-only` (re-check a snapshot on disk against its stored sha1sums).

Run as yourself, not root, or the snapshots land root-owned and nobody else can
prune them:

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$REFCACHE_HOST:/cache" ...
```

## cisTarget set sizes

Measured from upstream `Content-Length`, 2026-08. Pick with `CISTARGET_SETS`
(comma-separated); the default is everything, **≈72 GB**.

| Set | Size | Needed for |
|---|---:|---|
| `region_hg38` | 45.7 GB | SCENIC+ on human ATAC peaks |
| `region_mm10` | 24.2 GB | SCENIC+ on mouse ATAC peaks |
| `gene_hg38` | 1.15 GB | pySCENIC, human |
| `gene_mm10` | 0.90 GB | pySCENIC, mouse |
| `motif2tf` | 0.20 GB | essentially every cisTarget run |

The two region-based sets are the entire cost. If you are only running
pySCENIC on RNA, `gene_hg38,gene_mm10,motif2tf` is ~2.3 GB.

Feathers are verified against the `.sha1sum.txt` files published alongside
them. This is slow on a 33 GB ranking database, but a silently truncated
feather yields *plausible wrong* motif enrichments rather than an error, so
the refresh pays it once.

## Consuming from a container

`init-project.sh --refcache PATH` emits `REFCACHE_ROOT=/refcache` into the
compose service's `environment:` block (not into `.env`) and binds the cache
read-only at `/refcache`. In analysis code, resolve from the env var
— never hardcode a snapshot tag:

```python
import os, pathlib
ct = pathlib.Path(os.environ["REFCACHE_ROOT"], "cistarget", "current")
rankings = ct / "region_based" / "mm10_screen_v10_clust.regions_vs_motifs.rankings.feather"
motif2tf = ct / "motif2tf" / "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"
```

## Adding a source

Drop `sources/<name>.sh` defining `SOURCE_NAME`, `src_fetch <dir>` and
`src_verify <dir>` (optionally `src_manifest_extra <dir>`, `src_dry_run`, and
`SOURCE_SNAPSHOT_TAG`). Add a `# desc:` line so it shows up in `--list`.
Everything else — dated snapshot, `MANIFEST.json`, atomic flip, pruning — is
the driver's job.

Good candidates to fold in next: the ad-hoc `download_*.sh` scripts under the
data lake's `_scripts/`, which currently have no snapshotting or manifests.
