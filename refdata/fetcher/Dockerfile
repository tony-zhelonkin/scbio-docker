# refdata-fetcher — one small image that can run any source under refdata/sources/.
#
# Replaces the per-source coresh-updater image. Credentials are passed at
# `docker run` time and never baked in; the analysis dev container never sees
# them, which is the reason fetching is a separate image at all.
#
# Build (context is refdata/, NOT refdata/fetcher — the COPYs below reach up to
# refcache.sh and sources/):
#   docker build -t refdata-fetcher:latest -f refdata/fetcher/Dockerfile refdata
#
# Run:
#   docker run --rm -it \
#     -v /data2/users/shared/refcache:/cache \
#     -e SYNAPSE_AUTH_TOKEN=... \
#     refdata-fetcher:latest coresh
#
#   docker run --rm -it \
#     -v /data2/users/shared/refcache:/cache \
#     -e CISTARGET_SETS=gene_hg38,gene_mm10,motif2tf \
#     refdata-fetcher:latest cistarget
FROM python:3.12-slim

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      ca-certificates curl jq coreutils findutils \
 && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir 'synapseclient>=4.0.0'

WORKDIR /work
COPY refcache.sh /opt/refdata/refcache.sh
COPY sources/ /opt/refdata/sources/
RUN chmod +x /opt/refdata/refcache.sh

# Bind-mount the cache root here.
ENV CACHE_ROOT=/cache

ENTRYPOINT ["/opt/refdata/refcache.sh"]
CMD ["--list"]
