Contributing

This is a personal toolkit. If you want to improve something, great—open a PR or issue. Keep it small and practical.

Quick build check
- scripts/build.sh --tag scdock-r-dev:dev-local
- docker run --rm scdock-r-dev:dev-local bash -lc 'scripts/poststart_sanity.sh'

Layout
- docker/: Dockerfiles and installers
- scripts/: build/init/sanity helpers
- .devcontainer/: VS Code config
- docs/: longer notes
