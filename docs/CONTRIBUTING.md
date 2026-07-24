Contributing

This project was created as a personal toolkit. If you want this repo helpful, and wish to contribute, great — open a PR or issue. Keep it small and practical.

Quick build check (sanity script is mounted, not baked into the image)
- scripts/build.sh --tag scdock-r-dev:dev-local
- docker run --rm -v "$PWD:/repo" -w /repo scdock-r-dev:dev-local bash scripts/poststart_sanity.sh

Layout
- docker/: Dockerfiles and installers
- scripts/: build/init/sanity helpers
- .devcontainer/: VS Code config
- docs/: longer notes
