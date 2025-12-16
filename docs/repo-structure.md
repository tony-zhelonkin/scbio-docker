# Repository Structure

## scbio-docker

```
scbio-docker/
├── docker/
│   └── base/
│       └── Dockerfile         # Primary multi-stage build
├── .devcontainer/
│   ├── install_R_core.R       # Core R stack installer
│   ├── install_renv_project.R # renv init/restore wrapper
│   ├── install_httpgd.R       # CRAN-first + GitHub fallback
│   ├── install_quarto.sh      # Quarto installation
│   ├── create_layered_venv.sh # Runtime layered venv helper
│   └── scripts/
│       └── poststart_sanity.sh
├── docker/
│   └── requirements/
│       ├── base.txt           # Core Python stack
│       ├── squid.txt          # Spatial transcriptomics
│       ├── atac.txt           # scATAC-seq tools
│       └── comms.txt          # Cell communication
├── scripts/
│   └── build.sh               # Docker build wrapper
├── templates/
│   ├── base/                  # Standard project template
│   │   └── README.md
│   ├── docs/
│   │   ├── README.md.template
│   │   ├── tasks.md.template
│   │   ├── notes.md.template
│   │   └── .env.example
│   ├── config/
│   │   ├── config.R.template
│   │   ├── pipeline.yaml.template
│   │   └── color_config.R.template
│   └── .vscode/
│       └── settings.json
├── toolkits/
│   └── SciAgent-toolkit/      # Git submodule
├── docs/
│   ├── architecture.md
│   ├── build.md
│   ├── devops.md
│   └── ...
├── init-project.sh            # Project scaffolding script
├── CLAUDE.md                  # AI context for this repo
└── README.md                  # Overview + links
```

## Project Structure (After init-project.sh + setup-ai.sh)

```
my-project/
├── .devcontainer/
│   ├── devcontainer.json
│   ├── docker-compose.yml
│   ├── .env
│   └── scripts/
│       └── poststart_sanity.sh
├── .vscode/
│   └── settings.json
├── 00_data/
│   ├── raw/
│   ├── processed/
│   └── references/
├── 01_modules/
│   ├── RNAseq-toolkit/        # Git submodule
│   └── SciAgent-toolkit/      # Git submodule
├── 02_analysis/
│   ├── config/
│   │   ├── config.R
│   │   ├── pipeline.yaml
│   │   ├── color_config.R
│   │   └── analysis_config.yaml  # Created by setup-ai.sh
│   └── helpers/
├── 03_results/
│   ├── checkpoints/
│   ├── plots/
│   └── tables/
├── logs/
├── .claude/                   # Created by setup-ai.sh
│   ├── agents/
│   └── skills/
├── .mcp.json                  # Created by setup-ai.sh
├── CLAUDE.md                  # Created by setup-ai.sh
├── GEMINI.md                  # Created by setup-ai.sh
├── AGENTS.md                  # Created by setup-ai.sh
├── context.md                 # Created by setup-ai.sh
├── tasks.md
├── notes.md
├── README.md
└── .gitignore
```

## Separation of Concerns

| Repository | Responsibility |
|------------|----------------|
| **scbio-docker** | Docker images, container setup, project directory structure |
| **SciAgent-toolkit** | AI tools, MCP servers, agents/skills, methodology guidelines |

## Key Paths

- **Image build**: `scripts/build.sh` or `docker build -f docker/base/Dockerfile`
- **Project setup**: `./init-project.sh ~/projects/my-analysis`
- **AI setup**: `./01_modules/SciAgent-toolkit/scripts/setup-ai.sh` (inside container)
- **Guidelines**: `01_modules/SciAgent-toolkit/docs/guidelines/` (single source of truth)
