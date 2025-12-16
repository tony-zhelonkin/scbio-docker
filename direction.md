# Unified Implementation Plan: scbio-docker + SciAgent-toolkit Architecture

**Created**: 2025-12-16
**Status**: Approved for Implementation

---

## Executive Summary

**Strict separation of concerns:**

| Repository | Responsibility | AI Involvement |
|------------|----------------|----------------|
| **scbio-docker** | Docker images, container spin-up, project directory structure | NONE - only root CLAUDE.md/AGENTS.md for its own codebase |
| **SciAgent-toolkit** | ALL AI tooling: installation, MCP servers, profiles, templates, methodology guidelines | OWNS everything AI-related for analysis projects |

**Flow:**
```
1. scbio-docker/init-project.sh → Creates project structure (NO AI files)
2. User opens project in container (VS Code Dev Container)
3. SciAgent-toolkit/scripts/setup-ai.sh → Installs AI tools, creates CLAUDE.md/GEMINI.md/.claude/
4. User runs analysis with AI assistance
```

---

## Key Decisions (User Feedback Incorporated)

| Decision | Rationale |
|----------|-----------|
| `01_scripts/` → `01_modules/` | Contains reusable modules, not project-specific scripts |
| `analysis_config.yaml` at `02_analysis/config/` | Co-locate with other configs |
| setup-ai.sh creates `.claude/`, `.claude/agents/`, `.claude/skills/` | Claude Code requires these directories |
| Consolidate to `context.md`, `tasks.md`, `notes.md` | Minimal set, clear purposes |
| Base role only (for now) | Future roles system deferred |

---

## Target Architecture

### File Ownership

```
scbio-docker/                              # CONTAINER INFRASTRUCTURE ONLY
├── CLAUDE.md                              # FOR scbio-docker CODEBASE ONLY
├── AGENTS.md                              # FOR scbio-docker CODEBASE (optional)
│
├── docs/
│   ├── architecture.md                    # Docker architecture
│   ├── build.md                           # Build instructions
│   └── guidelines/                        # → DELETE (move to SciAgent-toolkit)
│
├── templates/
│   └── docs/
│       ├── CLAUDE.md.template             # → DELETE (SciAgent-toolkit owns)
│       ├── plan.md.template               # → DELETE (replaced by context.md)
│       ├── tasks.md.template              # KEEP
│       └── notes.md.template              # KEEP
│
├── toolkits/
│   └── SciAgent-toolkit/                  # Submodule
│
└── init-project.sh                        # MODIFY: 01_scripts→01_modules, use context.md

SciAgent-toolkit/                          # ALL AI INFRASTRUCTURE + METHODOLOGY
├── CLAUDE.md                              # FOR SciAgent-toolkit CODEBASE ONLY
├── AGENTS.md                              # FOR SciAgent-toolkit CODEBASE (optional)
│
├── docs/
│   ├── guidelines/                        # CANONICAL - methodology for AI agents
│   │   ├── README.md
│   │   ├── core_architecture.md
│   │   ├── data_processing.md
│   │   ├── gsea_analysis.md
│   │   ├── checkpoint_caching.md
│   │   ├── master_tables.md
│   │   ├── visualization.md
│   │   └── code_style.md
│   │
│   ├── ISSUES.md                          # Architecture issues
│   └── TROUBLESHOOTING.md                 # Common problems
│
├── templates/
│   ├── vendor/                            # Templates for ANALYSIS PROJECTS
│   │   ├── CLAUDE.md.template             # AI context for Claude
│   │   ├── GEMINI.md.template             # CREATE - AI context for Gemini
│   │   ├── AGENTS.md.template             # CREATE - Universal rules
│   │   ├── context.md.template            # CREATE - Scientific context
│   │   └── analysis_config.yaml.template  # Project parameters
│   │
│   ├── mcp-profiles/                      # MCP configurations
│   │   └── *.mcp.json                     # FIX: Issues 001, 002, 004, 005
│   │
│   └── gemini-profiles/                   # Gemini configurations
│
├── roles/                                 # Role definitions (YAML specs, NO agent/skill copies)
│   ├── base.yaml                          # Default role
│   ├── coding.yaml                        # Code-focused role (future)
│   └── research.yaml                      # Research-focused role (future)
│
├── agents/                                # SINGLE SOURCE OF TRUTH - all agent.md files
│   ├── bioinf-librarian.md
│   ├── rnaseq-methods-writer.md
│   └── figure-caption-generator.md
│
├── skills/                                # SINGLE SOURCE OF TRUTH - all skill.md files
│   └── [skills go here]
│
└── scripts/
    ├── setup-ai.sh                        # MODIFY: Create .claude/ directories
    ├── setup_mcp_infrastructure.sh
    ├── configure_mcp_servers.sh
    ├── switch-mcp-profile.sh              # FIX: API key substitution
    └── mcp_servers/
        └── setup_*.sh
```

### Project Structure (After Both Tools Run)

```
my-project/
├── .devcontainer/                         # From scbio-docker init-project.sh
│   ├── devcontainer.json
│   ├── docker-compose.yml
│   └── .env                               # API keys (gitignored)
│
├── 00_data/                               # From scbio-docker
│   ├── raw/
│   ├── processed/
│   └── references/
│
├── 01_modules/                            # RENAMED from 01_scripts/
│   ├── RNAseq-toolkit/                    # Submodule (optional)
│   └── SciAgent-toolkit/                  # Submodule
│       └── docs/guidelines/               # Referenced, not copied
│
├── 02_analysis/
│   ├── config/
│   │   ├── config.R                       # From scbio-docker
│   │   ├── pipeline.yaml                  # From scbio-docker
│   │   ├── color_config.R                 # From scbio-docker
│   │   └── analysis_config.yaml           # From SciAgent-toolkit setup-ai.sh
│   └── helpers/                           # Project-specific scripts
│
├── 03_results/
│   ├── checkpoints/
│   ├── plots/
│   └── tables/
│
├── .claude/                               # From SciAgent-toolkit setup-ai.sh
│   ├── agents/                            # Populated based on role
│   │   └── [symlinks or copies from role]
│   ├── skills/                            # Populated based on role
│   │   └── [symlinks or copies from role]
│   └── settings.local.json                # Generated by profile switcher
│
├── .gemini/                               # From SciAgent-toolkit switch-mcp-profile.sh
│   └── settings.json
│
├── .mcp.json                              # From SciAgent-toolkit setup-ai.sh
├── tooluniverse-env/                      # From SciAgent-toolkit setup-ai.sh
│
├── CLAUDE.md                              # From SciAgent-toolkit setup-ai.sh
├── GEMINI.md                              # From SciAgent-toolkit setup-ai.sh
├── AGENTS.md                              # From SciAgent-toolkit setup-ai.sh
│
├── context.md                             # Scientific context, broad goals
├── tasks.md                               # Concrete implementation steps
└── notes.md                               # Web research notes
```

---

## Implementation Phases

### Phase 1: Fix Critical MCP Issues

#### 1.1 Fix ISSUE-002: API Key Substitution

**File**: `SciAgent-toolkit/scripts/switch-mcp-profile.sh`
**Lines**: ~129-133

```python
# ADD after existing substitutions:
content = content.replace('\${GEMINI_API_KEY}', os.environ.get('GEMINI_API_KEY', ''))
content = content.replace('\${OPENAI_API_KEY}', os.environ.get('OPENAI_API_KEY', ''))
content = content.replace('\${CONTEXT7_API_KEY}', os.environ.get('CONTEXT7_API_KEY', ''))
```

#### 1.2 Evaluate PAL Invocation Strategy

Test wrapper vs direct uvx (see TROUBLESHOOTING.md). Document winner.

#### 1.3 Fix ISSUE-005: research-full.mcp.json

Add explicit `--include-tools` argument.

#### 1.4 Fix ISSUE-004: Profile Validation

Add validation function before applying profile.

---

### Phase 2: Enhance setup-ai.sh

#### 2.1 Activate Base Role (Creates .claude/ directories)

**Add to setup-ai.sh** (after MCP configuration):

```bash
# Step: Activate base role (creates .claude/agents/, .claude/skills/)
log_info "Activating base role..."

# activate-role.sh reads roles/base.yaml and symlinks
# specified agents/skills from agents/ and skills/ to .claude/
"${SCIAGENT_SCRIPTS}/activate-role.sh" base --project-dir "${PROJECT_DIR}"
```

This delegates to `activate-role.sh` which:
1. Reads `roles/base.yaml` to get list of agents and skills
2. Creates `.claude/agents/` and `.claude/skills/` directories
3. Symlinks specified files from `agents/` and `skills/` folders

#### 2.2 Create New Templates

**GEMINI.md.template**:
```markdown
# GEMINI.md - Gemini AI Context for {{PROJECT_ID}}

> **Thin wrapper.** Full methodology: `01_modules/SciAgent-toolkit/AGENTS.md`

---

## Project: {{PROJECT_ID}}

## Your Role: Research Focus

Large context window - ideal for:
- Literature synthesis
- Cross-referencing data sources
- Documentation analysis

## Guidelines

See: `01_modules/SciAgent-toolkit/docs/guidelines/`

## Context

See: `context.md`
```

**AGENTS.md.template**:
```markdown
# AGENTS.md - Universal AI Rules for {{PROJECT_ID}}

> Full guidelines: `01_modules/SciAgent-toolkit/docs/guidelines/`

## Critical Rules

1. Annotate genes BEFORE filtering
2. Use `filterByExpr()` not manual thresholds
3. Cache anything >1 minute
4. Colorblind-safe palettes

## Config

Import from: `02_analysis/config/analysis_config.yaml`
```

**context.md.template**:
```markdown
# Project Context: {{PROJECT_ID}}

**Created:** {{DATE}}

---

## Scientific Question

[What biological question are you investigating?]

## Datasets

[What data do you have? Sample sizes, conditions, etc.]

## Hypotheses

[What do you expect to find?]

## Analysis Goals

[What outputs do you need? DE lists, pathway analysis, figures?]

---

## Broad Direction

[High-level strategy and approach]
```

#### 2.3 Template Installation in setup-ai.sh

```bash
# Step: Create AI context files
log_info "Creating AI context files..."

VENDOR_TEMPLATES="${SCIAGENT_SCRIPTS}/../templates/vendor"
PROJECT_NAME=$(basename "$PROJECT_DIR")

substitute_template() {
    local template="$1"
    local output="$2"

    if [ -f "$template" ]; then
        cp "$template" "$output"
        sed -i "s|{{PROJECT_ID}}|${PROJECT_NAME}|g" "$output"
        sed -i "s|{{DATE}}|$(date +%Y-%m-%d)|g" "$output"
        log_ok "  Created $(basename "$output")"
    fi
}

# Only create if not exists (preserve user customizations)
[ ! -f "${PROJECT_DIR}/CLAUDE.md" ] && \
    substitute_template "${VENDOR_TEMPLATES}/CLAUDE.md.template" "${PROJECT_DIR}/CLAUDE.md"

[ ! -f "${PROJECT_DIR}/GEMINI.md" ] && \
    substitute_template "${VENDOR_TEMPLATES}/GEMINI.md.template" "${PROJECT_DIR}/GEMINI.md"

[ ! -f "${PROJECT_DIR}/AGENTS.md" ] && \
    substitute_template "${VENDOR_TEMPLATES}/AGENTS.md.template" "${PROJECT_DIR}/AGENTS.md"

[ ! -f "${PROJECT_DIR}/context.md" ] && \
    substitute_template "${VENDOR_TEMPLATES}/context.md.template" "${PROJECT_DIR}/context.md"

# analysis_config.yaml goes to 02_analysis/config/
mkdir -p "${PROJECT_DIR}/02_analysis/config"
[ ! -f "${PROJECT_DIR}/02_analysis/config/analysis_config.yaml" ] && \
    substitute_template "${VENDOR_TEMPLATES}/analysis_config.yaml.template" \
        "${PROJECT_DIR}/02_analysis/config/analysis_config.yaml"
```

---

### Phase 3: Update scbio-docker

#### 3.1 Rename 01_scripts → 01_modules

**File**: `scbio-docker/init-project.sh`

```bash
# Change line creating directories:
for dir in 00_data/raw 00_data/processed 00_data/references \
           01_modules 02_analysis/config 02_analysis/helpers \
           03_results/checkpoints 03_results/plots 03_results/tables logs; do
    mkdir -p "${PROJECT_DIR}/${dir}"
done
```

Update submodule paths throughout script.

#### 3.2 Consolidate Doc Templates

- Delete `templates/docs/CLAUDE.md.template` (SciAgent-toolkit owns)
- Delete `templates/docs/plan.md.template` (replaced by context.md)
- Keep `templates/docs/tasks.md.template`
- Keep `templates/docs/notes.md.template`

#### 3.3 Remove CLAUDE.md Creation from init-project.sh

Delete lines ~304-314 that create CLAUDE.md.

---

### Phase 4: Create Role System Structure

#### 4.1 Role Definition Format

**Create**: `SciAgent-toolkit/roles/base.yaml`

```yaml
# base.yaml - Default analysis role
name: base
description: Default bioinformatics analysis role
mcp_profile: coding

# Agents to activate (reference by filename without .md)
agents:
  - bioinf-librarian
  - rnaseq-methods-writer

# Skills to activate (reference by filename without .md)
skills: []
  # - checkpoint-caching  # Example for future
```

#### 4.2 Role Activation Script

**Create**: `SciAgent-toolkit/scripts/activate-role.sh`

```bash
#!/usr/bin/env bash
# activate-role.sh - Activates a role by populating .claude/ directories
#
# Usage: ./activate-role.sh [role-name] [--project-dir DIR]
#
# Reads role YAML, symlinks specified agents/skills to .claude/

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLKIT_DIR="$(dirname "$SCRIPT_DIR")"
ROLE="${1:-base}"
PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"

ROLE_FILE="${TOOLKIT_DIR}/roles/${ROLE}.yaml"
AGENTS_DIR="${TOOLKIT_DIR}/agents"
SKILLS_DIR="${TOOLKIT_DIR}/skills"

CLAUDE_DIR="${PROJECT_DIR}/.claude"

# Parse YAML (simple grep/sed for portability)
get_yaml_list() {
    local file="$1"
    local key="$2"
    awk "/^${key}:$/,/^[a-z]/" "$file" | grep "^\s*-" | sed 's/.*- //'
}

log_info() { echo -e "[INFO] $*"; }
log_ok() { echo -e "[OK] $*"; }

# Create directories
mkdir -p "${CLAUDE_DIR}/agents"
mkdir -p "${CLAUDE_DIR}/skills"

# Clear existing symlinks (for role switching)
find "${CLAUDE_DIR}/agents" -type l -delete 2>/dev/null || true
find "${CLAUDE_DIR}/skills" -type l -delete 2>/dev/null || true

log_info "Activating role: ${ROLE}"

# Symlink agents
for agent in $(get_yaml_list "$ROLE_FILE" "agents"); do
    src="${AGENTS_DIR}/${agent}.md"
    if [ -f "$src" ]; then
        ln -sf "$src" "${CLAUDE_DIR}/agents/${agent}.md"
        log_ok "  Agent: ${agent}"
    fi
done

# Symlink skills
for skill in $(get_yaml_list "$ROLE_FILE" "skills"); do
    src="${SKILLS_DIR}/${skill}.md"
    if [ -f "$src" ]; then
        ln -sf "$src" "${CLAUDE_DIR}/skills/${skill}.md"
        log_ok "  Skill: ${skill}"
    fi
done

log_ok "Role '${ROLE}' activated"
```

#### 4.3 Integration with setup-ai.sh

setup-ai.sh calls activate-role.sh after creating directories:

```bash
# In setup-ai.sh, after MCP configuration:

# Activate base role (creates .claude/agents/, .claude/skills/ with symlinks)
log_info "Activating base role..."
"${SCIAGENT_SCRIPTS}/activate-role.sh" base --project-dir "${PROJECT_DIR}"
```

#### 4.4 Single Source of Truth

```
SciAgent-toolkit/
├── agents/                    # ALL agents live here (single source)
│   ├── bioinf-librarian.md
│   └── rnaseq-methods-writer.md
│
├── skills/                    # ALL skills live here (single source)
│   └── [future skills]
│
├── roles/                     # YAML definitions only (no copies)
│   └── base.yaml              # References agents/skills by name
│
└── scripts/
    └── activate-role.sh       # Reads YAML, symlinks to .claude/

Project/.claude/               # Populated by activate-role.sh
├── agents/                    # Symlinks → SciAgent-toolkit/agents/
│   ├── bioinf-librarian.md → ...
│   └── rnaseq-methods-writer.md → ...
└── skills/                    # Symlinks → SciAgent-toolkit/skills/
    └── [active skills]
```

---

### Phase 5: Clean Up Duplication

#### 5.1 Delete from scbio-docker:
- `docs/guidelines/rnaseq-analysis-guidelines.md`
- `templates/docs/CLAUDE.md.template`
- `templates/docs/plan.md.template`

#### 5.2 Update References

Update any documentation pointing to old locations.

---

## Validation Checklist

### Phase 1 - Fix MCP Issues:
- [ ] `switch-mcp-profile.sh` substitutes API keys in .mcp.json
- [ ] PAL invocation approach tested and documented in TROUBLESHOOTING.md
- [ ] research-full.mcp.json has --include-tools
- [ ] Profile validation warns about missing dependencies

### Phase 2 - Enhance setup-ai.sh:
- [ ] setup-ai.sh calls `activate-role.sh base`
- [ ] `.claude/agents/` and `.claude/skills/` populated with symlinks
- [ ] setup-ai.sh creates CLAUDE.md, GEMINI.md, AGENTS.md, context.md
- [ ] analysis_config.yaml placed in `02_analysis/config/`

### Phase 3 - Update scbio-docker:
- [ ] init-project.sh creates `01_modules/` not `01_scripts/`
- [ ] init-project.sh does NOT create CLAUDE.md
- [ ] Old AI templates deleted from scbio-docker

### Phase 4 - Create Role System:
- [ ] `roles/base.yaml` exists with agents/skills list
- [ ] `activate-role.sh` script reads YAML, creates symlinks
- [ ] All agents live in single `agents/` folder
- [ ] All skills live in single `skills/` folder

### Phase 5 - Clean Up:
- [ ] `scbio-docker/docs/guidelines/` deleted
- [ ] `scbio-docker/templates/docs/CLAUDE.md.template` deleted

---

## Test Protocol

### Full Integration:
```bash
# 1. Create project with scbio-docker
cd scbio-docker
./init-project.sh /tmp/test-project basic-rna

# Verify structure
ls /tmp/test-project/01_modules/  # Should exist
ls /tmp/test-project/01_scripts/  # Should NOT exist
ls /tmp/test-project/CLAUDE.md    # Should NOT exist (yet)

# 2. Open in container (or simulate)
cd /tmp/test-project

# 3. Run AI setup
./01_modules/SciAgent-toolkit/scripts/setup-ai.sh

# 4. Verify AI files created
ls -la CLAUDE.md GEMINI.md AGENTS.md context.md

# 5. Verify role activation (symlinks to SciAgent-toolkit)
ls -la .claude/agents/
# Should show symlinks: bioinf-librarian.md -> ../../01_modules/SciAgent-toolkit/agents/...

ls -la .claude/skills/
# Should be empty for base role (no skills defined yet)

# 6. Verify config placement
ls -la 02_analysis/config/analysis_config.yaml

# 7. Verify MCP API key substitution
export GEMINI_API_KEY="test-key-12345"
./01_modules/SciAgent-toolkit/scripts/switch-mcp-profile.sh coding
grep "test-key-12345" .mcp.json  # Should find the actual key
```

### Test Role Switching (Future):
```bash
# Switch to a different role
./01_modules/SciAgent-toolkit/scripts/activate-role.sh research

# Verify .claude/agents/ and .claude/skills/ updated
ls -la .claude/agents/
```

---

## Summary

**Changes to scbio-docker:**
- Rename `01_scripts/` → `01_modules/`
- Delete AI templates (SciAgent-toolkit owns)
- Keep tasks.md, notes.md templates
- Remove CLAUDE.md creation from init-project.sh

**Changes to SciAgent-toolkit:**
- Fix MCP issues (API keys, validation, tool limits)
- Create `activate-role.sh` script that reads role YAML and symlinks agents/skills
- Integrate role activation into setup-ai.sh
- Add template copying (CLAUDE.md, GEMINI.md, AGENTS.md, context.md)
- Create analysis_config.yaml in 02_analysis/config/
- Create new templates: GEMINI.md.template, AGENTS.md.template, context.md.template
- Create `roles/base.yaml` role definition

**Single Sources of Truth:**
- **Agents**: `SciAgent-toolkit/agents/` (all agent.md files)
- **Skills**: `SciAgent-toolkit/skills/` (all skill.md files)
- **Roles**: `SciAgent-toolkit/roles/*.yaml` (YAML specs referencing agents/skills)
- **Guidelines**: `SciAgent-toolkit/docs/guidelines/`
- **API keys**: `.env` files (gitignored)

**Clear Ownership:**
- scbio-docker = container infrastructure, directory structure (NO AI)
- SciAgent-toolkit = ALL AI files, methodology, MCP configuration, roles
