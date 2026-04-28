# Template System Guide

This guide explains how the scbio-docker template system works and how to maintain it.

## Overview

The `init-project.sh` script initializes new bioinformatics projects using a **hybrid template approach**:
- **Static templates** for configuration files and documentation
- **Generated files** for devcontainer infrastructure (devcontainer.json, docker-compose.yml, .env)
- **Copied resources** for helper scripts and documentation

## Template Directory Structure

```
scbio-docker/
├── init-project.sh                  # Main initialization script
└── templates/
    ├── base/                        # Base template (default)
    │   └── README.md                # Template-specific docs
    ├── config/                      # Configuration templates
    │   ├── config.R.template
    │   ├── pipeline.yaml.template
    │   └── color_config.R.template
    ├── docs/                        # Documentation templates
    │   ├── README.md.template
    │   ├── tasks.md.template
    │   ├── notes.md.template
    │   └── .env.example
    ├── devcontainer/                # Devcontainer resources
    │   ├── MCP_AUTH_SETUP.md
    │   ├── PYTHON_VENV_GUIDE.md
    │   └── scripts/
    │       ├── setup_claude_mcp.sh
    │       └── source_env.sh
    └── .vscode/
        └── settings.json
```

## Template Substitution Variables

Templates use `{{VARIABLE}}` syntax for substitution:

| Variable | Description | Example |
|----------|-------------|---------|
| `{{PROJECT_NAME}}` | Project directory basename | DC_cancer |
| `{{PROJECT_PATH}}` | Full project path | /scratch/.../DC_cancer |
| `{{DATE}}` | Current date | 2026-02-13 |
| `{{TEMPLATE_TYPE}}` | Template name | base |
| `{{IMAGE_VERSION}}` | Docker image version | scdock-r-dev:v0.5.2 |
| `{{SPECIES}}` | Species name | Mus musculus |
| `{{SPECIES_DB}}` | Species database code | MM |
| `{{GENOME_BUILD}}` | Genome build | mm10 |
| `{{SCBIO_DOCKER_PATH}}` | Path to scbio-docker | /data1/.../scbio-docker |

## What Gets Created

### Generated Files (Programmatic)
These are **generated** by init-project.sh, NOT copied from templates:
- `.devcontainer/devcontainer.json` - VS Code devcontainer config
- `.devcontainer/docker-compose.yml` - Docker compose services
- `.devcontainer/.env` - Environment variables (UID, GID, paths, resource limits)
- `.gitignore` - Git ignore rules

**Why generated?** Allows dynamic customization:
- Project-specific paths and names
- User-specific UID/GID
- Data mount configurations
- Resource limit settings

### Copied and Substituted
These are **copied from templates/** with variable substitution:
- `02_analysis/config/config.R` ← `templates/config/config.R.template`
- `02_analysis/config/pipeline.yaml` ← `templates/config/pipeline.yaml.template`
- `02_analysis/config/color_config.R` ← `templates/config/color_config.R.template`
- `README.md` ← `templates/docs/README.md.template`
- `tasks.md` ← `templates/docs/tasks.md.template`
- `notes.md` ← `templates/docs/notes.md.template`
- `.env.example` ← `templates/docs/.env.example`
- `.vscode/settings.json` ← `templates/.vscode/settings.json`

### Copied As-Is
These are **copied without substitution**:
- `.devcontainer/MCP_AUTH_SETUP.md` ← `templates/devcontainer/MCP_AUTH_SETUP.md`
- `.devcontainer/PYTHON_VENV_GUIDE.md` ← `templates/devcontainer/PYTHON_VENV_GUIDE.md`
- `.devcontainer/scripts/setup_claude_mcp.sh` ← `templates/devcontainer/scripts/`
- `.devcontainer/scripts/source_env.sh` ← `templates/devcontainer/scripts/`
- `.devcontainer/scripts/poststart_sanity.sh` ← fallback from `scripts/` or generated

### Created Empty
Standard directory structure with `.gitkeep` files:
```
00_data/{raw,processed,references}
01_modules/
02_analysis/{config,helpers}
03_results/{checkpoints,plots,tables}
logs/
```

## AI Context Files (Separate Phase)

AI context files are **NOT created by init-project.sh**. They are created later by:
```bash
./01_modules/SciAgent-toolkit/scripts/setup-ai.sh
```

This creates:
- `CLAUDE.md`, `GEMINI.md`, `AGENTS.md`, `context.md`
- `02_analysis/config/analysis_config.yaml`
- `.mcp.json` (MCP server configuration)
- `.claude/agents/` and `.claude/skills/` (symlinks to toolkit)

**Why separate?** Allows projects without AI tooling and keeps AI infrastructure versioned in SciAgent-toolkit.

## Maintaining Templates

### Adding a New Template File

1. Create the template in appropriate subdirectory:
   ```bash
   cd /data1/users/antonz/pipeline/scbio-docker/templates
   vi config/new_config.yaml.template
   ```

2. Add variable placeholders (`{{PROJECT_NAME}}`, etc.)

3. Update `init-project.sh` to copy the template:
   ```bash
   if [ -f "${TEMPLATES_DIR}/config/new_config.yaml.template" ]; then
       cp "${TEMPLATES_DIR}/config/new_config.yaml.template" \
          "${PROJECT_DIR}/02_analysis/config/new_config.yaml"
       sed -i "s|{{PROJECT_NAME}}|${PROJECT_NAME}|g" \
          "${PROJECT_DIR}/02_analysis/config/new_config.yaml"
   fi
   ```

4. Test with a new project initialization

### Updating Existing Templates

1. **Never edit generated project files** - edit the template source
2. Update template in `templates/` directory
3. Commit changes to scbio-docker repo
4. **Existing projects are not affected** - templates only used at init time
5. Document breaking changes in CHANGELOG.md

### Template Best Practices

#### ✅ DO:
- Use clear, descriptive variable names
- Add comments explaining substituted values
- Include example/placeholder values
- Document required vs optional configuration
- Keep templates minimal - only essentials
- Version control all templates

#### ❌ DON'T:
- Include project-specific data or paths
- Commit deprecated content to templates
- Include API keys or credentials (use placeholders)
- Create nested `deprecated/` folders
- Copy entire project structures as templates

## Deprecation Strategy

### For Templates
- **Remove deprecated content** - templates should be clean and current
- **No `deprecated/` folders** in templates directory
- **Document breaking changes** in git history and CHANGELOG.md

### For Individual Projects
When a project evolves and has deprecated content:

1. **Create project-level `.archive/` directory**:
   ```bash
   mkdir -p .archive/{feature}_superseded_{YYYYMMDD}
   ```

2. **Move deprecated content** with descriptive naming:
   ```bash
   mv .devcontainer/deprecated/ .archive/devcontainer_superseded_20260213/
   ```

3. **Add README.md** to archived content explaining why deprecated

4. **Never nest** `deprecated/deprecated/` folders - flatten when archiving

5. **Do NOT include in new project templates** - only current code

### Example: Clean Project Structure
```
project/
├── .archive/                         # Old code (gitignored)
│   ├── peak_atlas_superseded_20260106/
│   └── b2_5_monolithic_superseded_20260202/
├── .devcontainer/                    # Only current config (no deprecated/)
├── 00_data/
├── 01_modules/
├── 02_analysis/
└── 03_results/
```

## Common Issues and Solutions

### Issue: `.env` has wrong WORKSPACE_FOLDER after merge

**Cause:** User manually copied `.devcontainer/` before running init-project, old .env was preserved.

**Fix (as of 2026-02-13):** init-project.sh now always regenerates .env with correct paths while preserving API keys.

**Manual fix for existing projects:**
```bash
cd .devcontainer
# Edit .env and change WORKSPACE_FOLDER to match project path
# Or use relative path: WORKSPACE_FOLDER=..
```

### Issue: Deprecated folder copied to new project

**Cause:** User manually copied from existing project instead of using init-project.sh.

**Fix:**
```bash
cd project/.devcontainer
rm -rf deprecated/
git add -u && git commit -m "Remove deprecated content"
```

**Prevention:** Always use `init-project.sh`, never manually copy entire `.devcontainer/` from another project.

### Issue: Template references old naming (01_scripts vs 01_modules)

**Cause:** Templates not updated during directory rename migration.

**Fix (completed 2026-02-13):** Updated `config.R.template` to use `DIR_MODULES` instead of `DIR_SCRIPTS`.

### Issue: Templates reference deprecated files (plan.md vs context.md)

**Cause:** AI context file naming changed but templates not updated.

**Fix (completed 2026-02-13):** Updated `README.md.template` to reference `context.md` instead of `plan.md`.

## Testing Changes

After modifying templates or init-project.sh:

1. **Create test project:**
   ```bash
   init-project.sh /tmp/test_project_$$
   ```

2. **Verify structure:**
   ```bash
   cd /tmp/test_project_$$
   tree -L 2 -a
   ```

3. **Check substitutions:**
   ```bash
   grep -r "{{" .  # Should be empty (all variables substituted)
   ```

4. **Verify .devcontainer:**
   ```bash
   cat .devcontainer/.env | grep WORKSPACE_FOLDER  # Should be ".."
   ls -la .devcontainer/  # No deprecated/ folder
   ```

5. **Clean up:**
   ```bash
   rm -rf /tmp/test_project_$$
   ```

## Version History

### 2026-02-13 - Template Modernization
- Fixed `config.R.template`: `DIR_SCRIPTS` → `DIR_MODULES`
- Fixed `README.md.template`: `plan.md` → `context.md`
- Added devcontainer documentation templates (MCP_AUTH_SETUP.md, PYTHON_VENV_GUIDE.md)
- Added devcontainer helper scripts to templates
- Fixed `.env` generation to always use correct WORKSPACE_FOLDER
- Implemented API key preservation when updating existing .env

### Previous Versions
- See git history: `git log --oneline init-project.sh templates/`

## References

- Main script: `init-project.sh`
- Template directory: `templates/`
- SciAgent-toolkit: `toolkits/SciAgent-toolkit/`
- Example projects: `/scratch/current/antonz/projects/`
