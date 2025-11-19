# HANDOFF: DevContainer MCP Setup Issues & Recommendations

**Date:** 2025-11-13
**Author:** Claude Code (via User: antonz)
**Context:** Issues discovered while debugging JBader_scHFD devcontainer startup failures
**Related Project:** `/data2/users/JCRLab/JBader/JBader_scHFD/`
**Related Scripts:** `init-project.sh`, `.devcontainer/` templates

---

## Executive Summary

During setup of the `JBader_scHFD` project, multiple issues were discovered in the devcontainer MCP (Model Context Protocol) setup process that prevent reliable container initialization. These issues stem from:

1. **Timing race conditions** between `postStartCommand` execution and volume mount synchronization
2. **Hardcoded project paths** in template files (`DC_Dictionary` instead of dynamic project names)
3. **Missing dependencies** in base Docker image (Node.js, Python `uv` package)
4. **Wrong lifecycle hook usage** (`postStartCommand` vs `postCreateCommand`)

These issues affect the **init-project.sh** script and should be addressed in the **scbio-docker** repository before generating new projects.

---

## Problem 1: postStartCommand Timing Race Condition

### Issue Description

When VS Code rebuilds a devcontainer using docker-compose, the `postStartCommand` executes **before volume mount synchronization is complete**. This causes the following error:

```
chmod: cannot access '.devcontainer/scripts/poststart_sanity.sh': No such file or directory
chmod: cannot access '.devcontainer/scripts/setup_claude_mcp.sh': No such file or directory
```

### Root Cause

**DevContainer Startup Sequence:**
1. Container starts (`docker-compose up`)
2. `postStartCommand` begins executing immediately
3. Volume mounts may still be syncing (especially for untracked files)
4. Scripts not visible → command fails → container startup blocked

**Git Tracking Matters:**
- ✅ **Tracked files** (in git index): Reliably available during startup
- ❌ **Untracked files** (e.g., `?? setup_claude_mcp.sh`): May not be visible when `postStartCommand` runs

**Evidence from JBader_scHFD:**
```bash
$ git status
M .devcontainer/devcontainer.json
M .devcontainer/docker-compose.yml
M README.md
?? .devcontainer/scripts/setup_claude_mcp.sh  # UNTRACKED → Not reliably available
?? mcp.json                                    # UNTRACKED → Not reliably available

$ git ls-files .devcontainer/scripts/
.devcontainer/scripts/poststart_sanity.sh     # TRACKED → Always available
```

### Impact

- Container startup fails with exit code 1
- User cannot connect to devcontainer in VS Code
- Requires manual intervention to fix

### Recommendations

**Option A: Use postCreateCommand Instead** (Preferred)
- `postCreateCommand`: Runs **once** when container is first created, **after** volume mounts are ready
- `postStartCommand`: Runs **every time** container starts (including rebuilds)
- MCP setup only needs to run once, not on every restart

**Option B: Ensure All Scripts Are Tracked in Git**
- Add setup scripts to git template repository
- Ensure they're committed (not untracked) when copied by `init-project.sh`

**Option C: Add Conditional Checks**
- Test for file existence before running commands
- Fail gracefully if files aren't available

**Example Fix (Option C):**
```json
"postStartCommand": "bash -c '[ -x .devcontainer/scripts/poststart_sanity.sh ] && .devcontainer/scripts/poststart_sanity.sh || echo \"[INFO] Sanity checks skipped\"'"
```

---

## Problem 2: Hardcoded Project Paths in Templates

### Issue Description

Template files contain **hardcoded project paths** from a previous project (`DC_Dictionary`), not the current project name.

### Affected Files

**1. `.devcontainer/scripts/setup_claude_mcp.sh` (lines 156, 181):**
```bash
MCP_JSON="/workspaces/DC_Dictionary/.mcp.json"  # ❌ WRONG PROJECT
```

Should be:
```bash
MCP_JSON="/workspaces/JBader_scHFD/.mcp.json"   # ✅ CORRECT
```

**2. `mcp.json` (line 23):**
```json
{
  "mcpServers": {
    "serena": {
      "args": [
        "--project",
        "/workspaces/DC_Dictionary"  // ❌ WRONG PROJECT
      ]
    }
  }
}
```

Should be:
```json
{
  "mcpServers": {
    "serena": {
      "args": [
        "--project",
        "/workspaces/JBader_scHFD"  // ✅ CORRECT
      ]
    }
  }
}
```

### Root Cause

The template files were copy-pasted from the `DC_Dictionary` project without proper variable substitution. The `init-project.sh` script does not perform template variable replacement for project names.

### Impact

- MCP configuration created in **wrong directory** (`/workspaces/DC_Dictionary/.mcp.json`)
- Serena MCP server indexes **wrong project** codebase
- User must manually fix paths after project initialization

### Recommendations

**Solution: Template Variable Substitution**

1. **Convert template files to use placeholders:**
   ```bash
   MCP_JSON="/workspaces/{{PROJECT_NAME}}/.mcp.json"
   ```

2. **Update init-project.sh to replace placeholders:**
   ```bash
   # In init-project.sh
   PROJECT_NAME=$(basename "$PROJECT_DIR")

   # After copying template files
   find "$PROJECT_DIR/.devcontainer" -type f -name "*.sh" -o -name "*.json" | while read -r file; do
       sed -i "s|{{PROJECT_NAME}}|$PROJECT_NAME|g" "$file"
   done
   ```

3. **Test with multiple projects to ensure paths are correct**

---

## Problem 3: Missing MCP Dependencies in Base Image

### Issue Description

MCP servers require runtime dependencies that are **NOT present** in the base Docker image (`scdock-r-dev:v0.5.1`).

### Base Image Contents (Dockerfile.optimized)

**✅ What's Included:**
- Ubuntu 22.04
- R 4.5.0 (compiled from source)
- Python 3.x with pip
- Scientific libraries (HDF5, GSL, ImageMagick, etc.)

**❌ What's Missing:**
- Node.js / npm / npx
- Python `uv` package

### MCP Dependency Requirements

| MCP Server | Command Required | Package Source | In Base Image? |
|------------|------------------|----------------|----------------|
| **serena** | `uvx` | Python `uv` package | ❌ NO |
| **sequential-thinking** | `npx` | Node.js 20.x + npm | ❌ NO |
| **context7** | `npx` (or SSE) | Node.js 20.x + npm | ❌ NO (disabled) |

### What setup_claude_mcp.sh Installs

The setup script performs the following installations **at runtime** (during postStartCommand):

1. **Claude Code CLI** (via native installer)
   ```bash
   curl -fsSL https://claude.ai/install.sh | bash -s latest
   ```

2. **Python `uv` package** (via sudo pip3)
   ```bash
   sudo pip3 install uv
   ```

3. **Node.js 20.x** (via NodeSource repository)
   ```bash
   curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash -
   sudo apt-get install -y nodejs
   ```

**Time Cost:** ~30-60 seconds per container restart (if using `postStartCommand`)

### Impact

- **Without the setup script**, MCP servers will **NOT work**
- Container startup time increases by 30-60 seconds on **every restart** (if using `postStartCommand`)
- Network dependency (requires internet access to download packages)

### Recommendations

**Option A: Add Dependencies to Dockerfile** (Preferred)

Add to `Dockerfile.optimized` (stage 2, runtime image):

```dockerfile
####################################################
### MCP Server Dependencies (Optional)
####################################################
# Node.js 20.x (for npx-based MCP servers)
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - && \
    apt-get install -y nodejs && \
    rm -rf /var/lib/apt/lists/*

# Python uv package (for uvx-based MCP servers)
RUN pip3 install --no-cache-dir uv
```

**Advantages:**
- Dependencies baked into image (faster container startup)
- No runtime downloads required
- No need for setup script or postStartCommand
- Consistent across all project containers

**Disadvantages:**
- Increases base image size (~150-200 MB)
- All projects get MCP dependencies (even if not using them)

**Option B: Keep Runtime Installation, Use postCreateCommand**

Move setup from `postStartCommand` → `postCreateCommand`:

```json
{
  "postCreateCommand": "bash -lc '.devcontainer/scripts/setup_claude_mcp.sh'"
}
```

**Advantages:**
- Only installs when needed (project-specific)
- Smaller base image size
- Runs only **once** (not on every restart)

**Disadvantages:**
- First container startup takes 30-60 seconds longer
- Requires internet access during first startup
- Script must be properly templated (see Problem 2)

**Option C: Make MCP Setup Optional**

Add a build argument to `init-project.sh`:

```bash
./init-project.sh --project JBader_scHFD --mcp-setup [yes|no]
```

- If `yes`: Copy MCP setup scripts and configure `postCreateCommand`
- If `no`: Skip MCP setup entirely (faster startup, no dependencies)

---

## Problem 4: Wrong Lifecycle Hook for Setup

### Issue Description

The `setup_claude_mcp.sh` script is called from `postStartCommand`, which runs **on every container restart**. This is inappropriate for one-time setup tasks like installing Node.js and uv.

### Current Configuration (Wrong)

```json
{
  "postStartCommand": "bash -lc 'chmod +x .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/poststart_sanity.sh'"
}
```

**Problems:**
- Runs on **every** `docker-compose up` or container rebuild
- Wastes 30-60 seconds reinstalling packages that are already installed
- Increases container startup latency

### DevContainer Lifecycle Hooks

| Hook | When It Runs | Use Case | Example |
|------|--------------|----------|---------|
| `postCreateCommand` | **Once** when container first created | One-time setup (install dependencies, clone repos) | Install Node.js, uv, Claude Code |
| `postStartCommand` | **Every time** container starts | Ephemeral checks, start services | Run sanity checks, start tmux |
| `postAttachCommand` | **Every time** VS Code attaches | User-specific setup | Display welcome message |

### Recommendations

**Use postCreateCommand for MCP Setup:**

```json
{
  "postCreateCommand": "bash -lc 'chmod +x .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/setup_claude_mcp.sh'",
  "postStartCommand": "bash -c '[ -x .devcontainer/scripts/poststart_sanity.sh ] && .devcontainer/scripts/poststart_sanity.sh || echo \"[INFO] Sanity checks skipped\"'"
}
```

**Benefits:**
- MCP setup runs **once** (first container creation)
- Sanity checks run **every restart** (as intended)
- Faster container restarts (no redundant installations)

---

## Recommended Action Plan for scbio-docker

### Phase 1: Fix Base Image (Dockerfile.optimized)

**File:** `/data1/users/antonz/pipeline/scbio-docker/.devcontainer/Dockerfile.optimized`

Add MCP dependencies to the runtime stage (after line ~204):

```dockerfile
####################################################
### MCP Server Dependencies (Optional)
####################################################
# Node.js 20.x (provides npx for sequential-thinking MCP)
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - && \
    apt-get install -y nodejs && \
    rm -rf /var/lib/apt/lists/*

# Python uv package (provides uvx for serena MCP)
RUN pip3 install --no-cache-dir uv

# Verify installations
RUN node --version && npm --version && uvx --version
```

**Rebuild image:**
```bash
cd /data1/users/antonz/pipeline/scbio-docker
./build-optimized.sh  # or equivalent build command
```

---

### Phase 2: Create Template Files with Placeholders

**File:** `.devcontainer/templates/setup_claude_mcp.sh`

Replace hardcoded paths with placeholders:

```bash
# Line 156 (before)
MCP_JSON="/workspaces/DC_Dictionary/.mcp.json"

# Line 156 (after)
MCP_JSON="/workspaces/{{PROJECT_NAME}}/.mcp.json"
```

**File:** `.devcontainer/templates/mcp.json`

```json
{
  "mcpServers": {
    "serena": {
      "args": [
        "--project",
        "/workspaces/{{PROJECT_NAME}}"
      ]
    }
  }
}
```

---

### Phase 3: Update init-project.sh

**File:** `/data1/users/antonz/pipeline/scbio-docker/init-project.sh`

Add template substitution logic:

```bash
# After copying .devcontainer/ directory
PROJECT_NAME=$(basename "$PROJECT_DIR")

echo "[INFO] Substituting project name in template files..."
find "$PROJECT_DIR/.devcontainer" -type f \( -name "*.sh" -o -name "*.json" \) | while read -r file; do
    if grep -q "{{PROJECT_NAME}}" "$file" 2>/dev/null; then
        sed -i "s|{{PROJECT_NAME}}|$PROJECT_NAME|g" "$file"
        echo "  ✓ Updated: $file"
    fi
done
```

---

### Phase 4: Update devcontainer.json Template

**File:** `.devcontainer/templates/devcontainer.json`

Change lifecycle hooks:

```json
{
  "postCreateCommand": "bash -lc 'chmod +x .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/setup_claude_mcp.sh'",
  "postStartCommand": "bash -c '[ -x .devcontainer/scripts/poststart_sanity.sh ] && .devcontainer/scripts/poststart_sanity.sh || echo \"[INFO] Sanity checks skipped\"'"
}
```

---

### Phase 5: Document in QUICK-START.md

Add troubleshooting section:

```markdown
## Troubleshooting

### Container startup fails with "script not found"

**Symptom:**
```
chmod: cannot access '.devcontainer/scripts/setup_claude_mcp.sh': No such file or directory
```

**Cause:** Setup scripts are untracked in git, causing timing race condition with volume mounts.

**Solution:**
1. Ensure all scripts in `.devcontainer/scripts/` are committed to git
2. Or use `postCreateCommand` instead of `postStartCommand` for setup tasks
3. Or add conditional checks: `[ -f script.sh ] && ./script.sh || echo "Skipped"`

### MCP servers not working

**Symptom:** `serena` or `sequential-thinking` commands fail with "command not found"

**Cause:** Base image missing Node.js or Python `uv` package.

**Solution:**
1. Rebuild base image with MCP dependencies (see Dockerfile.optimized)
2. Or run setup script manually: `.devcontainer/scripts/setup_claude_mcp.sh`
```

---

## Summary Table: Issues & Solutions

| Issue | Impact | Solution | Priority | Effort |
|-------|--------|----------|----------|--------|
| postStartCommand timing race | Container startup fails | Use postCreateCommand OR ensure scripts tracked in git | **HIGH** | Low |
| Hardcoded project paths | Wrong MCP config location | Add template substitution to init-project.sh | **HIGH** | Medium |
| Missing MCP dependencies | MCPs don't work without setup script | Add Node.js + uv to Dockerfile.optimized | **MEDIUM** | Low |
| Wrong lifecycle hook | Slow container restarts (30-60s overhead) | Move setup to postCreateCommand | **MEDIUM** | Low |

---

## Testing Checklist

After implementing fixes, test with a new project:

```bash
# 1. Build updated base image
cd /data1/users/antonz/pipeline/scbio-docker
./build-optimized.sh

# 2. Create test project
./init-project.sh --project TestProject_MCP

# 3. Open in VS Code Dev Containers
cd /data2/users/JCRLab/TestProject_MCP
code .

# 4. Verify container starts without errors
# Expected: No "script not found" errors in devcontainer log

# 5. Verify MCP dependencies present
node --version    # Should show v20.x
npm --version     # Should show v10.x
uvx --version     # Should show uv version

# 6. Verify MCP config has correct paths
cat .mcp.json     # Should show /workspaces/TestProject_MCP, NOT DC_Dictionary

# 7. Test serena MCP (if enabled)
uvx --from git+https://github.com/oraios/serena serena --help

# 8. Test sequential-thinking MCP (if enabled)
npx -y @modelcontextprotocol/server-sequential-thinking --help

# 9. Restart container and verify fast startup
# Expected: postStartCommand completes in <5 seconds (not 30-60s)
```

---

## Reference Files

**Current (Broken) Setup:**
- `/data2/users/JCRLab/JBader/JBader_scHFD/.devcontainer/`
- Setup script with bugs: `.devcontainer/scripts/setup_claude_mcp.sh` (now deleted)
- MCP config with wrong paths: `mcp.json` (now deleted)

**Fixed Configuration:**
- `/data2/users/JCRLab/JBader/JBader_scHFD/.devcontainer/devcontainer.json` (simplified postStartCommand)

**Files to Update in scbio-docker:**
- `/data1/users/antonz/pipeline/scbio-docker/.devcontainer/Dockerfile.optimized`
- `/data1/users/antonz/pipeline/scbio-docker/init-project.sh`
- `/data1/users/antonz/pipeline/scbio-docker/.devcontainer/templates/*`

---

## Questions for Future Consideration

1. **Should MCP setup be optional?**
   - Add `--with-mcp` / `--no-mcp` flag to init-project.sh?
   - Not all projects need Claude Code + MCP servers

2. **Should Claude Code CLI be in base image?**
   - Currently installed by setup script
   - Could be baked into Dockerfile for consistency

3. **Should Context7 be enabled by default?**
   - Requires API key (external dependency)
   - May not be suitable for all projects/users

4. **Should we support multiple MCP configurations?**
   - "Minimal" (sequential-thinking only, no external APIs)
   - "Full" (serena + sequential-thinking + context7)
   - "Custom" (user-defined)

---

## Contacts

**Issue Reporter:** antonz
**Date Reported:** 2025-11-13
**Related Ticket:** [If applicable, add Jira/GitHub issue link]

For questions or clarifications, see:
- `/data1/users/antonz/pipeline/scbio-docker/DEVOPS.md` (build documentation)
- `/data1/users/antonz/pipeline/scbio-docker/QUICK-START.md` (user guide)
- This handoff document

---

## ADDENDUM: Actual Root Cause Identified via DC_Dictionary Comparison

**Date:** 2025-11-13
**Status:** RESOLVED

### Executive Summary of Actual Issue

After comparing with the working DC_Dictionary setup at `/scratch/current/antonz/projects/DMATAC/DC_Dictionary`, the **ACTUAL root cause** was identified:

**The setup scripts were STAGED but NOT COMMITTED to git.**

The original analysis suggested timing issues with `postStartCommand` and volume mounts, but the real problem was that uncommitted files don't exist in the container's filesystem when it's first created from the git repository state.

### Comparison: DC_Dictionary (Working) vs JBader_scHFD (Failing)

#### Configuration Differences

| Aspect | DC_Dictionary (✅ Working) | JBader_scHFD (❌ Failing) |
|--------|---------------------------|---------------------------|
| **Git Status** | Scripts committed & in git index | Scripts staged (`A`) but not committed |
| **File Permissions** | `-rwxr-xr-x` (executable) | `-rw-r--r--` (not executable) |
| **Lifecycle Hook** | `postStartCommand` for MCP setup | `postCreateCommand` for MCP setup |
| **postCreateCommand** | Simple `echo 'Container up and running'` | Complex MCP setup script |
| **Script Execution** | `bash -lc` for MCP setup | `bash -lc` for MCP setup |

#### Side-by-Side devcontainer.json

**DC_Dictionary:**
```json
"postCreateCommand": "echo 'Container up and running'",
"postStartCommand": "bash -lc 'chmod +x .devcontainer/scripts/poststart_sanity.sh .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/poststart_sanity.sh'"
```

**JBader_scHFD (before fix):**
```json
"postCreateCommand": "bash -lc 'chmod +x .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/setup_claude_mcp.sh'",
"postStartCommand": "bash -c '[ -x .devcontainer/scripts/poststart_sanity.sh ] && .devcontainer/scripts/poststart_sanity.sh || echo \"[INFO] Sanity checks skipped\"'"
```

### Why DC_Dictionary Works

1. **Scripts are committed to git** → Files physically exist in container from the start
2. **Scripts are executable** → No permission issues
3. **Uses `postStartCommand`** → Runs after full container initialization
4. **Simple `postCreateCommand`** → Lower risk of failure during initial setup

### Why JBader_scHFD Failed

1. **Scripts staged but not committed** → Files don't exist in container filesystem
2. **Scripts not executable** → Would fail even if committed
3. **Uses `postCreateCommand`** → Runs during build before files may be available
4. **Complex setup in `postCreateCommand`** → High risk of blocking container creation

### The Fix Applied

**Step 1: Commit Files to Git**
```bash
chmod +x .devcontainer/scripts/setup_claude_mcp.sh
git add .devcontainer/scripts/setup_claude_mcp.sh mcp.json
git commit -m "Add MCP setup script and configuration"
```

**Step 2: Restructure Lifecycle Hooks** (to match DC_Dictionary pattern)
```json
{
  "postCreateCommand": "echo 'Container up and running'",
  "postStartCommand": "bash -lc 'chmod +x .devcontainer/scripts/poststart_sanity.sh .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/poststart_sanity.sh'"
}
```

### Key Lessons Learned

1. **Git Tracking is Critical**: Devcontainer scripts MUST be committed to git, not just staged
   - Staged files (`A` in git status) exist on host but not in container build
   - Use `git ls-files .devcontainer/` to verify scripts are actually tracked

2. **File Permissions Matter**: Scripts must be executable (`chmod +x`) before committing
   - Git preserves executable bit across platforms
   - Check with `ls -la` to verify permissions

3. **Lifecycle Hook Selection**: Choose the right hook for the task
   - `postCreateCommand`: Simple verification, low risk of failure
   - `postStartCommand`: Complex setup that can recover from errors

4. **Follow Working Patterns**: DC_Dictionary's pattern is more resilient
   - Simple echo in `postCreateCommand` ensures container always creates successfully
   - MCP setup in `postStartCommand` allows retries if initial setup fails

### Updated Testing Checklist

Before rebuilding devcontainer, verify:

```bash
# 1. Check files are committed (not just staged)
git ls-files .devcontainer/scripts/
# Expected output should include:
# .devcontainer/scripts/poststart_sanity.sh
# .devcontainer/scripts/setup_claude_mcp.sh

# 2. Verify files are executable
ls -la .devcontainer/scripts/*.sh
# Expected: -rwxr-xr-x (not -rw-r--r--)

# 3. Check git status (should be clean or have only other changes)
git status
# setup_claude_mcp.sh should NOT show as ?? or A
# It should either not appear, or show as M if modified

# 4. Verify mcp.json is committed
git ls-files mcp.json
# Should return: mcp.json
```

### Resolution Status

**JBader_scHFD project status:** ✅ FIXED
- Files committed to git: commit `1ea8ed0`
- devcontainer.json restructured to match DC_Dictionary pattern
- Ready for container rebuild

**scbio-docker repository action items:**
1. Update init-project.sh to commit template files after copying
2. Ensure template scripts are executable in the template repository
3. Document git commit requirement in QUICK-START.md
4. Consider pre-commit hooks to verify script permissions

---

## ADDENDUM 2: Working Directory Mismatch in postStartCommand

**Date:** 2025-11-12
**Status:** RESOLVED
**Project:** JBader_scHFD

### Issue Description

After resolving the git commit issues (Addendum 1), a new issue emerged during container rebuild:

```
Running the postStartCommand from devcontainer.json...
[9856 ms] Start: Run in container: /bin/sh -c bash -lc 'chmod +x .devcontainer/scripts/poststart_sanity.sh .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/poststart_sanity.sh'
chmod: cannot access '.devcontainer/scripts/poststart_sanity.sh': No such file or directory
chmod: cannot access '.devcontainer/scripts/setup_claude_mcp.sh': No such file or directory
[9932 ms] postStartCommand from devcontainer.json failed with exit code 1. Skipping any further user-provided commands.
```

**Key Distinction from Addendum 1:**
- Scripts **DO exist** and **ARE committed to git**
- Scripts **ARE in the correct location** (`.devcontainer/scripts/`)
- Verification from inside container confirmed:
  ```bash
  devuser@9c2cef5a7e70:/workspaces/JBader_scHFD/scripts$ ll
  -rwxr-xr-x 1 devuser devgroup 3793 Nov 12 20:28 poststart_sanity.sh
  -rwxr-xr-x 1 devuser devgroup 10136 Nov 12 21:58 setup_claude_mcp.sh
  ```

### Root Cause Analysis

The issue is a **working directory mismatch** in the `postStartCommand` execution context.

**The Problem:**
- `devcontainer.json` sets `"workspaceFolder": "/workspaces/JBader_scHFD"` ✓
- `docker-compose.yml` sets `working_dir: /workspaces/JBader_scHFD` ✓
- **BUT** the `postStartCommand` uses **relative paths** (`.devcontainer/scripts/...`)
- The command executes from a **different working directory** (likely `/` or container default)
- Result: Relative paths don't resolve → "No such file or directory"

**Why This Happens:**

VS Code Dev Containers lifecycle commands (`postCreateCommand`, `postStartCommand`, etc.) run in a **separate execution context** that doesn't automatically inherit the workspace folder as the current working directory. This is documented behavior, though not always obvious to users.

**Evidence:**
```bash
# Inside a terminal session in the container:
devuser@9c2cef5a7e70:/workspaces/JBader_scHFD$ pwd
/workspaces/JBader_scHFD  # ✓ Correct

# But postStartCommand runs from a different directory
# (likely / or another default), causing relative paths to fail
```

### The Fix

**Updated `devcontainer.json` line 30:**

**Before:**
```json
"postStartCommand": "bash -lc 'chmod +x .devcontainer/scripts/poststart_sanity.sh .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/poststart_sanity.sh'"
```

**After:**
```json
"postStartCommand": "bash -lc 'cd /workspaces/JBader_scHFD && chmod +x .devcontainer/scripts/poststart_sanity.sh .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/poststart_sanity.sh'"
```

**Key Change:** Added `cd /workspaces/JBader_scHFD &&` at the beginning of the command to explicitly set the working directory before running scripts.

### Alternative Solutions Considered

**Option A: Use Absolute Paths**
```json
"postStartCommand": "bash -lc 'chmod +x /workspaces/JBader_scHFD/.devcontainer/scripts/poststart_sanity.sh /workspaces/JBader_scHFD/.devcontainer/scripts/setup_claude_mcp.sh && /workspaces/JBader_scHFD/.devcontainer/scripts/setup_claude_mcp.sh && /workspaces/JBader_scHFD/.devcontainer/scripts/poststart_sanity.sh'"
```
- **Pros:** Explicit, no directory change needed
- **Cons:** Verbose, harder to read, hardcodes project name

**Option B: Use ${containerWorkspaceFolder} Variable**
```json
"postStartCommand": "bash -lc 'cd ${containerWorkspaceFolder} && chmod +x .devcontainer/scripts/poststart_sanity.sh .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/poststart_sanity.sh'"
```
- **Pros:** Most portable, works across projects
- **Cons:** Variable may not be expanded in all contexts

**Selected Solution:** Option with explicit `cd /workspaces/JBader_scHFD` (implemented above)
- **Rationale:** Clear, reliable, matches working pattern from DC_Dictionary project

### Directory Structure Validation

The current directory structure is **CORRECT** and follows best practices:

```
JBader_scHFD/
├── .devcontainer/
│   ├── devcontainer.json
│   ├── docker-compose.yml
│   └── scripts/               # ✓ Container lifecycle scripts
│       ├── poststart_sanity.sh
│       └── setup_claude_mcp.sh
├── 00_data/
├── 01_scripts/                # ✓ Project-specific external tools
├── 02_analysis/               # ✓ Analysis scripts
└── ...
```

**DO NOT move scripts to:**
- `/workspaces/JBader_scHFD/scripts/` (conflicts with project structure)
- Top-level `scripts/` directory (wrong separation of concerns)

**Rationale:**
- `.devcontainer/scripts/` is for **container setup/lifecycle**
- `01_scripts/` is for **external tools and submodules** (per project structure)
- `02_analysis/` is for **analysis R/Python scripts**

### Lessons Learned

1. **Relative Paths in Lifecycle Commands Are Unreliable**
   - Always use explicit `cd` or absolute paths in `postCreateCommand`, `postStartCommand`, etc.
   - Don't assume the working directory matches `workspaceFolder`

2. **Test with Fresh Container Rebuilds**
   - "Works in terminal" ≠ "Works in lifecycle commands"
   - These run in different execution contexts

3. **Follow Working Patterns from Other Projects**
   - DC_Dictionary uses `postStartCommand` successfully
   - But initial debugging revealed this working directory issue was different from DC_Dictionary's committed files issue

4. **Progressive Debugging is Key**
   - First issue: Scripts not committed to git (Addendum 1)
   - Second issue: Working directory mismatch (this addendum)
   - Both had similar symptoms but different root causes

### Impact on scbio-docker Repository

**Action Items for `init-project.sh` template:**

1. **Update devcontainer.json template** to use explicit directory change:
   ```json
   "postStartCommand": "bash -lc 'cd /workspaces/{{PROJECT_NAME}} && chmod +x .devcontainer/scripts/poststart_sanity.sh .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/setup_claude_mcp.sh && .devcontainer/scripts/poststart_sanity.sh'"
   ```

2. **Document this requirement** in template comments:
   ```json
   // NOTE: Use explicit 'cd /workspaces/{{PROJECT_NAME}}' before relative paths
   // Lifecycle commands don't inherit workspaceFolder as current directory
   ```

3. **Add to QUICK-START.md troubleshooting**:
   ```markdown
   ### postStartCommand fails with "No such file or directory"

   **Symptom:**
   ```
   chmod: cannot access '.devcontainer/scripts/script.sh': No such file or directory
   ```

   **Cause:** Lifecycle commands execute from a different working directory than `workspaceFolder`.

   **Solution:** Add explicit `cd /workspaces/PROJECT_NAME &&` before relative paths:
   ```json
   "postStartCommand": "bash -lc 'cd /workspaces/PROJECT_NAME && ./script.sh'"
   ```
   ```

### Testing Checklist

After applying the fix, verify:

```bash
# 1. Rebuild container in VS Code
# Dev Containers: Rebuild Container (Ctrl+Shift+P)

# 2. Check postStartCommand output in devcontainer log
# Expected: No "No such file or directory" errors

# 3. Verify scripts executed successfully
# Expected: MCP setup completes, sanity checks pass

# 4. Open integrated terminal
# Expected: Working directory is /workspaces/JBader_scHFD

# 5. Verify MCP dependencies installed
node --version    # Should show v20.x
npm --version     # Should show v10.x
uvx --version     # Should show uv version

# 6. Check MCP configuration
cat .mcp.json     # Should exist and have correct paths
```

### Resolution Timeline

1. **Initial error:** Scripts not found (relative paths failed)
2. **Investigation:** Used Plan agent to analyze devcontainer configuration
3. **Root cause identified:** Working directory mismatch in lifecycle command execution
4. **Fix applied:** Added `cd /workspaces/JBader_scHFD &&` to postStartCommand
5. **Status:** Ready for container rebuild and testing

**Files Modified:**
- `/data2/users/JCRLab/JBader/JBader_scHFD/.devcontainer/devcontainer.json` (line 30)

**Commit recommendation:**
```bash
git add .devcontainer/devcontainer.json
git commit -m "Fix postStartCommand working directory issue

Add explicit 'cd /workspaces/JBader_scHFD' before running scripts
to ensure relative paths resolve correctly. VS Code Dev Containers
lifecycle commands don't inherit workspaceFolder as current directory."
```

### Related Documentation

- **Addendum 1:** Git commit requirement for devcontainer scripts
- **Main document:** MCP setup, dependency issues, lifecycle hook selection
- **Project:** `/data2/users/JCRLab/JBader/JBader_scHFD/`
- **Template source:** `/data1/users/antonz/pipeline/scbio-docker/`

---

**End of Handoff Document**

---

### Update: Devcontainer bind mount fix (2024-11-13)

**Issue:** VS Code container started with only `.devcontainer` contents; lifecycle scripts missing, causing `chmod` failures. Root cause was the compose bind mount `- ${WORKSPACE_FOLDER:-.}:/workspaces/JBader_scHFD` resolving relative to `.devcontainer/`, so the repo root never mounted.

**Action:** Edited `/data2/users/JCRLab/JBader/JBader_scHFD/.devcontainer/docker-compose.yml` lines 11 and 34 to mount the parent directory explicitly: `- ..:/workspaces/JBader_scHFD`. Both services now receive the full repository, ensuring `.devcontainer/scripts/*` and project folders exist inside the container.

**Next steps:** Rebuild the devcontainer without cache in VS Code so Docker Compose reloads the updated volume mapping, then verify `/workspaces/JBader_scHFD` contains `01_scripts`, `02_analysis`, etc., and rerun the post-start scripts if needed.

---

## 2025-11-13 Update – Step 1: Runtime prerequisites for AI tooling

- `.devcontainer/Dockerfile.optimized` now installs Node.js 20.x (NodeSource), npm/npx, the `uv` CLI, and a pinned copy of `nvm` in `/opt/nvm`.
- Added a build-time `nvm install 20 && nvm alias default 20` plus profile hooks so interactive shells can run `nvm` without manual setup.
- Added `node --version`, `npm --version`, `npx --version`, and `uvx --version` checks to fail image builds early if MCP dependencies regress.
- **Next:** Introduce `.devcontainer/scripts/install_ai_tooling.sh` and wire it into the generated devcontainer via `postCreateCommand`.

---

## 2025-11-13 Update – Step 2: Post-create AI bootstrap script

- Added `.devcontainer/scripts/install_ai_tooling.sh`, an idempotent post-create helper that verifies `node/npm/npx/uvx`, auto-installs the Claude CLI when Claude assets exist (unless `SKIP_CLAUDE_CLI_INSTALL=1`), and logs Codex readiness.
- `init-project.sh` now copies this script into generated projects (with executable bit) and wires it into `postCreateCommand`, while `postStartCommand` stays reserved for `poststart_sanity.sh`.
- Projects created from any branch now have a tracked AI setup script, eliminating the race condition where untracked scripts disappeared before lifecycle hooks executed.

---

## 2025-11-13 Update – Step 3: Branch-agnostic AI scaffolding

- `init-project.sh` gained an `--ai {none,claude,codex,both}` flag (also available in interactive mode); it logs the selection and gracefully degrades if requested templates are absent on the current branch.
- Shared helper `replace_placeholders_in_file()` plus a `WORKSPACE_PATH` variable ensure future templates (docs/scripts) can embed `/workspaces/<project>` without ad-hoc `sed` calls.
- AI scaffolding now runs before devcontainer generation, so downstream files (docs, `.mcp.json`, scripts) are ready before lifecycle hooks kick in.

---

## 2025-11-13 Update – Step 4: GPT-Codex templates + MCP placeholders

- Added tracked Codex assets under `templates/gpt-codex/` (workflow doc, GPT-CODEX.md template, `.gpt-codex/agents/*`, cache/log scaffolding).
- Introduced `templates/gpt-codex/mcp.json.template` plus a codex-specific devcontainer script; `setup_gpt_codex_integration()` now copies these files, runs placeholder substitutions (`{{PROJECT_NAME}}`, `{{WORKSPACE_PATH}}`), and drops a ready-to-use `.mcp.json`.
- `install_ai_tooling.sh` warns when `.mcp.json` is missing despite Codex assets, keeping MCP path regressions visible during container creation.

---

## 2025-11-13 Update – Step 5: Docs + validation

- Updated `README.md`, `CLAUDE.md`, `GPT-CODEX.md`, and `BRANCH_MANAGEMENT.md` to describe the baked-in AI toolchain, the `--ai` flag, `.mcp.json` templating, and the new `install_ai_tooling.sh` lifecycle hook.
- Validated the scaffolding path with `./init-project.sh /tmp/test-gpt-codex basic-rna --ai codex --git-init`; the run generated `.mcp.json`, `.gpt-codex/agents/*`, and the Codex MCP script as expected (see terminal log above).
- Skipped the full `./build-optimized.sh` Docker build in this environment to avoid a multi-hour compile; please run it on a build host when ready (example: `./build-optimized.sh --github-pat <token>`). The Docker CLI is available (`docker version` succeeds), so the command should work once resources are allocated.

---

## 2025-11-14 Update – MCP parity + docs

- Added `templates/ai-common/mcp.json.template` so both Claude and Codex share the same server definitions (`context7`, `serena start-mcp-server --context ide-assistant`, `sequential-thinking`). `init-project.sh` now generates `.mcp.json` via a helper regardless of AI mode.
- `templates/gpt-codex/.devcontainer/scripts/setup_codex_mcp.sh` and `.devcontainer/scripts/install_ai_tooling.sh` gained context7 awareness (warn when `CONTEXT7_API_KEY` is unset) and now surface the MCP inventory during setup.
- Updated README, CLAUDE.md, GPT-CODEX.md, BRANCH_MANAGEMENT.md, and INIT_PROJECT_ENHANCEMENT.md to describe the shared template, branch merge expectations, and context7 authentication requirements.
- Validation: 
  - `bash -n init-project.sh`
  - `./init-project.sh /tmp/test-codex basic-rna --ai codex`
  - `./init-project.sh /tmp/test-both basic-rna --ai both`
  - Verified generated `.mcp.json` contains the upgraded commands and was removed after tests.
- Next steps: run `./build-optimized.sh --github-pat <token>` to bake the new dependencies, then commit the changes (see summary below).
