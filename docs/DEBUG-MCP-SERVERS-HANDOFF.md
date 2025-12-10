# MCP Servers Not Detected by Claude Code - Debug Handoff

**Date:** 2025-12-10
**Status:** UNRESOLVED
**Environment:** VS Code Dev Container with scbio-docker base image

---

## Problem Summary

After running `setup_mcp_infrastructure.sh`, Claude Code reports "No MCP servers configured" even though:
- `.mcp.json` file exists in project root
- MCP servers (sequential-thinking, tooluniverse) are installed and working
- Claude Code CLI is installed and functional

---

## Environment Details

- **Container:** scbio-docker v0.5.2 Dev Container
- **Project path:** `/workspaces/12868-EH`
- **SciAgent-toolkit path:** `/workspaces/12868-EH/01_scripts/SciAgent-toolkit` (submodule)
- **Claude Code version:** Latest (installed via `curl -fsSL https://claude.ai/install.sh`)
- **Node.js:** Available via nvm at `/opt/nvm/`

---

## Issues Discovered and Fixed

### Issue 1: npm EACCES Permission Denied

**Symptom:**
```
npm error code EACCES
npm error syscall mkdir
npm error path /opt/nvm/versions/node/v22.16.0/lib/node_modules/@anthropic-ai
```

**Root Cause:** Container's `/opt/nvm/` is read-only. npm trying to install global packages there.

**Fix Applied:** Added `configure_npm_prefix()` to `scripts/common.sh`:
```bash
configure_npm_prefix() {
    local npm_prefix="$HOME/.npm-global"
    mkdir -p "$npm_prefix/bin" "$npm_prefix/lib"
    npm config set prefix "$npm_prefix"
    export PATH="$npm_prefix/bin:$PATH"
}
```

**File:** `toolkits/SciAgent-toolkit/scripts/common.sh` (lines 103-136)

---

### Issue 2: Claude Code Installation Hanging

**Symptom:** Installation hangs indefinitely at "Installing Claude Code native build latest..."

**Root Cause:** Network issues or large download with no timeout.

**Fix Applied:** Added timeout wrapper in `scripts/install_claude.sh`:
```bash
INSTALL_TIMEOUT=300
timeout $INSTALL_TIMEOUT bash -c 'curl -fsSL https://claude.ai/install.sh | bash -s latest'
```

**File:** `toolkits/SciAgent-toolkit/scripts/install_claude.sh` (lines 50-56)

---

### Issue 3: `claude doctor` Hanging

**Symptom:** After successful install, `claude doctor` hangs indefinitely in container environment.

**Root Cause:** `claude doctor` performs checks that may hang in non-interactive/container environments.

**Fix Applied:** Removed `claude doctor` call entirely:
```bash
# Skip diagnostics - claude doctor can hang in container environments
log_ok "Skipping diagnostics (run 'claude doctor' manually if needed)"
```

**File:** `toolkits/SciAgent-toolkit/scripts/install_claude.sh` (lines 75-77)

---

### Issue 4: PATH Not Propagating from Subshells

**Symptom:** `claude` command not found after installation, even though install_claude.sh exports PATH.

**Root Cause:** Child bash scripts export PATH but it doesn't propagate to parent process.

**Fix Applied:** Added explicit PATH exports in orchestrator after each install:
```bash
# After Claude install:
export PATH="$HOME/.local/bin:$PATH"

# After Codex/Gemini install:
export PATH="$HOME/.npm-global/bin:$PATH"
```

**File:** `toolkits/SciAgent-toolkit/scripts/setup_mcp_infrastructure.sh` (lines 159-160, 179-180, 199-200)

---

### Issue 5: `.mcp.json` Created in Wrong Directory

**Symptom:** `.mcp.json` created in SciAgent-toolkit directory instead of user's project root.

**Root Cause:** `PROJECT_DIR` was set to `${SCRIPT_DIR}/..` instead of `${PWD}`.

**Fix Applied:** Changed PROJECT_DIR to use PWD:
```bash
# Use PWD as project root - the directory where user invoked the script
PROJECT_DIR="${PWD}"
```

**File:** `toolkits/SciAgent-toolkit/scripts/setup_mcp_infrastructure.sh` (lines 37-39)

---

### Issue 6: Invalid JSON in `.mcp.json`

**Symptom:** Claude Code can't parse `.mcp.json` - invalid JSON syntax.

**Example of invalid JSON found:**
```json
    }
,
    "tooluniverse": {
```

**Root Cause:** Bash heredoc JSON generation put commas on separate lines.

**Fix Applied:** Rewrote `configure_mcp_servers.sh` to use `claude mcp add` CLI instead of manual JSON generation. Falls back to Python-based JSON generation.

**File:** `toolkits/SciAgent-toolkit/scripts/configure_mcp_servers.sh` (complete rewrite)

---

### Issue 7: Configure Script Skipping Due to Existing File

**Symptom:** Script exits early when `.mcp.json` exists, even if it's invalid.

**Root Cause:** Script checks for file existence but not validity.

**Fix Applied:**
1. Added `--force` mode to delete and regenerate config
2. Orchestrator now passes `--force` by default

**Files:**
- `toolkits/SciAgent-toolkit/scripts/configure_mcp_servers.sh` (lines 70-80)
- `toolkits/SciAgent-toolkit/scripts/setup_mcp_infrastructure.sh` (line 306)

---

### Issue 8: Missing `.claude/settings.json`

**Symptom:** Even with valid `.mcp.json`, Claude Code doesn't enable the servers.

**Root Cause:** Claude Code discovers servers from `.mcp.json` but requires explicit enablement via `enabledMcpjsonServers` in `.claude/settings.json`.

**Fix Applied:** Script now creates `.claude/settings.json`:
```json
{
  "enabledMcpjsonServers": ["sequential-thinking", "tooluniverse"]
}
```

**File:** `toolkits/SciAgent-toolkit/scripts/configure_mcp_servers.sh` (lines 287-321)

---

### Issue 9: Missing `type` Field in Server Config

**Symptom:** Claude Code may not recognize servers without explicit `type` field.

**Root Cause:** Some Claude Code versions require `"type": "stdio"` in server config.

**Fix Applied:** Python fallback now includes type field:
```python
servers["sequential-thinking"] = {
    "type": "stdio",
    "command": "npx",
    "args": ["-y", "@modelcontextprotocol/server-sequential-thinking"]
}
```

**File:** `toolkits/SciAgent-toolkit/scripts/configure_mcp_servers.sh` (Python fallback section)

---

## Current State of Code

### Latest Commits

**SciAgent-toolkit** (commit `b7eece7`):
- `configure_mcp_servers.sh` - Uses `claude mcp add` CLI as primary approach
- `setup_mcp_infrastructure.sh` - Passes `--force` flag

**scbio-docker** (commit `0db2311`):
- Updated submodule reference

### Files Modified

1. `toolkits/SciAgent-toolkit/scripts/common.sh`
   - Added `configure_npm_prefix()` function

2. `toolkits/SciAgent-toolkit/scripts/install_claude.sh`
   - Added 5-minute timeout
   - Removed `claude doctor` call

3. `toolkits/SciAgent-toolkit/scripts/setup_mcp_infrastructure.sh`
   - Changed `PROJECT_DIR` to `${PWD}`
   - Added PATH propagation after CLI installs
   - Passes `--force` to configure script

4. `toolkits/SciAgent-toolkit/scripts/configure_mcp_servers.sh`
   - Complete rewrite to use `claude mcp add` CLI
   - Falls back to Python JSON generation
   - Creates `.claude/settings.json` for server enablement

---

## What Has Been Tested

1. **Syntax validation:** All scripts pass `bash -n` validation
2. **Submodule commits:** Successfully pushed to GitHub
3. **Container testing:** User reports issue still present after pulling latest

---

## What Has NOT Been Tested

1. **Fresh container with no prior config:** Need to test in clean environment
2. **`claude mcp add` command behavior:** May have different flags/behavior than expected
3. **`claude mcp list` output:** Need to verify servers appear after configuration
4. **Interactive `/mcp` command:** Need to verify servers show in Claude Code UI

---

## Remaining Investigation Areas

### 1. Verify `claude mcp add` CLI syntax

The script uses:
```bash
claude mcp add --transport stdio --scope project sequential-thinking -- \
    npx -y @modelcontextprotocol/server-sequential-thinking
```

**Need to verify:**
- Is `--transport stdio` correct? (vs `--type stdio`)
- Is `--scope project` correct? (vs `--scope local` or `-s project`)
- Does `--` separator work correctly?

### 2. Check if `claude mcp add` requires authentication

The CLI may require being logged in before adding servers. Script doesn't handle this.

### 3. Verify `.mcp.json` location expectations

Claude Code may look for `.mcp.json` in:
- Current working directory
- Project root (defined how?)
- `~/.claude/` global config

### 4. Check Claude Code settings file location

The script creates `.claude/settings.json` but Claude Code may expect:
- `.claude/settings.local.json`
- Different structure for `enabledMcpjsonServers`

### 5. MCP server initialization timing

Servers may need to be running before Claude Code can detect them.

---

## Suggested Next Steps

### Step 1: Manual Testing in Container

```bash
# 1. Clean slate
rm -f /workspaces/12868-EH/.mcp.json
rm -rf /workspaces/12868-EH/.claude

# 2. Test claude mcp add manually
cd /workspaces/12868-EH
claude mcp add --help  # Check exact syntax

# 3. Try adding one server manually
claude mcp add sequential-thinking -- npx -y @modelcontextprotocol/server-sequential-thinking

# 4. Check result
claude mcp list
cat .mcp.json
ls -la .claude/
```

### Step 2: Check Claude Code Documentation

- Verify exact syntax for `claude mcp add`
- Verify `.mcp.json` schema requirements
- Verify settings file structure for enabling servers

### Step 3: Debug Script Execution

```bash
# Run with verbose output
bash -x /workspaces/12868-EH/01_scripts/SciAgent-toolkit/scripts/configure_mcp_servers.sh --force 2>&1 | tee /tmp/mcp-config-debug.log
```

### Step 4: Check Claude Code Logs

```bash
# Claude Code may have debug logs
ls -la ~/.claude/
cat ~/.claude/*.log 2>/dev/null
```

---

## Reference: Expected Final State

### `.mcp.json` (project root)

```json
{
  "mcpServers": {
    "sequential-thinking": {
      "type": "stdio",
      "command": "npx",
      "args": ["-y", "@modelcontextprotocol/server-sequential-thinking"]
    },
    "tooluniverse": {
      "type": "stdio",
      "command": "uv",
      "args": ["--directory", "/workspaces/12868-EH/01_scripts/SciAgent-toolkit/scripts/tooluniverse-env", "run", "tooluniverse-mcp", "--exclude-tool-types", "PackageTool"]
    }
  }
}
```

### `.claude/settings.json` (project root)

```json
{
  "enabledMcpjsonServers": ["sequential-thinking", "tooluniverse"]
}
```

---

## Git Log of Fixes

```
b7eece7 fix: Use claude mcp add CLI for robust MCP configuration
dd050bb Fix: Update SciAgent-toolkit submodule for improved npm installations and MCP timeouts
... (earlier commits with partial fixes)
```

---

## Contact / Repository Links

- **scbio-docker:** https://github.com/tony-zhelonkin/scbio-docker (branch: dev)
- **SciAgent-toolkit:** https://github.com/tony-zhelonkin/SciAgent-toolkit (branch: main)

---

*This handoff document was created to facilitate debugging by the next agent or developer.*
