# Branch Management Vignette: dev + dev-claude-integration

**Purpose:** Guide for maintaining the two-branch architecture in scbio-docker
**Audience:** Repository maintainers
**Last Updated:** 2025-11-11

---

## Branch Architecture Overview

This repository uses a **stacked branch strategy** with graceful degradation:

```
main (production releases)
 │
 ├─ dev (core functionality)
 │   └─ dev-claude-integration (AI workflow extension)
```

### Branch Purposes

| Branch | Purpose | Contains | Users |
|--------|---------|----------|-------|
| `main` | Stable releases | Tagged versions (v0.5.1, v0.5.2, etc.) | Production users |
| `dev` | Active development | Latest features, Docker enhancements | All users |
| `dev-claude-integration` | AI workflow layer | dev + Claude Code integration files | Claude Code users |
| `dev-gpt-codex-integration` | AI workflow layer | dev + GPT-Codex integration files | Codex CLI users |

---

## Key Design Principle: Conditional Files

The `init-project.sh` script uses **branch-aware conditional logic**:

```bash
# Lines 224-270 in init-project.sh
if [ -d "${TEMPLATES_DIR}/claude" ]; then
    echo "Setting up Claude Code integration..."
fi

if [ -d "${TEMPLATES_DIR}/gpt-codex" ]; then
    echo "Setting up GPT-Codex integration..."
fi
```

**Result:**
- On `dev`: Script skips AI files (directories don’t exist)
- On `dev-claude-integration`: Script includes Claude files (`templates/claude/`)
- On `dev-gpt-codex-integration`: Script includes GPT-Codex files (`templates/gpt-codex/`)
- **Same script works on all branches** without modification

**New in v0.5.2:**
- `init-project.sh` exposes `--ai {none,claude,codex,both}` (and an interactive prompt). Selecting Codex/Claude copies the matching templates *if they exist* on the current branch and now also generates `.mcp.json` with the correct `/workspaces/<project>` path using `templates/ai-common/mcp.json.template`.
- `.devcontainer/scripts/install_ai_tooling.sh` ships with every scaffolded project and runs from VS Code’s `postCreateCommand` to install the Claude CLI (when applicable), validate Node.js 20.x + `uvx`, and warn if `context7` is enabled without `CONTEXT7_API_KEY`.
- When you change `templates/ai-common/mcp.json.template` or AI scripts, immediately merge dev → `dev-claude-integration` and dev → `dev-gpt-codex-integration` so both branches stay aligned on MCP behavior.

---

## When to Sync Branches

### Trigger: New commits on `dev`

Whenever you commit significant work to `dev`, sync it to both AI branches:

**Examples of sync-worthy changes:**
- Docker image enhancements (v0.5.2 package additions)
- `init-project.sh` improvements
- Template updates in `templates/docs/`
- Documentation updates (CLAUDE.md, GPT-CODEX.md, DEVOPS.md, README.md)
- Build script changes (`build-optimized.sh`)

**NOT sync-worthy:**
- Changes already in the AI branches
- Experimental/WIP commits
- Branch-specific files (`templates/claude/`, `templates/gpt-codex/`)

---

## How to Check Branch Status

### 1. View branch divergence

```bash
# How many commits is dev ahead of dev-claude-integration?
git log dev-claude-integration..dev --oneline

# Example output:
# 34cc7a2 Release v0.5.2 + consolidated QUICK-START.md
# dca19d9 Enhance init-project.sh with interactive configuration
```

If this shows commits, `dev` has changes that need to be merged.

```bash
# How many commits is dev-claude-integration ahead of dev?
git log dev..dev-claude-integration --oneline

# Example output:
# 5e0f14e Document init-project.sh enhancement in plan.md and tasks.md
# 0f89707 Add Claude Code integration extension (stacked on dev branch)
```

If this shows commits, the AI branch has unique work (Claude or GPT-Codex templates, docs, etc.).

### 2. View file differences

```bash
# What files differ between dev and an AI branch?
git diff dev dev-claude-integration --name-status
git diff dev dev-gpt-codex-integration --name-status

# Example output:
# M    init-project.sh
# A    templates/claude/CLAUDE.md.template
# A    templates/gpt-codex/GPT-CODEX.md.template
```

**Expected differences:**
- `templates/claude/` directory (only in dev-claude-integration)
- `templates/gpt-codex/` directory (only in dev-gpt-codex-integration)
- `init-project.sh` contains copy logic for whichever AI templates exist
- `plan.md`, `tasks.md`, AI READMEs may have branch-specific doc sections

**Unexpected differences = branches out of sync**

### 3. Visualize branch history

```bash
git log --oneline --graph --all --decorate -20
```

Look for:
- ✅ Merge commits connecting dev → dev-claude-integration
- ❌ Diverging parallel commits (indicates out of sync)

---

## Step-by-Step: Syncing Branches

### Scenario: You committed v0.5.2 work to `dev`, need to sync to `dev-claude-integration`

#### Step 1: Verify uncommitted work is committed

```bash
git checkout dev
git status

# If you see uncommitted changes, commit them:
git add <files>
git commit -m "Your commit message"
```

#### Step 2: Check what needs merging

```bash
git log dev-claude-integration..dev --oneline
```

If this shows commits, proceed to Step 3.

#### Step 3: Switch to target branch

```bash
git checkout dev-claude-integration
```

#### Step 4: Merge dev into dev-claude-integration

```bash
git merge dev -m "Merge dev into dev-claude-integration: <summary>

<Details about what's being merged>

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

**Expected result:**
```
Merge made by the 'ort' strategy.
 .devcontainer/Dockerfile.optimized |   8 +-
 CLAUDE.md                          |  33 ++-
 HANDOFF.md                         | 533 ++++++++++++++++++
 ...
```

#### Step 5: Verify merge success

```bash
# Check that Claude templates still exist
ls -la templates/claude/
# Should show: CLAUDE.md.template, WORKFLOW.md, .claude/

# Check that v0.5.2 changes are present
grep "v0.5.2" CLAUDE.md | head -3

# View commit graph
git log --oneline --graph --all --decorate -10
# Should show merge commit connecting branches
```

#### Step 6: Test both branches

See "Testing After Merge" section below.

---

### Scenario: You committed core work to `dev`, need to sync to `dev-gpt-codex-integration`

Steps mirror the Claude workflow—just swap branch names:

1. **Ensure `dev` is clean** (`git status`).
2. **Inspect divergence**  
   `git log dev-gpt-codex-integration..dev --oneline`
3. **Checkout AI branch**  
   `git checkout dev-gpt-codex-integration`
4. **Merge dev**  
   `git merge dev -m "Merge dev into dev-gpt-codex-integration: <summary>"`  
   Add a note referencing Codex/GPT work if useful.
5. **Verify GPT assets**  
   - `ls templates/gpt-codex/`  
   - Ensure `init-project.sh` copy block still exists (grep for `gpt-codex`).
6. **Run tests** – scaffold a sample project and confirm `GPT-CODEX.md`, `WORKFLOW-gpt-codex.md`, `.gpt-codex/` appear.

Then continue with the shared testing checklist.

## Handling Merge Conflicts

### Common conflict scenarios

1. **init-project.sh conflicts** (most likely)
   - Cause: Both branches modified the script
   - Resolution: Keep dev-claude-integration version (has Claude logic)
   - Verify: Check that Claude file copying section (lines 224-248) is intact

2. **CLAUDE.md conflicts**
   - Cause: Documentation updates on both branches
   - Resolution: Manual merge, keep both changes
   - Verify: Ensure v0.5.2 section exists AND Claude-specific notes remain

3. **Template conflicts** (unlikely)
   - Cause: Overlapping template changes
   - Resolution: Prefer AI branch for its templates (`templates/claude/` or `templates/gpt-codex/`), dev for `templates/docs/`

### Resolving conflicts

```bash
git merge dev
# Auto-merging init-project.sh
# CONFLICT (content): Merge conflict in init-project.sh

# Open conflicted files
code init-project.sh

# Resolve conflicts (keep dev-claude-integration version with Claude logic)
# Remove conflict markers: <<<<<<<, =======, >>>>>>>

# Stage resolved files
git add init-project.sh

# Complete merge
git commit
```

---

## Testing After Merge

### Test 1: Verify branches have correct files

**On dev branch:**
```bash
git checkout dev
ls templates/claude/ 2>/dev/null && echo "❌ ERROR: Claude templates exist on dev!" || echo "✅ OK: No Claude templates"
ls templates/gpt-codex/ 2>/dev/null && echo "❌ ERROR: GPT-Codex templates exist on dev!" || echo "✅ OK: No GPT-Codex templates"
```

**On dev-claude-integration branch:**
```bash
git checkout dev-claude-integration
ls templates/claude/ && echo "✅ OK: Claude templates exist" || echo "❌ ERROR: Claude templates missing!"
```

**On dev-gpt-codex-integration branch:**
```bash
git checkout dev-gpt-codex-integration
ls templates/gpt-codex/ && echo "✅ OK: GPT-Codex templates exist" || echo "❌ ERROR: GPT-Codex templates missing!"
```

### Test 2: Test init-project.sh on dev branch

```bash
git checkout dev
./init-project.sh /tmp/test-dev basic-rna --git-init

# Verify created files
ls /tmp/test-dev/CLAUDE.md 2>/dev/null && echo "❌ ERROR: CLAUDE.md should not exist" || echo "✅ OK"
ls /tmp/test-dev/GPT-CODEX.md 2>/dev/null && echo "❌ ERROR: GPT-CODEX.md should not exist" || echo "✅ OK"
ls /tmp/test-dev/plan.md || echo "❌ ERROR: plan.md missing"
ls /tmp/test-dev/tasks.md || echo "✅ OK"

# Cleanup
rm -rf /tmp/test-dev
```

### Test 3: Test init-project.sh on dev-claude-integration branch

```bash
git checkout dev-claude-integration
./init-project.sh /tmp/test-claude basic-rna --git-init

# Verify created files
ls /tmp/test-claude/CLAUDE.md || echo "❌ ERROR: CLAUDE.md missing"
ls /tmp/test-claude/WORKFLOW.md || echo "❌ ERROR: WORKFLOW.md missing"
ls /tmp/test-claude/.claude/agents/ || echo "❌ ERROR: .claude/agents missing"
ls /tmp/test-claude/plan.md || echo "✅ OK"

# Cleanup
rm -rf /tmp/test-claude
```

### Test 4: Test init-project.sh on dev-gpt-codex-integration branch

```bash
git checkout dev-gpt-codex-integration
./init-project.sh /tmp/test-gpt basic-rna --git-init

# Verify created files
ls /tmp/test-gpt/GPT-CODEX.md || echo "❌ ERROR: GPT-CODEX.md missing"
ls /tmp/test-gpt/WORKFLOW-gpt-codex.md || echo "❌ ERROR: WORKFLOW-gpt-codex.md missing"
ls /tmp/test-gpt/.gpt-codex/agents/ || echo "❌ ERROR: .gpt-codex/agents missing"

# Cleanup
rm -rf /tmp/test-gpt
```

### Test 5: Verify Docker build still works

```bash
git checkout dev
./build-optimized.sh --tag scdock-r-dev:v0.5.2-test

# Should complete without errors
# Image size should be ~20-25GB
```

---

## Common Pitfalls and Solutions

### Pitfall 1: Forgot to commit before merging

**Symptom:** Merge brings in uncommitted changes from dev
```bash
git status
# Shows modified files not related to merge
```

**Solution:**
```bash
git merge --abort  # Cancel merge
git checkout dev
git add <files>
git commit -m "Commit message"
git checkout dev-claude-integration
git merge dev  # Retry
```

### Pitfall 2: Accidentally modified Claude files on dev

**Symptom:** `templates/claude/` exists on dev branch
```bash
git checkout dev
ls templates/claude/  # This should NOT exist
```

**Solution:**
```bash
git checkout dev
git rm -rf templates/claude/
git commit -m "Remove Claude templates from dev branch (belong in dev-claude-integration)"
```

### Pitfall 3: init-project.sh not working after merge

**Symptom:** Script fails with "templates/claude not found" on dev-claude-integration

**Solution:**
```bash
git checkout dev-claude-integration
ls templates/claude/  # Should exist

# If missing:
git checkout origin/dev-claude-integration -- templates/claude/
git commit -m "Restore Claude templates"
```

### Pitfall 4: Pushed to wrong remote

**Symptom:** Accidentally pushed dev-claude-integration changes to main

**Solution:**
```bash
# Revert remote branch (dangerous, use with caution)
git push origin <correct-branch-name> --force-with-lease

# Better: Create new branch from correct state
git checkout main
git checkout -b main-backup  # Save current state
git reset --hard origin/main  # Reset to remote
git push origin main --force-with-lease
```

---

## Branch Release Workflow

### Releasing to main (production)

When `dev` is stable and ready for release:

```bash
# 1. Ensure dev is clean and tested
git checkout dev
git status  # Should be clean

# 2. Merge dev into main
git checkout main
git merge dev -m "Merge dev into main: Release v0.5.2"

# 3. Tag the release
git tag -a v0.5.2 -m "Release v0.5.2: Enhanced package set

New capabilities:
- Chromatin accessibility analysis (chromVAR, motifmatchr)
- Motif analysis (TFBSTools, JASPAR2022)
- Cell type annotation (SingleR, celldex)
- Runtime package compilation (libgsl-dev, libhdf5-dev)

Bug fixes:
- Fixed safe_install() meta-package detection"

# 4. Push to remote
git push origin main
git push origin v0.5.2

# 5. Sync dev-claude-integration with released changes
git checkout dev-claude-integration
git merge main  # Brings in any main-specific changes
```

### Hotfix workflow (emergency fix on main)

If critical bug found in production:

```bash
# 1. Create hotfix branch from main
git checkout main
git checkout -b hotfix/v0.5.2.1

# 2. Fix the bug
# ... make changes ...
git commit -m "Hotfix: Fix critical bug"

# 3. Merge into main
git checkout main
git merge hotfix/v0.5.2.1
git tag -a v0.5.2.1 -m "Hotfix: Description"
git push origin main v0.5.2.1

# 4. Merge into dev (so fix doesn't get lost)
git checkout dev
git merge hotfix/v0.5.2.1
git push origin dev

# 5. Merge into dev-claude-integration
git checkout dev-claude-integration
git merge dev
git push origin dev-claude-integration

# 6. Delete hotfix branch
git branch -d hotfix/v0.5.2.1
```

---

## Cheat Sheet: Quick Commands

### Check branch status
```bash
# What's on dev that's not on dev-claude-integration?
git log dev-claude-integration..dev --oneline

# What files differ?
git diff dev dev-claude-integration --name-status

# Visual history
git log --oneline --graph --all --decorate -20
```

### Sync branches (fast)
```bash
git checkout dev-claude-integration
git merge dev
git push origin dev-claude-integration
```

### Test init-project.sh (both branches)
```bash
# dev branch test
git checkout dev
./init-project.sh /tmp/test-dev basic-rna && ls /tmp/test-dev/

# dev-claude-integration test
git checkout dev-claude-integration
./init-project.sh /tmp/test-claude basic-rna && ls /tmp/test-claude/CLAUDE.md

# Cleanup
rm -rf /tmp/test-{dev,claude}
```

### Undo a merge (if caught immediately)
```bash
git merge --abort  # Before committing
git reset --hard HEAD~1  # After committing (dangerous!)
```

---

## FAQ

### Q: Why two AI branches instead of one with feature flags?

**A:** Separation of concerns:
- `dev` = Universal core (anyone can use)
- `dev-claude-integration` = Opt-in AI layer (Claude users only)
- `dev-gpt-codex-integration` = Opt-in AI layer (Codex users only)

This allows non-AI users to clone/use `dev` without seeing AI-specific files, and keeps the repo structure clean.

### Q: Can I work directly on dev-claude-integration / dev-gpt-codex-integration?

**A:** Yes, for branch-specific changes (`templates/claude/` or `templates/gpt-codex/`, agent stubs, AI docs). For core changes (Docker, R packages, init-project.sh core logic), work on `dev` first, then merge.

### Q: How often should I sync branches?

**A:** After every significant commit to `dev`. Daily if actively developing.

### Q: What if I forget to sync?

**A:** Branches will diverge. Not fatal, but makes merging harder. Use `git log` to identify divergence, then merge.

### Q: Can I delete dev-claude-integration and recreate it?

**A:** Technically yes, but not recommended. Better to merge regularly to preserve commit history.

### Q: Should I squash merges?

**A:** No. Keep merge commits (`--no-ff` is default). This preserves branch history and makes it clear where features were integrated.

---

## Summary Checklist

After every dev commit, run this checklist:

- [ ] Committed all changes on `dev`
- [ ] Checked `git log dev-claude-integration..dev` to see what needs merging
- [ ] Switched to `dev-claude-integration`
- [ ] Merged `dev` into `dev-claude-integration`
- [ ] Verified `templates/claude/` still exists on `dev-claude-integration`
- [ ] Verified `templates/claude/` does NOT exist on `dev`
- [ ] Checked `git log dev-gpt-codex-integration..dev`
- [ ] Switched to `dev-gpt-codex-integration`
- [ ] Merged `dev` into `dev-gpt-codex-integration`
- [ ] Verified `templates/gpt-codex/` still exists on `dev-gpt-codex-integration`
- [ ] Verified `templates/gpt-codex/` does NOT exist on `dev`
- [ ] Tested `init-project.sh` on both branches
- [ ] Pushed both branches to remote if needed

---

**Maintainer Notes:**
- This vignette should be updated when branch strategy changes
- Keep examples up-to-date with actual repo state
- Add new pitfalls as they're discovered

**Last Review:** 2025-11-11 (v0.5.2 release)
