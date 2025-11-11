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

---

## Key Design Principle: Conditional Files

The `init-project.sh` script uses **branch-aware conditional logic**:

```bash
# Lines 224-248 in init-project.sh
if [ -d "${TEMPLATES_DIR}/claude" ]; then
    echo "Setting up Claude Code integration..."
    # Copy CLAUDE.md, WORKFLOW.md, .claude/agents/
fi
```

**Result:**
- On `dev`: Script skips Claude files (directory doesn't exist)
- On `dev-claude-integration`: Script includes Claude files (directory exists)
- **Same script works on both branches** without modification

---

## When to Sync Branches

### Trigger: New commits on `dev`

Whenever you commit significant work to `dev`, sync it to `dev-claude-integration`:

**Examples of sync-worthy changes:**
- Docker image enhancements (v0.5.2 package additions)
- `init-project.sh` improvements
- Template updates in `templates/docs/`
- Documentation updates (CLAUDE.md, DEVOPS.md, README.md)
- Build script changes (`build-optimized.sh`)

**NOT sync-worthy:**
- Changes already in `dev-claude-integration`
- Experimental/WIP commits
- Branch-specific files (templates/claude/)

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

If this shows commits, `dev-claude-integration` has unique work (usually Claude integration files).

### 2. View file differences

```bash
# What files differ between branches?
git diff dev dev-claude-integration --name-status

# Example output:
# M    init-project.sh              (modified)
# A    templates/claude/CLAUDE.md.template  (added in dev-claude-integration)
# A    templates/claude/WORKFLOW.md         (added in dev-claude-integration)
```

**Expected differences:**
- `templates/claude/` directory (only in dev-claude-integration)
- `init-project.sh` has Claude file copying logic (lines 224-248)
- `plan.md`, `tasks.md` may have Claude-specific documentation

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
   - Resolution: Prefer dev-claude-integration for templates/claude/, dev for templates/docs/

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
```

**On dev-claude-integration branch:**
```bash
git checkout dev-claude-integration
ls templates/claude/ && echo "✅ OK: Claude templates exist" || echo "❌ ERROR: Claude templates missing!"
```

### Test 2: Test init-project.sh on dev branch

```bash
git checkout dev
./init-project.sh /tmp/test-dev basic-rna --git-init

# Verify created files
ls /tmp/test-dev/CLAUDE.md 2>/dev/null && echo "❌ ERROR: CLAUDE.md should not exist" || echo "✅ OK"
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

### Test 4: Verify Docker build still works

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

### Q: Why two branches instead of one with feature flags?

**A:** Separation of concerns:
- `dev` = Universal core (anyone can use)
- `dev-claude-integration` = Opt-in AI layer (Claude users only)

This allows non-Claude users to clone/use `dev` without seeing AI-specific files, and keeps the repo structure clean.

### Q: Can I work directly on dev-claude-integration?

**A:** Yes, for Claude-specific changes (templates/claude/, agent stubs). For core changes (Docker, R packages, init-project.sh core logic), work on `dev` first, then merge.

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
- [ ] Tested `init-project.sh` on both branches
- [ ] Pushed both branches to remote if needed

---

**Maintainer Notes:**
- This vignette should be updated when branch strategy changes
- Keep examples up-to-date with actual repo state
- Add new pitfalls as they're discovered

**Last Review:** 2025-11-11 (v0.5.2 release)
