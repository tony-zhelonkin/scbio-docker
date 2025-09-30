# Git Branch Recovery Vignette: Fixing Incorrect Merge Strategy

## Table of Contents
1. [The Situation](#the-situation)
2. [What Went Wrong](#what-went-wrong)
3. [Understanding Git Merge Strategies](#understanding-git-merge-strategies)
4. [The Recovery Process](#the-recovery-process)
5. [Detailed Step-by-Step Recovery](#detailed-step-by-step-recovery)
6. [Prevention and Best Practices](#prevention-and-best-practices)
7. [Quick Reference Commands](#quick-reference-commands)

---

## The Situation

### Initial State
I had three branches with different development histories:
- **`lt_dev`**: My active, modern development branch with the latest work (v0.4.0+)
- **`main`**: Outdated, hadn't been updated in a long time, went in a different direction (v0.2)
- **`dev`**: Also outdated, changed significantly since it diverged from `lt_dev` (v0.3)

### The Goal
I wanted to unify all branches with `lt_dev` taking precedence, essentially making `main` and `dev` match `lt_dev` exactly while preserving Git history for auditing purposes.

### The Problem
I chose **Option 2: Merge with "ours" strategy** thinking it would adopt `lt_dev`'s content. However, I misunderstood the strategy direction, and it did the opposite - it kept the OLD content from `main` and `dev`, overwriting my modern `lt_dev` work.

---

## What Went Wrong

### The Incorrect Commands

```bash
# Step 1: Merged lt_dev INTO main using "ours" strategy
git checkout main
git merge -s ours lt_dev -m "Merge lt_dev into main, adopting lt_dev's direction"
# ❌ WRONG: This kept main's OLD content

# Step 2: Merged main INTO lt_dev (fast-forward)
git checkout lt_dev
git merge main
# ❌ WRONG: This made lt_dev adopt main's old content

# Step 3: Same mistake repeated for dev
git checkout dev
git merge -s ours lt_dev -m "Merge lt_dev into dev, adopting lt_dev's direction"
git checkout lt_dev
git merge dev
# ❌ WRONG: This further corrupted lt_dev with dev's old content
```

### Understanding the Mistake

The `-s ours` strategy means:
- "Create a merge commit that connects both branches"
- "But keep OUR (current branch's) content, ignoring THEIRS (incoming branch's) content"

So when I ran:
```bash
git checkout main
git merge -s ours lt_dev
```

Git interpreted this as:
- Current branch: `main` (the OLD branch)
- Incoming branch: `lt_dev` (the NEW branch)
- Strategy: Keep `main`'s content, ignore `lt_dev`'s changes
- Result: `main` stayed old, `lt_dev` got corrupted when merged back

**The cardinal rule I violated**: When using `-s ours`, the current branch's content WINS. I was on the OLD branches, so the old content won.

---

## Understanding Git Merge Strategies

### Common Merge Strategies

1. **Default merge (no strategy flag)**
   ```bash
   git merge branch-name
   ```
   - Git tries to intelligently merge changes
   - Creates merge conflicts if there are overlapping changes
   - Best for: Normal development where both branches have valid changes

2. **`-s ours` strategy**
   ```bash
   git merge -s ours other-branch
   ```
   - Keeps CURRENT branch's content entirely
   - Ignores ALL changes from the other branch
   - Creates a merge commit for history
   - Best for: When you want to mark two branches as merged in history but discard the incoming changes
   - **Dangerous if misunderstood!**

3. **`-X ours` option (different from `-s ours`!)**
   ```bash
   git merge -X ours other-branch
   ```
   - Performs a normal merge
   - On conflicts, automatically chooses CURRENT branch's version
   - Still merges non-conflicting changes from both sides
   - Best for: When you want to merge but automatically resolve conflicts in your favor

4. **`-X theirs` option**
   ```bash
   git merge -X theirs other-branch
   ```
   - Performs a normal merge
   - On conflicts, automatically chooses the INCOMING branch's version
   - Best for: When you want to merge but automatically resolve conflicts in their favor

### What I Should Have Used

For my situation (making `main` and `dev` match `lt_dev` completely), the correct approach is:

**Option A: Hard reset (simplest)**
```bash
git checkout main
git reset --hard lt_dev
git push --force origin main
```

**Option B: Merge from the right direction**
```bash
# Merge main INTO lt_dev, preferring lt_dev's changes
git checkout lt_dev
git merge -X theirs main -m "Adopt lt_dev's direction"
# Then reset main to match
git checkout main
git reset --hard lt_dev
```

---

## The Recovery Process

### Overview

The recovery relied on Git's **reflog** - a safety net that tracks EVERY move of HEAD (branch pointer), even if commits become unreachable through normal means.

Git's reflog keeps a history for ~30-90 days by default, which means you can almost always undo recent mistakes.

### Key Concepts

1. **Reflog**: A per-branch log of where that branch's HEAD has been
2. **SHA/Commit hash**: A unique identifier for each commit (e.g., `82862bb`)
3. **HEAD**: A pointer to the current commit you're on
4. **Reset**: Moving a branch pointer to a different commit

---

## Detailed Step-by-Step Recovery

### Phase 1: Identifying the Problem

When I switched back to `main` after the merge, I saw:
```bash
git checkout main
ls -la  # Old directory structure!
```

The files looked old, not like my modern `lt_dev`. This was the first clue something went very wrong.

### Phase 2: Understanding What Happened

I checked the branch state:
```bash
git branch
#   dev
# * lt_dev
#   main
```

Then checked `lt_dev`:
```bash
git checkout lt_dev
ls -la  # Also old structure!
```

Both `lt_dev` AND the other branches had the OLD code. The merge had overwritten my modern work.

### Phase 3: Finding the Good Commit

This is the CRITICAL step. I used `git reflog` to find the last good state of `lt_dev`:

```bash
git reflog lt_dev -10
```

Output:
```
c9c8b9e lt_dev@{0}: merge dev: Fast-forward
c0ac519 lt_dev@{1}: merge main: Fast-forward
82862bb lt_dev@{2}: commit: README update          ← FOUND IT!
36f9058 lt_dev@{3}: commit: Pre v0.4.1 build state: update logic, README
528f0ae lt_dev@{4}: commit: Successfull build
0bf05cc lt_dev@{5}: commit: Poststart logic modified
2b26683 lt_dev@{6}: commit: Separating pyenv envs to baase, atac, squid
e74240a lt_dev@{7}: commit: Two-image approach
9848c06 lt_dev@{8}: commit: Cleaning up repo, changing course
4accd54 lt_dev@{9}: branch: Created from HEAD
```

**How to read this:**
- `lt_dev@{0}` is the current (broken) state
- `lt_dev@{1}` is the previous state (also broken)
- `lt_dev@{2}` is TWO moves ago - before any merges started

The entries show:
- `@{0}` and `@{1}`: "Fast-forward" merges that corrupted the branch
- `@{2}`: "commit: README update" - a regular commit, meaning this was normal development
- The commit hash is `82862bb`

**Why `82862bb` was the right commit:**
1. It's the last entry that says "commit" (not "merge")
2. It's right before the problematic merges started
3. The commit message "README update" was recognizable as recent work

### Phase 4: Verifying the Commit

Before resetting everything, I could have verified this commit was correct:

```bash
# View the commit details
git show 82862bb

# Or checkout in detached HEAD state to inspect
git checkout 82862bb
ls -la  # Check if files look right
cat README.md  # Verify content
git log --oneline -5  # Check recent history

# Return to branch
git checkout lt_dev
```

In a real scenario, you'd want to verify before proceeding. In our case, I was confident based on the reflog.

### Phase 5: The Recovery - Resetting lt_dev

```bash
# Make sure I'm on lt_dev
git checkout lt_dev

# Reset lt_dev to the good commit
git reset --hard 82862bb
```

Output:
```
HEAD is now at 82862bb README update
```

**What `git reset --hard` does:**
- `--hard`: Resets three things:
  1. The branch pointer (where `lt_dev` points)
  2. The staging area (index)
  3. The working directory (actual files)
- Moves `lt_dev` to point at `82862bb`
- Discards all changes after that commit
- **DANGEROUS**: Can lose uncommitted work (but that's what we wanted here)

After this, I verified:
```bash
ls -la  # Modern structure restored!
git log --oneline -5  # Check commit history
```

### Phase 6: Recovering main and dev

Now I needed to reset `main` and `dev` to their pre-merge states, then properly unify them.

**Finding main's good state:**
```bash
git reflog main -5
```

Output:
```
c0ac519 main@{0}: merge lt_dev: Merge made by the 'ours' strategy.
18b1b33 main@{1}: clone: from github.com:tony-zhelonkin/scbio-docker.git
```

Here, `main@{1}` (commit `18b1b33`) was the last state before our bad merge.

**Resetting main:**
```bash
git checkout main
git reset --hard 18b1b33
```

**Finding dev's good state:**
```bash
git reflog dev -5
```

Output:
```
82862bb dev@{0}: reset: moving to 82862bb    ← Wait, this was from our previous reset attempt
c9c8b9e dev@{1}: merge lt_dev: Merge made by the 'ours' strategy.
4accd54 dev@{2}: commit: Add R packages      ← FOUND IT!
d9c3f4d dev@{3}: commit: Fixing install_R_packages.R bug
541c67b dev@{4}: commit: Dockerfile user group bug
```

Here, `dev@{2}` (commit `4accd54`) was the last clean state.

**Resetting dev:**
```bash
git checkout dev
git reset --hard 4accd54
```

### Phase 7: The Correct Unification

Now all three branches were restored to their pre-merge states. To unify them CORRECTLY with `lt_dev` taking precedence:

```bash
# Option 1: Simple hard reset (my choice)
git checkout main
git reset --hard 82862bb  # Make main identical to lt_dev

git checkout dev
git reset --hard 82862bb  # Make dev identical to lt_dev

# Verify all branches point to the same commit
git log --oneline --all --graph -10
```

Output shows all branches at the same commit:
```
* 82862bb README update          ← All three branches here
* 36f9058 Pre v0.4.1 build state: update logic, README
* 528f0ae Successfull build
...
```

### Phase 8: Pushing to Remote (requires SSH access)

```bash
# Push all three branches, forcing the remote to accept the reset
git push --force origin main dev lt_dev
```

**Note**: `--force` is necessary because we're rewriting history. The remote has the "broken" commits, and we're replacing them with the "fixed" state.

**WARNING about `--force`:**
- Use with caution on shared branches
- Can cause problems for collaborators who pulled the broken state
- In this case, it's safe because I was the only developer
- Alternative: `--force-with-lease` is safer for shared repos

### Phase 9: Cleanup (optional)

Since `main` is now the primary branch and `lt_dev` served its purpose:

```bash
# Switch to main for future work
git checkout main

# Delete lt_dev locally (optional)
git branch -d lt_dev

# Delete lt_dev remotely (optional)
git push origin --delete lt_dev
```

---

## Prevention and Best Practices

### Before Merging

1. **Understand your merge strategy**
   ```bash
   # If you want to keep EVERYTHING from branch-to-merge:
   git checkout old-branch
   git reset --hard new-branch

   # If you want to merge with new-branch winning conflicts:
   git checkout target-branch
   git merge -X theirs source-branch
   ```

2. **Test in a backup branch first**
   ```bash
   # Create a test branch
   git checkout -b test-merge

   # Try your merge
   git merge -s ours other-branch

   # Inspect result
   ls -la
   git log --oneline -5

   # If good, repeat on real branch; if bad, delete test branch
   git checkout main
   git branch -D test-merge
   ```

3. **Check what will happen**
   ```bash
   # See what commits are different
   git log main..lt_dev --oneline
   git log lt_dev..main --oneline

   # See file differences
   git diff main lt_dev
   ```

### During a Merge

1. **Don't force-push immediately**
   - Test locally first
   - Verify files look correct
   - Check git log looks right
   - Only then push

2. **Stop if something looks wrong**
   ```bash
   # If merge completed but looks wrong:
   git merge --abort  # Only works if merge is in progress

   # Or if merge is complete:
   git reset --hard HEAD~1  # Go back one commit
   ```

### After a Mistake

1. **Don't panic - reflog saves everything**
   ```bash
   git reflog
   ```

2. **Find the last good commit systematically**
   ```bash
   # For current branch
   git reflog -20

   # For a specific branch
   git reflog branch-name -20

   # Look for the last "commit" entry before "merge" entries
   ```

3. **Verify before resetting**
   ```bash
   # Check out the commit temporarily
   git checkout <commit-hash>

   # Inspect
   ls -la
   cat important-file.txt

   # Return to branch
   git checkout branch-name
   ```

4. **Reset carefully**
   ```bash
   # Soft reset: moves branch pointer, keeps changes staged
   git reset --soft <commit-hash>

   # Mixed reset: moves pointer, unstages changes, keeps files
   git reset --mixed <commit-hash>

   # Hard reset: moves pointer, discards everything
   git reset --hard <commit-hash>
   ```

### Communication

If working with a team:

1. **Before force-pushing**, notify team members
2. **Provide instructions** for them to reset their local branches:
   ```bash
   git fetch origin
   git reset --hard origin/main
   ```
3. **Document what happened** (like this vignette!)

---

## Quick Reference Commands

### Finding the Problem

```bash
# Check current branch and status
git status
git branch

# View recent history
git log --oneline -20
git log --oneline --all --graph -20

# See all branches and their commits
git log --oneline --all --decorate
```

### Using Reflog to Find Good Commits

```bash
# View reflog for current branch
git reflog

# View reflog for specific branch
git reflog branch-name -20

# View detailed reflog
git reflog show --all

# Find when you were on a specific commit
git reflog | grep "commit-hash"
```

### Inspecting a Commit Before Resetting

```bash
# Checkout commit without changing branch
git checkout <commit-hash>

# Or view commit details
git show <commit-hash>
git show <commit-hash>:path/to/file.txt

# Compare with current state
git diff HEAD <commit-hash>

# Return to branch
git checkout branch-name
```

### Resetting Branches

```bash
# Reset current branch to specific commit (DESTRUCTIVE)
git reset --hard <commit-hash>

# Reset to a certain number of commits back
git reset --hard HEAD~3  # Go back 3 commits

# Reset branch to match another branch exactly
git checkout branch-to-reset
git reset --hard branch-to-match

# Reset branch to match remote
git fetch origin
git reset --hard origin/branch-name
```

### Fixing Multiple Branches

```bash
# Typical recovery scenario
git checkout branch1
git reset --hard <good-commit-hash>

git checkout branch2
git reset --hard <good-commit-hash>

git checkout branch3
git reset --hard <good-commit-hash>

# Verify all branches are at expected commit
git log --oneline --all --graph -10
```

### Pushing Fixed Branches

```bash
# Force push single branch (DESTRUCTIVE to remote)
git push --force origin branch-name

# Force push multiple branches
git push --force origin branch1 branch2 branch3

# Safer force push (fails if remote changed since last fetch)
git push --force-with-lease origin branch-name
```

### Creating Backups Before Risky Operations

```bash
# Create backup branch
git branch backup-main main
git branch backup-dev dev
git branch backup-lt_dev lt_dev

# If disaster strikes, restore from backup
git checkout main
git reset --hard backup-main

# Delete backups when done
git branch -D backup-main backup-dev backup-lt_dev
```

---

## Common Scenarios and Solutions

### Scenario 1: Made wrong merge, haven't pushed yet

```bash
# Go back one commit
git reset --hard HEAD~1

# Or find exact commit in reflog
git reflog
git reset --hard HEAD@{5}  # Example: 5 moves ago
```

### Scenario 2: Made wrong merge, already pushed

```bash
# Fix locally first
git reset --hard <good-commit>

# Force push (coordinate with team!)
git push --force origin branch-name
```

### Scenario 3: Merged wrong branch direction

```bash
# Find pre-merge state in reflog
git reflog -20

# Reset to before merge
git reset --hard <commit-before-merge>

# Merge in correct direction
git merge correct-branch
```

### Scenario 4: Lost commits after reset

```bash
# Commits are in reflog for ~30 days
git reflog -50

# Find lost commit
git show <lost-commit-hash>

# Restore it
git cherry-pick <lost-commit-hash>
# Or
git reset --hard <lost-commit-hash>
```

### Scenario 5: Want to unify branches preserving history

```bash
# Option A: Merge with preference
git checkout target-branch
git merge -X theirs source-branch

# Option B: Rebase (rewrites history)
git checkout feature-branch
git rebase main

# Option C: Hard reset (simplest, loses history)
git checkout target-branch
git reset --hard source-branch
```

---

## Glossary

- **HEAD**: Pointer to your current commit
- **Branch**: A movable pointer to a commit
- **Reflog**: Log of all HEAD movements (per branch)
- **SHA/Hash**: Unique identifier for a commit (e.g., `82862bb`)
- **Fast-forward**: Moving a branch pointer forward when no divergence exists
- **Merge commit**: A commit with two parents (merges two branches)
- **Detached HEAD**: Checking out a commit directly (not a branch)
- **Reset**: Moving a branch pointer to a different commit
- **Hard reset**: Reset + discard all changes in working directory
- **Soft reset**: Reset but keep changes staged
- **Mixed reset**: Reset but keep changes in working directory (unstaged)
- **Force push**: Overwrite remote history (dangerous!)
- **Cherry-pick**: Apply a specific commit to current branch

---

## Conclusion

Git's reflog is your safety net. Even if you make catastrophic mistakes with merges, resets, or rebases, you can almost always recover using reflog to find the last good commit.

The key lessons from this incident:

1. **Understand merge strategies** before using them (`-s ours` vs `-X ours`)
2. **Test on a backup branch** first for risky operations
3. **Check results** before force-pushing
4. **Use `git reflog`** to find good commits when things go wrong
5. **Verify commits** before resetting (`git show`, `git checkout`)
6. **Force-push carefully** and coordinate with team

Remember: In Git, you can almost always undo your mistakes if you act before garbage collection runs (~30-90 days). The reflog is there to save you!

---

## Additional Resources

- [Git documentation on reset](https://git-scm.com/docs/git-reset)
- [Git documentation on reflog](https://git-scm.com/docs/git-reflog)
- [Git merge strategies](https://git-scm.com/docs/merge-strategies)
- [Atlassian Git reflog tutorial](https://www.atlassian.com/git/tutorials/rewriting-history/git-reflog)

---

**Document created**: 2025-09-30
**Author**: Recovery procedure documented by Claude Code
**Purpose**: Reference guide for Git branch recovery using reflog
