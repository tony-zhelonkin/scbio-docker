# Branches (simple guide)

This is a personal toolkit. Branching is kept simple.

Branches
- dev: main working branch.
- dev-claude-integration: optional AI/Claude extras layered on top of dev.
- main: stable snapshots when needed.

Syncing dev → dev-claude-integration
```bash
git checkout dev-claude-integration
git merge dev
git push origin dev-claude-integration
```

Quick checks
```bash
# Commits on dev not yet in dev-claude-integration
git log dev-claude-integration..dev --oneline | head

# Files that differ
git diff --name-status dev dev-claude-integration | head
```

Releasing
```bash
git checkout main
git merge dev -m "Release <tag>"
git tag -a <tag> -m "<note>"
git push origin main <tag>
```

That’s it—keep it lightweight.

