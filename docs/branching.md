# Branches (simple guide)

Created as a personal toolkit to support working with remote computational resources in interactive VS Code session, have ended up with minimal branching.

Active branches
- dev: working branch.
- main: stable snapshots/releases.

Archived AI branches
- Previously used: dev-claude-integration and dev-gpt-codex-integration.
- Now archived as tags: `archived/dev-claude-integration`, `archived/dev-gpt-codex-integration`.
- Agent/MCP setup will probably live in SciAgent‑toolkit: https://github.com/tony-zhelonkin/SciAgent-toolkit

Quick checks
```bash
# Commits on dev not yet in main
git log main..dev --oneline | head

# Files that differ
git diff --name-status main dev | head
```

Release (dev → main)
```bash
git checkout main
git merge dev -m "Release <tag>"
git tag -a <tag> -m "<note>"
git push origin main <tag>
```

Try to keep it lightweight.
