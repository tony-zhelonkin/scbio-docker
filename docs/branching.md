# Branches (simple guide)

Created as a personal toolkit for working with remote computational resources in interactive VS Code sessions, this repo has ended up with minimal branching.

## Active branches

| Branch | Role |
|--------|------|
| `main` | Default branch; stable snapshots/releases. |
| `dev`  | Working branch. |

## Archived AI branches

Earlier AI-integration work lived on `dev-claude-integration` and `dev-gpt-codex-integration`. These are now kept only as tags under the `archive/` prefix:

```bash
git tag -l 'archive/*'
# archive/dev-claude-integration
# archive/dev-gpt-codex-integration
# archive/dev-pre-rewrite
# archive/dev-restructure
# archive/main-pre-rewrite

# Inspect an archived branch
git checkout archive/dev-claude-integration
```

Agent/MCP setup now lives in SciAgent-toolkit: https://github.com/tony-zhelonkin/SciAgent-toolkit

## Quick checks

```bash
# Commits on dev not yet in main
git log main..dev --oneline | head

# Files that differ
git diff --name-status main dev | head
```

## Release (dev → main)

```bash
git checkout main
git merge dev -m "Release <tag>"
git tag -a <tag> -m "<note>"
git push origin main <tag>
```

Try to keep it lightweight.
