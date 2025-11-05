# Handoff Writer Agent

**Status:** Stub (Not yet implemented)
**Purpose:** Auto-generate handoff.md at end of Claude Code session
**Trigger:** User command `/handoff` or automatic at session end

---

## Overview

The Handoff Writer agent automates the creation and updating of `handoff.md` files to ensure smooth context transfer between:
- Claude Code instances (across sessions)
- Container switches (dev-core ↔ dev-archr)
- Collaborators (passing work to team members)
- Work resumption (after breaks or days away)

---

## Intended Functionality

### Inputs

1. **tasks.md** - Parse status changes since last handoff
   - Identify newly completed substeps (✅)
   - Identify blocked substeps (⛔) and their reasons
   - Identify currently in-progress substep (🔄)

2. **Filesystem scan** - Detect new outputs
   - New checkpoints in `03_results/checkpoints/`
   - New plots in `03_results/plots/`
   - New log files in `logs/`

3. **Session chat history** - Extract key information
   - Errors encountered and solutions attempted
   - Parameter choices and rationale
   - Technical notes and quirks discovered

4. **plan.md** - Check for new Known Issues entries

### Processing

1. **Compare** current tasks.md with previous handoff.md
2. **Identify** what changed (completed, blocked, new files)
3. **Extract** relevant technical details from chat
4. **Determine** next immediate steps based on tasks.md

### Outputs

Generate/update `handoff.md` with structure:

```markdown
# Handoff Note

**Date:** YYYY-MM-DD HH:MM
**From:** Claude Code Instance
**To:** Next instance / Collaborator
**Current Stage:** S#.# — {name}

## Context
[2-3 sentence summary]

## Progress Since Last Handoff
- [x] S#.# — {substep name}
- [ ] S#.# — {in progress, incomplete}

## Current Status
**Working on:** S#.# — {active substep}
**Blocked by:** {issue or "None"}

## Next Steps
1. {First immediate action}
2. {Second action}

## File Locations
**Outputs Created This Session:**
- Checkpoint: {path}
- Plot: {path}

## Technical Notes
- Parameter choices: {rationale}
- Known issues: {any quirks}
```

---

## Usage Examples

### Manual Trigger

```bash
# In Claude Code interface
/handoff

# Agent scans tasks.md, filesystem, chat
# Writes handoff.md
# Displays: "✓ handoff.md updated"
```

### Automatic Trigger (Future)

```json
// In .claude/config.json
{
  "agents": {
    "handoff-writer": {
      "enabled": true,
      "triggers": [
        "session_end",
        "container_switch"
      ]
    }
  }
}
```

---

## Implementation Notes

**Dependencies:**
- Read access to `tasks.md`, `plan.md`
- Filesystem read access (checkpoints, logs)
- Chat history parsing

**Challenges:**
- Determining what's "relevant" from long chat histories
- Handling cases where tasks.md wasn't updated properly
- Avoiding overwrite of manual edits to handoff.md

**Proposed solution:**
- Use timestamped handoff files: `handoff_YYYYMMDD_HHMMSS.md`
- Keep symbolic link `handoff.md` → latest
- User can manually merge if needed

---

## Testing Criteria

Agent is successful when:
- [ ] Handoff captures all completed substeps
- [ ] Next steps are actionable and specific
- [ ] File paths are correct and exist
- [ ] Technical notes include relevant parameter choices
- [ ] Blockers are clearly documented with context

---

## Future Enhancements

1. **Smart summarization** - Use LLM to extract key technical details from chat
2. **Diff mode** - Show what changed since last handoff
3. **Integration with stage-reviewer** - Include review results in handoff
4. **Slack/Teams integration** - Post handoff summary to team channel

---

**Status:** Awaiting implementation
**Priority:** Medium (Manual handoff works, but automation would save time)
**See also:** `WORKFLOW.md` section on "handoff.md" usage
