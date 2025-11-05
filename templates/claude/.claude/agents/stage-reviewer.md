# Stage Reviewer Agent

**Status:** Stub (Not yet implemented)
**Purpose:** Assess whether a stage is ready for completion (✅)
**Trigger:** User command `/review-stage S#` or before marking stage complete

---

## Overview

The Stage Reviewer agent acts as a quality gate to ensure stages are truly complete before moving forward. It prevents:
- Incomplete implementations marked as "done"
- Skipped testing steps
- Missing checkpoints or outputs
- Deviation from plan.md strategy

---

## Intended Functionality

### Inputs

1. **tasks.md** - Current stage definition
   - All substeps for the stage
   - Testing criteria for each substep
   - Expected checkpoints

2. **Filesystem** - Verify outputs exist
   - Checkpoints: `03_results/checkpoints/S#_*.rds`
   - Plots: `03_results/plots/S#/`
   - Tables: `03_results/tables/S#/` (if applicable)

3. **Log files** - Check for errors
   - Parse `logs/S#_*.log` for errors, warnings
   - Identify unhandled exceptions

4. **plan.md** - Ensure alignment with strategy
   - Check if stage achieves stated goals
   - Verify no "Known Issues" block this stage

5. **Code** - Scan analysis scripts
   - Ensure testing code is present (not just analysis)
   - Check for `saveRDS()` / `write_h5ad()` calls

### Processing

1. **Completeness check:**
   - All substeps marked ✅?
   - All expected checkpoints exist?
   - All expected plots generated?

2. **Quality check:**
   - Log files show no critical errors?
   - Testing sections in scripts executed?
   - Visual QC plots created?

3. **Alignment check:**
   - Stage goals from plan.md achieved?
   - Results make biological sense?
   - No blockers preventing next stage?

### Outputs

Generate review report with recommendation:

```markdown
# Stage Review: S# — {Stage Name}

**Reviewed:** YYYY-MM-DD HH:MM
**Status:** [GREEN LIGHT / REFINE / BLOCKED]

## Completeness (X/Y)
- [x] All substeps completed
- [ ] Missing checkpoint: S#.#_name.rds
- [x] All plots generated

## Quality (X/Y)
- [x] No critical errors in logs
- [x] Testing code executed
- [ ] Visual QC plot missing for {aspect}

## Alignment
- [x] Stage goals achieved
- [x] Results biologically plausible
- [x] No blockers for next stage

## Recommendation

**GREEN LIGHT**: Stage S# is complete. Safe to mark all substeps ✅ and proceed to S#next.

**REFINE**: Address these issues before proceeding:
1. {Issue 1}
2. {Issue 2}

**BLOCKED**: Cannot proceed due to:
- {Blocker description}
- Suggested action: {What to do}

## Next Actions

[Specific steps to complete the stage or move forward]
```

---

## Usage Examples

### Manual Review

```bash
# In Claude Code interface
/review-stage S1

# Agent scans tasks.md, filesystem, logs
# Generates review report
# Recommendation: GREEN LIGHT / REFINE / BLOCKED
```

### Automated Review (Future)

```bash
# In tasks.md, user marks stage as "Review requested"
## Stage 1: QC and Filtering 🔄 (Review)

# Agent automatically runs review
# Updates tasks.md with review results
```

---

## Review Criteria

### GREEN LIGHT (✅ Ready to proceed)

- All substeps completed
- All checkpoints exist and loadable
- All expected plots/tables created
- No critical errors in logs
- Testing sections executed
- Results align with plan.md goals

### REFINE (⚠️ Needs work)

- 1-2 substeps incomplete or failed tests
- Minor issues in output quality
- Visual QC suggests potential problems
- Non-critical warnings in logs

### BLOCKED (⛔ Cannot proceed)

- Critical substep failed
- Missing essential checkpoints
- Errors block next stage
- Results don't align with expected biology
- Known Issue in plan.md unresolved

---

## Implementation Notes

**Dependencies:**
- Read access to `tasks.md`, `plan.md`
- Filesystem read/validate access
- Log file parsing
- Code scanning (basic AST analysis)

**Challenges:**
- Determining "biological plausibility" requires domain knowledge
- Log parsing is fragile (different tools, formats)
- Defining "complete" is sometimes subjective

**Proposed solution:**
- Focus on objective checks (files exist, no errors)
- User confirms subjective aspects ("Do results look reasonable?")
- Agent highlights anomalies for user review

---

## Testing Criteria

Agent is successful when:
- [ ] Correctly identifies missing checkpoints
- [ ] Detects errors in log files
- [ ] Recognizes incomplete substeps
- [ ] Provides actionable recommendations
- [ ] No false positives (flagging complete stages as incomplete)
- [ ] No false negatives (green-lighting incomplete stages)

---

## Integration Points

### With tasks.md

Agent updates tasks.md directly:

```markdown
## Stage 1: QC and Filtering ✅

**Review Results (YYYY-MM-DD):**
- Status: GREEN LIGHT
- All checks passed
- Ready for Stage 2
```

### With handoff-writer

Stage review results included in handoff.md:

```markdown
## Stage Review Summary
- S1: GREEN LIGHT (2025-01-15)
- S2: In progress
```

---

## Future Enhancements

1. **Checkpoint validation** - Load and inspect checkpoints, verify structure
2. **Plot analysis** - Use computer vision to detect issues (e.g., empty plots, failed renderings)
3. **Statistical checks** - Verify metrics within expected ranges (e.g., "20k-30k cells")
4. **Integration with plan.md** - Auto-update "Known Issues" if stage reveals problems
5. **Learning mode** - Track common issues across projects, warn proactively

---

## Example Review Session

```
User: /review-stage S1

Agent:
Scanning S1 substeps...
✓ S1.1 complete (checkpoint exists)
✓ S1.2 complete (checkpoint exists)
✓ S1.3 complete (checkpoint exists)
✓ S1.4 complete (checkpoint exists)

Checking outputs...
✓ 4 checkpoints found
✓ 3 plots in 03_results/plots/S1/
⚠️  Expected plot missing: S1_tss_vs_nucleosome.pdf

Checking logs...
✓ No critical errors
⚠️  Warning: "High nucleosome signal in 15% of cells"

Checking alignment with plan.md...
✓ Stage goal: "Remove low-quality cells" → ACHIEVED
✓ Next stage (S2) has no dependencies blocked

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

RECOMMENDATION: REFINE

Address before proceeding:
1. Generate missing S1_tss_vs_nucleosome.pdf plot
2. Review cells with high nucleosome signal (consider stricter filtering)

Once addressed, run /review-stage S1 again.
```

---

**Status:** Awaiting implementation
**Priority:** High (Prevents premature stage completion, catches issues early)
**See also:** `WORKFLOW.md` section on "Testing Requirements"
