# Claude Code Workflow for Bioinformatics Projects

> "If I had an hour to vibe code a solution, I'd spend the first 50 minutes on setting the damn context."
> — [DeanOnDelivery](https://www.reddit.com/user/DeanOnDelivery/)

---

## Core Philosophy

**Context > Speed**

Bioinformatics analyses are complex, multi-stage processes. The most critical skill is **context management**, not raw coding speed.

**Key principle:** Never dump the whole project at once. This **always** leads to buggy code when:
- Context becomes too long (Claude loses track)
- Multiple codebases are mixed (confusion about file structure)
- Testing is skipped (bugs accumulate)

**Instead:** Work **phase-by-phase**, clearing context between stages.

---

## 🗂️ Ad-Hoc Documentation System

Create these files as needed to manage context:

### 1. **plan.md** (Strategic Compass)

**Purpose:** High-level strategy, scientific narrative, known issues
**Length:** ~2000-3000 words
**When to create:** Beginning of project
**When to update:** After major milestones or discovering new issues

**Contents:**
- Scientific question and hypotheses
- Known issues and proposed solutions
- Data inventory and file paths
- Analysis narrative (Levels 1-4: QC → Annotation → Differential → Functional)

**How to use with Claude:**
- In plan mode, have Claude iterate over your rough plan to refine the global picture
- Ask Claude: "Review this plan.md and suggest improvements to the analysis strategy"
- Use Claude to identify gaps in your approach

---

### 2. **tasks.md** (Execution Tracker)

**Purpose:** Concrete stages, substeps, testing protocols
**Length:** ~2000-2500 words (grows over time)
**When to create:** After finalizing plan.md
**When to update:** Real-time during execution

**Contents:**
- Numbered stages (S1, S2, S3...)
- Substeps with clear deliverables
- Testing criteria for each substep
- Status markers (⏸️ pending, 🔄 in progress, ✅ complete, ⛔ blocked)

**How to use with Claude:**
- Ask Claude to create definitive tasks from your plan.md
- Example: "Based on plan.md, create detailed tasks.md with testable substeps"
- Mark each substep as ✅ in tasks.md before moving to the next
- **Critical:** Only ONE substep should be 🔄 at a time

---

### 3. **notes.md** (Research Findings)

**Purpose:** External resources, troubleshooting notes, discoveries
**Length:** Unlimited (append-only)
**When to create:** As needed
**When to update:** After internet searches, when discovering solutions

**Contents:**
- URLs to relevant papers, Stack Overflow, documentation
- Explanations of complex concepts
- Troubleshooting steps that worked
- Benchmarking results

**How to use with Claude:**
- Have Claude write findings to notes.md when it does internet searches
- Example: "Search for scATAC-seq peak calling best practices and document in notes.md"
- Reference notes.md when encountering similar issues later

---

### 4. **handoff.md** (Session Continuity)

**Purpose:** Context transfer between Claude Code instances
**Length:** ~500-1000 words per handoff
**When to create:** End of work session, before container switch, when passing to collaborator
**When to update:** Every session boundary

**Contents:**
- What was accomplished this session
- What's blocked (and why)
- Next immediate steps
- File paths for inputs/outputs created
- Technical notes (parameter choices, quirks)

**How to use with Claude:**
- Ask Claude to document what it accomplished: "Update handoff.md with this session's progress"
- Read handoff.md first when resuming work
- Especially critical for container switches (e.g., dev-core → dev-archr)

---

## 🚀 Recommended Workflow

### Phase 1: Planning (Plan Mode)

**Goal:** Define the strategy without writing code yet

1. **Run** `/init` (or similar) to start planning mode
2. **Write** a rough draft of plan.md with:
   - Scientific question
   - Expected data
   - High-level analysis stages
3. **Ask Claude** to iterate: "Review this plan and refine the analysis strategy"
4. **Iterate** 2-3 rounds until plan is solid
5. **Ask Claude** to create tasks.md: "Create detailed tasks.md based on this plan, with testing steps"
6. **Review** tasks.md for completeness

**Time investment:** 30-60 minutes
**Benefit:** Prevents wasted time on wrong approaches

---

### Phase 2: Execution (Phase-by-Phase)

**Goal:** Implement one stage at a time, thoroughly

#### Before Each Stage:

1. **Clear the chat** (fresh context)
2. **Claude reads:**
   - CLAUDE.md (auto-loaded)
   - plan.md (for context)
   - tasks.md (for current stage)
   - handoff.md (if resuming)
3. **Mark stage as 🔄** in tasks.md

#### During Each Stage:

4. **Execute substeps** one at a time
5. **Test after each substep** (visual QC, metrics validation)
6. **Save checkpoints** liberally
7. **Update tasks.md** status in real-time
8. **If stuck:** Mark ⛔, document in plan.md "Known Issues", create handoff.md

#### After Each Stage:

9. **Mark all substeps as ✅** in tasks.md
10. **Clear the chat again** (remove stale context)
11. **Document progress** in handoff.md
12. **Move to next stage**

**Why this works:**
- Fresh context for each stage (no context overflow)
- Systematic testing prevents bug accumulation
- Clear progress tracking
- Easy to resume after breaks

---

## 🧠 Tips for Claude Code Sessions

### 1. **Start fresh frequently**

❌ **Wrong:** Keep adding to 100-message conversation
✅ **Right:** Clear chat every 1-2 stages

**Benefit:** Claude maintains focus, doesn't conflate earlier context with current work

---

### 2. **Use double ESC to go back in chat**

If Claude misunderstands something, **don't try to convince it otherwise**.

❌ **Wrong:** "No, you're wrong, I meant X not Y"
✅ **Right:** Double ESC, rephrase from a different angle

**Benefit:** Avoids argumentative loops, saves time

---

### 3. **Write context down, don't rely on message history**

❌ **Wrong:** "Remember earlier when I said X?"
✅ **Right:** Write X in plan.md/notes.md, reference it

**Benefit:** Context persists across sessions, no memory drift

---

### 4. **Be very precise about what to execute**

❌ **Wrong:** "Implement Stages 1-3"
✅ **Right:** "Execute Stage 1, Substep 1.1: Load raw data and filter by QC metrics"

**Benefit:** Claude focuses on one thing, does it well

---

### 5. **Mark each phase as done before moving on**

❌ **Wrong:** "Almost done with Stage 1, let's start Stage 2"
✅ **Right:** Complete Stage 1, test, mark ✅, clear chat, then start Stage 2

**Benefit:** Ensures every feature is finished and tested

---

## 🏗️ Building Complex Projects Cleanly

**Protocol for large projects:**

1. **Week 1:** Planning only
   - Write plan.md with all stages
   - Create tasks.md with detailed substeps
   - Review with collaborators/advisors

2. **Week 2-N:** Execute stage-by-stage
   - 1 stage per day (or per session)
   - Clear chat between stages
   - Update handoff.md at end of each day

3. **Final week:** Integration and polishing
   - Cross-stage QC
   - Figure generation
   - Documentation

**This approach works for:**
- Multi-omics integration
- Large-scale differential analysis
- Complex workflows (10+ stages)
- Collaborative projects

---

## 📊 Example Session Structure

### Session 1: Stage 1 Implementation

```
[Clear chat]
User: "Read plan.md and tasks.md. Execute Stage 1: QC and filtering."

[Claude implements S1.1-S1.3]

User: "All tests passed. Mark Stage 1 as complete in tasks.md."

[Claude updates tasks.md]

User: "Update handoff.md with progress."

[Claude writes handoff]

[End session]
```

---

### Session 2: Stage 2 Implementation

```
[Clear chat - fresh context]
User: "Read CLAUDE.md, plan.md, tasks.md, and handoff.md. Start Stage 2: Clustering."

[Claude resumes from S1 checkpoints]

User: "S2.1 complete. Move to S2.2."

[Continue...]
```

---

## 🤖 Future: Claude Code Agents

**(Under development in this branch)**

These agents will automate parts of the workflow:

### 1. **handoff-writer**

**Purpose:** Auto-generate handoff.md at end of session

**Trigger:** User runs `/handoff` command

**Actions:**
- Scans tasks.md for ✅ ⛔ changes
- Identifies new checkpoints created
- Documents any errors encountered
- Writes structured handoff.md

---

### 2. **stage-reviewer**

**Purpose:** Assess whether stage is ready for ✅

**Trigger:** User runs `/review-stage S1`

**Actions:**
- Checks if all substeps are complete
- Verifies checkpoints exist
- Reviews test results
- Gives recommendation: "Refine", "Test more", or "Green light"

---

### 3. **context-optimizer** (future)

**Purpose:** Identify when to clear chat

**Trigger:** Automatic or user runs `/check-context`

**Actions:**
- Estimates current context length
- Identifies stale context
- Suggests when to start fresh session

---

## ✅ Success Criteria

**You know the workflow is working when:**

- [ ] Projects complete in predictable time
- [ ] Code is tested at every step
- [ ] You can resume after weeks away
- [ ] Collaborators can pick up your work
- [ ] Bugs are caught early, not late
- [ ] Documentation stays synchronized

---

## 🚨 Warning Signs You're Doing It Wrong

**Red flags:**

- ❌ Chat history >50 messages
- ❌ Multiple stages implemented but none tested
- ❌ "Almost done" for >3 days
- ❌ Tasks.md out of sync with reality
- ❌ Can't remember what file does what
- ❌ Starting over because earlier code was buggy

**If you see these:** Stop, clear chat, go back to last ✅ checkpoint, resume phase-by-phase

---

## 📚 Further Reading

- **plan.md** in this project - Example of strategic planning
- **tasks.md** in this project - Example of execution tracking
- **CLAUDE.md** in this project - Lightweight context file design

---

**End of WORKFLOW.md**

**Remember:** Context > Speed. Spend 50 minutes on context, 10 on coding.
