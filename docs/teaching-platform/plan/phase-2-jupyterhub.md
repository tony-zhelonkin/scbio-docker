# Phase 2: The "Thin Slice" Environment

**Goal:** Build the complete `scdock-teaching:v1.0` image and validate the full toolstack for a single user, ensuring security and complete dependency coverage.

## Tasks
1.  **Image Layering & Hardening (scdock-teaching:v1.0):**
    - Integrate `Positron Server` and `jupyter-positron-server`.
    - Install missing case study R dependencies (`maplet`, `webr`, `colorRamps`, `ggbreak`, `ggforce`, `extrafont`, `moonBook`, `openxlsx`, `writexl`, `xlsx`).
    - Remove `sudo` access from the `devuser` to prevent privilege escalation.
    - Pre-install the AI agent globally (`npm install -g @mariozechner/pi-coding-agent`).
2.  **Toolstack Validation:**
    - Launch the container manually using strict security flags (`--cap-drop=ALL`, `--security-opt=no-new-privileges`).
    - Verify Positron is accessible via the browser.
    - Run the `pi` agent, use `/login` to configure OpenRouter, and test a "Hello World" task: "Summarize the data in /data/research".
3.  **Proxy Refinement:**
    - Ensure `jupyter-server-proxy` correctly routes Positron and RStudio.

## Success Criteria
- `scdock-teaching:v1.0` builds without errors and contains no `sudo` privileges.
- A single user can log in, open Positron, and securely launch `pi` with their own OpenRouter key.
- The terminal inside Positron has all dependencies available.
- Container security restrictions (`cap-drop=ALL`) do not break IDE functionality.
