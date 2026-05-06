# Phase 4: Day Zero Provisioning

**Goal:** Automate the setup of student environments so they are ready for analysis upon first login.

## Tasks
1.  **Home Directory Pre-seeding:**
    - Script the creation of 12 student home directories on the NVMe.
    - Pre-clone the workshop repository and example scripts into each home.
    - Set correct ownership (`chown 1000:1000`).
2.  **Mount Finalization:**
    - Configure the `ro` mount for `/data/research` in `jupyterhub_config.py`.
    - Ensure student `scratch` directories are mapped to the NVMe.
3.  **API Key Injection:**
    - Configure the Hub to inject `GEMINI_API_KEY` and `CLAUDE_API_KEY` into each container's environment.

## Success Criteria
- A new student logs in and finds the class repo already in their home folder.
- The student can run the example R script immediately without path errors.
- AI agents are pre-authenticated via the injected environment variables.
