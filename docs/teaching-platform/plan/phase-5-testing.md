# Phase 5: Stress Testing & Hardening

**Goal:** Validate the system's stability under the full class load and finalize the deployment.

## Tasks
1.  **Concurrency Simulation:**
    - Use a script to log in 12 concurrent "robot" users.
    - Trigger a memory-intensive Seurat analysis in each container.
    - Monitor host telemetry (CPU/RAM/VRAM) to ensure no system-wide instability.
2.  **GPU Contention Test:**
    - Run 3-4 concurrent GPU-accelerated tasks (e.g., scVI training) to observe performance degradation.
3.  **Final Hardening:**
    - Disable unnecessary Jupyter extensions.
    - Set up a basic `systemd` service for JupyterHub to ensure it restarts on host reboot.

## Success Criteria
- 12 users active simultaneously with 40GB RAM limits enforced.
- Host RAM usage stays within the 512GB budget (max ~480GB for students).
- Documentation of "Common Troubleshooting" steps for students is completed.
