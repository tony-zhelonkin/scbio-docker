# Implementation Plan: Metabolism 2026 Teaching Platform

This directory contains the step-by-step implementation plan for the Dell Precision 7960 teaching platform. The plan follows a **"deep thin slice"** DevOps methodology: we build the full stack (Data -> Hub -> Image -> IDE -> AI) for a single user first to validate the architecture before scaling to the full class.

## Phases Index

1.  **[Phase 1: Foundations & Data Skeleton](phase-1-image.md):** Host preparation, Zenodo data acquisition, and skeleton configuration.
2.  **[Phase 2: The "Thin Slice" Environment](phase-2-jupyterhub.md):** Building the complete Docker image (Positron + AI Agents) and validating the toolchain.
3.  **[Phase 3: JupyterHub Orchestration](phase-3-storage.md):** Implementing the central hub, resource limits (40GB/6C), and shared GPU access.
4.  **[Phase 4: Day Zero Provisioning](phase-4-integration.md):** Automating student home directory seeding and read-only data mounting.
5.  **[Phase 5: Stress Testing & Hardening](phase-5-testing.md):** Simulating 12 concurrent users and validating performance under load.
