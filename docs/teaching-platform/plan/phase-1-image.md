# Phase 1: Foundations & Data Skeleton

**Goal:** Establish the physical storage structure and acquire the workshop dataset.

## Tasks
1.  **Host Directory Structure:**
    - Create `/mnt/nvme/students/` (Home persistence).
    - Create `/mnt/nvme/scratch/` (Fast temporary IO).
    - Create `/mnt/hdd/research_data/metabolism_v0.3.4/` (Archival storage).
2. **Data Acquisition:**
    - Download Zenodo record `7874228` (~2GB).
    - Verify checksums for `data_for_scripts_v0.3.4.tar.gz` and `pancancer_metabolomics_v.0.3.4.tar.gz`.
    - Extract and structure data for read-only mounting.
3.  **Skeleton Config:**
    - Create a minimal `docker-compose.yml` to test mounting logic and UIDs.
    - Validate that a test container can write to NVMe but only read from HDD.

## Success Criteria
- Data is present on the host and verified.
- Directory permissions allow the `1000:1000` user (student) to write to home/scratch.
- A manual `docker run` test confirms the mount logic.
