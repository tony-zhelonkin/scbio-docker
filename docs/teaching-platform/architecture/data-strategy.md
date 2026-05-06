# Data Strategy & Persistence

Handling multi-omics data for 12 students requires a clear separation between shared reference data and student-generated results. While the workstation features a large 16TB HDD for archival storage, the specific workshop dataset is a high-density ~2GB collection.

## 1. Storage Layout on Workstation

| Path | Hardware | Use Case | Permissions |
| :--- | :--- | :--- | :--- |
| `/mnt/hdd/research_data` | 16TB HDD | **Workshop Dataset (v0.3.4, ~2GB)**: Raw Zenodo files, Seurat/SCE objects | Read-Only (Student), RW (Admin) |
| `/mnt/nvme/student_homes` | 1TB NVMe | Student scripts, `.Rhistory`, small figures | Read-Write (Student) |
| `/mnt/nvme/scratch` | 1TB NVMe | Temporary R/Python cache, temp files | Read-Write (Student) |

## 2. Shared Reference Data (Read-Only)

To prevent students from accidentally modifying or deleting the workshop dataset, it must be mounted with the `:ro` flag in Docker. The dataset comprises three main files from Zenodo (Record: 7874228, DOI: 10.5281/zenodo.7150252):
- `data_for_scripts_v0.3.4.tar.gz` (1.4 GB)
- `pancancer_metabolomics_v.0.3.4.tar.gz` (408 MB)
- `supplementary_dataset_v0.3.4.zip` (8.1 MB)

**JupyterHub Config Example:**
```python
c.DockerSpawner.volumes = {
    '/mnt/hdd/research_data/metabolism_v0.3.4': {
        'bind': '/data/research',
        'mode': 'ro'
    }
}
```

## 3. Persistent Student Work

Student work should persist even if the container is destroyed. We use named Docker volumes or host-mapped directories.

**Option A: Named Volumes (Recommended for isolation)**
```python
c.DockerSpawner.volumes = {
    'jupyterhub-user-{username}': '/home/devuser'
}
```

**Option B: Host-mapped (Recommended for easy admin access)**
```python
c.DockerSpawner.volumes = {
    '/mnt/nvme/student_homes/{username}': '/home/devuser'
}
```

## 4. Pre-loading the Environment

For a 2.5h class, you don't want students to run `git clone`. The `scbio-docker` image (or the persistent home directory) should be pre-populated.

**Teacher's Pre-Class Script:**
```bash
# For each student, pre-clone the class repos
for i in {01..12}; do
  mkdir -p /mnt/nvme/student_homes/student$i
  # Copy or clone the project template
  cp -r /workspaces/project /mnt/nvme/student_homes/student$i/
  chown -R 1000:1000 /mnt/nvme/student_homes/student$i
done
```

## 6. Capacity, Scaling & Operational Notes

### Resource Partitioning (Single-Host Orchestration)

The Dell Precision 7960 (96 Cores / 512GB RAM) is partitioned to ensure stability even if students run recursive or "leaky" code.

| Profile | CPU Limit | RAM Limit | Max Concurrent Users | Total Consumption |
| :--- | :--- | :--- | :--- | :--- |
| **Standard Class** | 6 Cores | 40 GB | 12 Students | 72 Cores / 480 GB RAM |
| **Scaling (Max)** | 4 Cores | 20 GB | 22 Students | 88 Cores / 440 GB RAM |
| **Instructor** | 12 Cores | 64 GB | 1 Admin | 12 Cores / 64 GB RAM |

### Scaling Strategy
If the class size increases beyond 12:
1.  **Memory Throttling**: The `c.DockerSpawner.mem_limit` in `jupyterhub_config.py` must be reduced (e.g., to `20G`) to prevent the host from hitting OOM (Out of Memory) and swapping to the HDD, which would degrade performance for everyone.
2.  **CPU Oversubscription**: 96 physical cores can handle 20+ students even with a 6-core limit because students rarely hit 100% utilization simultaneously. Docker will distribute cycles efficiently.

### Storage Efficiency (Docker Layering)
*   **Zero Duplication**: The `scdock-teaching:v1.0` image resides as a single set of read-only layers on the NVMe drive.
*   **Shared Data**: The `/data/research` mount (HDD) is a single physical copy shared by all 12-20 containers via pointer-based mounting.
*   **Total Disk Footprint**: ~20GB (Image) + ~2GB (Shared Data) + `(N students * ~500MB home dir)`. The environment remains lightweight regardless of class size.

### Authentication & Access
*   **Entrance**: `http://<HOST_IP>:8000`
*   **ID/Pass**: Users log in with any unique username; the `DummyAuthenticator` allows the instructor to set a single global password for the session.
*   **Persistence**: Containers are ephemeral, but the `/home/devuser` directory is mapped to the host NVMe. If a container crashes, the student simply logs back in to find their work exactly where they left it.
