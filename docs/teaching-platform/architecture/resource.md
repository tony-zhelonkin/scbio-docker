# Resource Management & Partitioning

The Dell Precision 7960 is a powerhouse, but multi-omics analysis (especially single-cell) can be extremely memory-intensive. Proper partitioning ensures one student doesn't crash the entire machine.

## 1. Hardware Budget (Assumes 12 Students)

| Resource | Total | Per Student (Fair Share) | Strategy |
| :--- | :--- | :--- | :--- |
| **CPU Threads** | 72 | 6 | Soft limit (burst allowed) |
| **RAM** | 512GB | ~42GB | Hard limit (OOM protection) |
| **GPU VRAM** | 32GB | ~2.6GB | Shared (Time-sliced) |
| **Scratch Disk** | 1TB NVMe | ~80GB | Quota or shared scratch |

## 2. Docker/JupyterHub Configuration

JupyterHub allows setting resource limits in `jupyterhub_config.py`.

### Suggested Limits
```python
# c.Spawner.cpu_limit: Maximum cores a student can use
c.Spawner.cpu_limit = 8 

# c.Spawner.cpu_guarantee: Minimum cores reserved for a student
c.Spawner.cpu_guarantee = 2

# c.Spawner.mem_limit: Maximum RAM a student can use
# 40GB is extremely generous and handles almost all single-cell tasks.
c.Spawner.mem_limit = '40G'

# c.Spawner.mem_guarantee: Minimum RAM reserved
c.Spawner.mem_guarantee = '16G'
```

## 3. GPU Management (NVIDIA RTX 5000)

Since the RTX 5000 Ada supports Compute Capability 8.9 but **not** Multi-Instance GPU (MIG) which is reserved for A/H-series, we have two options for the GPU:

1.  **Shared Mode (Default):** All students see the same 32GB VRAM. If one student runs a large model, others might see "Out of Memory". This is usually fine for a 2.5h class where GPU usage is intermittent (e.g., training a small scVI model).
2.  **Fractional Assignment:** Using the [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/index.html) and Docker, you can limit VRAM visible to containers, though it's less strictly enforced than MIG.

## 4. Storage Strategy

- **System & Home (NVMe):** Faster IO for student scripts and temporary objects.
- **Large Data (HDDs):** The workshop dataset (~2GB) and any larger archival data reside on the 16TB HDD and should be mounted **Read-Only** to all containers.
    - Example: `/mnt/hdd/research_data/metabolism_v0.3.4` -> `/data/research:ro`
- **Cleanup:** Implement a script to prune containers and non-persistent volumes after the class.
