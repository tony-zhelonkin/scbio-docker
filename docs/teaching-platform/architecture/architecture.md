# Architecture Overview: Docker vs. VMs

For a teaching platform on a single high-end workstation (Dell Precision 7960), the choice of virtualization layer significantly impacts performance, scalability, and ease of management.

## 1. Comparison: Why Docker wins for this use case

| Feature | Virtual Machines (VMs) | Docker Containers (Recommended) |
| :--- | :--- | :--- |
| **Overhead** | High (Guest OS per student) | Minimal (Shared Host Kernel) |
| **Startup Time** | Minutes | Seconds |
| **RAM Usage** | Fixed allocation (Rigid) | Dynamic sharing (Flexible) |
| **GPU Access** | Complex (PCIe Passthrough/vGPU) | Native (NVIDIA Container Toolkit) |
| **Storage** | Virtual Disks (Opaque) | Bind Mounts (Transparent/Shared) |
| **Scale** | 10-15 students (RAM bound) | 12 students (Highly Performant) |

**Conclusion:** A VM-based architecture would waste ~60-80GB of RAM just running 12 copies of an OS. Docker allows nearly all 512GB of RAM to be used for actual data analysis, giving each student a massive ~40GB buffer.

## 2. Proposed Architecture: JupyterHub Stack

The "Lean" architecture relies on **JupyterHub** as the orchestration layer. Despite the name, it is not limited to Jupyter Notebooks; it serves as a gateway to any web-based IDE.

### Workflow Components
1.  **Gateway (JupyterHub):** Handles authentication (PAM, OAuth, or simple dummy auth for workshops) and user management.
2.  **Orchestrator (DockerSpawner):** When a student logs in, JupyterHub instructs the Docker engine to start a new container from the `scbio-docker` image.
3.  **User Environment (`scbio-docker`):**
    - **RStudio:** Accessed via `jupyter-rsession-proxy`.
    - **VS Code:** Accessed via `jupyter-server-proxy` + `code-server`.
    - **JupyterLab:** Native support.
4.  **Hardware Interface:**
    - **CPU/RAM:** Managed via Docker cgroups (limits/reservations).
    - **GPU:** Exposed via `--gpus all` (shared compute).

### Logical Diagram
```text
[ Students (Browser) ]
      |
      v
[ Port 8000: JupyterHub ]
      |
      +-- [ Container: Student 1 (scbio-docker) ] --> RStudio / VS Code
      +-- [ Container: Student 2 (scbio-docker) ] --> RStudio / VS Code
      +-- [ Container: Student N (scbio-docker) ] --> RStudio / VS Code
      |
      v
[ Shared Host Resources ]
(36C/72T CPU, 512GB RAM, RTX 5000, 16TB HDD)
```

## 3. Benefits for the Teaching Platform
- **Zero Install:** Students only need a web browser.
- **Identical Environments:** Every student has the exact same package versions and data paths.
- **Efficiency:** The system only consumes resources for active sessions.
- **Persistent Work:** Student home directories are mounted from the host, so their work survives container restarts.
