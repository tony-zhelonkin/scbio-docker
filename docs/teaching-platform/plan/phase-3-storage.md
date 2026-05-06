# Phase 3: JupyterHub Orchestration

**Goal:** Implement the management layer to handle multiple students and enforce resource constraints.

## Tasks
1.  **Hub Installation:**
    - Install `jupyterhub` and `dockerspawner` on the host.
    - Set up the `DummyAuthenticator` with the global password `metabolism2026`.
2.  **Spawner Configuration:**
    - Link `DockerSpawner` to the `scdock-teaching:v1.0` image.
    - Implement resource limits: 6 CPU cores and 40GB RAM per container.
    - **Security Hardening**: Configure `DockerSpawner` with `cap_drop=["ALL"]` and `extra_host_config={"security_opt": ["no-new-privileges"]}`.
    - Enable NVIDIA GPU access via `extra_host_config`.
3.  **Security & Networking:**
    - Configure the Hub to listen on Port 8000.
    - Ensure the Hub can communicate with spawned containers via the Docker bridge.

## Success Criteria
- Navigating to Port 8000 prompts for login.
- Successful login spawns the student container with the correct 40GB/6C limits.
- `nvidia-smi` inside the spawned container shows the RTX 5000.
- `id` and `sudo` checks inside the container confirm restricted privileges.
