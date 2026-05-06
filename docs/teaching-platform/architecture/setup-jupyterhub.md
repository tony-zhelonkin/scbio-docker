# Implementation Guide: Setting up JupyterHub

This guide outlines the steps to install and configure JupyterHub on the Dell Precision 7960 to serve the `scbio-docker` image to students.

## 1. Prerequisites

- **Ubuntu 22.04/24.04** installed on the workstation.
- **Docker Engine** installed and running.
- **NVIDIA Container Toolkit** (if using the RTX 5000).
- **Node.js & NPM** (for the configurable-http-proxy).

## 2. Installation

```bash
# Install proxy
npm install -g configurable-http-proxy

# Install JupyterHub and DockerSpawner
pip install jupyterhub dockerspawner
```

## 3. JupyterHub Configuration (`jupyterhub_config.py`)

Create a configuration file to link JupyterHub with Docker.

```python
import os

c = get_config()

# 1. Use DockerSpawner to launch containers
c.JupyterHub.spawner_class = 'dockerspawner.DockerSpawner'

# 2. Point to the scbio-docker image
c.DockerSpawner.image = 'scdock-r-dev:v0.5.2'

# 3. Networking: JupyterHub must be reachable by containers
import netifaces
docker_ip = netifaces.ifaddresses('docker0')[netifaces.AF_INET][0]['addr']
c.JupyterHub.hub_ip = docker_ip

# 4. User Environment: Mount home directories and data
notebook_dir = os.environ.get('DOCKER_NOTEBOOK_DIR') or '/home/devuser'
c.DockerSpawner.notebook_dir = notebook_dir
c.DockerSpawner.volumes = {
    'jupyterhub-user-{username}': notebook_dir,
    '/mnt/data/research': {'bind': '/workspaces/project/00_Data', 'mode': 'ro'}
}

# 5. Resource Limits (from resource.md)
c.Spawner.mem_limit = '16G'
c.Spawner.cpu_limit = 4

# 6. Authentication (Workshop Mode: Dummy Auth)
# WARNING: Only for trusted internal networks!
c.JupyterHub.authenticator_class = 'dummy'
c.DummyAuthenticator.password = "metabolism2026"
```

## 4. Enabling RStudio & VS Code

To make RStudio and VS Code available in the browser, the `scbio-docker` image needs the following packages (usually installed via `pip` inside the container):

- `jupyter-rsession-proxy`: Routes RStudio Server through JupyterHub.
- `jupyter-server-proxy`: Routes any web service.
- `code-server`: For VS Code in the browser.

### RStudio Server Note
The `scbio-docker` image should have `rstudio-server` installed. If it doesn't, add it to the `Dockerfile`:
```dockerfile
RUN apt-get update && apt-get install -y rstudio-server
```

## 5. Starting the Hub

```bash
jupyterhub -f jupyterhub_config.py
```

Students can then navigate to `http://<workstation-ip>:8000`, log in with any username, and their personal Docker container will spin up automatically.
