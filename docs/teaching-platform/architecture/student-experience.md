# Student Experience Workflow

The primary goal of the "Lean" architecture is to minimize friction for students with little computational background. This document describes the user journey.

## 1. Access & Login

1.  **URL:** Students receive a local URL (e.g., `http://192.168.1.50:8000`).
2.  **Login:** Using the "Dummy Authenticator," students enter a unique username (e.g., `student01`) and a shared workshop password.
3.  **Spawning:** Upon clicking "Start Server," JupyterHub pulls (or uses local) `scbio-docker` and starts a container. This takes < 10 seconds.

## 2. Choosing an IDE

Once the container starts, students see the JupyterLab interface. From the "Launcher" tab, they can choose:

- **RStudio:** Opens a full RStudio Server session in a new browser tab.
- **VS Code:** Opens `code-server` (VS Code in the browser) for a more developer-focused experience.
- **Terminal:** A `bash` terminal inside their container.
- **Notebooks:** Standard Jupyter notebooks with R/Python kernels.

## 3. Analysis Workflow

- **Data Access:** All required datasets (Zenodo, pancancer-metabolomics) are pre-mounted at `/workspaces/project/00_Data`. Students can see them but cannot delete them.
- **Code Execution:** Students open the pre-cloned repositories (e.g., `CAMP-shiny-app`) and run scripts using the provided snippets.
- **Plotting:** 
    - In **RStudio**, plots appear in the "Plots" pane.
    - In **VS Code**, `httpgd` routes plots to a browser window.
    - In **Jupyter**, plots are inline.

## 4. Saving Work

- Any files created in `/home/devuser` (or the configured `notebook_dir`) are saved to a persistent Docker volume on the workstation's NVMe drive.
- If a student's browser crashes or they log out, their container keeps running (for a configurable timeout), and their work is preserved.

## 5. Session Termination

- At the end of the class, students log out. 
- The teacher can use the JupyterHub Admin panel to "Stop All Servers," freeing up all hardware resources immediately.
