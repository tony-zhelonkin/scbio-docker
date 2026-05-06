# IDE & AI Agent Integration

This document outlines how to integrate modern IDEs (Positron) and AI agents into the `scbio-docker` teaching platform.

## 1. Supporting Positron

[Positron](https://github.com/posit-dev/positron) is the next-generation data science IDE from Posit. It is built on top of Code OSS (VS Code) but specialized for R and Python.

### How to run Positron in the browser via JupyterHub:
You can use `jupyter-positron-server`, which leverages `jupyter-server-proxy` to tunnel Positron through the student's Jupyter session.

**Implementation Steps:**
1.  **Request Teaching License:** Contact `academic-licenses@posit.co` to get a free teaching license for the web-based Positron Server.
2.  **Modify Dockerfile:**
    - Install the proxy: `pip install jupyter-positron-server`
    - Download and extract the Positron Server binary to `/opt/positron`.
    - Add `/opt/positron/bin` to the `PATH`.
    - Place the `license.lic` at `/opt/positron/resources/activation/linux/x64/license.lic` (or set the environment variable).
3.  **Result:** Students will see a "Positron" icon in their JupyterLab Launcher. Clicking it opens a full IDE in a new tab.

---

## 2. Running AI Agents in the Terminal

The teaching platform is fully compatible with CLI-based AI agents (like Gemini CLI or Claude Code). Students can interact with them directly in the terminal within their IDE (Positron, VS Code, or RStudio).

### Prerequisites (Pre-installed in `scbio-docker`):
- **Node.js 20 LTS:** For running JS/TS-based agents and MCP servers.
- **`uv` / `uvx`:** For fast Python tool execution.
- **Python `toml`:** For configuration management.

### Setup for Students (Runtime):
Because AI agents often require project-specific configuration (like `.mcp.json` with absolute paths), we recommend a runtime setup script.

1.  **Initialize Project:** Use `./init-project.sh --with-submodules` to attach the `SciAgent-toolkit`.
2.  **Student Command:** Inside the container terminal, students run:
    ```bash
    ./01_modules/SciAgent-toolkit/scripts/setup-ai.sh --minimal
    ```
    *This takes ~2-3 minutes and installs the CLI agents into the container.*

3.  **Usage:**
    Students can then type `gemini` or `claude` in the terminal to start an agentic session.
    - **Example:** `"Gemini, help me write a maplet pipe to filter these 100 metabolites."`

### Security Note (API Keys):
For a classroom setting, **do not** bake API keys into the Docker image.
- **Recommended:** Provide students with a temporary API key to paste into their terminal, or use JupyterHub's environment variable feature to inject a key into every student container at startup.

```python
# jupyterhub_config.py
c.DockerSpawner.environment = {
    'GEMINI_API_KEY': 'your-teaching-key-here'
}
```

---

## 3. Recommended IDE Choice: Positron vs. VS Code

| Feature | Positron (Server) | VS Code (code-server) |
| :--- | :--- | :--- |
| **R Support** | Native (Ark Kernel), Data Explorer | Extension-based (httpgd) |
| **Python Support** | First-class | Industry standard |
| **AI Integration** | Terminal-based agents | Terminal + Extensions |
| **Licensing** | Academic License Required | Open Source (MIT) |

**Recommendation:** If you can obtain the academic license, **Positron** provides a "batteries-included" experience that PhD students (often familiar with RStudio) will find more intuitive than raw VS Code, while still supporting the terminal-based AI agents they need.
