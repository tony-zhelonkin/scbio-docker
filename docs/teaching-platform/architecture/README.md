# Teaching Platform Documentation

This folder contains the architectural design and setup instructions for transforming the `scbio-docker` environment into a multi-user teaching platform hosted on high-end workstation hardware.

## Project Goal
To provide a seamless, browser-based computational environment for 12 students (typically Metabolism PhD students) to perform multi-omics analysis without local installation overhead.

## Hardware Context: Dell Precision 7960
- **CPU:** Intel Xeon w9-3475X (36C/72T)
- **RAM:** 512GB DDR5
- **GPU:** NVIDIA RTX 5000 Ada (32GB)
- **Storage:** 1TB NVMe (System/Scratch), 16TB HDD (Data)

## Documentation Index

1.  **[Architecture Overview](architecture.md):** The high-level design choosing Docker-based orchestration over VMs.
2.  **[Resource Management](resource.md):** How to partition the 512GB RAM and 72 threads across students.
3.  **[IDE & AI Agent Integration](ide-and-ai.md):** Running Positron and CLI-based AI agents.
4.  **[JupyterHub Setup](setup-jupyterhub.md):** Implementation guide for the central management layer.
5.  **[Student Experience](student-experience.md):** Workflow for students from login to analysis.
6.  **[Data Persistence](data-strategy.md):** Handling shared datasets and student work.

## Implementation Roadmap
The project is divided into 5 distinct phases for implementation.
See the **[Implementation Plan](../plan/README.md)** for the step-by-step coding sessions.

## Quick Recommendation
**Use JupyterHub with DockerSpawner.**
It provides the leanest, most scalable architecture for a single-node high-performance workstation. It allows students to use RStudio, Jupyter, or VS Code through a single web portal while sharing the underlying hardware efficiently.
