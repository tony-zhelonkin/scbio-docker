# Bioinformatics Docker Container 

This is a Docker image recipe for versatile Bioinformatics tools, enclosed in a Docker container. 

# Current development 

Current development successfull builld is `Dockerfile.dev`, which includes R with single-cell R libraries, radian for effective, and a minimal python venv setup. The docker is developed with the idea to work with it under Visual Studio Code Dev Container functionality, using remote server computational resources. 

# Change log 

## v0.1 build (current)
Successfull build with functional R setup and various R single-cell bioinformatics libraries. 

### Build problems (planning to resolve in the next build)

##### Problem #1
tidyverse didn\`t get installed. 
1. Need to add packages during build 

```dockerfile
RUN apt-get update && apt-get install -y \
    libfribidi-dev \
    libharfbuzz-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    pkg-config
```
2. Need to add these CRAN libraries to install `tidyverse`:
* textshaping
* ragg
****

##### Problem #2
R process run on a single core. 
**Solution**
```R
# Load parallel processing packages
library(future)
library(future.apply)

# Set up parallel processing
# You can adjust the number of cores (12 in your case)
plan("multicore", workers = 12)

# Enable parallel processing for Seurat operations
options(future.globals.maxSize = 8000 * 1024^2) # memory limit to 8GB for 32GB RAM system
options(future.globals.maxSize = 32000 * 1024^2) # Set to 32GB for 128GB RAM system
options(future.globals.maxSize = 24000 * 1024^2) # Set to 24GB for 62RAM system
```

****

# Future direction 

I plan to expand the functionality of the build to include: 
* various NGS tools for alignment, QC and preprocessing of sequencing data (from my training [RNAseq_pipelineDock](https://github.com/tony-zhelonkin/RNAseq_pipelineDock))
* Python single-cell tools 
* JB Genome Browser 