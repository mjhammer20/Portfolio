Variant Calling Analysis - Docker container

This directory contains a minimal Docker container with a reproducible environment for running variant calling analysis, including common bioinformatics tools using a conda environment. Additionally, a snakemake workflow for running variant calling analysis on a single sample is included.

Files:
- Dockerfile - builds an image with a conda env named `vca` containing common tools
- environment.yml - conda environment specification
- entrypoint.sh - activates the conda env then runs the provided command
- docker-compose.yml - convenience file for running the container with the current directory mounted
- .dockerignore - files to exclude from the build context

Quick commands (zsh):

# Build the Docker image (from project root)
docker build -t vca_image:latest .

# Run an interactive container with the project mounted
docker run --rm -it -v "$PWD":/work -w /work vca_image:latest

# Or with docker-compose
docker-compose up --build

Notes:
- The image uses Miniconda to create the `vca` environment from `environment.yml`.
- After the container starts, common bioinformatics tools like `bwa`, `gatk`, `samtools`, and `snakemake` will be available.

