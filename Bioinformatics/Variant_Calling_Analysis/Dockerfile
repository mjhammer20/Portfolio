# Dockerfile for Variant Calling Analysis (VCA) project
# Uses Miniconda to create a reproducible environment with common bioinformatics tools
FROM continuumio/miniconda3:latest

WORKDIR /work

# Copy conda environment and create it
COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && conda clean -afy

# Ensure the conda env binaries are on PATH
ENV PATH=/opt/conda/envs/vca/bin:$PATH

# Add a simple entrypoint that activates the env and runs commands
COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["bash"]
