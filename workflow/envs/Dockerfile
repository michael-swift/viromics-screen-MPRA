# Use a base image with Miniconda
FROM continuumio/miniconda3

# Install Mamba
RUN conda install mamba -n base -c conda-forge

# Copy the environment YAML file into the Docker image
COPY bioinfo.yaml /tmp/bioinfo.yaml

# Create the Conda environment from the YAML file using Mamba
RUN mamba env create -f /tmp/bioinfo.yaml

# Activate the environment and ensure it's activated by default
RUN echo "source activate bioinfo" > ~/.bashrc
ENV PATH /opt/conda/envs/bioinfo/bin:$PATH

# Set the entrypoint to bash
CMD ["/bin/bash"]