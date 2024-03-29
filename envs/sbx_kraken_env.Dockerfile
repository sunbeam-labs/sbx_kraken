FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/sbx_kraken_env

COPY envs/sbx_kraken_env.yml ./

# Install environment
RUN mamba env create --file sbx_kraken_env.yml --name sbx_kraken_env

ENV PATH="/opt/conda/envs/sbx_kraken_env/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "-n", "sbx_kraken_env", "/bin/bash", "-c"]

# Run
CMD "bash"