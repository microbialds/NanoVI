# containers/Dockerfile — NanoVI Python environment
# Builds the container used by vi_stats, vi_inference, vi_output,
# vi_database, vi_collapse, and vi_combine processes.

FROM python:3.11.17-slim-trixie

LABEL maintainer="NanoVI developers"
LABEL description="Python environment for NanoVI variational inference pipeline"
LABEL version="1.1.0"

# procps provides the ps that Nextflow uses to collect task metrics
RUN apt-get update && \
    apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies (exact pins; every package ships a wheel for 3.11)
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt && \
    rm /tmp/requirements.txt

# Set working directory
WORKDIR /data

# Default command
CMD ["python3"]
