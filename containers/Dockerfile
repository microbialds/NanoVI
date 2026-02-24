# containers/Dockerfile — NanoVI Python environment
# Builds the container used by vi_stats, vi_inference, vi_output,
# vi_database, vi_collapse, and vi_combine processes.

FROM python:3.8-slim

LABEL maintainer="NanoVI developers"
LABEL description="Python environment for NanoVI variational inference pipeline"
LABEL version="1.0.0"

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        procps \
        gcc \
        libc6-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt && \
    rm /tmp/requirements.txt

# Set working directory
WORKDIR /data

# Default command
CMD ["python3"]
