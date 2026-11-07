# containers/report/Dockerfile — NanoVI run report environment (REPORT process, label 'report')
# Installs containers/report/environment.yml into micromamba's base environment. Build from this folder:
#   docker build -t ghcr.io/microbialds/nanovi-report:<version> containers/report/

FROM mambaorg/micromamba:2.10.0-debian13-slim@sha256:e1476182cd40eb1f7c0cfae4adec9bfd78f3d7f05c61368306abfdc8939c1198

LABEL maintainer="NanoVI developers"
LABEL description="Run report environment for NanoVI (Typst, Jinja, matplotlib, Altair, Inter)"
LABEL version="1.2.0"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes && \
    rm /tmp/environment.yml

# Nextflow runs the task script with bash, not through the image's entrypoint, so the
# environment's programs go on PATH directly
ENV PATH=/opt/conda/bin:$PATH

WORKDIR /data
CMD ["python3"]
