"""NanoVI version: the one place the Python code keeps it.

Kept free of imports so that CI can read it without the pipeline's dependencies.
It must match manifest.version in nextflow.config, the top CHANGELOG release,
the Dockerfile LABEL and the image tag in conf/base.config (CI checks this).
"""
__version__ = "1.1.1"
