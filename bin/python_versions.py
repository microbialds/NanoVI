#!/usr/bin/env python3
"""Print an nf-core style versions.yml block for a Python step.

Usage: python_versions.py "<task.process>" > versions.yml
"""
import platform
import sys
from importlib.metadata import PackageNotFoundError, version

PACKAGES = ["numpy", "pandas", "scipy", "pysam", "biopython"]


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: python_versions.py <process name>")
    # Quoted, so that YAML keeps versions such as 1.80 as strings
    lines = [f'"{sys.argv[1]}":', f'    python: "{platform.python_version()}"']
    for package in PACKAGES:
        try:
            lines.append(f'    {package}: "{version(package)}"')
        except PackageNotFoundError:
            lines.append(f"    {package}: not installed")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
