#!/usr/bin/env bash
# bioconda build script for Strainify
set -euo pipefail

SHARE="$PREFIX/share/strainify"
mkdir -p "$SHARE" "$PREFIX/bin"

# Copy pipeline assets
cp -r main.nf nextflow.config conf modules subworkflows bin "$SHARE/"

# Install the CLI launcher
install -m 0755 strainify "$PREFIX/bin/strainify"
