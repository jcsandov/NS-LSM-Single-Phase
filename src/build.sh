#!/bin/bash

# Load necessary modules
echo "Loading necessary modules..."
ml purge
ml iimpi/2018.04

# Source Spack environment
echo "Loading Spack environment..."
. ~/spack/share/spack/setup-env.sh

# Ensure Spack recognizes the Intel compiler
echo "Ensuring Spack recognizes Intel compiler..."
spack compiler find

# Load OpenBLAS built with the Intel 18.0.5 compiler
echo "Loading OpenBLAS..."
if ! spack load openblas %intel@18.0.5; then
  echo "Error: Failed to load OpenBLAS module."
  exit 1
fi

# Get OpenBLAS install directory
OPENBLAS_DIR=$(spack location -i 'openblas%intel@18.0.5')
if [ -z "$OPENBLAS_DIR" ]; then
  echo "Error: OpenBLAS directory not found."
  exit 1
fi

# Set LD_LIBRARY_PATH if necessary
export LD_LIBRARY_PATH=$OPENBLAS_DIR/lib:$LD_LIBRARY_PATH

# Run make
echo "Running make..."
make "$@"
