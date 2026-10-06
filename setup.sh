#!/bin/bash
set -e  # Exit on error

# Install packages
brew install cmake hdf5 tclap superlu

# Variables
ARMADILLO_VERSION="15.0.1"
ARMADILLO_DIR="armadillo-code"
INSTALL_PREFIX="$HOME/opt/armadillo"
BUILD_DIR="build"

# Clone or update repo
if [ -d "$ARMADILLO_DIR" ]; then
    echo "Armadillo source exists."
else
    echo "Cloning Armadillo $ARMADILLO_VERSION source..."
    git clone --branch "$ARMADILLO_VERSION" --depth 1 https://gitlab.com/conradsnicta/armadillo-code.git "$ARMADILLO_DIR"
fi

# Create build directory
cd "$ARMADILLO_DIR"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Configure and build
cmake .. -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" -DARMA_USE_SUPERLU=ON
make -j$(nproc)
make install

echo "Armadillo with SuperLU installed in $INSTALL_PREFIX"

