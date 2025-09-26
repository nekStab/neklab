#!/bin/bash

# Welcome message
echo "LightKrylov Setup Script"
echo "--------------------"
echo "This script will perform the following actions upon your confirmation:"
echo "1. Remove existing LightKrylov directory (if found)."
echo "2. Clone the LightKrylov repository."
echo "3. Build and export LightKrylov."
echo "4. Run the LightKrylov unit-tests."
echo "--------------------"

# Detect the operating system
OS=$(uname)
should_clone="no"

# Install dependencies
read -p "Do you want to install/update necessary packages? (y/n): " confirm
if [ "$confirm" == "y" ] || [ "$confirm" == "Y" ]; then
    echo "Installing dependencies..."
    if [ "$OS" == "Linux" ]; then
        conda config --add channels conda-forge
        conda install fpm
    elif [ "$OS" == "Darwin" ]; then
        brew tap fortran-lang/fortran
        brew update
        brew install fpm
    else
        echo "Unsupported operating system. Exiting."
        exit 1
    fi
else
    echo "Skipping dependencies installation."
fi

# Check for existing Nek5000 directory
if [ -d "LightKrylov" ]; then
    read -p "Directory LightKrylov exists. Do you want to remove it? (y/n): " confirm
    if [ "$confirm" == "y" ] || [ "$confirm" == "Y" ]; then
        echo "Removing LightKrylov directory..."
        rm -rf LightKrylov
        should_clone="yes"
    else
        echo "LightKrylov directory will be retained."
    fi
else
    should_clone="yes"
fi

# Clone LightKrylov repository
if [ "$should_clone" == "yes" ]; then
    read -p "Do you want to clone the LightKrylov repository? (y/n): " confirm
    if [ "$confirm" == "y" ] || [ "$confirm" == "Y" ]; then
        echo "Cloning LightKrylov repository..."
        git clone https://github.com/nekStab/LightKrylov.git
        cd LightKrylov
        git checkout main
    else
        echo "Skipping cloning."
    fi
else
    cd LightKrylov
fi




############################################
######                                 #####
######     INSTALLING LIGHT KRYLOV     #####
######                                 #####
############################################

# Automatically clean-up previous builds.
fpm clean --all

# Sets the compilation options depending on the detected compiler.
if command -v mpiifort >/dev/null 2>&1; then
    FPM_FFLAGS="-Ofast -xHost -g -traceback"
    FPM_FC="mpiifort"
else
    FPM_FFLAGS="-march=native -O3 -funroll-loops -DMPI"
    FPM_FC="mpifort"
fi

export FPM_CC
export FPM_FFLAGS
export FPM_FC

# Install LightKrylov.
fpm install

# Run the unit tests.
read -p "Do you want to run the tests for LightKrylov? (y/n): " confirm
if [ "$confirm" == "y" ] || [ "$confirm" == "Y" ]; then
    fpm test
fi

# Finalize.
echo "LightKrylov setup complete."
