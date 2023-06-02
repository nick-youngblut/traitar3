#!/bin/bash

# OS install
sudo apt-get update -y
## tools
sudo apt-get install -y tree screen curl python3.8

# dependencies
#mamba install python=3.8 bioconda::hmmer scipy pandas numpy matplotlib
micromamba env create -f conda_env.yml -y

# package install
python -m pip install --upgrade pip
pip install pytest pytest-dependency pytest-console-scripts
pip install .
