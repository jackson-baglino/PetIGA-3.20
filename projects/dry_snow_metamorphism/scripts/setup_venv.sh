#!/bin/bash

# Specify Python version (Modify as needed)
PYTHON_VERSION="python3.9"

# Check if the venv/ directory already exists
if [ ! -d "venv_DSM" ]; then
    echo "Creating virtual environment..."
    $PYTHON_VERSION -m venv venv_DSM
else
    echo "Virtual environment already exists. Skipping creation."
fi

# Activate virtual environment
source venv_DSM/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install dependencies
pip install -r requirements.txt

echo "✅ Virtual environment setup complete!"