#!/bin/bash

# Specify Python version (Modify as needed)
PYTHON_VERSION="python3.11"

# Check if the venv/ directory already exists
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    $PYTHON_VERSION -m venv venv
else
    echo "Virtual environment already exists. Skipping creation."
fi

# Activate virtual environment
source venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install dependencies
pip install -r requirements.txt

echo "✅ Virtual environment setup complete!"