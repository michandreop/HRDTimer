#!/bin/bash

set -e

echo "🔧 Creating conda environment from environment.yml..."
conda env create -f environment.yml

echo " --- Activating environment 'hrdtimer_env'..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hrdtimer_env

echo " -- Installing MuSiCal dependencies..."
conda install -y numpy scipy scikit-learn matplotlib pandas seaborn

echo " -- Cloning and installing MuSiCal from Park Lab GitHub..."
if [ ! -d "MuSiCal" ]; then
    git clone https://github.com/parklab/MuSiCal.git
fi

cd MuSiCal
pip install .
cd ..

echo " ------------- All setup complete! HRDTimer env is ready."
