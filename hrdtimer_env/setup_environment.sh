#!/bin/bash

set -e

echo "🔧 Creating conda environment from environment.yml..."
conda env create -f environment.yml

echo " --- Activating environment 'hrdtimer_env'..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hrdtimer_env

echo " -- Cloning and installing MuSiCal from Park Lab GitHub..."
if [ ! -d "MuSiCal" ]; then
    git clone https://github.com/parklab/MuSiCal.git
fi

cd MuSiCal
pip install .
cd ..

echo " -- Installing SigProfilerMatrixGenerator reference genomes..."

python - <<EOF
from SigProfilerMatrixGenerator import install as genInstall

# Install GRCh37
genInstall.install('GRCh37', rsync=False, bash=True)

# Install GRCh38
genInstall.install('GRCh38', rsync=False, bash=True)
EOF

echo " ------------- All setup complete! HRDTimer env is ready."
