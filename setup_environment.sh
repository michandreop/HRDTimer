#!/bin/bash

set -euo pipefail

# Usage:
#   bash setup_environment.sh                 # creates a named env 'hrdtimer_env'
#   bash setup_environment.sh /path/to/env    # creates the env at that prefix
#
# The prefix form is handy on HPC where you want the env on a large shared
# filesystem (e.g. /scratch/myproject) instead of your home directory.
#
# Optional environment variable overrides (set before running):
#   MAMBA_ROOT_PREFIX   where micromamba stores its package cache / envs
#                       default: $HOME/micromamba
#   MUSICAL_DIR         where the MuSiCal source is cloned
#                       default: $HOME/tools/MuSiCal
#   MUSICAL_REF         git tag / branch of MuSiCal to install
#                       default: v1.0.0

# Resolve the repo root from this script's location so it works from any cwd.
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"

ENV_PREFIX="${1:-}"

# Keep micromamba's root (and its package cache) off a quota-limited filesystem
# if you need to. Override before running if you want a different location.
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"

# Permanent location for the MuSiCal source clone. An editable (-e) install
# points here forever, so it must NOT be a temp dir. Override if you like.
MUSICAL_DIR="${MUSICAL_DIR:-$HOME/tools/MuSiCal}"
MUSICAL_REF="${MUSICAL_REF:-v1.0.0}"

# --- Locate the micromamba binary -------------------------------------------
# In interactive shells, `micromamba` is usually only a bash *function* defined
# in ~/.bashrc, which does not exist inside `bash script.sh`. So we resolve the
# real executable here. Order of preference:
#   1) $MAMBA_EXE  (exported by `micromamba shell init`)
#   2) a binary on PATH
#   3) the conventional location under $MAMBA_ROOT_PREFIX/bin
if [ -n "${MAMBA_EXE:-}" ] && [ -x "${MAMBA_EXE:-}" ]; then
    MM="$MAMBA_EXE"
elif command -v micromamba >/dev/null 2>&1 && [ -x "$(command -v micromamba)" ]; then
    MM="$(command -v micromamba)"
elif [ -x "$MAMBA_ROOT_PREFIX/bin/micromamba" ]; then
    MM="$MAMBA_ROOT_PREFIX/bin/micromamba"
else
    echo "ERROR: could not find the micromamba binary." >&2
    echo "Install it first (see README.md for instructions), then re-run this script." >&2
    exit 1
fi
echo "Using micromamba binary: $MM"

if [ -n "$ENV_PREFIX" ]; then
    echo "Creating micromamba environment at prefix: $ENV_PREFIX"
    # Strip the 'name:' line so it doesn't clash with -p/--prefix.
    tmp_yml="$(mktemp --suffix=.yml)"
    grep -v '^name:' environment.yml > "$tmp_yml"
    "$MM" create -y -p "$ENV_PREFIX" -f "$tmp_yml"
    rm -f "$tmp_yml"
    ENV_REF=(-p "$ENV_PREFIX")
else
    echo "Creating micromamba environment 'hrdtimer_env' from environment.yml..."
    "$MM" create -y -n hrdtimer_env -f environment.yml
    ENV_REF=(-n hrdtimer_env)
fi

# From here on we use `micromamba run` instead of `micromamba activate`.
# `run` executes each command inside the env directly, which is far more
# robust under `set -euo pipefail` than activating in a non-interactive shell.

echo " -- Installing MuSiCal ($MUSICAL_REF) from a source clone (editable)..."
# We clone the repo and install with -e instead of `pip install git+...`,
# because MuSiCal's packaging does not declare its musical/data/*.csv files,
# so a wheel install silently drops them (e.g. COSMIC_v3p2_SBS_WGS.csv).
# An editable install resolves against the clone's source tree, which has them.
if [ ! -d "$MUSICAL_DIR/.git" ]; then
    mkdir -p "$(dirname "$MUSICAL_DIR")"
    git clone --branch "$MUSICAL_REF" https://github.com/parklab/MuSiCal.git "$MUSICAL_DIR"
else
    echo "    (clone already exists at $MUSICAL_DIR; fetching $MUSICAL_REF)"
    git -C "$MUSICAL_DIR" fetch --tags --quiet
    git -C "$MUSICAL_DIR" checkout --quiet "$MUSICAL_REF"
fi
"$MM" run "${ENV_REF[@]}" pip install -e "$MUSICAL_DIR"

# Sanity check that the data files actually made it into the import path.
"$MM" run "${ENV_REF[@]}" python - <<'PYEOF'
import os, musical
p = os.path.join(os.path.dirname(musical.__file__), "data", "COSMIC_v3p2_SBS_WGS.csv")
assert os.path.exists(p), f"MuSiCal data file missing: {p}"
print("    OK: MuSiCal data files present ->", p)
PYEOF

echo " -- Installing the HRDTimer package (editable)..."
"$MM" run "${ENV_REF[@]}" pip install -e "$REPO_ROOT"

echo " -- Registering a Jupyter kernel (so the env shows up in VS Code / JupyterLab)..."
"$MM" run "${ENV_REF[@]}" pip install ipykernel
"$MM" run "${ENV_REF[@]}" python -m ipykernel install --user \
    --name hrdtimer_env --display-name "Python (hrdtimer_env)"

echo " -- Installing SigProfilerMatrixGenerator reference genomes..."
"$MM" run "${ENV_REF[@]}" python - <<'PYEOF'
from SigProfilerMatrixGenerator import install as genInstall

# Install GRCh37
genInstall.install('GRCh37', rsync=False, bash=True)

# Install GRCh38
genInstall.install('GRCh38', rsync=False, bash=True)
PYEOF

echo " ------------- All setup complete! HRDTimer env is ready."
echo ""
echo "To activate the environment in future sessions, run:"
if [ -n "$ENV_PREFIX" ]; then
    echo "    micromamba activate $ENV_PREFIX"
else
    echo "    micromamba activate hrdtimer_env"
fi
