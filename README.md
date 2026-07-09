# HRDTimer

A computational toolkit for estimating the timing of homologous recombination deficiency (HRD) in whole-genome–duplicated tumors.

# Installation

HRDTimer uses [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) to manage its environment. micromamba is a fast, standalone conda-compatible package manager with no base environment overhead — it is recommended over full Anaconda/Miniconda installs.

## 1. Install micromamba (if you don't have it)

Run the official one-liner installer:

```bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

Follow the prompts. When asked for an installation prefix, the default (`~/micromamba`) is fine for most users. On HPC clusters where your home directory has a storage quota, point it to a larger filesystem instead:

```bash
# Example for an HPC user — adjust the path to your own project space
MAMBA_ROOT_PREFIX=/scratch/myproject/micromamba "${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

After the installer finishes, open a new shell (or `source ~/.bashrc`) so the `micromamba` command is available.

## 2. Install the HRDTimer environment

Clone the repo and run the setup script:

```bash
git clone https://github.com/parklab/HRDTimer.git
cd HRDTimer
bash setup_environment.sh
```

This creates the `hrdtimer_env` micromamba environment and installs all dependencies (including MuSiCal and the reference genomes) **and** the `hrdtimer` package itself.

**HPC / shared filesystem variant** — if you want the environment stored on a large shared filesystem instead of your home directory, pass the desired path as the first argument:

```bash
bash setup_environment.sh /scratch/myproject/hrdtimer_env
```

After the script finishes, activate the environment with:

```bash
micromamba activate hrdtimer_env
# or, for the prefix form:
micromamba activate /scratch/myproject/hrdtimer_env
```

You can then `import hrdtimer` from anywhere, or run the notebooks in `notebooks/`.

### Optional overrides

The script respects the following environment variables if you need to customize paths:

| Variable | Default | Purpose |
|---|---|---|
| `MAMBA_ROOT_PREFIX` | `~/micromamba` | micromamba package cache / named-env root |
| `MUSICAL_DIR` | `~/tools/MuSiCal` | where the MuSiCal source is cloned |
| `MUSICAL_REF` | `v1.0.0` | MuSiCal git tag / branch to install |

Example:

```bash
MUSICAL_DIR=/scratch/myproject/MuSiCal bash setup_environment.sh
```

# Usage

After setting up and activating the environment (`micromamba activate hrdtimer_env`), HRDTimer can be run through the provided notebooks. A complete example workflow is available in the [example_workflow.ipynb](notebooks/example_workflow.ipynb) notebook.

> **Note:** One of HRDTimer's preprocessing steps involves running [MutationTimeR](https://github.com/gerstung-lab/MutationTimeR) on the VCFs of interest. For simplicity, the tutorial uses preprocessed example VCFs, so this step is already applied in the example data.  
> An example of how to implement MutationTimeR yourself can be found in the script: `/scripts/run_timer.R`.

The `example_workflow.ipynb` demonstrates typical HRDTimer analyses using the `example_data/` VCF files, including:

- **Preprocessing:** Assigning mutation-specific signatures and timing probabilities to individual mutations
- **Bootstrapping:** Estimating confidence in mutational signature analysis
- **Timing analyses:** Inferring the timing of Whole-Genome Duplication (WGD) and the onset of Homologous Recombination Deficiency (HRD) on a molecular time scale
- **Age translation:** Converting molecular timing to human age using inferred SBS1 acceleration models
