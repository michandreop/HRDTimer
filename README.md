# HRDTimer

A computational toolkit for estimating the timing of homologous recombination deficiency (HRD) in whole-genome–duplicated tumors.

# Installation

You can set up the HRDTimer environment using `conda`. You can install `conda` through [Anaconda](https://www.anaconda.com/docs/main).
The full installation typically takes less than a minute on a standard laptop. The Conda environment and all required dependencies can be installed by running the provided setup script:

```bash
cd HRDTimer/hrdtimer_env
bash setup_environment.sh
```

> **Note:** HRDTimer is not packaged as a Python module — after activating the environment, simply open and run the notebooks in the `notebooks/` directory using Jupyter Notebook or VS Code.

# Usage

After setting up and activating the environment, HRDTimer can be run through the provided notebooks. A complete example workflow is available in the [example_workflow.ipynb](notebooks/example_workflow.ipynb) notebook.

> **Note:** One of HRDTimer's preprocessing steps involves running [MutationTimeR](https://github.com/gerstung-lab/MutationTimeR) on the VCFs of interest. For simplicity, the tutorial uses preprocessed example VCFs, so this step is already applied in the example data.  
> An example of how to implement MutationTimeR yourself can be found in the script: `/scripts/run_timer.R`.

The `example_workflow.ipynb` demonstrates typical HRDTimer analyses using the `example_data/` VCF files, including:

- **Preprocessing:** Assigning mutation-specific signatures and timing probabilities to individual mutations
- **Bootstrapping:** Estimating confidence in mutational signature analysis
- **Timing analyses:** Inferring the timing of Whole-Genome Duplication (WGD) and the onset of Homologous Recombination Deficiency (HRD) on a molecular time scale
- **Age translation:** Converting molecular timing to human age using inferred SBS1 acceleration models
