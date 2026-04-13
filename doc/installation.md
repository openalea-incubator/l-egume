# Installation

You must use conda environment : <https://docs.conda.io/en/latest/index.html>

## Users

### Create a new environment with legume installed in there

```bash

mamba create -n legume -c openalea3 -c conda-forge  openalea.legume
mamba activate legume
```

Install legume in a existing environment

```bash
mamba install -c openalea3 -c conda-forge openalea.legume
```

### (Optional) Test your installation

```bash
mamba install -c conda-forge pytest
git clone https://github.com/openalea/legume.git
cd legume/test; pytest
```

## Developers

### Install From source

```bash
# Install dependency with conda
mamba env create -n phm -f conda/environment.yml
mamba activate legume

# Clone legume and install
git clone https://github.com/openalea/legume.git
cd legume
pip install .

# (Optional) Test your installation
cd test; pytest
```
