# Asymmetric_Immiscibility_Simulations

[![Actions Status][actions-badge]][actions-link]
[![Documentation Status][rtd-badge]][rtd-link]

[![PyPI version][pypi-version]][pypi-link]
[![Conda-Forge][conda-badge]][conda-link]
[![PyPI platforms][pypi-platforms]][pypi-link]

[![GitHub Discussion][github-discussions-badge]][github-discussions-link]

<!-- SPHINX-START -->

![GitPic](https://github.com/SoftLivingMatter/Asymmetric_Immiscibility_Simulations/assets/68254269/9e1da8c1-3fbf-42c7-a747-d9f13ccc619d)


This is a repository of Python scripts associated with the manuscript
"Asymmetric oligomerization state and sequence patterning can tune multiphase
condensate miscibility" published on Bioxriv on March 12, 2023, and all
subsequent versions of this manuscript. It contains simulation scripts for MD
software [HOOMD-Blue 2.9.7](https://hoomd-blue.readthedocs.io/en/v2.9.7/) with
the plugin [azplugins](https://github.com/mphowardlab/azplugins/) to perform
direct coexistence NPAT simulations for estimating the relative miscibility of
model disordered proteins with oligomerization effects. 

## Installation
Because the current dependencies are not available on pypi or conda, we have to
build azplugins and hoomd from source before installing this package. 

**In the MyDella terminal:**

1. Create and load the environment (recommend putting 1st and 3rd line in a shell script):
```bash
module load anaconda3/2022.5
conda create --name ais python=3.9
conda activate ais
```

2. Within the Asymmetric_Immiscibility_Simulations directory, install
the necessary dependencies.  Installing this package provides the dependencies
to build hoomd, which is necessary for azplugins
```bash
pip install gsd==2.1.2
pip install -e .  # install ".[test]" for development dependencies

# in any directory
wget https://github.com/mphowardlab/azplugins/archive/refs/tags/v0.12.0.tar.gz -O azplugins-v0.12.0.tar.gz
tar -xzf azplugins-v0.12.0.tar.gz
wget https://github.com/glotzerlab/hoomd-blue/releases/download/v2.9.7/hoomd-v2.9.7.tar.gz
tar -xzf hoomd-v2.9.7.tar.gz
rm azplugins-v0.12.0.tar.gz hoomd-v2.9.7.tar.gz 

cd hoomd-v2.9.7/hoomd
ln -sr ../../azplugins-0.12.0/azplugins azplugins 
mkdir ../build && cd ../build
cmake ../ -DCMAKE_INSTALL_PREFIX=`python3 -c "import site; print(site.getsitepackages()[0])"`

make install -j4  # or more cores
```

3. Test
Once building finishes, check if import hoomd and azplugins work -
```bash
python -c 'import hoomd ; from hoomd import azplugins'
```
The above should produce no output any directory if hoomd is installed correctly.

## Example Usage

1. In folder InitialConfig, the script `GenInitConfig.py` creates an initial
   configuration by initializing a system of 147 KE1x3 star polymers and 441
   KE7 polymers in a cubic box, taking as input a pre-equilibrated configuration
   of a single KE1x3 polymer in `start_ke1_3arm.gsd`. An initial configuration
   file named `start_ke1.gsd` is created. 
2. This initial configuration file is then utilized by the script
   `SlabResize.py` in folder `NVT` to compress the cubic box to a size
   20nm<sup>3</sup> at constant temperature T=250K, following which the
   z-dimension of the simulation box is extended to 120nm by unwrapping the
   coordinates, to produce a configuration file `box2slab_extend_250.gsd`.
3. A direct coexistence NVT run is then performed by the script
   `SlabEquilibrate.py` at T=250K to equilibrate the system and produce
   coexisting dense and dilute phases. This produces a trajectory file named
   `restart_tmp1_250.gsd`.
4. The equilibrated simulation trajectory from the NVT run is then utilized to
   perform a direct coexistence NPAT simulation using script
   `SlabNPATProduction.py`. The final production simulation trajectory of the NPAT
   can then be analyzed to investigate the miscibility of model disordered protein
   sequences and its dependence on oligomerization state and sequence identity. 
   
## Acknowledgments

Code for setting up the HPS interaction model and NVT simulations was adapted
from [slab_builder](https://github.com/Roshan-M-Regy/slab_builder/tree/79283702a304556b46e53eeaede0f6a706299a86).
Regy, R. M.; Zheng, W.; Mittal, J. Theory of biological phase separation,
Liquid-Liquid Phase Coexistence and Membraneless Organelles. in Liquid-Liquid
Phase Coexistence and Membraneless Organelles (ed. Keating, C. D.) (2020).

<!-- prettier-ignore-start -->
[actions-badge]:            https://github.com/SoftLivingMatter/asymmetric-immiscibility-simulations/workflows/CI/badge.svg
[actions-link]:             https://github.com/SoftLivingMatter/asymmetric-immiscibility-simulations/actions
[conda-badge]:              https://img.shields.io/conda/vn/conda-forge/asymmetric-immiscibility-simulations
[conda-link]:               https://github.com/conda-forge/asymmetric-immiscibility-simulations-feedstock
[github-discussions-badge]: https://img.shields.io/static/v1?label=Discussions&message=Ask&color=blue&logo=github
[github-discussions-link]:  https://github.com/SoftLivingMatter/asymmetric-immiscibility-simulations/discussions
[pypi-link]:                https://pypi.org/project/asymmetric-immiscibility-simulations/
[pypi-platforms]:           https://img.shields.io/pypi/pyversions/asymmetric-immiscibility-simulations
[pypi-version]:             https://img.shields.io/pypi/v/asymmetric-immiscibility-simulations
[rtd-badge]:                https://readthedocs.org/projects/asymmetric-immiscibility-simulations/badge/?version=latest
[rtd-link]:                 https://asymmetric-immiscibility-simulations.readthedocs.io/en/latest/?badge=latest

<!-- prettier-ignore-end -->
