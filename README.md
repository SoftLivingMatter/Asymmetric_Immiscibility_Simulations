# Asymmetric_Immiscibility_Simulations

![GitPic](https://github.com/SoftLivingMatter/Asymmetric_Immiscibility_Simulations/assets/68254269/9e1da8c1-3fbf-42c7-a747-d9f13ccc619d)


This is a repository of Python scripts associated with the manuscript "Asymmetric oligomerization state and sequence patterning can tune multiphase condensate miscibility" published on Bioxriv on March 12, 2023, and all subsequent versions of this manuscript. It contains simulation scripts for MD software [HOOMD-Blue 2.9.7](https://hoomd-blue.readthedocs.io/en/v2.9.7/) with the plugin [azplugins](https://github.com/mphowardlab/azplugins/) to perform direct coexistence NPAT simulations for estimating the relative miscibility of model disordered proteins with oligomerization effects. 

## Example Usage

1. In folder InitialConfig, the script ``` GenInitConfig.py ``` creates an initial configuration by initializing a system of 147 KE1x3 star polymers and 441 KE7 polymers in a cubic box, taking as input a pre-equilibrated configuration of a single KE1x3 polymer in ``` start_ke1_3arm.gsd ```. An initial configuration file named ``` start_ke1.gsd ``` is created. 
2. This intial configuration file is then utilized to compress the cubic box to a size 20nm<sup>3</sup>, following which the z-dimension of the simulation box is extended to 120nm by unwrapping the coordinates, to produce a configuration file ``` box2slab_extend_250.gsd ```. 

## Acknowledgments

Code for setting up HPS model NVT simulations was adapted from [slab_builder](https://github.com/Roshan-M-Regy/slab_builder/tree/79283702a304556b46e53eeaede0f6a706299a86).
Regy, R. M.; Zheng, W.; Mittal, J. Theory of biological phase separation, Liquid-Liquid Phase Coexistence and Membraneless Organelles. in Liquid-Liquid Phase Coexistence and Membraneless Organelles (ed. Keating, C. D.) (2020).
