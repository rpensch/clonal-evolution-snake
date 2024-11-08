# Config

## Input file format

The variant files used as input for PyClone-VI should follow the requirements described [here](https://github.com/Roth-Lab/pyclone-vi). 

## Define input files

Define input files to be used in the pipeline with the config parameter `input_files`. The tsv-file should include two columns named `sample` and `input_file`. 

`sample` has to be an id that describes the individual that the sample or the multiple samples (e.g. samples from the same patient but multiple timepoints) are assigned to. 

The column `input_file` should include the path to the input file created according to PyClone-VI specifications. 

## Parameters

Set the parameters required for running the pipeline in the config file `config/config.yaml`. 

- Pyclone-VI parameters
    - `pyclone-c`: Refers to the original Pyclone-VI parameter `-c` that sets the number of clusters to use while fitting. 
    - `pyclone-d`: Refers to the original Pyclone-VI parameter `-d` that defines the probability density model. Options are `beta-binomial` and `binomial`. 
    - `pyclone-g`: Refers to the original Pyclone-VI parameter `-g` that sets the number of grid points used for approximating the posterior distribution. 
    - `pyclone-r`: Refers to the original Pyclone-VI parameter `-r` that sets the number of random restarts of variational inference.

- Pyclone-VI filtering parameters
    - `min_cluster_size`: Defines the minimum required number of mutations in a cluster to pass filtering. 
    - `min_founder_size`: Defines the minimum size of the founding clone as a fraction of the total number of mutations (e.g. 0.10)
    - `min_cap`: Defines the minimum median cluster assignment probability (CAP) of a cluster to pass filtering (e.g. 0.5)
    - `adjust_ccf`: The cellular prevalence of sublcones can be ajusted downwards so that no subclone can have a higher cellular prevalence than the founder clone in any sample. Set to `True` to apply the adjustment or `False` to keep the original cellular prevalence values.