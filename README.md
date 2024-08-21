# Clonal evolution

This snakemake workflow infers clonal population structure from somatic single-nucleotide and indel variants identified by whole-genome sequencing of single or matched tumor samples from the same patient. Deconvolution is performed with [PyClone-VI](https://github.com/Roth-Lab/pyclone-vi) and clonal ordering with [Clonevol](https://github.com/hdng/clonevol). 

## Parameters

Set the parameters required for running the pipeline in the config file `config/config.yaml`. 

- Pyclone-VI parameters
    - `pyclone-c`: Refers to the original Pyclone-VI parameter `-c` that sets the number of clusters to use while fitting. 
    - `pyclone-d`: Refers to the original Pyclone-VI parameter `-d` that defines the probability density model. Options are `beta-binomial` and `binomial`. 
    - `pyclone-g`: Refers to the original Pyclone-VI parameter `-g` that sets the number of grid points used for approximating the posterior distribution. 
    - `pyclone-r`: Refers to the original Pyclone-VI parameter `-r` that sets the number of random restarts of variational inference.

- Pyclone-VI filtering parameters
    - `min_cluster_size`: Defines the minimum required number of mutations in a cluster to pass filtering. 

- ClonEvol parameters:
    - `clonevol_model`: Sets the ClonEvol cancer initiation model. Options are `polyclonal` and `monoclonal`. 

## Running on Uppmax

1. Load the conda module 

```
module load conda
source conda_init.sh
```

2. Activate snakemake environment

The Uppmax conda and snakemake modules are incompatible, so a conda environment with a snakemake installation. 

```
conda activate workflow/envs/snakemake
```

3. Run the pipeline
 
```
snakemake --cores 2 --sdm conda --conda-frontend conda
```

## To dos:

- [ ] Add Pyclone-VI seed parameter to config
- [ ] Sample_ids cannot include '-' for clonevols
- [ ] Add Filtering for multi-sample analysis
