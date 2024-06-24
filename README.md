# Deep mutational scanning of ZIKV NS5 (RdRp) protein

Experiments were performed by Blake Richardson and Matt Evans. Analysis was performed by Jesse Bloom, David Bacsik, and Caroline Kikawa

## Overview

The Evan's lab performed DMS on the ZIKV NS5 (RdRp) protein using a [tiled subamplicon approach](https://jbloomlab.github.io/dms_tools2/bcsubamp.html). A mutation's effect on ZIKV growth was assessed by comparing the passaged virus library to the original plasmid stock. 

See [results/summary/](results/summary/) for a markdown summary of the results for each tile over the genome (e.g., [results/summary/dms_tile_1_analysis.md](results/summary/dms_tile_1_analysis.md), etc).

See [this file](results/all_tiles/alltiles_host_adaptation.csv) for the results from all tiles in a single table. The mutational scan was performed in two cell lines: Huh75 (human) and C636 (mosquito).

## Running the analysis

Activate the *ZIKV_DMS_NS5_EvansLab* [conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment if it exists. Otherwise, create it from the [environment.yml](environment.yml) file:

```bash
conda env create -f environment.yml
```

And activate the environment:

```bash
conda activate ZIKV_DMS_NS5_EvansLab
```

The analysis is run by the [snakemake](https://snakemake.readthedocs.io/) pipeline detailed in [Snakefile](Snakefile). This pipeline runs the Jupyter notebook [dms_tile_analysis.py.ipynb](dms_tile_analysis.ipynb) for each tile using the information specified in the [config.yml](config.yml) file, generating a markdown summary for each tile.

Run the pipeline with following command:

```bash
snakemake
```

If you've got access to the Hutch's cluster, run the bash script:

```bash
sbatch run_Snakemake.bash
```

## Input data

The input data are in [./data/](data):

 - `./data/tile_*_amplicon.fasta`: amplicons for each tile of the barcoded-subamplicon sequencing.

 - `./data/tile_*_subamplicon_alignspecs.txt`: the alignment specs for the [barcoded subamplicon sequencing](https://jbloomlab.github.io/dms_tools2/bcsubamp.html) for each amplicon.

 - `./data/tile_*_samplelist.csv`: all the samples that we sequenced and the locations of the associated deep-sequencing data for each amplicon.

 - [./data/6WCZ.pdb](data/6WCZ.pdb) is the [6WCZ](https://www.rcsb.org/structure/6wcz) PDB file of ZIKV NS5 bound to human STAT2.

 - [./data/NS5_STAT2_joined.pdb](data/NS5_STAT2_joined.pdb) is a PDB file provided by Matt Evans that manually combines several other relevant PDBs.
