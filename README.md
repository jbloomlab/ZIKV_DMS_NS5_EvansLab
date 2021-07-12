# Deep mutational scanning of ZIKV NS5 protein
Experiments by Blake Richardson and Matt Evans.
Analysis by Jesse Bloom.

## Results
For a summary of the results, see [results/summary/](results/summary/), which has Markdown summaries for the analysis of each tile (e.g., [results/summary/dms_tile_1_analysis.md](results/summary/dms_tile_1_analysis.md), etc).

Other results are placed in [./results/](results), although not all files are tracked in the GitHub repo.

## Running analysis
First activate the *ZIKV_DMS_NS5_EvansLab* [conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment for the analysis.
If you have not already created this environment, build it from [environment.yml](ZIKV_DMS_NS5_EvansLab) with:

    conda env create -f environment.yml

Then activate the environment with:

    conda activate ZIKV_DMS_NS5_EvansLab

The analysis is run by the [snakemake](https://snakemake.readthedocs.io/) pipeline in [Snakefile](Snakefile).
Essentially, this pipeline runs the Jupyter notebook [dms_tile_analysis.ipynb](dms_tile_analysis.ipynb) for each deep mutational scanning tile, with the tile information specified in [config.yml](config.yml).
To run the pipeline using 36 jobs, use the command:

    snakemake -j 36

Add the `--keep-incomplete` flag if you don't want to delete results on an error.
To run on the Hutch cluster using `slurm`, do:

    sbatch -c 36 run_Snakemake.bash


## Input data
The input data are in [./data/](data):

 - `./data/tile_*_amplicon.fasta`: amplicons for each tile of the barcoded-subamplicon sequencing.

 - `./data/tile_*_subamplicon_alignspecs.txt`: the alignment specs for the [barcoded subamplicon sequencing](https://jbloomlab.github.io/dms_tools2/bcsubamp.html) for each amplicon.

 - `./data/tile_*_samplelist.csv`: all the samples that we sequenced and the locations of the associated deep-sequencing data for each amplicon.

 - [./data/6WCZ.pdb](data/6WCZ.pdb) is the [6WCZ](https://www.rcsb.org/structure/6wcz) PDB file of ZIKV NS5 bound to human STAT2.

## Visualizing results with `dms-view`
To visualize results with `dms-view`, first go to the `dms-view` homepage at [https://dms-view.github.io/](https://dms-view.github.io/).
Then you need to navigate to the protein structure and data you want to show, and paste those in the data box (upper left) and protein structure box (upper right).
Importantly, you need to enter the `Raw` values.
Specifically:
 1. Go to [https://github.com/jbloomlab/ZIKV_DMS_NS5_EvansLab/blob/main/results/dms-view/data_all_tiles.csv](https://github.com/jbloomlab/ZIKV_DMS_NS5_EvansLab/blob/main/results/dms-view/data_all_tiles.csv), click on `Raw`, and then enter the resulting link in the upper-right box for a CSV on the `dms-view` page.
 2. Go to [https://github.com/jbloomlab/ZIKV_DMS_NS5_EvansLab/blob/main/data/6WCZ.pdb](https://github.com/jbloomlab/ZIKV_DMS_NS5_EvansLab/blob/main/data/6WCZ.pdb), click on the `Raw` button, and enter the resulting link on the upper-right CSV box on `dms-view` for the protein structure.
 3. Go to [https://github.com/jbloomlab/ZIKV_DMS_NS5_EvansLab/blob/main/data/dms-view_manifest.md](https://github.com/jbloomlab/ZIKV_DMS_NS5_EvansLab/blob/main/data/dms-view_manifest.md), click on the `Raw` button, and enter the resulting link in the lower box for the manifest.
The resulting `dms-view` link can be copied and shared.
However, as long as this repository remains "private" the link will expire and we can't make a fully permanent one.
