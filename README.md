# Deep mutational scanning of ZIKV NS5 protein
Experiments by Blake Richardson and Matt Evans.
Analysis by Jesse Bloom.

## Results
For a summary of the results, see [results/summary/dms_analysis.md](results/summary/dms_analysis.md), which is the Markdown summary of running the Jupyter notebook [dms_analysis.ipynb](dms_analysis.ipynb).

Other results are placed in [./results/](results), although not all files are tracked in the GitHub repo.

## Running analysis
First activate the [conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment for the analysis.
If you are using the *BloomLab* software on the Fred Hutch computing cluster, you can do this just with:

    conda activate ZIKV_DMS

Otherwise, first build the *ZIKV_DMS* [conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment from the [environment_pinned.yml](environment_pinned.yml) or [environment_unpinned.yml](environment_unpinned.yml) file (depending on whether you want fully pinned or unpinned versions), then activate it as above.

After you have activated the *ZIKV_DMS* [conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment, simply run the Python Jupyter notebook [dms_analysis.ipynb](dms_analysis.ipynb).

To run the notebook automatically and build the HTML summary linked to above, simply run the bash script [run_nbs.bash](run_nbs.bash).
On the Hutch cluster, you will first want to grab a node with 16 cores before doing this.

## Input data
The input data are in [./data/](data):

 - [./data/tile_1_amplicon.fasta](data/tile_1_amplicon.fasta): amplicon for the first tile of the barcoded-subamplicon sequencing.

 - [./data/subamplicon_alignspecs.txt](data/subamplicon_alignspecs.txt): the alignment specs for the [barcoded subamplicon sequencing](https://jbloomlab.github.io/dms_tools2/bcsubamp.html).

 - [./data/samplelist.csv](data/samplelist.csv): all the samples that we sequenced and the locations of the associated deep-sequencing data.

 - [./data/6WCZ.pdb](data/6WCZ.pdb) is the [6WCZ](https://www.rcsb.org/structure/6wcz) PDB file of ZIKV NS5 bound to human STAT2.
