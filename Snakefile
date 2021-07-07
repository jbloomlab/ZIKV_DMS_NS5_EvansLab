"""``snakemake`` pipeline that runs analysis."""

configfile: 'config.yml'

rule all:
    input:
        expand("results/{tile}",
               tile=config['tiles']),
        expand("results/notebooks/dms_{tile}_analysis.ipynb",
               tile=config['tiles'])

rule dms_tile_analysis:
    """Analyze DMS data for a tile."""
    input:
        amplicon="data/{tile}_amplicon.fasta",
        alignspecs="data/{tile}_subamplicon_alignspecs.txt",
        samplelist="data/{tile}_samplelist.csv",
    output: resultsdir=directory("results/{tile}")
    params: errpre=lambda wc: config['tiles'][wc.tile]['errpre']
    threads: config['max_cpus']
    conda: 'environment.yml'
    log: notebook='results/notebooks/dms_{tile}_analysis.ipynb'
    notebook: 'dms_tile_analysis.py.ipynb'
