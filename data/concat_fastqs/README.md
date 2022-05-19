# concat_fastqs
This directory concatenates multiple FASTQ files into a single FASTQ file. We do this because some samples were sequenced multiple times, and we want to pool the sequencing reads before aligning.

## Input
The list of files to be concatenated is in the `input_fastq_list.csv` file. Files are specified with sample info, sequencing date, and a path to the Read 1 FASTQ file.

## Output
The concatentated FASTQ files are stored in this directory in the `concat_fastq_out` subdirectory. They are named for the sample.