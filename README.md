# Charge tRNA-Seq
This repository provides code and examples to process charge tRNA-Seq data as described in our [manuscript](www.somelink.com).

(example)[projects/example/process_data.ipynb]

### Input data
This method takes in raw fastq paired-end reads.
To save disc space bzip2 compression is used and thus an input requirement.
Commonly, fastq files are gzipped but this can be changed using `bzip2`.
To queue files for parallel processing from gzipped to bzipped the following command can be used:  
`ls *.gz | parallel "gunzip -c {} | bzip2 > {.}.bz2"`


### Dependencies
The following needs to be installed and in the enviroment:
* AdapterRemoval v2 ([link](https://adapterremoval.readthedocs.io))
    * AdapterRemoval v3 may be working but this has not been tested
* SWIPE ([link](https://github.com/torognes/swipe))
* makeblastdb
    * Only required if making a new tRNA database
    * Install the command line BLAST tools ([link](https://www.ncbi.nlm.nih.gov/books/NBK569861/#intro_Installation))
* ImageMagick
    * Only required to render and display pdf pages in the Jupyter notebook example


The following Python packages and required and can be install with `pip` or `conda` commands:
* pandas
* Biopython
* numpy
    * Some newer functions are used so update
* scipy
* seaborn
* matplotlib
* mpire
* jellyfish
* json_stream
* logomaker
* wand
    * Only required to render and display pdf pages in the Jupyter notebook example






Only tested on Linux and MacOS.








