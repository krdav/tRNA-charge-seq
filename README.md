# Charge tRNA-Seq
This repository provides code and examples to process charge tRNA-Seq data as described in our [manuscript](https://www.biorxiv.org/content/10.1101/2023.07.31.551363v1).
After installing the dependencies (see below) the best way to test the code is to run the example provided [here](projects/example/process_data.ipynb).
This is a minimal example of processing, going from raw reads to data analysis plots and shows how to use many of the functionalities provided.
To process your own samples, copy the example folder, rename it and use it as a boilerplate to fill in your own sample list and change the processing notebook to perform the processing/plotting you want.

Only tested on Linux and MacOS, probably does not work on Windows.



### Input data
The method takes in raw paired-end reads in fastq format.
To save disc space bzip2 compression is used and thus is an input requirement.
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


The following Python packages are required and can be install with `pip` or `conda` commands:
* jupyterlab
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

We recommend using an [Anaconda](https://www.anaconda.com/download) Python install which already has many of the above packages installed by default.


