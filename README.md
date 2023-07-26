# Charge tRNA-Seq
This repository provides code and examples to process charge tRNA-Seq data as described in our [manuscript](www.somelink.com).

### Input data
This method takes in raw fastq paired-end reads.
To save disc space bzip2 compression is used and thus an input requirement.
Commonly, fastq files are gzipped but this can be changed using `bzip2`.
To queue files for parallel processing from gzipped to bzipped the following command can be used:  
`ls *.gz | parallel "gunzip -c {} | bzip2 > {.}.bz2"`


### Dependencies
The following needs to be installed and in the enviroment:
* AdapterRemoval v2 ([link](https://adapterremoval.readthedocs.io))
*  Have not been tested with adapter removal v3
* Swipe (https://github.com/torognes/swipe)
* makeblastdb (required to make masked tRNA database, so not super important)
* ImageMagick

Dependencies (Python packages):
* pandas
* Biopython
* numpy
* scipy
* seaborn
* matplotlib
* mpire
* jellyfish
* json_stream
* logomaker
* wand
Using some newer functions in numpy, so be aware to update python packages.



For makeblastdb install the command line BLAST tools:
https://www.ncbi.nlm.nih.gov/books/NBK569861/#intro_Installation


Only tested on Linux and MacOS.








