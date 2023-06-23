# tRNA-charge-seq

To save disc space bzip2 compression is used and thus required as input.
To change common gzip compression files can be processed using the following:  
`ls *.gz | parallel "gunzip -c {} | bzip2 > {.}.bz2"`



Dependencies (installed and in Env):
* AdapterRemoval v2 (https://adapterremoval.readthedocs.io)
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






### Future dev
Make unit test with a few subsampled input files gather some of the final output and verify checksum.
- Set seed for adaper removal 







