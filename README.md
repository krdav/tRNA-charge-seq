# tRNA-charge-seq

To save disc space bzip2 compression is used and thus required as input.
To change common gzip compression files can be processed using the following:  
`ls *.gz | parallel "gunzip -c {} | bzip2 > {.}.bz2"`



Dependencies (installed and in Env):
* AdapterRemoval (https://adapterremoval.readthedocs.io)
* Swipe (https://github.com/torognes/swipe)

Dependencies (Python packages):
* pandas
* Biopython
* numpy
* seaborn
* matplotlib



Using some newer functions in numpy, so be aware to update python packages.
