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



### Future dev
Make unit test with a few subsampled input files gather some of the final output and verify checksum.
- Set seed for adaper removal 


Could de-duplicate the reads before running SWIPE and just store a tag with the count
- Could do this by running "find shared" to find the seperate out shared sequences
- Keep an account of how many of these shared each file had (write to JSON)
- Run the shared file for alignment
- Include this share alignment and the accounting JSON in the stats collection (would have to add a count column to the stats.csv.bz2 files)
- The goal will be to make a dictionary like: dict[seq][sample_id] = count
- Could extract all fastq.bz2 files and write their sequences to a single file
- Then sort this file and loop through it. If the last line was equal to the current then dict[seq] = dict()
- Dump these sequences as the "shared sequences" file
- When preparing files for SWIPE converting from fastq to fasta, exclude sequences found in dict[seq] and add the count.




