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
* SWIPE — see [SWIPE aligner](#swipe-aligner) below
    * A faster fork is bundled in [`swipe/`](swipe/) and is recommended
    * The official SWIPE ([link](https://github.com/torognes/swipe)) also works
    * For Apple silicon use Rosetta2 to compile and run
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
    * Does not work for Apple silicon (unknown bug). Turn off using `stream=False` during alignment and stats collection
* logomaker
* wand
    * Only required to render and display pdf pages in the Jupyter notebook example
* natsort


We recommend using an [Anaconda](https://www.anaconda.com/download) Python install which already has many of the above packages installed by default.


### SWIPE aligner
The heavy lifting in the pipeline is the read-to-tRNA alignment, performed by
[SWIPE](https://github.com/torognes/swipe) (Smith-Waterman with inter-sequence
SIMD parallelisation).

A lightly modified fork is bundled in [`swipe/`](swipe/) and is **recommended**.
This workload aligns many short reads against a small tRNA database, and the
pipeline keeps only the best-scoring hit(s) for each read. The fork is tuned for
exactly this case and adds:
* `--best_only` — report only the best-scoring hit(s) per read, skipping the
  Smith-Waterman traceback and output for the lower-scoring hits that the
  pipeline would discard anyway;
* `--outfmt 10` — a compact, tab-separated output the pipeline reads directly,
  without an intermediate XML reformatting step;
* an inline single-thread search path, avoiding per-query thread overhead (for
  this workload more threads are slower, so parallelism is better applied across
  samples via the `n_jobs` argument of `run_parallel`).

On the example data (≈42k reads) this roughly halves the alignment time
(≈29 s → ≈16 s) and shrinks the intermediate output by ~90 % (≈207 MB → ≈16 MB).
Results are **identical** to the official SWIPE — the fork only avoids work the
pipeline does not use.

#### Installing the bundled fork
```
cd swipe
make
```
This builds a `swipe` executable in the `swipe/` directory. Put it on your
`PATH`, e.g.:
```
sudo cp swipe /usr/local/bin/swipe        # or copy/symlink anywhere on $PATH
```
On Apple silicon the bundled `Makefile` builds an x86_64 binary that runs under
Rosetta2; no changes are needed.

#### Using the official SWIPE instead
The official SWIPE works too. The pipeline auto-detects which `swipe` is on your
`PATH`: with the fork it uses the fast `--best_only --outfmt 10` path, and with
the official SWIPE it falls back to the standard XML output. The results are the
same either way. Detection can be overridden with the `swipe_compact` argument
of `SWIPE_align` (`'auto'` by default; `True` forces the fork path, `False`
forces the compatible XML path).


