import sys, os, shutil, bz2, random, resource, warnings, subprocess, copy, re, glob, json
from subprocess import Popen, PIPE, STDOUT
import xml.etree.ElementTree as ET
from pathlib import Path
from Bio import Seq, SeqIO, SearchIO, SeqRecord, bgzf
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd
import numpy as np
from mpire import WorkerPool
import jellyfish


### Plotting imports ###
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import logomaker as lm
palette = list(mcolors.TABLEAU_COLORS.keys())
sns.set_theme(style="ticks", palette="muted")
sns.set_context("talk")
%matplotlib inline



def freq2ratio(freq):
    return(freq / (1 - freq))

def ratio2freq(ratio):
    return(ratio / (ratio + 1))


class TRNA_plot:
    pass

