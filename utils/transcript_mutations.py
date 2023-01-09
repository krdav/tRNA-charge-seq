import sys, os, shutil, bz2, random, resource, warnings, subprocess, copy, re, glob, json, contextlib
from subprocess import Popen, PIPE, STDOUT
import xml.etree.ElementTree as ET
from pathlib import Path
from Bio import Seq, SeqIO, SearchIO, SeqRecord, bgzf
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import Bio.Data.CodonTable
import pandas as pd
import numpy as np
from mpire import WorkerPool
import jellyfish
from Bio import Align


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




class TM_analysis:
    '''
    
    '''
    def __init__(self, dir_dict, sample_df, pull_default=False):
        pass






