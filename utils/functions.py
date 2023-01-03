import sys, os, subprocess, copy, shutil, re, glob, bz2, json
import xml.etree.ElementTree as ET
from pathlib import Path
from Bio import Seq, SeqIO, SearchIO, SeqRecord
import pandas as pd

def myfun():

    return(0)


def indices(lst, element):
    result = []
    offset = -1
    while True:
        try:
            offset = lst.index(element, offset+1)
        except ValueError:
            return result
        result.append(offset)


def fast_fasta_count(filename):
    '''See: https://stackoverflow.com/a/9631635'''
    def blocks(files, size=65536):
        while True:
            b = files.read(size)
            if not b: break
            yield b

    with open(filename, "r", encoding="utf-8", errors='ignore') as f:
        return(sum(bl.count(">") for bl in blocks(f)))



def fast_fastq_count_bz(filename):
    '''See: https://stackoverflow.com/a/9631635'''
    def blocks(files, size=65536):
        while True:
            b = files.read(size)
            if not b: break
            yield b

    with bz2.open(filename, 'rt', encoding="utf-8", errors='ignore') as f:
        return(sum(bl.count("@") for bl in blocks(f)))














