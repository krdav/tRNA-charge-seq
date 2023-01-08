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


def sample_df_to_dict(sample_df):
    sample_dict = dict()
    for index, row in sample_df.iterrows():
        sample_dict[row['sample_name_unique']] = {
        'sample_name': row['sample_name'],
        'replicate': row['replicate'],
        'barcode': row['barcode'],
        'species': row['species']
        }
    return(sample_dict)

DNAcompRNA = {a: b for a, b in zip('ATGC', 'UACG')}
def anticodon2codon(anticodon):
    codon = ''.join([DNAcompRNA[b] for b in anticodon[::-1]])
    return(codon)

def read_tRNAdb_info(tRNA_database):
    # Read the tRNA database to find the length of each sequence:
    tRNA_data = dict()
    for species in tRNA_database:
        for record in SeqIO.parse(tRNA_database[species], "fasta"):
            tRNA_data[record.id] = dict()
            tRNA_data[record.id]['len'] = len(record.seq)
            tRNA_data[record.id]['codon'] = anticodon2codon(record.id.split('-')[2])
            tRNA_data[record.id]['anticodon'] = record.id.split('-')[2]
            tRNA_data[record.id]['amino_acid'] = record.id.split('-')[1]
    return(tRNA_data)

def index_to_sample_df(sample_df, index_df):
    # Read index sequences into dict:
    index_dict = dict()
    for t, i, s in zip(index_df['type'].values, index_df['id'].values, index_df['sequence'].values):
        if t not in index_dict:
            index_dict[t] = dict()
        index_dict[t][i] = s

    # Add index sequences to dataframe:
    sample_df['P5_index_seq'] = [index_dict['P5_index'][i] for i in sample_df['P5_index'].values]
    sample_df['P7_index_seq'] = [index_dict['P7_index'][i] for i in sample_df['P7_index'].values]
    sample_df['barcode_seq'] = [index_dict['barcode'][i] for i in sample_df['barcode'].values]
    return(sample_df)

def downsample_raw_input(sample_df, inp_file_df, NBdir, data_dir, seq_dir, downsample_absolute=1e5, downsample_fold=False, overwrite=True):
    '''
    This functions provides a way of downsampling the input files
    before read processing. This enables the user to test the entire
    sample processing pipeline at a 100 or 1000 fold reduced time.
    '''
    inp_file_df = copy.deepcopy(inp_file_df)
    # Check files exists before starting:
    for index, row in inp_file_df.iterrows():
        fnam_mate1 = '{}/{}/{}/{}'.format(NBdir, data_dir, seq_dir, row['fastq_mate1_filename'])
        fnam_mate2 = '{}/{}/{}/{}'.format(NBdir, data_dir, seq_dir, row['fastq_mate2_filename'])
        assert(os.path.exists(fnam_mate1))
        assert(os.path.exists(fnam_mate2))

    # Create a file tag/extension to indicate
    # that files are downsampled:
    if downsample_absolute:
        DS_ext = '_DSA-{}k'.format(round(downsample_absolute // 1000))
    elif downsample_fold:
        DS_ext = '_DSF-{}'.format(round(downsample_fold))
    DS_dir = seq_dir + DS_ext

    # Create folder for files:
    DS_dir_abs = '{}/{}/{}'.format(NBdir, data_dir, DS_dir)
    try:
        os.mkdir(DS_dir_abs)
    except:
        if overwrite:
            shutil.rmtree(DS_dir_abs)
            os.mkdir(DS_dir_abs)
        else:
            raise Exception('Folder exists and overwrite set to false: {}'.format(DS_dir_abs))

    # Do the downsampling:
    fnam_mate1_lst = list()
    fnam_mate2_lst = list()
    for index, row in inp_file_df.iterrows():
        fnam_mate1 = row['fastq_mate1_filename']
        fnam_mate2 = row['fastq_mate2_filename']
        fnam_mate1_in = '{}/{}/{}/{}'.format(NBdir, data_dir, seq_dir, fnam_mate1)
        fnam_mate2_in = '{}/{}/{}/{}'.format(NBdir, data_dir, seq_dir, fnam_mate2)
        mate1_in = bz2.open(fnam_mate1_in, "rt")
        mate2_in = bz2.open(fnam_mate2_in, "rt")

        fnam_mate1_DS = '.'.join(fnam_mate1.split('.')[0:-2]) + DS_ext + '.' + '.'.join(fnam_mate1.split('.')[-2:])
        fnam_mate1_lst.append(fnam_mate1_DS)
        fnam_mate2_DS = '.'.join(fnam_mate2.split('.')[0:-2]) + DS_ext + '.' + '.'.join(fnam_mate2.split('.')[-2:])
        fnam_mate2_lst.append(fnam_mate2_DS)
        fnam_mate1_out = '{}/{}/{}/{}'.format(NBdir, data_dir, DS_dir, fnam_mate1_DS)
        fnam_mate2_out = '{}/{}/{}/{}'.format(NBdir, data_dir, DS_dir, fnam_mate2_DS)
        mate1_out = bz2.open(fnam_mate1_out, "wt")
        mate2_out = bz2.open(fnam_mate2_out, "wt")

        Nsampled = 0
        for P1, P2, in zip(FastqGeneralIterator(mate1_in), FastqGeneralIterator(mate2_in)):
            if downsample_absolute:
                if P1[0].split(' ')[0] == P2[0].split(' ')[0]:
                    mate1_out.write("@{}\n{}\n+\n{}\n".format(*P1))
                    mate2_out.write("@{}\n{}\n+\n{}\n".format(*P2))
                    Nsampled += 1
                if Nsampled >= downsample_absolute:
                    break
            elif downsample_fold:
                if P1[0].split(' ')[0] == P2[0].split(' ')[0] and random.randint(1, downsample_fold) == 1:
                    mate1_out.write("@{}\n{}\n+\n{}\n".format(*P1))
                    mate2_out.write("@{}\n{}\n+\n{}\n".format(*P2))
        mate1_in.close()
        mate2_in.close()
        mate1_out.close()
        mate2_out.close()

    # Rename the mate 1/2 filenames:
    inp_file_df['fastq_mate1_filename_new'] = fnam_mate1_lst
    inp_file_df['fastq_mate2_filename_new'] = fnam_mate2_lst
    cols = ['fastq_mate1_filename', 'fastq_mate2_filename', 'fastq_mate1_filename_new', 'fastq_mate2_filename_new']
    # Merge:
    sample_df = sample_df.reset_index() # Preserve original order
    sample_df = sample_df.merge(inp_file_df[cols], on=['fastq_mate1_filename', 'fastq_mate2_filename']).sort_values(by=['index']).drop(columns=cols[0:2]+['index']).rename(columns={cols[2]: cols[0], cols[3]: cols[1]})
    sample_df = sample_df.reset_index(drop=True)
    inp_file_df = inp_file_df.drop(columns=cols[0:2]).rename(columns={cols[2]: cols[0], cols[3]: cols[1]})
    return(sample_df, inp_file_df, DS_dir)

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








