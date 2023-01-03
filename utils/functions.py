import sys, os, subprocess, copy, shutil, re, glob, bz2, json, random
from subprocess import Popen, PIPE, STDOUT
import xml.etree.ElementTree as ET
from pathlib import Path
from Bio import Seq, SeqIO, SearchIO, SeqRecord, bgzf
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd
import numpy as np
from mpire import WorkerPool



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


def downsample_raw_input(sample_df, NBdir, data_folder, seq_folder, downsample_absolute=1e5, downsample_fold=False, overwrite=True):
    '''
    This functions provides a way of downsampling the input files
    before read processing. This enables the user to test the entire
    sample processing pipeline at a 100 or 1000 fold reduced time.
    '''
    # Get filenames from the sample information:
    inp_file_df = sample_df[['fastq_mate1_filename', 'fastq_mate2_filename', 'P5_index', 'P7_index', 'P5_index_seq', 'P7_index_seq']].copy().drop_duplicates()

    # Check files exists before starting:
    for index, row in inp_file_df.iterrows():
        fnam_mate1 = '{}/{}/{}/{}'.format(NBdir, data_folder, seq_folder, row['fastq_mate1_filename'])
        fnam_mate2 = '{}/{}/{}/{}'.format(NBdir, data_folder, seq_folder, row['fastq_mate2_filename'])
        assert(os.path.exists(fnam_mate1))
        assert(os.path.exists(fnam_mate2))

    # Create a file tag/extension to indicate
    # that files are downsampled:
    if downsample_absolute:
        DS_ext = '_DSA-{}k'.format(round(downsample_absolute // 1000))
    elif downsample_fold:
        DS_ext = '_DSF-{}'.format(round(downsample_fold))
    DS_folder = seq_folder + DS_ext

    # Create folder for files:
    DS_folder_abs = '{}/{}/{}'.format(NBdir, data_folder, DS_folder)
    try:
        os.mkdir(DS_folder_abs)
    except:
        if overwrite:
            shutil.rmtree(DS_folder_abs)
            os.mkdir(DS_folder_abs)
        else:
            raise Exception('Folder exists and overwrite set to false: {}'.format(DS_folder_abs))

    # Do the downsampling:
    fnam_mate1_lst = list()
    fnam_mate2_lst = list()
    for index, row in inp_file_df.iterrows():
        fnam_mate1 = row['fastq_mate1_filename']
        fnam_mate2 = row['fastq_mate2_filename']
        fnam_mate1_in = '{}/{}/{}/{}'.format(NBdir, data_folder, seq_folder, fnam_mate1)
        fnam_mate2_in = '{}/{}/{}/{}'.format(NBdir, data_folder, seq_folder, fnam_mate2)
        mate1_in = bz2.open(fnam_mate1_in, "rt")
        mate2_in = bz2.open(fnam_mate2_in, "rt")

        fnam_mate1_DS = '.'.join(fnam_mate1.split('.')[0:-2]) + DS_ext + '.' + '.'.join(fnam_mate1.split('.')[-2:])
        fnam_mate1_lst.append(fnam_mate1_DS)
        fnam_mate2_DS = '.'.join(fnam_mate2.split('.')[0:-2]) + DS_ext + '.' + '.'.join(fnam_mate2.split('.')[-2:])
        fnam_mate2_lst.append(fnam_mate2_DS)
        fnam_mate1_out = '{}/{}/{}/{}'.format(NBdir, data_folder, DS_folder, fnam_mate1_DS)
        fnam_mate2_out = '{}/{}/{}/{}'.format(NBdir, data_folder, DS_folder, fnam_mate2_DS)
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
    sample_df = sample_df.merge(inp_file_df[cols], on=['fastq_mate1_filename', 'fastq_mate2_filename']).drop(columns=cols[0:2]).rename(columns={cols[2]: cols[0], cols[3]: cols[1]})
    
    return(sample_df, DS_folder)




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







class AR_merge():
    '''
    This class is used to merge the paired end reads using AdapterRomoval.
    '''
    def __init__(self, sample_df, NBdir, data_folder, seq_folder, AdapterRemoval_dir, MIN_READ_LEN, overwrite=True, AR_threads=4):
        # Input:
        self.NBdir, self.data_folder, self.seq_folder, self.AdapterRemoval_dir, self.MIN_READ_LEN = NBdir, data_folder, seq_folder, AdapterRemoval_dir, MIN_READ_LEN
        # AdapterRomoval input:
        self.adapter1_tmp = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC<P7_index>ATCTCGTATGCCGTCTTCTGCTTG'
        self.adapter2_tmp = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT<P5_index>GTGTAGATCTCGGTGGTCGCCGTATCATT'
        self.AR_cmd_tmp = ["AdapterRemoval", "--bzip2", "--preserve5p", "--collapse", "--minalignmentlength", "10", "--threads", str(AR_threads)]

        # Get filenames from the sample information:
        self.inp_file_df = sample_df[['fastq_mate1_filename', 'fastq_mate2_filename', 'P5_index', 'P7_index', 'P5_index_seq', 'P7_index_seq']].copy().drop_duplicates()

        # Check files exists before starting:
        for _, row in self.inp_file_df.iterrows():
            self.fnam_mate1 = '{}/{}/{}/{}'.format(NBdir, data_folder, seq_folder, row['fastq_mate1_filename'])
            self.fnam_mate2 = '{}/{}/{}/{}'.format(NBdir, data_folder, seq_folder, row['fastq_mate2_filename'])
            assert(os.path.exists(self.fnam_mate1))
            assert(os.path.exists(self.fnam_mate2))

        # Create folder for files:
        self.AdapterRemoval_dir_abs = '{}/{}/{}'.format(NBdir, data_folder, AdapterRemoval_dir)
        try:
            os.mkdir(self.AdapterRemoval_dir_abs)
        except:
            if overwrite:
                shutil.rmtree(self.AdapterRemoval_dir_abs)
                os.mkdir(self.AdapterRemoval_dir_abs)
            else:
                raise Exception('Folder exists and overwrite set to false: {}'.format(self.AdapterRemoval_dir))

    def run_parallel(self, n_jobs=4):
        os.chdir(self.AdapterRemoval_dir_abs)
        try:
            data = list(self.inp_file_df.iterrows())
            with WorkerPool(n_jobs=n_jobs) as pool:
                results = pool.map(self.__start_AR, data)
            self.__collect_stats()
            os.chdir(self.NBdir)
        except Exception as err:
            os.chdir(self.NBdir)
            raise err
    
    def run_serial(self):
        os.chdir(self.AdapterRemoval_dir_abs)
        try:
            results = [self.__start_AR(index, row) for index, row in self.inp_file_df.iterrows()]
            self.__collect_stats()
            os.chdir(self.NBdir)
        except Exception as err:
            os.chdir(self.NBdir)
            raise err

    def __start_AR(self, index, row):
        AR_cmd = self.AR_cmd_tmp.copy()
        basename = '{}-{}'.format(row['P5_index'], row['P7_index'])
        adapter1 = self.adapter1_tmp.replace('<P7_index>', row['P7_index_seq'])
        adapter2 = self.adapter2_tmp.replace('<P5_index>', row['P5_index_seq'])

        AR_cmd.extend(['--minlength', str(self.MIN_READ_LEN)])
        AR_cmd.extend(['--adapter1', adapter1])
        AR_cmd.extend(['--adapter2', adapter2])
        AR_cmd.extend(['--basename', basename])
        AR_cmd.extend(['--file1', '{}'.format(self.fnam_mate1)])
        AR_cmd.extend(['--file2', '{}'.format(self.fnam_mate2)])

        # Run AdapterRemoval, collect log in "file":
        with Popen(AR_cmd, stdout=PIPE, stderr=STDOUT) as p, open('{}_logfile.txt'.format(basename), 'a') as file:
            file.write('Starting subprocess with command:')
            file.write(str(AR_cmd))
            file.write('\n')
            for line in p.stdout:
                file.write(line.decode('utf-8'))
            file.write('\n****** DONE ******\n\n\n')
        return(1)

    def __collect_stats(self):
        N_pairs = list()
        N_merged = list()
        for _, row in self.inp_file_df.iterrows():
            basename = '{}-{}'.format(row['P5_index'], row['P7_index'])
            with open('{}.settings'.format(basename), 'r') as fh:
                for line in fh:
                    if 'Total number of read pairs:' in line:
                        N_pairs.append(int(line.split(':')[1][1:]))
                    if 'Number of full-length collapsed pairs:' in line:
                        N_merged.append(int(line.split(':')[1][1:]))
        # Write stats:
        self.inp_file_df['N_pairs'] = N_pairs
        self.inp_file_df['N_merged'] = N_merged
        self.inp_file_df['percent_successfully_merged'] = self.inp_file_df['N_merged'].values / self.inp_file_df['N_pairs'].values *100
        self.inp_file_df.to_excel('merge_stats.xlsx')







