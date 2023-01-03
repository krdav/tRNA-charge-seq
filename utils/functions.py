import sys, os, subprocess, copy, shutil, re, glob, bz2, json, random, resource, warnings
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


def downsample_raw_input(sample_df, inp_file_df, NBdir, data_folder, seq_folder, downsample_absolute=1e5, downsample_fold=False, overwrite=True):
    '''
    This functions provides a way of downsampling the input files
    before read processing. This enables the user to test the entire
    sample processing pipeline at a 100 or 1000 fold reduced time.
    '''
    inp_file_df = copy.deepcopy(inp_file_df)
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
    inp_file_df = inp_file_df.drop(columns=cols[0:2]).rename(columns={cols[2]: cols[0], cols[3]: cols[1]})
    return(sample_df, inp_file_df, DS_folder)




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
    def __init__(self, sample_df, inp_file_df, NBdir, data_folder, seq_folder, AdapterRemoval_dir, MIN_READ_LEN, AR_threads=4):
        # Input:
        self.inp_file_df, self.NBdir, self.data_folder, self.seq_folder, self.AdapterRemoval_dir, self.MIN_READ_LEN = inp_file_df, NBdir, data_folder, seq_folder, AdapterRemoval_dir, MIN_READ_LEN
        # AdapterRomoval input:
        self.adapter1_tmp = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC<P7_index>ATCTCGTATGCCGTCTTCTGCTTG'
        self.adapter2_tmp = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT<P5_index>GTGTAGATCTCGGTGGTCGCCGTATCATT'
        self.AR_cmd_tmp = ["AdapterRemoval", "--bzip2", "--preserve5p", "--collapse", "--minalignmentlength", "10", "--threads", str(AR_threads)]

        # Check files exists before starting:
        for _, row in self.inp_file_df.iterrows():
            self.fnam_mate1 = '{}/{}/{}/{}'.format(NBdir, data_folder, seq_folder, row['fastq_mate1_filename'])
            self.fnam_mate2 = '{}/{}/{}/{}'.format(NBdir, data_folder, seq_folder, row['fastq_mate2_filename'])
            assert(os.path.exists(self.fnam_mate1))
            assert(os.path.exists(self.fnam_mate2))

    def make_dir(self, overwrite=True):
        # Create folder for files:
        self.AdapterRemoval_dir_abs = '{}/{}/{}'.format(self.NBdir, self.data_folder, self.AdapterRemoval_dir)
        try:
            os.mkdir(self.AdapterRemoval_dir_abs)
        except:
            if overwrite:
                shutil.rmtree(self.AdapterRemoval_dir_abs)
                os.mkdir(self.AdapterRemoval_dir_abs)
            else:
                raise Exception('Folder exists and overwrite set to false: {}'.format(self.AdapterRemoval_dir))
        return(self.AdapterRemoval_dir_abs)

    def run_parallel(self, n_jobs=4):
        os.chdir(self.AdapterRemoval_dir_abs)
        try:
            data = list(self.inp_file_df.iterrows())
            with WorkerPool(n_jobs=n_jobs) as pool:
                results = pool.map(self.__start_AR, data)
            self.__collect_stats()
            os.chdir(self.NBdir)
            return(self.inp_file_df)
        except Exception as err:
            os.chdir(self.NBdir)
            raise err
    
    def run_serial(self):
        os.chdir(self.AdapterRemoval_dir_abs)
        try:
            results = [self.__start_AR(index, row) for index, row in self.inp_file_df.iterrows()]
            self.__collect_stats()
            os.chdir(self.NBdir)
            return(self.inp_file_df)
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



class BC_split():
    '''
    This class is used to split fastq files based on barcodes.
    '''
    def __init__(self, sample_df, inp_file_df, NBdir, data_folder, AdapterRemoval_dir_abs, BC_dir):
        # Input:
        self.sample_df, self.inp_file_df, self.NBdir, self.data_folder, self.AdapterRemoval_dir_abs, self.BC_dir = sample_df, inp_file_df, NBdir, data_folder, AdapterRemoval_dir_abs, BC_dir

        # Check files exists before starting:
        for _, row in self.inp_file_df.iterrows(): # Pull out each merged fastq file
            basename = '{}-{}'.format(row['P5_index'], row['P7_index'])
            merged_fastq_fn = '{}/{}.collapsed.bz2'.format(AdapterRemoval_dir_abs, basename)
            assert(os.path.exists(merged_fastq_fn))

    def make_dir(self, overwrite=True):
        # Create folder for files:
        self.BC_dir_abs = '{}/{}/{}'.format(self.NBdir, self.data_folder, self.BC_dir)
        try:
            os.mkdir(self.BC_dir_abs)
        except:
            if overwrite:
                shutil.rmtree(self.BC_dir_abs)
                os.mkdir(self.BC_dir_abs)
            else:
                raise Exception('Folder exists and overwrite set to false: {}'.format(self.BC_dir_abs))
        return(self.BC_dir_abs)

    def run_parallel(self, n_jobs=4):
        # Must check that not too many file handles are opened at the same time:
        max_fh = resource.getrlimit(resource.RLIMIT_NOFILE)[0]
        fh_div = len(self.sample_df)*3 // max_fh
        try:
            assert(fh_div == 0)
        except:
            raise('Large number of samples, could cause problems with too many open filehandles under parallel processing. Either switch to serial processing, or implement parallel processing in smaller chunks.')

        # Run parallel:
        data = list(self.inp_file_df.iterrows())
        with WorkerPool(n_jobs=n_jobs) as pool:
            results = pool.map(self.__split_file, data)
        self.__collect_stats(results)
        return(self.sample_df, self.inp_file_df)

    def run_serial(self):
        results = [self.__split_file(index, row) for index, row in self.inp_file_df.iterrows()]
        self.__collect_stats(results)
        return(self.sample_df, self.inp_file_df)
    
    def __split_file(self, index, row):
        basename = '{}-{}'.format(row['P5_index'], row['P7_index'])
        merged_fastq_fn = '{}/{}.collapsed.bz2'.format(self.AdapterRemoval_dir_abs, basename)
        # List the barcodes and associated sample names:
        mask = (self.sample_df['P5_index'] == row['P5_index']) & (self.sample_df['P7_index'] == row['P7_index'])
        bc_fh = [(k, v, bz2.open('{}/{}.fastq.bz2'.format(self.BC_dir_abs, v), "wt")) for k, v in zip(self.sample_df[mask]['barcode_seq'].values, self.sample_df[mask]['sample_name_unique'].values)]
        unmapped_fh = bz2.open('{}/{}_unmapped.fastq.bz2'.format(self.BC_dir_abs, basename), "wt")

        # Collect stats:
        Nmapped = 0
        Nunmapped = 0
        Ncc = {k:0 for k in self.sample_df[mask]['sample_name_unique'].values}
        Ncca = {k:0 for k in self.sample_df[mask]['sample_name_unique'].values}
        Ntot = {k:0 for k in self.sample_df[mask]['sample_name_unique'].values}

        # Iterate over each record in the fastq file:
        with bz2.open(merged_fastq_fn, "rt") as input_fh:
            for title, seq, qual in FastqGeneralIterator(input_fh):
                # Search for barcodes and write to barcode specific file:
                found = False
                for bc, sample_name, fh in bc_fh:
                    if all(l1==l2 for l1, l2 in zip(seq[-len(bc):], bc) if l2 != 'N'):
                        found = True
                        # Add barcode sequence to title:
                        title = title + ':' + seq[-len(bc):]
                        fh.write("@{}\n{}\n+\n{}\n".format(title, seq[:-len(bc)], qual[:-len(bc)]))
                        Nmapped += 1
                        Ntot[sample_name] += 1
                        # Count if CC, CCA or not:
                        if seq[-(len(bc)+2):-len(bc)] == 'CC':
                            Ncc[sample_name] += 1
                        elif seq[-(len(bc)+3):-len(bc)] == 'CCA':
                            Ncca[sample_name] += 1
                        break
                if not found:
                    Nunmapped += 1
                    unmapped_fh.write("@{}\n{}\n+\n{}\n".format(title, seq, qual))
        for _, _, fh in bc_fh:
            fh.close()
        unmapped_fh.close()
        return(basename, Nmapped, Nunmapped, Ncc, Ncca, Ntot)

    def __collect_stats(self, results):
        # Unfold the results output:
        Ncc_union = dict()
        Ncca_union = dict()
        Ntot_union = dict()
        Nmapped_dict = dict()
        Nunmapped_dict = dict()
        for res in results:
            basename, Nmapped, Nunmapped, Ncc, Ncca, Ntot = res
            Nmapped_dict[basename] = Nmapped
            Nunmapped_dict[basename] = Nunmapped
            Ncc_union.update(Ncc)
            Ncca_union.update(Ncca)
            Ntot_union.update(Ntot)

        # Sort mapped/unmapped:
        Nmapped_list = list()
        Nunmapped_list = list()
        for _, row in self.inp_file_df.iterrows():
            basename = '{}-{}'.format(row['P5_index'], row['P7_index'])
            Nmapped_list.append(Nmapped_dict[basename])
            Nunmapped_list.append(Nunmapped_dict[basename])

        # Add stats to input file info dataframe:
        self.inp_file_df['N_BC-mapped'] = Nmapped_list
        self.inp_file_df['N_BC-unmapped'] = Nunmapped_list
        self.inp_file_df['N_sum-check'] = self.inp_file_df['N_BC-mapped'] + self.inp_file_df['N_BC-unmapped']
        self.inp_file_df['percent_BC-mapped'] = self.inp_file_df['N_BC-mapped'].values / self.inp_file_df['N_merged'].values *100
        # Dump stats as Excel file:
        self.inp_file_df.to_excel('{}/index-pair_stats.xlsx'.format(self.BC_dir_abs))

        # Add stats to sample info dataframe:
        self.sample_df['N_total'] = [Ntot_union[sn] for sn in self.sample_df['sample_name_unique']]
        self.sample_df['N_CC'] = [Ncc_union[sn] for sn in self.sample_df['sample_name_unique']]
        self.sample_df['N_CCA'] = [Ncca_union[sn] for sn in self.sample_df['sample_name_unique']]
        self.sample_df['N_CCA+CC'] = self.sample_df['N_CCA'].values + self.sample_df['N_CC'].values
        self.sample_df['CCA+CC_percent_total'] = self.sample_df['N_CCA+CC'].values / self.sample_df['N_total'].values *100
        self.sample_df['percent_CCA'] = self.sample_df['N_CCA'].values / self.sample_df['N_CCA+CC'].values *100
        # Dump stats as Excel file:
        self.sample_df.to_excel('{}/sample_stats.xlsx'.format(self.BC_dir_abs))






