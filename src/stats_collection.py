import sys, os, shutil, bz2, random, resource, warnings, subprocess, copy, re, glob, json, gc
from subprocess import Popen, PIPE, STDOUT
import json_stream
import xml.etree.ElementTree as ET
from pathlib import Path
from Bio import Seq, SeqIO, SearchIO, SeqRecord, bgzf
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd
import numpy as np
from mpire import WorkerPool
import jellyfish



class STATS_collection:
    '''
    This class is used to collect statistics from the
    alignment results.
    '''
    def __init__(self, dir_dict, tRNA_data, sample_df, common_seqs=None, ignore_common_count=False):
        self.stats_csv_header = ['readID', 'common_seq', 'sample_name_unique', 'sample_name', 'replicate', 'barcode', 'tRNA_annotation', 'align_score', 'unique_annotation', 'tRNA_annotation_len', 'align_5p_idx', 'align_3p_idx', 'align_5p_nt', 'align_3p_nt', 'codon', 'anticodon', 'amino_acid', '5p_cover', '3p_cover', '5p_non-temp', '3p_non-temp', '5p_UMI', '3p_BC', 'count']
        self.stats_agg_cols = ['sample_name_unique', 'sample_name', 'replicate', 'barcode', 'tRNA_annotation', 'align_score', 'unique_annotation', 'tRNA_annotation_len', 'align_5p_idx', 'align_3p_idx', 'align_5p_nt', 'align_3p_nt', '5p_cover', '3p_cover', '5p_non-temp', '3p_non-temp', 'codon', 'anticodon', 'amino_acid', 'count']
        self.stats_agg_cols2 = ['sample_name_unique', 'sample_name', 'replicate', 'barcode', 'tRNA_annotation', 'tRNA_annotation_len', 'unique_annotation', 'align_3p_nt', 'codon', 'anticodon', 'amino_acid', 'count']

        # Input:
        self.tRNA_data, self.sample_df = tRNA_data, sample_df
        self.dir_dict = dir_dict
        self.common_seqs_fnam = common_seqs
        self.common_seqs_info = dict() # Name to sequence for common sequence

        # Check files exists before starting:
        self.align_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['align_dir'])
        self.UMI_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['UMI_dir'])
        for _, row in self.sample_df.iterrows():
            SWres_fnam = '{}/{}_SWalign.json.bz2'.format(self.align_dir_abs, row['sample_name_unique'])
            assert(os.path.exists(SWres_fnam))
            trimmed_fn = '{}/{}_UMI-trimmed.fastq.bz2'.format(self.UMI_dir_abs, row['sample_name_unique'])
            assert(os.path.exists(trimmed_fn))
            common_obs_fn = '{}/{}_common-seq-obs.json'.format(self.align_dir_abs, row['sample_name_unique'])
            if not self.common_seqs_fnam is None:
                assert(os.path.exists(common_obs_fn))
            elif self.common_seqs_fnam is None and ignore_common_count is False and os.path.exists(common_obs_fn):
                raise Exception('Found common sequence counts for {}\n'
                                'See: {}\n'
                                'But common sequences are not specified. Either specify common sequences '
                                'or explicitly set ignore_common_count=True'.format(row['sample_name_unique'], common_obs_fn))

        # Make name to sequence dictionary for common sequences:
        if not self.common_seqs_fnam is None:
                # We can only allow one species if using common sequences.
                # Multiple species would require running the alignment on common sequences
                # several times, defeating the purpose, but also making the code much
                # more complicated.
                sp_set = set(self.sample_df['species'].values)
                if len(sp_set) > 1:
                    raise Exception('Only one species allowed in sample sheet when using common sequences.')
                self.common_seqs_sp = list(sp_set)[0]

                print('Using common sequences...')
                assert(self.common_seqs_fnam[-4:] == '.bz2')
                with bz2.open(self.common_seqs_fnam, "rt") as input_fh:
                    for ridx, record in enumerate(SeqIO.parse(input_fh, "fasta")):
                        assert(ridx == int(record.id))
                        self.common_seqs_info[record.id] = str(record.seq)

    def make_dir(self, overwrite=True):
        # Create folder for files:
        self.stats_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['stats_dir'])
        try:
            os.mkdir(self.stats_dir_abs)
        except:
            if overwrite:
                shutil.rmtree(self.stats_dir_abs)
                os.mkdir(self.stats_dir_abs)
            else:
                raise Exception('Folder exists and overwrite set to false: {}'.format(self.stats_dir_abs))
        return(self.stats_dir_abs)

    def run_parallel(self, n_jobs=4, verbose=True):
        self.verbose = verbose
        if self.verbose:
            print('Collecting stats from:', end='')
        # Run parallel:
        data = list(self.sample_df.iterrows())
        with WorkerPool(n_jobs=n_jobs) as pool:
            results = pool.map(self.__collect_stats, data)
        self.__concat_stats(results)
        return(self.concat_df)

    def run_serial(self, verbose=True):
        self.verbose = verbose
        if self.verbose:
            print('Collecting stats from:', end='')
        results = [self.__collect_stats(index, row) for index, row in self.sample_df.iterrows()]
        self.__concat_stats(results)
        return(self.concat_df)
    
    def __collect_stats(self, index, row):
        if self.verbose:
            print('  {}'.format(row['sample_name_unique']), end='')

        # Write stats to file:
        stats_fnam = '{}/{}_stats.csv.bz2'.format(self.stats_dir_abs, row['sample_name_unique'])
        stats_agg_fnam = '{}/{}_stats_aggregate.csv'.format(self.stats_dir_abs, row['sample_name_unique'])
        with bz2.open(stats_fnam, 'wt') as stats_fh:
            # Print header to stats CSV file:
            print(','.join(self.stats_csv_header), file=stats_fh)
            self.__read_non_common(row, stats_fh)
            if not self.common_seqs_fnam is None:
                self.__read_common(row, stats_fh)

        # Filter data and aggregate to count charged/uncharged tRNAs
        # Read stats from stats CSV file:
        with bz2.open(stats_fnam, 'rt') as stats_fh:
            # Use "keep_default_na=False" to read an empty string
            # as an empty string and not as NaN:
            stat_df = pd.read_csv(stats_fh, keep_default_na=False)

        # Aggregate dataframe and write as CSV file:
        # stat_df['count'] = np.zeros(len(stat_df))  # dummy for groupby count
        agg_df = stat_df.groupby(self.stats_agg_cols[:-1], as_index=False).agg({"count": "sum"})
        agg_df.to_csv(stats_agg_fnam, header=True, index=False)

        return(stats_agg_fnam)

    def __read_non_common(self, row, stats_fh):
        # Extract info from UMI processed reads:
        trimmed_fn = '{}/{}_UMI-trimmed.fastq.bz2'.format(self.UMI_dir_abs, row['sample_name_unique'])
        UMI_info = dict()
        with bz2.open(trimmed_fn, 'rt') as fh_bz:
            for UMIread in SeqIO.parse(fh_bz, "fastq"):
                # The last two strings are the adapter sequence and the UMI:
                _3p_bc, _5p_umi = UMIread.description.split()[-1].split(':')[-2:]
                seq = str(UMIread.seq)
                UMI_info[UMIread.id] = {
                    '_5p_umi': _5p_umi,
                    'seq': seq
                }

        # Open the alignment results:
        SWres_fnam = '{}/{}_SWalign.json.bz2'.format(self.align_dir_abs, row['sample_name_unique'])
        with bz2.open(SWres_fnam, 'rt', encoding="utf-8") as SWres_fh:
            # Parse JSON data as a stream (saves memory),
            # i.e. as a transient dict-like object
            SWres = json_stream.load(SWres_fh)
            # Loop through each read in the alignment results:
            for readID, align_dict in SWres.persistent().items():
                common_seq = False
                # Skip reads that were not aligned:
                if not align_dict['aligned']:
                    continue

                # Collect all the information:
                sample_name_unique = row['sample_name_unique']
                sample_name = row['sample_name']
                replicate = row['replicate']
                barcode = row['barcode']
                tRNA_annotation = align_dict['name']
                tRNA_annotation_first = tRNA_annotation.split('@')[0]
                align_score = align_dict['score']
                unique_annotation = '@' not in tRNA_annotation
                tRNA_annotation_len = self.tRNA_data[tRNA_annotation_first]['len']
                align_5p_idx, align_3p_idx = align_dict['dpos'][0]
                align_5p_nt = align_dict['qseq'][0]
                align_3p_nt = align_dict['qseq'][-1]

                # Move index for reads with 3' cleaved A:
                if align_3p_idx == (tRNA_annotation_len - 1) and align_3p_nt == 'C':
                    align_3p_idx += 1
                codon = self.tRNA_data[tRNA_annotation_first]['codon']
                anticodon = self.tRNA_data[tRNA_annotation_first]['anticodon']
                amino_acid = self.tRNA_data[tRNA_annotation_first]['amino_acid']
                _5p_cover = align_5p_idx == 1
                _3p_cover = align_3p_idx == tRNA_annotation_len

                # Extract non-template bases from UMI processed reads:
                try:
                    readUMI = UMI_info[readID]
                except KeyError:
                    raise Exception('Read ID ({}) not found among UMI trimmed sequences. Did any of the fastq headers change such that there is a mismatch between headers in the alignment json and those in the trimmed UMIs?'.format(readID))
                qpos = align_dict['qpos'][0]
                _5p_non_temp = readUMI['seq'][0:(qpos[0]-1)]
                _3p_non_temp = readUMI['seq'][qpos[1]:]
                _5p_umi = readUMI['_5p_umi']
                _3p_bc = row['barcode_seq']
                # For "non-common" sequences multiple reads
                # have not been collapsed:
                count = 1

                # Print line to output csv file:
                line_lst = [readID, common_seq, sample_name_unique, sample_name, replicate, \
                            barcode, tRNA_annotation, align_score, unique_annotation, \
                            tRNA_annotation_len, align_5p_idx, align_3p_idx, align_5p_nt, \
                            align_3p_nt, codon, anticodon, amino_acid, _5p_cover, _3p_cover, \
                            _5p_non_temp, _3p_non_temp, _5p_umi, _3p_bc, count]
                csv_line = ','.join(map(str, line_lst))
                print(csv_line, file=stats_fh)

        # Free memory from taken by "UMI_info":
        del UMI_info
        gc.collect()

    def __read_common(self, row, stats_fh):
        # Read common sequences observations for this sample:
        common_obs_fn = '{}/{}_common-seq-obs.json'.format(self.align_dir_abs, row['sample_name_unique'])
        with open(common_obs_fn, 'r') as fh_in:
            common_obs = json.load(fh_in)

        # Open the alignment results:
        SWres_fnam = '{}/{}_SWalign.json.bz2'.format(self.align_dir_abs, 'common-seqs')
        with bz2.open(SWres_fnam, 'rt', encoding="utf-8") as SWres_fh:
            # Parse JSON data as a stream (saves memory),
            # i.e. as a transient dict-like object
            SWres = json_stream.load(SWres_fh)
            # Loop through each read in the alignment results:
            for readID, align_dict in SWres.persistent().items():
                common_seq = True
                readID_int = int(readID)
                # Skip reads that were not aligned:
                if not align_dict['aligned']:
                    continue
                # And skip if this sample did not have this common sequence:
                elif common_obs[readID_int] == 0:
                    continue

                # Collect all the information:
                sample_name_unique = row['sample_name_unique']
                sample_name = row['sample_name']
                replicate = row['replicate']
                barcode = row['barcode']
                tRNA_annotation = align_dict['name']
                tRNA_annotation_first = tRNA_annotation.split('@')[0]
                align_score = align_dict['score']
                unique_annotation = '@' not in tRNA_annotation
                tRNA_annotation_len = self.tRNA_data[tRNA_annotation_first]['len']
                align_5p_idx, align_3p_idx = align_dict['dpos'][0]
                align_5p_nt = align_dict['qseq'][0]
                align_3p_nt = align_dict['qseq'][-1]

                # Move index for reads with 3' cleaved A:
                if align_3p_idx == (tRNA_annotation_len - 1) and align_3p_nt == 'C':
                    align_3p_idx += 1
                codon = self.tRNA_data[tRNA_annotation_first]['codon']
                anticodon = self.tRNA_data[tRNA_annotation_first]['anticodon']
                amino_acid = self.tRNA_data[tRNA_annotation_first]['amino_acid']
                _5p_cover = align_5p_idx == 1
                _3p_cover = align_3p_idx == tRNA_annotation_len

                # Extract non-template bases from common reads:
                try:
                    seq = self.common_seqs_info[readID]
                except KeyError:
                    raise Exception('Read ID ({}) not found among UMI trimmed sequences. '
                                    'Did any of the fastq headers change such that there is '
                                    'a mismatch between headers in the alignment json and those '
                                    'in the trimmed UMIs?'.format(readID))
                qpos = align_dict['qpos'][0]
                _5p_non_temp = seq[0:(qpos[0]-1)]
                _3p_non_temp = seq[qpos[1]:]
                _5p_umi = ''  # UMI information is lost when using common sequences
                _3p_bc = row['barcode_seq']
                # For common sequences the add the read count:
                count = int(common_obs[readID_int])

                # Print line to output csv file:
                line_lst = [readID, common_seq, sample_name_unique, sample_name, replicate, \
                            barcode, tRNA_annotation, align_score, unique_annotation, \
                            tRNA_annotation_len, align_5p_idx, align_3p_idx, align_5p_nt, \
                            align_3p_nt, codon, anticodon, amino_acid, _5p_cover, _3p_cover, \
                            _5p_non_temp, _3p_non_temp, _5p_umi, _3p_bc, count]
                csv_line = ','.join(map(str, line_lst))
                print(csv_line, file=stats_fh)

    def __concat_stats(self, csv_paths):
        # Concatenate all the aggregated stats csv files:
        stats_agg_fnam = '{}/ALL_stats_aggregate.csv'.format(self.stats_dir_abs)
        with open(stats_agg_fnam, 'w') as fh_out:
            # Grap header:
            with open(csv_paths[0], 'r') as fh_in:
                print(fh_in.readline(), file=fh_out)            
            for path in csv_paths:
                with open(path, 'r') as fh_in:
                    next(fh_in) # burn the header
                    for line in fh_in:
                        print(line, file=fh_out)

        # Apply some filtering and aggregate again
        # to get a smaller output dataframe:
        stats_agg2_fnam = '{}/ALL_stats_aggregate_filtered.csv'.format(self.stats_dir_abs)
        # Read and store as dataframe:
        with open(stats_agg_fnam, 'r') as stats_fh:
            concat_df = pd.read_csv(stats_fh, keep_default_na=False)

        # For charge to be determined the 3' must be covered
        # and have no 3' non-template bases:
        row_mask = (concat_df['3p_cover']) & (concat_df['3p_non-temp'] == '')
        concat_df = concat_df.loc[row_mask, self.stats_agg_cols2]
        self.concat_df = concat_df.groupby(self.stats_agg_cols2[:-1], as_index=False).agg({"count": "sum"})
        self.concat_df.to_csv(stats_agg2_fnam, header=True, index=False)








