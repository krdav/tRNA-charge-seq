import sys, os, shutil, bz2, copy, contextlib, gc
import pickle
from natsort import natsorted
from subprocess import Popen, PIPE, STDOUT
from Bio import SeqIO
from Bio import Align
import pandas as pd
import numpy as np
from mpire import WorkerPool
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

### Plotting imports ###
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import logomaker as lm
palette = list(mcolors.TABLEAU_COLORS.keys())
sns.set_theme(style="ticks", palette="muted")
sns.set_context("talk")



class TM_analysis:
    '''
    This class is used generate statistics over the observed
    transcript mutations i.e. the mismatches and gaps in the
    alignment between a read and its tRNA transcript.

    Keyword arguments:
    common_seqs -- bzip2 compressed fasta file of commonly observed sequences to avoid duplicated alignments (default None)
    ignore_common_count -- Ignore common count even if X_common-seq-obs.json filename exists (default False)
    pull_default -- Attempt to load sample_df from default path (default False)
    check_exists -- Check if required files exist before starting (default True)
    overwrite_dir -- Overwrite old transcript analysis folder if any exists (default False)
    use_UMIcount -- Use UMI counts instead of read counts (default True)
    verbose -- Verbose printing collection progress (default True)
    '''
    def __init__(self, dir_dict, sample_df, tRNA_database, \
                 pull_default=False, common_seqs=None, ignore_common_count=False, \
                 overwrite_dir=False, verbose=True, check_exists=True, use_UMIcount=True):
        self.stats_csv_header = ['readID', 'common_seq', 'sample_name_unique', \
                                 'sample_name', 'replicate', 'barcode', 'tRNA_annotation', \
                                 'align_score', 'unique_annotation', 'tRNA_annotation_len', \
                                 'align_5p_idx', 'align_3p_idx', 'align_5p_nt', 'align_3p_nt', \
                                 'codon', 'anticodon', 'amino_acid', '5p_cover', '3p_cover', \
                                 '5p_non-temp', '3p_non-temp', '5p_UMI', '3p_BC', 'count']
        self.stats_csv_header_type = [str, bool, str, str, int, str, str, int, str, int, \
                                      int, int, str, str, str, str, str, bool, bool, \
                                      str, str, str, str, int]
        self.stats_csv_header_td = {nam:tp for nam, tp in zip(self.stats_csv_header, self.stats_csv_header_type)}

        # Input:
        self.sample_df, self.tRNA_database = sample_df, tRNA_database
        self.dir_dict = dir_dict
        self.common_seqs_fnam = common_seqs
        self.common_seqs_dict = dict() # map common sequences to index
        self.char_str = 'ACGTUN-'
        self.char_list = [c for c in self.char_str]
        self.char_dict = {c: i for i, c in enumerate(self.char_str)}
        self.tr_muts_masked = None
        self.use_UMIcount = use_UMIcount
        if self.use_UMIcount:
            self.count_col = 'UMIcount'
        else:
            self.count_col = 'count'
        
        self.stats_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['stats_dir'])
        self.align_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['align_dir'])
        self.UMI_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['UMI_dir'])
        # Attempt to load sample_df from default path:
        if pull_default:
            sample_df_path = '{}/{}'.format(self.align_dir_abs, 'sample_stats.xlsx')
            if os.path.isfile(sample_df_path):
                self.sample_df = pd.read_excel(sample_df_path, index_col=0)
            else:
                raise Exception('Default "sample_df" could not be found: {}'.format(sample_df_path))

        # Check files exists before starting:
        if check_exists:
            for sp in tRNA_database:
                assert(os.path.exists(tRNA_database[sp]))
            for _, row in self.sample_df.iterrows():
                stats_fnam = '{}/{}_stats.csv.bz2'.format(self.stats_dir_abs, row['sample_name_unique'])
                assert(os.path.exists(stats_fnam))
                common_obs_fn = '{}/{}_common-seq-obs.json'.format(self.align_dir_abs, row['sample_name_unique'])
                if not self.common_seqs_fnam is None:
                    assert(os.path.exists(common_obs_fn))
                elif self.common_seqs_fnam is None and ignore_common_count is False and os.path.exists(common_obs_fn):
                    raise Exception('Found common sequence counts for {}\n'
                                    'See: {}\n'
                                    'But common sequences are not specified. Either specify common sequences '
                                    'or explicitly set ignore_common_count=True'.format(row['sample_name_unique'], common_obs_fn))

        # Dictionary to store mutation info for each transcript:
        self.tr_muts = dict()
        # Make a dictionary template to fill out
        # for each sample that is read:
        self.tr_muts_tmp = dict()
        self.longest_tRNA = 0 # length of the longest tRNA
        # Read the tRNA transcripts:
        for species in tRNA_database:
            self.tr_muts_tmp[species] = dict()
            for record in SeqIO.parse(tRNA_database[species], "fasta"):
                seq_len = len(record.seq)
                if seq_len > self.longest_tRNA:
                    self.longest_tRNA = seq_len
                self.tr_muts_tmp[species][record.id] = dict()
                self.tr_muts_tmp[species][record.id]['seq'] = str(record.seq)
                self.tr_muts_tmp[species][record.id]['seq_len'] = seq_len
                # Position specific count matrix:
                self.tr_muts_tmp[species][record.id]['PSCM'] = np.zeros((seq_len, len(self.char_list)))
                self.tr_muts_tmp[species][record.id]['RTstops'] = np.zeros(seq_len)
                self.tr_muts_tmp[species][record.id]['mut_freq'] = np.zeros(seq_len)
                self.tr_muts_tmp[species][record.id]['gap_freq'] = np.zeros(seq_len)

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

            if verbose:
                print('Using common sequences...')
            assert(os.path.exists(self.common_seqs_fnam))
            assert(self.common_seqs_fnam[-4:] == '.bz2')
            with bz2.open(self.common_seqs_fnam, "rt") as input_fh:
                for ridx, record in enumerate(SeqIO.parse(input_fh, "fasta")):
                    seq = str(record.seq)
                    assert(ridx == int(record.id))
                    assert(not seq in self.common_seqs_dict)
                    self.common_seqs_dict[seq] = ridx

        # Make output folder:
        self._make_dir(overwrite_dir)

    def _make_dir(self, overwrite):
        # Create folder for files:
        self.TM_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['TM_dir'])
        try:
            os.mkdir(self.TM_dir_abs)
        except:
            if overwrite:
                shutil.rmtree(self.TM_dir_abs)
                os.mkdir(self.TM_dir_abs)
            else:
                print('Folder exists and overwrite set to false... Doing nothing.')

    def find_muts(self, unique_anno=True, match_score=1, mismatch_score=-2, \
                  open_gap_score=-3, extend_gap_score=-2, n_jobs=4, verbose=True, \
                  sample_list=None, max_5p_non_temp=10):
        '''
        Find mutations, gaps and RT stops in the input samples.
        For each read in the sample a new alignment is generated
        to the transcripts previously defined as annotations.

        Keyword arguments:
        unique_anno -- Only use reads with a unique annotation (default True)
        match_score -- Match score in the alignment between the read and its reference transcript (default 1)
        mismatch_score -- Mismatch penalty in the alignment between the read and its reference transcript (default -1)
        open_gap_score -- Open gap penalty in the alignment between the read and its reference transcript (default -2)
        extend_gap_score -- Extend gap penalty in the alignment between the read and its reference transcript (default -1)
        n_jobs -- Number of subprocesses to spawn in parallel (default 4)
        verbose -- Print status on samples processed (default True)
        sample_list -- Unique sample names of the samples to be processed. If None, all samples are processed (default None)
        max_5p_non_temp -- (default 10)
        '''
        if verbose:
            print('Collecting stats from:', end='')
        if sample_list is None:
            sample_list = [row['sample_name_unique'] for _, row in self.sample_df.iterrows()]
        
        # Find mutations in the transcripts for each file:
        data = list()
        for idx, row in self.sample_df.iterrows():
            if row['sample_name_unique'] in sample_list:
                data.append((idx, row, unique_anno, max_5p_non_temp, \
                             match_score, mismatch_score, open_gap_score, \
                             extend_gap_score, verbose))
        with WorkerPool(n_jobs=n_jobs) as pool:
            results = pool.map(self._collect_transcript_muts, data)
        # Fill out the transcript mutations per sample:
        for unam_res in results:
            unam, res = unam_res
            self.tr_muts[unam] = res
    
    def pickle_muts_write(self, pickle_name='saved_muts.pickle'):
        '''
        Pickle transcript mutations found using the find_muts method.
        '''
        if len(self.tr_muts) == 0:
            print('No mutations to pickle. First run the find_muts method.')
            return()
        pickle_fnam = '{}/{}'.format(self.TM_dir_abs, pickle_name)
        with open(pickle_fnam, 'wb') as fh:
            pickle.dump(self.tr_muts, fh)

    def pickle_muts_read(self, pickle_name='saved_muts.pickle'):
        '''
        Read previously generated transcript mutations from pickle file.

        '''
        pickle_fnam = '{}/{}'.format(self.TM_dir_abs, pickle_name)
        if not os.path.exists(pickle_fnam):
            print('No filename found here: {}'.format(pickle_fnam))
            return()
        with open(pickle_fnam, 'rb') as fh:
            self.tr_muts = pickle.load(fh)

    def _collect_transcript_muts(self, index, row, unique_anno, max_5p_non_temp, \
                                 match_score, mismatch_score, open_gap_score, \
                                 extend_gap_score, verbose):
        '''
        Find mutations, gaps and RT stop in a single sample.
        First, define the read sequences to work on, then
        align these to their respective reference and
        finally extract the mutation, gap and RT stop info.
        '''
        if verbose:
            print('  {}'.format(row['sample_name_unique']), end='')

        species = row['species']
        # Dictionary to store mutation info for each transcript:
        tr_muts_sp = copy.deepcopy(self.tr_muts_tmp)
        
        # Initiate the aligner:
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.match_score = match_score
        aligner.mismatch_score = mismatch_score
        aligner.open_gap_score = open_gap_score
        aligner.extend_gap_score = extend_gap_score

        # Deduplicate and count the input reads for alignment:
        trimmed_fn = '{}/{}_UMI-trimmed.fastq.bz2'.format(self.UMI_dir_abs, row['sample_name_unique'])
        dedup_seq_count = dict()
        with bz2.open(trimmed_fn, 'rt') as fh_bz:
            for UMIread in SeqIO.parse(fh_bz, "fastq"):
                seq = str(UMIread.seq)
                if seq in dedup_seq_count:
                    dedup_seq_count[seq]['count'] += 1
                else:
                    dedup_seq_count[seq] = dict()
                    dedup_seq_count[seq]['count'] = 1
                    # One id is enough, the other will have the same alignment:
                    dedup_seq_count[seq]['id'] = UMIread.id

        # Write deduplicated sequences to temporary file.
        # This is to save memory:
        dedup_tmp = '{}/tmp_{}.txt'.format(self.TM_dir_abs, row['sample_name_unique'])
        with open(dedup_tmp, 'w') as fh_dedup_out:
            for seq in dedup_seq_count:
                print('{}\t{}\t{}'.format(seq, \
                                          dedup_seq_count[seq]['id'], \
                                          dedup_seq_count[seq]['count']), \
                      file=fh_dedup_out)
        dedup_seq_count = None
        del dedup_seq_count
        gc.collect()

        # Read the stats file to get the old alignment annotations:
        stats_fnam = '{}/{}_stats.csv.bz2'.format(self.stats_dir_abs, row['sample_name_unique'])
        with bz2.open(stats_fnam, 'rt', encoding="utf-8") as stats_fh:
            sample_stats = pd.read_csv(stats_fh, keep_default_na=False, dtype=self.stats_csv_header_td)
        # Mask sequences with undesirable features:
        row_mask = (sample_stats['3p_cover']) & (sample_stats['3p_non-temp'] == '') & \
                   (sample_stats['5p_non-temp'].apply(len) <= max_5p_non_temp) & \
                   ((sample_stats['align_3p_nt'] == 'A') | (sample_stats['align_3p_nt'] == 'C'))
        sample_stats = sample_stats[row_mask]
        # Collect the sequence count annotation(s):
        ID2anno = dict()
        for rid, tan, count in zip(sample_stats['readID'].values, \
                                   sample_stats['tRNA_annotation'].values, \
                                   sample_stats[self.count_col].values):
            ID2anno[rid] = {'count': count, 'anno': tan.split('@')}
        sample_stats = None
        del sample_stats
        gc.collect()

        with open(dedup_tmp, 'r') as fh_dedup_in:
            for line in fh_dedup_in:
                seq, seq_id, dedup_count = line.split('\t')

                # Check if common sequence, then readID has changed:
                if not self.common_seqs_fnam is None and seq in self.common_seqs_dict:
                    readID = str(self.common_seqs_dict[seq])
                else:
                    readID = seq_id

                # Get list of annotations and the sequence count.
                # Notice, the count can be either raw read count
                # or UMI corrected counts:
                if readID in ID2anno:
                    anno_list = ID2anno[readID]['anno']
                    seq_count = ID2anno[readID]['count']
                else:
                    # Skip unaligned reads:
                    continue
                # Skip if multiple annotations found but unique requested:
                if unique_anno and len(anno_list) > 1:
                    continue

                # Generate all alignments first to enable a total count:
                alignments_anno = list()
                for anno in anno_list:
                    # Generate alignments:
                    target = tr_muts_sp[species][anno]['seq']
                    alignments = aligner.align(target, seq)
                    alignments_anno.append(alignments)

                # If multiple alignments with the same score these should be weighted
                # so one read contributes with one observation:
                weight = 1.0 / sum(len(algn) for algn in alignments_anno) * seq_count
                for anno, alignments in zip(anno_list, alignments_anno):
                    # Generate alignments:
                    target = tr_muts_sp[species][anno]['seq']
                    for alignment in alignments:
                        # Extract the alignment coordinates:
                        t_cor, q_cor = alignment.aligned
                        if (t_cor[-1, -1] + 1) < tr_muts_sp[species][anno]['seq_len']:
                            # Sometimes the change in alignment reward/penalties between 
                            # SWIPE and what is specified for transcript mutation analysis
                            # makes the alignment shift. This is particularly the case
                            # when a large gap is opened in the middle of the read and
                            # SWIPE aligns the right part to the 3p end of the reference
                            # and the above alignment aligns the left part to some other
                            # place in the reference.
                            # These are very rare events (probably less the 1/100,000),
                            # so we simply skip them:
                            continue
                        # Initiate the character observation matrix.
                        # Note that the number of observations is determined by the
                        # weight variable, hence the number could be less than 1:
                        count_mat = np.zeros((len(target), len(self.char_list)))
                        # Find gaps:
                        for i in range(1, len(t_cor)):
                            for j in range(t_cor[i][0] - t_cor[i-1][1]):
                                count_mat[t_cor[i-1][1] + j, self.char_dict['-']] = weight

                        # Find mismatches/mutations:
                        mut_idx = list()
                        for tran, qran in zip(t_cor, q_cor):
                            for ti, qi in zip(range(*tran), range(*qran)):
                                count_mat[ti, self.char_dict[seq[qi]]] = weight
                        tr_muts_sp[species][anno]['PSCM'] += count_mat

        os.remove(dedup_tmp)
        # Convert the count matrix to dataframe,
        # calculate mutation/gap/RTstop frequencies and return:
        for anno in tr_muts_sp[species]:
            tr_muts_sp[species][anno]['PSCM'] = pd.DataFrame(tr_muts_sp[species][anno]['PSCM'], columns=self.char_list)
            if tr_muts_sp[species][anno]['PSCM'].max().max() > 0:
                # Mutation frequencies:
                freq_mut = self._calc_mut_freq(tr_muts_sp, anno, species, gap_only=False, min_count_show=0)
                tr_muts_sp[species][anno]['mut_freq'] = freq_mut
                # Gap frequencies:
                freq_gap = self._calc_mut_freq(tr_muts_sp, anno, species, gap_only=True, min_count_show=0)
                tr_muts_sp[species][anno]['gap_freq'] = freq_gap
                # Fix ends (if CC-end, due to oxidation/cleavage, turn into CCA-end):
                end_obs = tr_muts_sp[species][anno]['PSCM']['C'].values[-2]
                end_ar = np.zeros(len(self.char_list))
                A_idx = self.char_list.index('A')
                end_ar[A_idx] = end_obs
                seq_len = tr_muts_sp[species][anno]['seq_len']
                tr_muts_sp[species][anno]['PSCM'].loc[seq_len-1, :] = end_ar
                # RT stops:
                RTstops_arr = self._calc_RTstops(tr_muts_sp, anno, species)
                tr_muts_sp[species][anno]['RTstops'] = RTstops_arr

        return((row['sample_name_unique'], tr_muts_sp))
 
    def _combine_tr_muts(self, sample_list, freq_avg_weighted=True):
        '''
        Combine mutation data from multiple samples.
        '''
        sample_list_cp = copy.deepcopy(sample_list)
        tr_muts_combi = copy.deepcopy(self.tr_muts_tmp)

        avg_count = 0 # for cumulative average
        for unam, sp_muts in self.tr_muts.items():
            # Skip if sample name not requested:
            if not unam in sample_list:
                continue
            # Pop from sample list to track that all
            # sample names have been found:
            sample_list_cp.pop(sample_list_cp.index(unam))
            
            # Combine mutation count matrices:
            species = list(sp_muts.keys())[0]
            for anno in sp_muts[species]:
                # Skip if no observations:
                if sp_muts[species][anno]['PSCM'].max().max() == 0:
                    continue
                # Fill the dictionary for all samples:
                if tr_muts_combi[species][anno]['PSCM'].max().max() == 0:
                    tr_muts_combi[species][anno]['PSCM'] = sp_muts[species][anno]['PSCM'].copy()
                    tr_muts_combi[species][anno]['mut_freq'] = sp_muts[species][anno]['mut_freq']
                    tr_muts_combi[species][anno]['gap_freq'] = sp_muts[species][anno]['gap_freq']
                    tr_muts_combi[species][anno]['RTstops'] = sp_muts[species][anno]['RTstops']
                else:
                    tr_muts_combi[species][anno]['PSCM'] += sp_muts[species][anno]['PSCM'].copy()
                    # Take the cumulative average for the frequencies:
                    tr_muts_combi[species][anno]['mut_freq'] = (sp_muts[species][anno]['mut_freq'] + avg_count*tr_muts_combi[species][anno]['mut_freq']) / (avg_count+1)
                    tr_muts_combi[species][anno]['gap_freq'] = (sp_muts[species][anno]['gap_freq'] + avg_count*tr_muts_combi[species][anno]['gap_freq']) / (avg_count+1)
                    tr_muts_combi[species][anno]['RTstops'] = (sp_muts[species][anno]['RTstops'] + avg_count*tr_muts_combi[species][anno]['RTstops']) / (avg_count+1)
            avg_count += 1

            # Get the frequency average weighted by total observations:
            if freq_avg_weighted:
                for anno in sp_muts[species]:
                    # Skip if no observations:
                    if sp_muts[species][anno]['PSCM'].max().max() == 0:
                        continue
                    tr_muts_combi[species][anno]['mut_freq'] = self._calc_mut_freq(tr_muts_combi, anno, species, gap_only=False, min_count_show=0)
                    tr_muts_combi[species][anno]['gap_freq'] = self._calc_mut_freq(tr_muts_combi, anno, species, gap_only=True, min_count_show=0)
                    tr_muts_combi[species][anno]['RTstops'] = self._calc_RTstops(tr_muts_combi, anno, species)

        if len(sample_list_cp) > 0:
            print('Following samples could not be found and therefore not combined: {}'.format(str(sample_list_cp)))
        return(tr_muts_combi)

    def write_transcript_mut(self, sample_list=None, species='human', \
                             csv_name='tr-mut_matrix', \
                             data_type='mut', \
                             mask_min_count=1, \
                             right_align=True):
        '''
        Export transcript mutations (mutations, gaps or RTstops)
        as csv file.

        Keyword arguments:
        sample_list -- Unique sample names of the samples to be used. If None, all samples are used (default None)
        species -- Species to plot (default 'human')
        data_type -- Data type to export. Choose between 'mut', 'gap' or 'RTstops' (default 'mut')
        mask_min_count -- Minimum observations, otherwise data is masked (default 1)
        right_align -- Align data to the right side of the row using the longest transcript to determine the number of columns (default True)
        '''

        if data_type not in ['mut', 'gap', 'RTstops']:
            print('The data_type input variable needs to be either \'mut\', \'gap\' or \'RTstops\'. Found {}'.format(data_type))

        # Write both mutation frequencies and positional counts:
        muts_csv_fnam_abs = '{}/{}_frequencies.csv'.format(self.TM_dir_abs, csv_name)
        counts_csv_fnam_abs = '{}/{}_counts.csv'.format(self.TM_dir_abs, csv_name)
        muts_fh = open(muts_csv_fnam_abs, 'w')
        counts_fh = open(counts_csv_fnam_abs, 'w')
        # Sequence positions are indexed from 1:
        seq_pos_str = ','.join(['P'+str(i+1) for i in range(self.longest_tRNA)])
        print('sample_name_unique,tRNA_annotation,tRNA_len,{}'.format(seq_pos_str), file=muts_fh)
        print('sample_name_unique,tRNA_annotation,tRNA_len,{}'.format(seq_pos_str), file=counts_fh)

        # Keep track of the samples requested
        # to warn if any was not found:
        # Get the mutations combined for the requested samples:
        if sample_list is None:
            sample_list = list(self.tr_muts.keys())
        sample_list_cp = copy.deepcopy(sample_list)
        for unam, sp_muts in self.tr_muts.items():
            # Skip if sample name not requested:
            if not unam in sample_list_cp:
                continue
            # Pop from sample list to track that all
            # sample names have been found:
            sample_list_cp.pop(sample_list_cp.index(unam))
            for anno in sp_muts[species]:
                freq_mut_print = np.zeros(self.longest_tRNA, dtype=float)
                counts_all_print = np.zeros(self.longest_tRNA, dtype=int)
                tRNA_len = sp_muts[species][anno]['seq_len']
                if data_type == 'gap':
                    freq_mut = sp_muts[species][anno]['gap_freq']
                elif data_type == 'mut':
                    freq_mut = sp_muts[species][anno]['mut_freq']
                else:
                    freq_mut = sp_muts[species][anno]['RTstops']

                # Enforce a minimum number of observations:
                counts_all = sp_muts[species][anno]['PSCM'].sum(1)
                min_count_mask = counts_all >= mask_min_count
                if min_count_mask.sum() == 0:
                    continue
                freq_mut[~min_count_mask] = 0

                # Left or right align the mutation frequencies
                # compared to the longest tRNA:
                if right_align:
                    freq_mut_print[-len(freq_mut):] = freq_mut
                    counts_all_print[-len(counts_all):] = counts_all
                else:
                    freq_mut_print[:len(freq_mut)] = freq_mut
                    counts_all_print[:len(counts_all)] = counts_all

                print('{},{},{},{}'.format(unam, anno, tRNA_len, ','.join(map(str, freq_mut_print))), file=muts_fh)
                print('{},{},{},{}'.format(unam, anno, tRNA_len, ','.join(map(str, counts_all_print))), file=counts_fh)
        muts_fh.close()
        counts_fh.close()
        if len(sample_list_cp) > 0:
            print('Did not find all the samples requested: {}'.format(str(sample_list_cp)))

    def plot_transcript_logo(self, topN=30, species='human', plot_name='tr-mut_logos', \
                             mito=False, sample_list=None):
        '''
        Simple logo plot of the nucleotide/gap observations broken 
        down by tRNA transcript.
        Will merge data from all samples in the provided sample list.
        
        Keyword arguments:
        topN -- Number of transcripts to plot. Taken from a list ordered by the number of observations (default 30)
        species -- Species to plot (default 'human')
        plot_name -- Plot filename (default 'tr-mut_logos')
        sample_list -- List of samples to combine for plotting. If None, all samples with mutation data are used (default None)
        '''

        # Get the mutations combined for the requested samples:
        if sample_list is None:
            sample_list = list(self.tr_muts.keys())
        tr_muts_combi = self._combine_tr_muts(sample_list)

        # Sort according to observations:
        anno_sorted = self._sort_anno(tr_muts_combi, species, mito=mito)
        if topN > len(anno_sorted):
            topN = len(anno_sorted)

        # Print each plot to the same PDF file:
        logo_fnam = '{}/{}.pdf'.format(self.TM_dir_abs, plot_name)
        with PdfPages(logo_fnam) as pp:
            for i in range(topN):
                anno, obs = anno_sorted[i]
                # Stop when no more observations:
                if obs == 0:
                    break
                # logomaker prints a warning when encountering "N"
                # characters, we don't want that:
                with contextlib.redirect_stdout(None):
                    logo_plot = lm.Logo(tr_muts_combi[species][anno]['PSCM'], color_scheme='classic');
                logo_plot.ax.set_title(anno, fontsize=15)
                logo_plot.ax.set_xlabel("5' to 3'")
                logo_plot.ax.set_ylabel("Count");
                logo_plot.fig.tight_layout()
                pp.savefig(logo_plot.fig, bbox_inches='tight')
                plt.close(logo_plot.fig)

    def plot_transcript_cov(self, topN=50, species='human', plot_name='tr-cov_matrix', \
                            png_dpi=False, no_plot_return=False, mito=False, \
                            sort_rows=True, sample_list=None, RTstops=False, \
                            min_obs=100):
        '''
        Plot the coverage of each transcript, with the additional option
        of plotting the RT PCR fall off.
        Will merge data from all samples in the provided sample list.
        Returns the mutation matrix as a pandas dataframe, the figure axes
        and the transcript annotations sorted.

        Keyword arguments:
        topN -- Number of transcripts to plot. Taken from a list ordered by the number of observations (default 50)
        species -- Species to plot (default 'human')
        plot_name -- Plot filename (default 'tr-cov_matrix')
        png_dpi -- DPI for plot as PNG. If False, only PDF is made (default False)
        no_plot_return -- Do not return plot (default False)
        mito -- Only plot mitochondrial transcripts (default False)
        RTstops -- Plot RT fall off instead of coverage (default False)
        min_obs -- Minimum observations of a position to be shown. Only valid for plotting RT stops (default 100)
        sort_rows -- Sort the plotted rows using hierarchical clustering (default True)
        sample_list -- List of samples to combine for plotting. If None, all samples with mutation data are used (default None)
        '''
        # Get the mutations combined for the requested samples:
        if sample_list is None:
            sample_list = list(self.tr_muts.keys())
        tr_muts_combi = self._combine_tr_muts(sample_list)

        # Sort according to observations:
        anno_sorted = self._sort_anno(tr_muts_combi, species, mito=mito)
        if topN > len(anno_sorted):
            topN = len(anno_sorted)
        
        # Find the observations and insert them into a matrix:
        obs_mat = np.zeros((topN, self.longest_tRNA))
        anno_topN = [anno_sorted[i][0] for i in range(topN)]
        for i in range(topN):
            anno = anno_topN[i]
            seq_len = tr_muts_combi[species][anno]['seq_len']
            if mito and 'mito' not in anno:
                continue
            # Count all observations for a given positions:
            counts_all = tr_muts_combi[species][anno]['PSCM'].sum(1).values

            # RT stops are pre-calculated:
            if RTstops:
                RTstops_arr = tr_muts_combi[species][anno]['RTstops']
                zero_mask = counts_all < min_obs
                RTstops_arr[zero_mask] = 0
                obs_mat[i, -seq_len:] = RTstops_arr
            else:
                # Insert from right side so CCA is always indexed
                # at the highest index.
                # Normalize so coverage at 3p is 100:
                obs_mat[i, -seq_len:] = 100 * counts_all / counts_all[-1]

        # Transform mutation matrix to dataframe:
        seq_idx_str = list(map(str, range(obs_mat.shape[1])))
        obs_mat_df = pd.DataFrame(obs_mat, columns=seq_idx_str, index=anno_topN)

        # Visually, it looks better when rows have been sorted
        # according to some sort of similarity.
        # This can be done by clustering:
        if sort_rows and len(obs_mat) > 2:
            dist_mat = pdist(obs_mat)
            Z = linkage(dist_mat, 'ward', optimal_ordering=True)
            dn = dendrogram(Z, no_plot=True)
            sorted_anno = [anno_topN[i] for i in dn['leaves']]
            obs_mat_df_plot = obs_mat_df.loc[sorted_anno, :]
        elif type(sort_rows) == list:
            sorted_anno = [sanno for sanno in sort_rows if sanno in anno_topN]
            mut_mat_df_plot = mut_mat_df.loc[sorted_anno, :]
        else:
            sorted_anno = anno_topN
            mut_mat_df_plot = mut_mat_df.loc[sorted_anno, :]

        # Adjust the plot size and generate it:
        y_scaler = 15/60
        y_len = y_scaler*len(obs_mat_df_plot)
        if y_len < 4:
            y_len = 4
        fig, ax = plt.subplots(1, 1, figsize=(25, y_len))
        g1 = sns.heatmap(obs_mat_df_plot, yticklabels=True, xticklabels=True, vmin=0, ax=ax);
        ax.yaxis.set_tick_params(labelsize=11)
        ax.xaxis.set_tick_params(labelsize=11)
        ax.set_xlabel("5' to 3' (filled in from right to left)", size=16)
        
        # Write to PDF file:
        plot_fnam = '{}/{}.pdf'.format(self.TM_dir_abs, plot_name)
        fig.tight_layout()
        fig.savefig(plot_fnam, bbox_inches='tight')
        if not png_dpi is False:
            fig.savefig(plot_fnam[:-4]+'.png', bbox_inches='tight', dpi=png_dpi)

        if no_plot_return:
            plt.close(fig)
            return((obs_mat_df, None, sorted_anno))
        else:
            return((obs_mat_df, fig, sorted_anno))

    def plot_transcript_mut(self, topN=50, species='human', plot_name='tr-mut_matrix', \
                            png_dpi=False, no_plot_return=False, mito=False, gap_only=False, \
                            sort_rows=True, min_count_show=100, sample_list=None,
                            freq_avg_weighted=True):
        '''
        Plot a heatmap of the transcript mutations in one
        or multiple samples averaged.
        Will merge data from all samples in the provided sample list.
        Returns the mutation matrix as a pandas dataframe, the figure axes
        and the transcript annotations sorted.

        Keyword arguments:
        topN -- Number of transcripts to plot. Taken from a list ordered by the number of observations (default 50)
        species -- Species to plot (default 'human')
        plot_name -- Plot filename (default 'tr-mut_matrix')
        png_dpi -- DPI for plot as PNG. If False, only PDF is made (default False)
        no_plot_return -- Do not return plot (default False)
        mito -- Only plot mitochondrial transcripts (default False)
        gap_only -- Only plot gap frequency (default False)
        sort_rows -- Sort the plotted rows using hierarchical clustering (default True)
        min_count_show -- Minimum observations of a position to be shown (default 100)
        sample_list -- List of samples to combine for plotting. If None, all samples with mutation data are used (default None)
        freq_avg_weighted -- Calculate the mutation frequency averaged across samples, weighted by the number observations (default True)
        '''
        # Get the mutations combined for the requested samples:
        if sample_list is None:
            sample_list = list(self.tr_muts.keys())
        tr_muts_combi = self._combine_tr_muts(sample_list, freq_avg_weighted=freq_avg_weighted)

        # Sort according to observations:
        anno_sorted = self._sort_anno(tr_muts_combi, species, mito=mito)
        if topN > len(anno_sorted):
            topN = len(anno_sorted)
        # Find the mutations and insert them into a matrix:
        mut_mat = np.zeros((topN, self.longest_tRNA))
        anno_topN = [anno_sorted[i][0] for i in range(topN)]
        anno_skip = set() # Skip these because of lacking observations
        for i in range(topN):
            anno = anno_topN[i]
            if gap_only:
                freq_mut = tr_muts_combi[species][anno]['gap_freq']
            else:
                freq_mut = tr_muts_combi[species][anno]['mut_freq']

            # Enforce a minimum number of observations:
            counts_all = tr_muts_combi[species][anno]['PSCM'].sum(1)
            min_count_mask = counts_all >= min_count_show
            if min_count_mask.sum() == 0 or freq_mut is None:
                anno_skip.add(anno)
            else:
                freq_mut[~min_count_mask] = 0
                # Insert from right side so CCA is always indexed
                # at the highest index:
                mut_mat[i, -len(freq_mut):] = freq_mut

        # Filter out the entries with no positions above min_count_show:
        anno_idx_nozero = [aidx for aidx, anno in enumerate(anno_topN) if not anno in anno_skip]
        anno_topN_nozero = [anno for anno in anno_topN if not anno in anno_skip]
        if len(anno_idx_nozero) == 0:
            print('Nothing left after filtering')
            return(1)
        mut_mat = mut_mat[anno_idx_nozero, :]

        # Transform mutation matrix to dataframe:
        seq_idx_str = list(map(str, range(mut_mat.shape[1])))
        mut_mat_df = pd.DataFrame(mut_mat, columns=seq_idx_str, index=anno_topN_nozero)

        # Visually, it looks better when rows have been sorted
        # according to some sort of similarity.
        # This can be done by clustering:
        if sort_rows is True and len(mut_mat) > 2:
            dist_mat = pdist(mut_mat)
            Z = linkage(dist_mat, 'ward', optimal_ordering=True)
            dn = dendrogram(Z, no_plot=True)
            sorted_anno = [anno_topN_nozero[i] for i in dn['leaves']]
            mut_mat_df_plot = mut_mat_df.loc[sorted_anno, :].copy()
        elif type(sort_rows) == list:
            sorted_anno = [sanno for sanno in sort_rows if sanno in anno_topN_nozero]
            mut_mat_df_plot = mut_mat_df.loc[sorted_anno, :].copy()
        else:
            sorted_anno = anno_topN_nozero
            mut_mat_df_plot = mut_mat_df.loc[sorted_anno, :].copy()

        # Adjust the plot size and generate it:
        y_scaler = 15/60
        y_len = y_scaler*len(mut_mat_df_plot)
        if y_len < 4:
            y_len = 4
        fig, ax = plt.subplots(1, 1, figsize=(25, y_len))
        g1 = sns.heatmap(mut_mat_df_plot, yticklabels=True, xticklabels=True, vmin=0, ax=ax);
        ax.yaxis.set_tick_params(labelsize=11)
        ax.xaxis.set_tick_params(labelsize=11)
        ax.set_xlabel("5' to 3' (filled in from right to left)", size=16)
        
        # Write to PDF file:
        plot_fnam = '{}/{}.pdf'.format(self.TM_dir_abs, plot_name)
        fig.tight_layout()
        fig.savefig(plot_fnam, bbox_inches='tight')
        if not png_dpi is False:
            fig.savefig(plot_fnam[:-4]+'.png', bbox_inches='tight', dpi=png_dpi)

        if no_plot_return:
            plt.close(fig)
            return((mut_mat_df, None, sorted_anno))
        else:
            return((mut_mat_df, fig, sorted_anno))

    def plot_transcript_mut_compare(self, species='human', \
                                    plot_name='tr-mut_matrix_comp', \
                                    no_plot_return=False, \
                                    mito=False, data_type='mut', \
                                    min_count_show=200, \
                                    sample_pairs=None, sample_unique_pairs=None, \
                                    sample_pairs_col='sample_name', \
                                    tr_compare_list_inp=None, \
                                    anno_substring_compare=None, \
                                    sample_list_exl=None, bc_list_exl=None,
                                    freq_avg_weighted=True, \
                                    topN=10, topN_select='max_diff'):
        '''
        Plot a heatmap of the pairwise comparison between the
        transcript mutation, gap or RT stop data of two samples.
        Will merge data from all samples in the provided sample list,
        except those excluded using the sample or barcode exclusion lists.

        Keyword arguments:
        sample_pairs -- Sample pairs to compare using the sample_df column defined by the sample_pairs_col input variable to merge replicates (default None)
        sample_pairs_col -- (default 'sample_name')
        sample_unique_pairs -- Sample pairs to compare using the 'sample_name_unique' column (default None)
        sample_list_exl -- List of unique sample names, using the sample_name_unique column, to exclude from any comparison (default None)
        bc_list_exl -- List of barcodes to exclude from any comparison (default None)
        anno_substring_compare -- Substring that must be in the annotation of the transcripts to compare (default None)
        tr_compare_list_inp -- List of transcript names to compare. If None, using the topN parameter (default None)
        topN -- Number of transcripts to plot. Taken from a list ordered by the difference between the sample pairs (default 10)
        topN_select -- The metric to calculate the difference between samples. Choose between 'max_diff', 'mean_diff' or 'mean_square' (default 'max_diff')
        species -- Species to plot (default 'human')
        plot_name -- Plot filename (default 'tr-mut_matrix_comp')
        no_plot_return -- Do not return plot (default False)
        mito -- Only plot mitochondrial transcripts (default False)
        data_type -- Data type to export. Choose between 'mut', 'gap' or 'RTstops' (default 'mut')
        min_count_show -- Minimum observations of a position to be shown (default 200)
        freq_avg_weighted -- Calculate the mutation frequency averaged across samples, weighted by the number observations (default True)
        '''

        if data_type not in ['mut', 'gap', 'RTstops']:
            print('The data_type input variable needs to be either \'mut\', \'gap\' or \'RTstops\'. Found {}'.format(data_type))

        # Handle if input is unique sample names
        # or sample names with replicates to merge:
        if not sample_pairs is None:
            assert(len(sample_pairs[0]) == len(sample_pairs[1]))
            # Exclude samples/barcodes:
            mask_exl = np.array([True]*len(self.sample_df))
            if not sample_list_exl is None:
                for snam in sample_list_exl:
                    mask_exl &= (self.sample_df['sample_name_unique'] != snam)
            if not bc_list_exl is None:
                for bc in bc_list_exl:
                    mask_exl &= (self.sample_df['barcode'] != bc)

            # Convert sample names to lists of unique sample names:
            pair_list = [[], []]
            name_list = [[], []]
            for sp1, sp2 in zip(*sample_pairs):
                mask1 = mask_exl & (self.sample_df[sample_pairs_col] == sp1)
                mask2 = mask_exl & (self.sample_df[sample_pairs_col] == sp2)
                if mask1.sum() > 0 and mask2.sum() > 0:
                    pair_list[0].append(list(self.sample_df.loc[mask1, 'sample_name_unique'].values))
                    pair_list[1].append(list(self.sample_df.loc[mask2, 'sample_name_unique'].values))
                    name_list[0].append(sp1)
                    name_list[1].append(sp2)
        elif not sample_unique_pairs is None:
            assert(len(sample_unique_pairs[0]) == len(sample_unique_pairs[1]))
            pair_list = [[], []]
            name_list = [[], []]
            for sp1, sp2 in zip(*sample_unique_pairs):
                pair_list[0].append([sp1])
                pair_list[1].append([sp2])
                name_list[0].append(sp1)
                name_list[1].append(sp2)
        else:
            raise Exception('Neither sample_pairs nor sample_unique_pairs lists were provided.')

        # Print each plot to the same PDF file:
        plot_fnam = '{}/{}.pdf'.format(self.TM_dir_abs, plot_name)
        with PdfPages(plot_fnam) as pp:
            # Process each pair:
            for name_idx, sp1l_sp2l in enumerate(zip(*pair_list)):
                sp1l, sp2l = sp1l_sp2l
                tr_muts_combi_s1 = self._combine_tr_muts(sp1l, freq_avg_weighted=freq_avg_weighted)
                tr_muts_combi_s2 = self._combine_tr_muts(sp2l, freq_avg_weighted=freq_avg_weighted)

                # Find the transcripts to compare:
                if not anno_substring_compare is None:
                    tr_compare_list = [anno for anno in tr_muts_combi_s1[species] if anno_substring_compare in anno]
                    # Sort list of annotations:
                    tr_compare_list = self._sort_freq_diff(tr_muts_combi_s1, \
                                                           tr_muts_combi_s2, \
                                                           species, min_count_show, \
                                                           data_type, mito, topN_select, \
                                                           list_to_sort=tr_compare_list)
                elif tr_compare_list_inp is None:
                    # Sort annotations according to largest distance between samples:
                    topN_anno = self._sort_freq_diff(tr_muts_combi_s1, tr_muts_combi_s2, \
                                                     species, min_count_show, data_type, \
                                                     mito, topN_select)
                    tr_compare_list = topN_anno[:topN]
                else:
                    tr_compare_list = tr_compare_list_inp

                # Get data for plotting.
                # The mutation matrix is going to contain the
                # transcript mutations for each sample,
                # one after the other, separated by an empty row:
                mut_mat = np.zeros((2, self.longest_tRNA))
                index_names = list()
                anno_idx = 0
                for anno in tr_compare_list:
                    if not anno in tr_muts_combi_s1[species] or not anno in tr_muts_combi_s2[species]:
                        continue
                    freq_mut_s1, min_count_mask_s1 = self._get_mut_freq_filted(tr_muts_combi_s1, species, anno, min_count_show, data_type)
                    freq_mut_s2, min_count_mask_s2 = self._get_mut_freq_filted(tr_muts_combi_s2, species, anno, min_count_show, data_type)
                    if freq_mut_s1 is None or freq_mut_s2 is None:
                        continue
                    # Positional coverage has to be fulfilled in both samples:
                    min_count_mask_s12 = min_count_mask_s1 & min_count_mask_s2
                    freq_mut_s1[~min_count_mask_s12] = 0
                    freq_mut_s2[~min_count_mask_s12] = 0

                    # Add an extra transcript to the mutation matrix:
                    if anno_idx > 0:
                        # The row that is spacing between transcript comparisons,
                        # should be filled with nan to map to white in the cmap:
                        space_row = np.ones((1, self.longest_tRNA)) * np.nan
                        new_rows = np.zeros((2, self.longest_tRNA))
                        mut_mat = np.vstack((mut_mat, space_row, new_rows))
                        index_names.append('')

                    # Add the transcript mutation frequencies from each
                    # sample, right aligned:
                    mut_mat[-2, -len(freq_mut_s1):] = freq_mut_s1
                    mut_mat[-1, -len(freq_mut_s2):] = freq_mut_s2
                    anno_short = anno.split('_')[-1]
                    index_names.append(name_list[0][name_idx]+'_'+anno_short)
                    index_names.append(name_list[1][name_idx]+'_'+anno_short)
                    anno_idx += 1
                # Skip if nothing was added:
                if anno_idx == 0:
                    continue

                # Transform mutation matrix to dataframe:
                seq_idx_str = list(map(str, range(mut_mat.shape[1])))
                mut_mat_df = pd.DataFrame(mut_mat, columns=seq_idx_str, index=index_names)

                # Adjust the plot size and generate it:
                y_scaler = 15/60
                y_len = y_scaler*len(mut_mat_df)
                if y_len < 4:
                    y_len = 4
                fig, ax = plt.subplots(1, 1, figsize=(25, y_len))
                g1 = sns.heatmap(mut_mat_df, yticklabels=True, xticklabels=True, vmin=0, ax=ax);
                ax.yaxis.set_tick_params(labelsize=11)
                ax.xaxis.set_tick_params(labelsize=11)
                ax.set_xlabel("5' to 3' (filled in from right to left)", size=16)
                plot_title = 'Sample comparison {} vs. {}'.format(name_list[0][name_idx], name_list[1][name_idx])
                ax.set_title(plot_title, fontsize=15)

                # Write to PDF file:
                fig.tight_layout()
                pp.savefig(fig, bbox_inches='tight')
                if no_plot_return:
                    plt.close(fig)

    def _sort_freq_diff(self, tr_muts_combi_s1, tr_muts_combi_s2, species, \
                        min_count_show, data_type, mito, topN_select, \
                        list_to_sort=None):
        '''
        Sort tRNA annotations according to the maximal difference
        in mutation or gap frequency observed between two samples.
        '''
        # Default, sort all annotations:
        if list_to_sort is None:
            anno_loop_list = tr_muts_combi_s1[species].keys()
        else:
            anno_loop_list = list_to_sort

        topN_anno = list()
        for anno in anno_loop_list:
            if not anno in tr_muts_combi_s1[species] and not anno in tr_muts_combi_s2[species]:
                continue
            # If mito is specified, skip non-mito annotations:
            if mito and 'mito' not in anno:
                continue
            freq_mut_s1, min_count_mask_s1 = self._get_mut_freq_filted(tr_muts_combi_s1, species, anno, min_count_show, data_type)
            freq_mut_s2, min_count_mask_s2 = self._get_mut_freq_filted(tr_muts_combi_s2, species, anno, min_count_show, data_type)
            if freq_mut_s1 is None or freq_mut_s2 is None:
                continue
            # Positional coverage has to be fulfilled in both samples:
            min_count_mask_s12 = min_count_mask_s1 & min_count_mask_s2
            freq_mut_s1[~min_count_mask_s12] = 0
            freq_mut_s2[~min_count_mask_s12] = 0

            # Find the requested distance between samples:
            if topN_select == 'max_diff':
                metric = np.max(np.abs(freq_mut_s1 - freq_mut_s2))
            elif topN_select == 'mean_diff':
                metric = np.mean(np.abs(freq_mut_s1 - freq_mut_s2))
            elif topN_select == 'mean_square':
                metric = np.mean(np.square(freq_mut_s1 - freq_mut_s2))
            else:
                raise Exception('"topN_select" input not understood: {}. Choose from "max_diff", "mean_diff", "mean_square".'.format(topN_select))
            topN_anno.append((anno, metric))

        # Sort annotations according to largest distance between samples:
        topN_anno = [tup[0] for tup in sorted(topN_anno, key=lambda x: x[1], reverse=True)]
        return(topN_anno)

    def plot_transcript_mut_cluster(self, species='human', \
                                    plot_name='tr-mut_matrix_clust', \
                                    mito=False, data_type='mut', \
                                    min_count_show=200, \
                                    anno_substring_incl=None, \
                                    sample_list_incl=None, \
                                    sample_list_exl=None, bc_list_exl=None,
                                    dist_metric='euclidean', \
                                    linkage_method='ward', \
                                    plot_compact=False, \
                                    right_align=True, \
                                    vmax=None):
        '''
        Cluster mutation, gap or RT stop data from 
        a set of samples for each transcript annotation requested.
        Then plot the clustered data as a heatmap
        for each transcript annotation.

        Keyword arguments:
        sample_list_incl -- List of unique sample names, using the "sample_name_unique" column, to include in clustering (default None)
        sample_list_exl -- List of unique sample names, using the "sample_name_unique" column, to exclude from clustering (default None)
        bc_list_exl -- List of barcodes to exclude from clustering (default None)
        anno_substring_incl -- Substring that must be in the annotation of the transcripts to cluster (default None)
        species -- Species to plot (default 'human')
        plot_name -- Plot filename (default 'tr-mut_matrix_clust')
        mito -- Only plot mitochondrial transcripts (default False)
        data_type -- Data type to cluster. Choose between 'mut', 'gap' or 'RTstops' (default 'mut')
        min_count_show -- Minimum observations of a position to be shown (default 200)
        dist_metric -- The metric to calculate the distance between samples. Passed directly to scipy.spatial.distance.pdist (default 'euclidean')
        linkage_method -- The clustering method used. Passed directly to scipy.cluster.hierarchy.linkage (default 'ward')
        plot_compact -- Make rows on heatmap compact (default False)
        right_align -- Right align data (default True)
        vmax -- Max value of heatmap range. Passed directly to seaborn.heatmap (default None)
        '''

        if data_type not in ['mut', 'gap', 'RTstops']:
            print('The data_type input variable needs to be either \'mut\', \'gap\' or \'RTstops\'. Found {}'.format(data_type))

        # If sample inclusion list is not provided
        # include all samples, minus those in the
        # exclusion lists:
        if sample_list_incl is None:
            # Exclude samples/barcodes:
            mask_incl = np.array([True]*len(self.sample_df))
            if not sample_list_exl is None:
                for snam in sample_list_exl:
                    mask_incl &= (self.sample_df['sample_name_unique'] != snam)
            if not bc_list_exl is None:
                for bc in bc_list_exl:
                    mask_incl &= (self.sample_df['barcode'] != bc)
        else:
            mask_incl = np.array([False]*len(self.sample_df))
            sample_list_incl_cp = copy.deepcopy(sample_list_incl)
            for snam in sample_list_incl:
                mask_found = (self.sample_df['sample_name_unique'] == snam)
                mask_incl |= mask_found
                if mask_found.sum() > 0:
                    sample_list_incl_cp.pop(sample_list_incl_cp.index(snam))
            if len(sample_list_incl_cp) > 0:
                print('Warning: {} of the provided samples were not found. {}'.format(len(sample_list_incl_cp), sample_list_incl_cp))

        # All the samples the cluster:
        sample_list = list(self.sample_df.loc[mask_incl, 'sample_name_unique'].values)
        tr_muts_combi = self._combine_tr_muts(sample_list, freq_avg_weighted=False)
        # Sort according to observations:
        anno_sorted = self._sort_anno(tr_muts_combi, species, mito=mito)
        anno_sorted = [anno_sorted[i][0] for i in range(len(anno_sorted))]
        # Get list of annotations to plot:
        if type(anno_substring_incl) == str:
            anno_list = [anno for anno in anno_sorted if anno_substring_incl in anno]
        elif type(anno_substring_incl) == list:
            anno_list = list()
            for substr in anno_substring_incl:
                anno_set = set(anno_list)
                anno_list.extend([anno for anno in anno_sorted if substr in anno and anno not in anno_set])
        else:
            anno_list = anno_sorted
            if not anno_substring_incl is None:
                print('Did not understand "anno_substring_incl" provided. Provide either list or str. Continuing with all annotations.')

        # Print each plot to the same PDF file:
        plot_fnam = '{}/{}.pdf'.format(self.TM_dir_abs, plot_name)
        with PdfPages(plot_fnam) as pp:
            # Process each transcript annotation:
            for anno in anno_list:
                # Find the mutations and insert them into a matrix:
                mut_mat = np.zeros((len(sample_list), self.longest_tRNA))
                tr_muts_snam = self._combine_tr_muts([sample_list[0]], freq_avg_weighted=False)
                _, master_min_count_mask = self._get_mut_freq_filted(tr_muts_snam, species, anno, \
                                                                     min_count_show, data_type)
                for row_i, snam in enumerate(sample_list):
                    # Get mutations for each sample:
                    tr_muts_snam = self._combine_tr_muts([snam], freq_avg_weighted=False)
                    freq_mut, min_count_mask = self._get_mut_freq_filted(tr_muts_snam, species, anno, \
                                                                         min_count_show, data_type)
                    # Skip annotation if sample is filtered:
                    if freq_mut is None:
                        print('Skipping annotation: {}\nNot enough counts for sample: {}'.format(anno, snam))
                        break
                    # Add mutations for matrix:
                    if right_align:
                        mut_mat[row_i, -len(freq_mut):] = freq_mut
                    else:
                        mut_mat[row_i, :len(freq_mut)] = freq_mut
                    master_min_count_mask &= min_count_mask
                
                # Do not attempt to plot filtered data:
                if freq_mut is None:
                    continue
                if master_min_count_mask.sum() == 0:
                    print('Skipping annotation: {}\nNot enough counts.'.format(anno))
                    continue

                # Mark masked positions by "X":
                X_mat = np.empty_like(mut_mat, dtype=str)
                X_mat[:, :] = 'x'
                unmask = np.array(['' if m else 'x' for m in master_min_count_mask])
                # Streamline min_count_mask for all samples:
                for row_i, snam in enumerate(sample_list):
                    if right_align:
                        mut_mat[row_i, -len(freq_mut):] *= master_min_count_mask
                        X_mat[row_i, -len(freq_mut):] = unmask
                    else:
                        mut_mat[row_i, :len(freq_mut)] *= master_min_count_mask
                        X_mat[row_i, :len(freq_mut)] = unmask
            
                # Transform mutation matrix to dataframe:
                seq_idx_str = list(map(str, range(mut_mat.shape[1])))
                mut_mat_df = pd.DataFrame(mut_mat, columns=seq_idx_str, index=sample_list)
                X_mat_df = pd.DataFrame(X_mat, columns=seq_idx_str, index=sample_list)

                # Cluster samples:
                dist_mat = pdist(mut_mat, metric=dist_metric)
                Z = linkage(dist_mat, method=linkage_method, optimal_ordering=True)
                dn = dendrogram(Z, no_plot=True)
                sorted_samples = [sample_list[i] for i in dn['leaves']]
                mut_mat_df_plot = mut_mat_df.loc[sorted_samples, :].copy()

                # Adjust the plot size and generate it:
                if plot_compact:
                    y_scaler = 10/60
                    y_len = y_scaler*len(mut_mat_df_plot)
                    if y_len < 2.7:
                        y_len = 2.7
                else:
                    y_scaler = 15/60
                    y_len = y_scaler*len(mut_mat_df_plot)
                    if y_len < 4:
                        y_len = 4

                # Generate figure:
                fig, ax = plt.subplots(1, 1, figsize=(25, y_len))
                g1 = sns.heatmap(mut_mat_df_plot, yticklabels=True, xticklabels=True, \
                                 vmin=0, vmax=vmax, ax=ax, annot=X_mat_df, fmt='1');
                ax.yaxis.set_tick_params(labelsize=11)
                ax.xaxis.set_tick_params(labelsize=11)
                if right_align:
                    ax.set_xlabel("5' to 3' (filled in from right to left)", size=16)
                else:
                    ax.set_xlabel("5' to 3' (filled in from left to right)", size=16)
                plot_title = 'Clustering on {}.'.format(anno)
                ax.set_title(plot_title, fontsize=15)

                # Write to PDF file:
                fig.tight_layout()
                pp.savefig(fig, bbox_inches='tight')
                plt.close(fig)

    def plot_transcript_mut_pos(self, tr_pos=34, species='human', \
                                idx_start=1, \
                                plot_name='tr-mut_pos', \
                                mito=False, \
                                min_count_show=200, \
                                anno_substring_incl=None, \
                                sample_list_incl=None, \
                                sample_rep_col='sample_name', \
                                sample_list_exl=None, bc_list_exl=None, \
                                xlabel='', xlabel_rot=0):
        '''
        Make barplots showing the frequency/percentage data for
        mutation, gap and RT stop on a single position from 
        a set of samples for each transcript annotation requested.


        Keyword arguments:

        tr_pos -- Transcript position, counted from left to right, to plot (default 34)
        idx_start -- "tr_pos" indexing (default 1)
        sample_list_incl -- List of unique sample names, using the "sample_name_unique" column, to include (default None)
        sample_list_exl -- List of unique sample names, using the "sample_name_unique" column, to exclude (default None)
        bc_list_exl -- List of barcodes to exclude (default None)
        anno_substring_incl -- Substring that must be in the annotation of the transcripts to plot (default None)
        species -- Species to plot (default 'human')
        plot_name -- Plot filename (default 'tr-mut_pos')
        mito -- Only plot mitochondrial transcripts (default False)
        min_count_show -- Minimum observations of a position to be shown (default 200)
        sample_rep_col -- Column in "sample_df" to use for defining sample replicates (default 'sample_name')
        xlabel -- Label on x-axis (default '')
        xlabel_rot -- Rotation of the transcript annotation used as xtick label (default 0)
        '''

        # Map transcript position to normal zero indexing:
        idx_pos = (tr_pos - idx_start)

        data_return = dict()  # Pick up and return data
        # If sample inclusion list is not provided
        # include all samples, minus those in the
        # exclusion lists:
        if sample_list_incl is None:
            # Exclude samples/barcodes:
            mask_incl = np.array([True]*len(self.sample_df))
            if not sample_list_exl is None:
                for snam in sample_list_exl:
                    mask_incl &= (self.sample_df['sample_name_unique'] != snam)
            if not bc_list_exl is None:
                for bc in bc_list_exl:
                    mask_incl &= (self.sample_df['barcode'] != bc)
        else:
            mask_incl = np.array([False]*len(self.sample_df))
            sample_list_incl_cp = copy.deepcopy(sample_list_incl)
            for snam in sample_list_incl:
                mask_found = (self.sample_df['sample_name_unique'] == snam)
                mask_incl |= mask_found
                if mask_found.sum() > 0:
                    sample_list_incl_cp.pop(sample_list_incl_cp.index(snam))
            if len(sample_list_incl_cp) > 0:
                print('Warning: {} of the provided samples were not found. {}'.format(len(sample_list_incl_cp), sample_list_incl_cp))

        # All the samples the cluster:
        sample_list = list(self.sample_df.loc[mask_incl, 'sample_name_unique'].values)
        tr_muts_combi = self._combine_tr_muts(sample_list, freq_avg_weighted=False)
        # Sort according to observations:
        anno_sorted = self._sort_anno(tr_muts_combi, species, mito=mito)
        anno_sorted = [anno_sorted[i][0] for i in range(len(anno_sorted))]
        # Get list of annotations to plot:
        if type(anno_substring_incl) == str:
            anno_list = [anno for anno in anno_sorted if anno_substring_incl in anno]
        elif type(anno_substring_incl) == list:
            anno_list = list()
            for substr in anno_substring_incl:
                anno_set = set(anno_list)
                anno_list.extend([anno for anno in anno_sorted if substr in anno and anno not in anno_set])
        else:
            anno_list = anno_sorted
            if not anno_substring_incl is None:
                print('Did not understand "anno_substring_incl" provided. Provide either list or str. Continuing with all annotations.')

        sample_list_mask = self.sample_df['sample_name_unique'].isin(sample_list)
        Nplot_groups = len(set(self.sample_df[sample_list_mask][sample_rep_col]))
        # Print each plot to the same PDF file:
        plot_fnam = '{}/{}.pdf'.format(self.TM_dir_abs, plot_name)
        with PdfPages(plot_fnam) as pp:
            # Process each transcript annotation:
            for anno in anno_list:
                plot_width = Nplot_groups * 4.5
                fig = plt.figure(figsize=(plot_width, 5))
                gs = fig.add_gridspec(1, 3)
                ax1 = fig.add_subplot(gs[:, 0])
                ax2 = fig.add_subplot(gs[:, 1])
                ax3 = fig.add_subplot(gs[:, 2])
                print_plot = True

                # Plot all three types of data:
                for data_type, ax in zip(['mut', 'gap', 'RTstops'], [ax1, ax2, ax3]):
                    # Find the mutations and insert them into a matrix:
                    tr_muts_snam = self._combine_tr_muts([sample_list[0]], freq_avg_weighted=False)
                    _, master_min_count_mask = self._get_mut_freq_filted(tr_muts_snam, species, anno, \
                                                                         min_count_show, data_type)
                    # Skip annotation if sample is filtered:
                    if master_min_count_mask is None:
                        print('Skipping annotation: {}\nNot enough counts for sample: {}'.format(anno, sample_list[0]))
                        print_plot = False
                        break
                    mut_mat = np.zeros((len(sample_list), len(master_min_count_mask)))
                    # Collect data for each sample:
                    for row_i, snam in enumerate(sample_list):
                        # Get mutations for each sample:
                        tr_muts_snam = self._combine_tr_muts([snam], freq_avg_weighted=False)
                        freq_mut, min_count_mask = self._get_mut_freq_filted(tr_muts_snam, species, anno, \
                                                                             min_count_show, data_type)
                        # Skip annotation if sample is filtered:
                        if freq_mut is None:
                            print('Skipping annotation: {}\nNot enough counts for sample: {}'.format(anno, snam))
                            break
                        # Add mutations for matrix:
                        mut_mat[row_i, :] = freq_mut
                        master_min_count_mask &= min_count_mask
                    
                    # Do not attempt to plot filtered data:
                    if freq_mut is None:
                        print_plot = False
                        continue
                    if not master_min_count_mask[idx_pos]:
                        print('Skipping annotation: {}\nNot enough counts.'.format(anno))
                        print_plot = False
                        continue

                    # Transform mutation matrix to dataframe:
                    mut_pos_df = pd.DataFrame(zip(sample_list, mut_mat[:, idx_pos]), \
                                              columns=['sample_name_unique', 'mut_freq'])
                    mut_pos_df = mut_pos_df.merge(self.sample_df, on='sample_name_unique', suffixes=('', '_2'))
                    assert(sample_rep_col in mut_pos_df.columns)
                    try:
                        data_return[anno][data_type] = mut_pos_df.copy()
                    except KeyError:
                        data_return[anno] = {mt: {} for mt in ['mut', 'gap', 'RTstops']}
                        data_return[anno][data_type] = mut_pos_df.copy()
                    
                    # Plot:
                    x_order = natsorted(set(mut_pos_df[sample_rep_col]))
                    g1 = sns.barplot(ax=ax, data=mut_pos_df, x=sample_rep_col, y='mut_freq', \
                                     capsize=0.1, edgecolor=".2", linewidth=2, alpha=0.8, \
                                     order=x_order)
                    g2 = sns.swarmplot(ax=ax, data=mut_pos_df, x=sample_rep_col, y='mut_freq', \
                                       color='grey', alpha=0.7, edgecolor='black', dodge=True, \
                                       linewidth=0.8, size=6, marker='X', warn_thresh=1, order=x_order)

                    g1.set(xlabel=xlabel)
                    g1.set_xticklabels(g1.get_xticklabels(), rotation=xlabel_rot)
                    ymin, ymax = g1.get_ylim()
                    if ymin < 0:
                        ymin = 0
                    if data_type == 'mut':
                        g1.set_title('Mutation frequency')
                        g1.set(ylabel='Frequency')
                    elif data_type == 'gap':
                        g1.set_title('Gap frequency')
                        g1.set(ylabel='Frequency')
                    else:
                        g1.set_title('RT stop')
                        g1.set(ylabel='Fall off (%)')
                        ymax = 100
                    g1.set(ylim=(ymin, ymax))

                plot_title = 'Mutation data for {}, position {}.'.format(anno, tr_pos)
                fig.suptitle(plot_title, fontsize=15)

                if print_plot:
                    # Write to PDF file:
                    fig.tight_layout()
                    pp.savefig(fig, bbox_inches='tight')
                plt.close(fig)

        return(data_return)

    def _get_mut_freq_filted(self, tr_muts_combi, species, anno, \
                             min_count_show, data_type):
        '''
        Filter the mutation or gap frequency array
        by a minimum number of observations.
        '''
        if data_type == 'gap':
            freq_mut = tr_muts_combi[species][anno]['gap_freq']
        elif data_type == 'mut':
            freq_mut = tr_muts_combi[species][anno]['mut_freq']
        else:
            freq_mut = tr_muts_combi[species][anno]['RTstops']
        # Enforce a minimum number of observations:
        counts_all = tr_muts_combi[species][anno]['PSCM'].sum(1)
        min_count_mask = counts_all >= min_count_show
        if min_count_mask.sum() == 0:
            return(None, None)
        freq_mut[~min_count_mask] = 0
        return(freq_mut, min_count_mask)

    def _sort_anno(self, tr_muts_combi, species, mito=False):
        '''
        Sort annotations (and species) according to the number
        of observations at the 3p end.
        '''
        anno2obs = dict()
        for anno in tr_muts_combi[species]:
            # If mito is specified, skip non-mito annotations:
            if mito and 'mito' not in anno:
                continue
            tRNA_len = tr_muts_combi[species][anno]['seq_len']
            if tr_muts_combi[species][anno]['PSCM'].max().max() > 0:
                # The all count in the second to last row i.e. C in CCA:
                obs = tr_muts_combi[species][anno]['PSCM'].sum(1).values[-2]
                anno2obs[anno] = obs
        anno_sorted = sorted(anno2obs.items(), key=lambda x: x[1], reverse=True)
        return(anno_sorted)

    def _calc_mut_freq(self, tr_muts_combi, anno, species, \
                       gap_only, min_count_show):
        '''
        Calculate the frequency of all mismatches (including gaps)
        or the frequency of gaps alone. Also allow for masking
        based on a minimum number of observations.
        '''
        tr_len = tr_muts_combi[species][anno]['seq_len']
        tr_seq = tr_muts_combi[species][anno]['seq']

        # Calculate the mutation (or gap) frequency per sequence position #
        # First count all observations for a given positions:
        counts_all = tr_muts_combi[species][anno]['PSCM'].sum(1)
        # Enforce a minimum number of observations:
        min_count_mask = counts_all >= min_count_show
        if min_count_mask.sum() == 0:
            return(None)
        counts_all.loc[~min_count_mask] = 0
        
        if gap_only:
            # Make a mask to count all the gaps:
            count_mask = np.zeros((tr_len, len(self.char_list)))
            tr_char_idx = [self.char_dict['-'] for _ in range(tr_len)]
            for seq_i, char_i in enumerate(tr_char_idx):
                count_mask[seq_i, char_i] = 1
            counts_gap = (tr_muts_combi[species][anno]['PSCM'] * count_mask).sum(1)
            # np.divide setting 0/0 to 0 i.e. no observations, no gap
            freq_mut = np.divide(counts_gap, counts_all, out=np.zeros_like(counts_gap), where=counts_all!=0)
        else:
            # Make a mask to count all the non-mutant characters:
            count_mask = np.zeros((tr_len, len(self.char_list)))
            tr_char_idx = [self.char_dict[nt] for nt in tr_seq]
            for seq_i, char_i in enumerate(tr_char_idx):
                count_mask[seq_i, char_i] = 1
            counts_cor = (tr_muts_combi[species][anno]['PSCM'] * count_mask).sum(1)
            # np.divide setting 0/0 to 1, i.e. no observation, no mutation
            freq_mut = 1 - np.divide(counts_cor, counts_all, out=np.ones_like(counts_cor), where=counts_all!=0)
        return(freq_mut.values)

    def _calc_RTstops(self, tr_muts_combi, anno, species):
        '''
        Calculate the RT stops, also referred to as "termination signal"
        in Wang et al. 2021.
        Calculated as: 100 * (cov(i+1) - cov(i)) / cov(i+1)
        '''
        # Count all observations for a given positions:
        counts_all = tr_muts_combi[species][anno]['PSCM'].sum(1).values
        counts_all_p1 = np.roll(counts_all, -1)
        counts_all_p1[-1] = counts_all[-1]
        # Calculate RTstops (similar to Wang et al. 2021):
        RTstops_arr = 100 * np.divide((counts_all_p1 - counts_all), counts_all_p1, \
                                      out=np.zeros_like(counts_all, dtype=np.float64), \
                                      where=counts_all_p1!=0)
        return(RTstops_arr)

    def mask_tRNA_database(self, min_mut_freq=0.84, min_pos_count=1000, min_tr_count=3000, \
                           frac_max_score=0.95, match_score=1, mismatch_score=-1, \
                           open_gap_score=-2, extend_gap_score=-2, sample_list=None,
                           freq_avg_weighted=True):
        '''
        Mask tRNA sequences based on the mutation frequency
        observed in the input samples.
        In order to extend the masking to transcripts with few
        or no observations, all the transcripts are aligned
        and those highly similar (defined by the frac_max_score input variable)
        inherit the maskings from transcripts with sufficient observations.

        Keyword arguments:
        min_mut_freq -- Minimum mutation frequency to be masked (default 0.84)
        min_pos_count -- Minimum observations of a position to be masked (default 1000)
        min_tr_count -- Minimum observations of a transcript to be masked (default 3000)
        frac_max_score -- The minimum fraction of the maximal alignment score to share masking between two transcripts (default 0.95)
        match_score -- Match score in alignment between two transcript (default 1)
        mismatch_score -- Mismatch penalty in alignment between two transcript (default -1)
        open_gap_score -- Open gap penalty in alignment between two transcript (default -2)
        extend_gap_score -- Extend gap penalty in alignment between two transcript (default -2)
        sample_list -- List of samples to combine for sequence masking. If None, all samples with mutation data are used (default None)
        freq_avg_weighted -- Calculate the mutation frequency averaged across samples, weighted by the number observations (default True)
        '''
        # Get the mutations combined for the requested samples:
        if sample_list is None:
            sample_list = list(self.tr_muts.keys())
        tr_muts_combi = self._combine_tr_muts(sample_list, freq_avg_weighted=freq_avg_weighted)

        # Find and store masked sequence #
        for species in tr_muts_combi:
            for anno in tr_muts_combi[species]:
                # Get the mutation frequency per sequence position:
                freq_mut = tr_muts_combi[species][anno]['mut_freq']
                # Enforce a minimum number of observations:
                counts_all = tr_muts_combi[species][anno]['PSCM'].sum(1)
                min_count_mask = counts_all >= min_pos_count
                if min_count_mask.sum() == 0 or freq_mut is None or counts_all.values[-2] < min_tr_count:
                    continue
                freq_mut[~min_count_mask] = 0

                # Find all positions that should be masked:
                mask_idx = set((freq_mut >= min_mut_freq).nonzero()[0])
                if len(mask_idx) == 0:
                    continue
                # Make masked sequence:
                tr_seq = tr_muts_combi[species][anno]['seq']
                seq_masked = ['N' if seq_i in mask_idx else char for seq_i, char in enumerate(tr_seq)]
                tr_muts_combi[species][anno]['seq_masked'] = ''.join(seq_masked)

        # Expand masked sequence based on similar sequences #
        # Initiate the aligner:
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = match_score
        aligner.mismatch_score = mismatch_score
        aligner.open_gap_score = open_gap_score
        aligner.extend_gap_score = extend_gap_score

        for species in tr_muts_combi:
            for anno_acpt in tr_muts_combi[species]: # acceptor
                seq_acpt = tr_muts_combi[species][anno_acpt]['seq']
                seq_acpt_cov = tr_muts_combi[species][anno_acpt]['PSCM'].sum(1)
                anticodon_acpt = anno_acpt.split('-')[2] # anticodon
                tr_muts_combi[species][anno_acpt]['seq_masked_exp'] = ''

                # Add masked characters from highly similar sequences:
                for anno_dnr in tr_muts_combi[species]: # donor
                    # Skip sequences with different anticodon:
                    if anticodon_acpt != anno_dnr.split('-')[2]:
                        continue
                    # Skip sequences with no masked characters to contribute:
                    elif not 'seq_masked' in tr_muts_combi[species][anno_dnr]:
                        continue

                    seq_dnr = tr_muts_combi[species][anno_dnr]['seq']
                    max_score = match_score * max([len(seq_acpt), len(seq_dnr)])
                    min_score = max_score * frac_max_score

                    # Perform alignment and take out one:
                    alignments = aligner.align(seq_acpt, seq_dnr)
                    alignment = alignments[0]
                    # Acceptor = target (t)
                    # Donor = query (q)
                    t_cor, q_cor = alignment.aligned

                    # Following enforces end-to-end coverage of the alignment
                    # on both seq_acpt and seq_dnr:
                    if t_cor[0][0] != 0 or q_cor[0][0] != 0 or t_cor[-1][1] != len(seq_acpt) or q_cor[-1][1] != len(seq_dnr):
                        continue
                    # Also enforce minimum alignment score:
                    elif alignment.score < min_score:
                        continue

                    # Concatenate the masked seq_dnr that aligns to the seq_acpt:
                    q_masked = tr_muts_combi[species][anno_dnr]['seq_masked']
                    q_masked_trans = list()
                    q_trans = list()
                    for i in range(0, len(q_cor)):
                        q_masked_trans.append(q_masked[q_cor[i][0]:q_cor[i][1]])
                        q_trans.append(seq_dnr[q_cor[i][0]:q_cor[i][1]])
                        # Insert gaps when these appear in the seq_dnr
                        # i.e. characters in the seq_acpt is unaligned:
                        if i < (len(q_cor)-1):
                            for j in range(t_cor[i+1][0] - t_cor[i][1]):
                                q_masked_trans.append('-')
                                q_trans.append('-')
                    q_masked_trans = ''.join(q_masked_trans)
                    q_trans = ''.join(q_trans)
                    assert(len(q_masked_trans) == len(seq_acpt))
                    assert(len(q_trans) == len(seq_acpt))

                    # Update masked sequence:
                    if anno_acpt == anno_dnr:
                        # Move the _exp regardless of positional coverage:
                        tr_muts_combi[species][anno_acpt]['seq_masked_exp'] = q_masked_trans
                    else:
                        if tr_muts_combi[species][anno_acpt]['seq_masked_exp'] == '':
                            t_masked = seq_acpt
                        else:
                            t_masked = tr_muts_combi[species][anno_acpt]['seq_masked_exp']

                        # If "N" is found in the masked seq_dnr then pick it out
                        # otherwise stick to acceptor sequence (seq_acpt):
                        # Also, require the acceptor of masking to be the same nucleotide
                        # as the donor (i.e. seq_dnr) [tc == qc]
                        # Finally, require low coverage to avoid forcing a mask
                        # when there is coverage to suggest otherwise.
                        seq_masked_exp_lst = list()
                        for cov, tc, qc, qm in zip(seq_acpt_cov, t_masked, q_trans, q_masked_trans):
                            if cov < min_pos_count and tc == qc and qm == 'N':
                                seq_masked_exp_lst.append(qm)
                            else:
                                seq_masked_exp_lst.append(tc)
                        seq_masked_exp = ''.join(seq_masked_exp_lst)
                        tr_muts_combi[species][anno_acpt]['seq_masked_exp'] = seq_masked_exp

        # Collect statistics on the masked sequences:
        self.mask_stats = {sp: {} for sp in tr_muts_combi}
        for species in tr_muts_combi:
            anno_list = list(tr_muts_combi[species].keys())
            mask_count = np.array([tr_muts_combi[species][anno]['seq_masked_exp'].count('N') for anno in anno_list])
            
            self.mask_stats[species]['mask_count'] = mask_count
            self.mask_stats[species]['mask_sum'] = sum(mask_count)
            self.mask_stats[species]['mask_mean'] = np.mean(mask_count)
            self.mask_stats[species]['mask_frac'] = (mask_count > 0).sum() / len(mask_count)
        
        # Add the masked sequences as an object attribute:
        self.tr_muts_masked = tr_muts_combi

    def write_masked_tRNA_database(self, out_dir='tRNA_database_masked'):
        '''
        Write the masked tRNA database with folder structure,
        sequences and BLAST index such that it can be used
        directly for a new alignment run.
        '''
        if self.tr_muts_masked is None:
            print('No masked sequences were found. Run the "mask_tRNA_database" function to generate masked sequences.')
            return(1)

        # Make out_dir folder:
        out_dir_abs = '{}/{}'.format(self.TM_dir_abs, out_dir)
        try:
            os.mkdir(out_dir_abs)
        except:
            shutil.rmtree(out_dir_abs)
            os.mkdir(out_dir_abs)

        tRNA_database_masked = dict()
        for species in self.tr_muts_masked:
            # Make out_dir species folder:
            out_dir_sp_abs = '{}/{}'.format(out_dir_abs, species)
            try:
                os.mkdir(out_dir_sp_abs)
            except:
                shutil.rmtree(out_dir_sp_abs)
                os.mkdir(out_dir_sp_abs)
            
            # Write the masked sequences as fasta:
            out_fnam = '{}-tRNAs.fa'.format(species)
            out_fnam_abs = '{}/{}'.format(out_dir_sp_abs, out_fnam)
            tRNA_database_masked[species] = out_fnam_abs
            with open(out_fnam_abs, 'w') as fh:
                for anno in self.tr_muts_masked[species]:
                    if self.tr_muts_masked[species][anno]['seq_masked_exp'] == '':
                        masked_seq = self.tr_muts_masked[species][anno]['seq']
                    else:
                        masked_seq = self.tr_muts_masked[species][anno]['seq_masked_exp']
                    fh.write('>{}\n{}\n'.format(anno, masked_seq))
            
            # Make BLAST DB on fasta:
            self._makeblastdb(out_dir_sp_abs, out_fnam)

        return(tRNA_database_masked)
    
    def _makeblastdb(self, file_dir, fnam):
        '''
        Build BLAST database (required for SWIPE).
        Notice, it is built with sequence type = protein
        This is to enable a custom scoring matrix for SWIPE.
        '''
        try:
            os.chdir(file_dir)
            cmd_str = 'makeblastdb -dbtype prot -in {} -blastdb_version 4'.format(fnam)
            cmd = cmd_str.split()
            log_fn = 'makeblastdb_logfile.txt'
            with Popen(cmd, stdout=PIPE, stderr=STDOUT) as p, open(log_fn, 'a') as file:
                file.write('Starting subprocess with command:')
                file.write(str(cmd))
                file.write('\n')
                for line in p.stdout:
                    file.write(line.decode('utf-8'))
                file.write('\n****** DONE ******\n\n\n')
            
            os.chdir(self.dir_dict['NBdir'])
        except Exception as err:
            os.chdir(self.dir_dict['NBdir'])
            raise err



