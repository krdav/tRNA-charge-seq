import sys, os, shutil, bz2, copy, contextlib
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
    This class is used generate statistics over the observed
    transcript mutations i.e. the mismatches and gaps in the
    alignment between a read and its tRNA transcript.
    '''
    def __init__(self, dir_dict, sample_df, tRNA_database, pull_default=False, common_seqs=None, ignore_common_count=False):
        # Input:
        self.sample_df, self.tRNA_database = sample_df, tRNA_database
        self.dir_dict = dir_dict
        self.common_seqs_fnam = common_seqs
        self.common_seqs_dict = dict() # map common sequences to index
        self.char_str = 'ACGTUN-'
        self.char_list = [c for c in self.char_str]
        self.char_dict = {c: i for i, c in enumerate(self.char_str)}
        
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
        # Read the tRNA transcripts:
        for species in tRNA_database:
            self.tr_muts[species] = dict()
            for record in SeqIO.parse(tRNA_database[species], "fasta"):
                self.tr_muts[species][record.id] = dict()
                self.tr_muts[species][record.id]['seq'] = str(record.seq)
                self.tr_muts[species][record.id]['seq_len'] = len(record.seq)
                # Position specific count matrix:
                self.tr_muts[species][record.id]['PSCM'] = list()

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
            assert(os.path.exists(self.common_seqs_fnam))
            assert(self.common_seqs_fnam[-4:] == '.bz2')
            with bz2.open(self.common_seqs_fnam, "rt") as input_fh:
                for ridx, record in enumerate(SeqIO.parse(input_fh, "fasta")):
                    seq = str(record.seq)
                    assert(ridx == int(record.id))
                    assert(not seq in self.common_seqs_dict)
                    self.common_seqs_dict[seq] = ridx

    def make_dir(self, overwrite=True):
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
        return(self.TM_dir_abs)

    def find_muts(self, unique_anno=True, match_score=1, mismatch_score=-1, open_gap_score=-2, extend_gap_score=-1, n_jobs=4, verbose=True, sample_list=None):
        self.verbose = verbose
        if self.verbose:
            print('Collecting stats from:', end='')
        if sample_list is None:
            sample_list = [row['sample_name_unique'] for _, row in self.sample_df.iterrows()]
        
        self.unique_anno = unique_anno
        self.match_score, self.mismatch_score, self.open_gap_score, self.extend_gap_score = match_score, mismatch_score, open_gap_score, extend_gap_score
        
        # Find mutations in the transcripts for each file:
        data = [(idx, row) for idx, row in self.sample_df.iterrows() if row['sample_name_unique'] in sample_list]
        with WorkerPool(n_jobs=n_jobs) as pool:
            results = pool.map(self.__collect_transcript_muts, data)
        self.__fill_tr_muts(results)

    def __collect_transcript_muts(self, index, row):
        if self.verbose:
            print('  {}'.format(row['sample_name_unique']), end='')

        species = row['species']
        # Dictionary to store mutation info for each transcript:
        tr_muts_sp = copy.deepcopy(self.tr_muts)
        
        # Initiate the aligner:
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.match_score = self.match_score
        aligner.mismatch_score = self.mismatch_score
        aligner.open_gap_score = self.open_gap_score
        aligner.extend_gap_score = self.extend_gap_score

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

        # Read the stats file to get the old alignment annotations:
        stats_fnam = '{}/{}_stats.csv.bz2'.format(self.stats_dir_abs, row['sample_name_unique'])
        with bz2.open(stats_fnam, 'rt', encoding="utf-8") as stats_fh:
            sample_stats = pd.read_csv(stats_fh, keep_default_na=False)
        ID2idx = {readID: i for i, readID in enumerate(sample_stats['readID'])}

        # To keep each read having weight 1:
        weight_dict = {anno: [] for anno in tr_muts_sp[species]}
        char_list = list() # 
        for seq in dedup_seq_count:
            # Check if common sequence, then readID has changed:
            if not self.common_seqs_fnam is None and seq in self.common_seqs_dict:
                readID = str(self.common_seqs_dict[seq])
            else:
                readID = dedup_seq_count[seq]['id']

            # Get the index in the alignment statistics:
            if readID in ID2idx:
                idx = ID2idx[readID]
            else:
                # Skip unaligned reads:
                continue
            stats_row = sample_stats.loc[idx, :]
            total_count = dedup_seq_count[seq]['count'] * stats_row['count']
            if not stats_row['3p_cover']:
                # Skip reads that do not have 3p coverage:
                continue
            # Get list of annotations:
            anno_list = stats_row['tRNA_annotation'].split('@')
            # Skip if multiple annotations found but unique requested:
            if self.unique_anno and len(anno_list) > 1:
                continue

            for anno in anno_list:
                # Generate alignments:
                target = tr_muts_sp[species][anno]['seq']
                alignments = aligner.align(target, seq)
                # If multiple alignments with the same score these should be weighted
                # so one read contributes with one observation:
                weight = 1.0 / len(alignments) * total_count
                for alignment in alignments:
                    weight_dict[anno].append(weight)
                    # Initiate the character array:
                    char_ar = np.empty(len(target), dtype='<U1')
                    # Extract the alignment coordinates:
                    t_cor, q_cor = alignment.aligned
                    if (t_cor[-1, -1] + 1) < tr_muts_sp[species][anno]['seq_len']:
                        raise Exception('Sequence should be 3p aligned but appears not to be...')

                    # Find gaps:
                    for i in range(1, len(t_cor)):
                        for j in range(t_cor[i][0] - t_cor[i-1][1]):
                            char_ar[t_cor[i-1][1] + j] = '-'

                    # Find mismatches/mutations:
                    mut_idx = list()
                    for tran, qran in zip(t_cor, q_cor):
                        for ti, qi in zip(range(*tran), range(*qran)):
                            char_ar[ti] = seq[qi]
                    tr_muts_sp[species][anno]['PSCM'].append(char_ar)

        # Convert the character arrays to observations:
        tr_muts_sp[species] = self.__count_char_matrix(tr_muts_sp[species], weight_dict)
        return(tr_muts_sp)

    def __count_char_matrix(self, tr_muts_sp, weight_dict):
        for anno in tr_muts_sp:
            # Skip if no observations:
            if len(tr_muts_sp[anno]['PSCM']) == 0:
                continue
            char_matrix = np.array(tr_muts_sp[anno]['PSCM'])
            # Count matrix (len tRNA, number char):
            char_count = np.zeros((char_matrix.shape[1], len(self.char_list)))
            # For each position in the transcript:
            for pos in range(char_matrix.shape[1]):
                # Sum the character observations:
                for char_i, char in enumerate(self.char_list):
                    char_count[pos, char_i] = sum((char_matrix[:, pos] == char)*weight_dict[anno])

            # Turn to dataframe:
            count_df = pd.DataFrame(char_count, columns=self.char_list)
            # col_mask = count_df.sum() != 0
            # count_df = count_df.loc[:, col_mask].copy()
            tr_muts_sp[anno]['PSCM'] = count_df
        return(tr_muts_sp)
    
    def __fill_tr_muts(self, results):
        for res in results:
            species = list(res.keys())[0]
            for anno in res[species]:
                # Skip if no observations:
                if len(res[species][anno]['PSCM']) == 0:
                    continue
                # Fill the dictionary for all samples:
                if len(self.tr_muts[species][anno]['PSCM']) == 0:
                    self.tr_muts[species][anno]['PSCM'] = res[species][anno]['PSCM']
                else:
                    self.tr_muts[species][anno]['PSCM'] += res[species][anno]['PSCM']

    def fix_end(self):
        for species in self.tr_muts:
            for anno in self.tr_muts[species]:
                PSCM_len = len(self.tr_muts[species][anno]['PSCM'])
                if PSCM_len > 0:
                    end_obs = self.tr_muts[species][anno]['PSCM']['C'].values[-2]
                    end_ar = np.zeros(len(self.char_list))
                    A_idx = self.char_list.index('A')
                    end_ar[A_idx] = end_obs
                    self.tr_muts[species][anno]['PSCM'].loc[PSCM_len-1, :] = end_ar

    def plot_transcript_logo(self, topN=30, species='human', plot_name='tr-mut_logos'):
        # Sort according to observations:
        anno2obs = dict()
        for anno in self.tr_muts[species]:
            if len(self.tr_muts[species][anno]['PSCM']) > 0:
                # The all count in the second to last row i.e. C in CCA:
                obs = self.tr_muts[species][anno]['PSCM'].sum(1).values[-2]
                anno2obs[anno] = obs
        anno_sorted = sorted(anno2obs.items(), key=lambda x: x[1], reverse=True)
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
                    logo_plot = lm.Logo(self.tr_muts[species][anno]['PSCM'], color_scheme='classic');
                logo_plot.ax.set_title(anno, fontsize=15)
                logo_plot.ax.set_xlabel("5' to 3'")
                logo_plot.ax.set_ylabel("Count");
                logo_plot.fig.tight_layout()
                pp.savefig(logo_plot.fig, bbox_inches='tight')
                plt.close(logo_plot.fig)

    def plot_transcript_cov(self, topN=50, species='human', plot_name='tr-cov_matrix', png_dpi=False, no_plot_return=False, mito=False, sort_rows=True):
        # Sort according to observations:
        anno2obs = dict()
        longest_tRNA = 0
        for anno in self.tr_muts[species]:
            # If mito is specified, skip non-mito annotations:
            if mito and 'mito' not in anno:
                continue
            tRNA_len = len(self.tr_muts[species][anno]['PSCM'])
            if tRNA_len > 0:
                # The all count in the second to last row i.e. C in CCA:
                obs = self.tr_muts[species][anno]['PSCM'].sum(1).values[-2]
                anno2obs[anno] = obs
                if tRNA_len > longest_tRNA:
                    longest_tRNA = tRNA_len
        anno_sorted = sorted(anno2obs.items(), key=lambda x: x[1], reverse=True)
        if topN > len(anno_sorted):
            topN = len(anno_sorted)
        
        # Find the observations and insert them into a matrix:
        obs_mat = np.zeros((topN, longest_tRNA))
        anno_topN = [anno_sorted[i][0] for i in range(topN)]
        for i in range(topN):
            anno = anno_topN[i]
            # Count all observations for a given positions:
            counts_all = self.tr_muts[species][anno]['PSCM'].sum(1)
            # Insert from right side so CCA is always indexed
            # at the highest index:
            obs_mat[i, -len(counts_all):] = counts_all

        # The read count for each transcript is weighted so they all start of with
        # the same coverage at the 3p and normalized so it sums to 100.
        # Weigh each transcript:
        obs_mat = obs_mat.T / obs_mat[:, -1]
        # Normalize so coverage at 3p is 100:
        obs_mat = 100 * obs_mat.T

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

    def plot_transcript_mut(self, topN=50, species='human', plot_name='tr-mut_matrix', png_dpi=False, no_plot_return=False, mito=False, gap_only=False, sort_rows=True, min_count_show=10):
        # Sort according to observations:
        anno2obs = dict()
        longest_tRNA = 0
        for anno in self.tr_muts[species]:
            # If mito is specified, skip non-mito annotations:
            if mito and 'mito' not in anno:
                continue
            tRNA_len = len(self.tr_muts[species][anno]['PSCM'])
            if tRNA_len > 0:
                # The all count in the second to last row i.e. C in CCA:
                obs = self.tr_muts[species][anno]['PSCM'].sum(1).values[-2]
                anno2obs[anno] = obs
                if tRNA_len > longest_tRNA:
                    longest_tRNA = tRNA_len
        anno_sorted = sorted(anno2obs.items(), key=lambda x: x[1], reverse=True)
        if topN > len(anno_sorted):
            topN = len(anno_sorted)
        # Find the mutations and insert them into a matrix:
        mut_mat = np.zeros((topN, longest_tRNA))
        anno_topN = [anno_sorted[i][0] for i in range(topN)]
        anno_skip = set() # Skip these because of lacking observations
        for i in range(topN):
            anno = anno_topN[i]
            tr_len = len(self.tr_muts[species][anno]['PSCM'])
            tr_seq = self.tr_muts[species][anno]['seq']

            # Calculate the mutation (or gap) frequency per sequence position #
            # First count all observations for a given positions:
            counts_all = self.tr_muts[species][anno]['PSCM'].sum(1)
            # Enforce a minimum number of observations:
            min_count_mask = counts_all >= min_count_show
            if min_count_mask.sum() == 0:
                anno_skip.add(anno)
                continue
            counts_all.loc[~min_count_mask] = 0
            
            if gap_only:
                # Make a mask to count all the gaps:
                count_mask = np.zeros((tr_len, len(self.char_list)))
                tr_char_idx = [self.char_dict['-'] for _ in range(tr_len)]
                for seq_i, char_i in enumerate(tr_char_idx):
                    count_mask[seq_i, char_i] = 1
                counts_gap = (self.tr_muts[species][anno]['PSCM'] * count_mask).sum(1)
                # np.divide setting 0/0 to 0 i.e. no observations, no gap
                freq_mut = np.divide(counts_gap, counts_all, out=np.zeros_like(counts_gap), where=counts_all!=0)
            else:
                # Make a mask to count all the non-mutant characters:
                count_mask = np.zeros((tr_len, len(self.char_list)))
                tr_char_idx = [self.char_dict[nt] for nt in tr_seq]
                for seq_i, char_i in enumerate(tr_char_idx):
                    count_mask[seq_i, char_i] = 1
                counts_cor = (self.tr_muts[species][anno]['PSCM'] * count_mask).sum(1)
                # np.divide setting 0/0 to 1, i.e. no observation, no mutation
                freq_mut = 1 - np.divide(counts_cor, counts_all, out=np.ones_like(counts_cor), where=counts_all!=0)

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
            mut_mat_df_plot = mut_mat_df.loc[sorted_anno, :]
        elif type(sort_rows) == list:
            sorted_anno = [sanno for sanno in sort_rows if sanno in anno_topN_nozero]
            mut_mat_df_plot = mut_mat_df.loc[sorted_anno, :]
        else:
            sorted_anno = anno_topN_nozero
            mut_mat_df_plot = mut_mat_df.loc[sorted_anno, :]

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
                

    def mask_tRNA_database(self, min_mut_freq=0.1, min_pos_count=10, min_tr_count=30, frac_max_score=0.90, match_score=1, mismatch_score=-1, open_gap_score=-2, extend_gap_score=-2):

        # Find and store masked sequence #
        for species in self.tr_muts:
            for anno in self.tr_muts[species]:
                tr_len = len(self.tr_muts[species][anno]['PSCM'])
                tr_seq = self.tr_muts[species][anno]['seq']
                if tr_len == 0:  # Skip no observations
                    continue
                
                # Calculate the mutation frequency per sequence position #
                # First, count all observations for a given positions:
                counts_all = self.tr_muts[species][anno]['PSCM'].sum(1)
                # Enforce a minimum number of observations:
                min_count_mask = counts_all >= min_pos_count
                if min_count_mask.sum() == 0 or counts_all.values[-2] < min_tr_count:
                    continue
                counts_all.loc[~min_count_mask] = 0
            
                # Make a mask to count all the non-mutant characters:
                count_mask = np.zeros((tr_len, len(self.char_list)))
                tr_char_idx = [self.char_dict[nt] for nt in tr_seq]
                for seq_i, char_i in enumerate(tr_char_idx):
                    count_mask[seq_i, char_i] = 1
                counts_cor = (self.tr_muts[species][anno]['PSCM'] * count_mask).sum(1)
                # np.divide setting 0/0 to 1, i.e. no observation, no mutation
                freq_mut = 1 - np.divide(counts_cor, counts_all, out=np.ones_like(counts_cor), where=counts_all!=0)
                
                # Find all positions that should be masked:
                mask_idx = set((freq_mut.values >= min_mut_freq).nonzero()[0])
                if len(mask_idx) == 0:
                    continue
                # Make masked sequence:
                seq_masked = ['N' if seq_i in mask_idx else char for seq_i, char in enumerate(tr_seq)]
                self.tr_muts[species][anno]['seq_masked'] = ''.join(seq_masked)

        # Expand masked sequence based on similar sequences #
        # Initiate the aligner:
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.match_score = match_score
        aligner.mismatch_score = mismatch_score
        aligner.open_gap_score = open_gap_score
        aligner.extend_gap_score = extend_gap_score

        for species in self.tr_muts:
            for anno1 in self.tr_muts[species]:
                target = self.tr_muts[species][anno1]['seq']
                ac1 = anno1.split('-')[2] # anticodon
                self.tr_muts[species][anno1]['seq_masked_exp'] = ''

                # Add masked characters from highly similar sequences:
                for anno2 in self.tr_muts[species]:
                    # Skip sequences with different anticodon:
                    if ac1 != anno2.split('-')[2]:
                        continue
                    # Skip sequences with no masked characters to contribute:
                    elif not 'seq_masked' in self.tr_muts[species][anno2]:
                        continue

                    query = self.tr_muts[species][anno2]['seq']
                    max_score = match_score * max([len(target), len(query)])
                    min_score = max_score * frac_max_score

                    # Perform alignment and take out one:
                    alignments = aligner.align(target, query)
                    alignment = alignments[0]
                    t_cor, q_cor = alignment.aligned

                    # Following enforces end-to-end coverage of the alignment
                    # on both target and query:
                    if t_cor[0][0] != 0 or q_cor[0][0] != 0 or t_cor[-1][1] != len(target) or q_cor[-1][1] != len(query):
                        continue
                    # Also enforce minimum alignment score:
                    elif alignment.score < min_score:
                        continue

                    # Concatenate the masked query that aligns to the target:
                    q_masked = self.tr_muts[species][anno2]['seq_masked']
                    q_masked_trans = list()
                    for i in range(0, len(q_cor)):
                        q_chunk = q_masked[q_cor[i][0]:q_cor[i][1]]
                        q_masked_trans.append(q_chunk)
                        # Insert gaps when these appear in the query
                        # i.e. characters in the target is unaligned:
                        if i < (len(q_cor)-1):
                            for j in range(t_cor[i+1][0] - t_cor[i][1]):
                                q_masked_trans.append('-')
                    q_masked_trans = ''.join(q_masked_trans)
                    assert(len(q_masked_trans) == len(target))

                    # Update masked sequence:
                    if self.tr_muts[species][anno1]['seq_masked_exp'] == '':
                        # If "N" is found in the masked query then pick it out
                        # otherwise stick to the target sequence:
                        seq_masked_exp = ''.join([qc if qc == 'N' else tc for tc, qc in zip(target, q_masked_trans)])
                        self.tr_muts[species][anno1]['seq_masked_exp'] = seq_masked_exp
                    else:
                        # Expand the masked character:
                        t_masked = self.tr_muts[species][anno1]['seq_masked_exp']
                        seq_masked_exp = ''.join([qc if qc == 'N' else tc for tc, qc in zip(t_masked, q_masked_trans)])
                        self.tr_muts[species][anno1]['seq_masked_exp'] = seq_masked_exp


        # Collect statistics on the masked sequences:
        self.mask_stats = {sp: {} for sp in self.tr_muts}
        for species in self.tr_muts:
            anno_list = list(self.tr_muts[species].keys())
            mask_count = np.array([self.tr_muts[species][anno]['seq_masked_exp'].count('N') for anno in anno_list])
            
            self.mask_stats[species]['mask_count'] = mask_count
            self.mask_stats[species]['mask_sum'] = sum(mask_count)
            self.mask_stats[species]['mask_mean'] = np.mean(mask_count)
            self.mask_stats[species]['mask_frac'] = (mask_count > 0).sum() / len(mask_count)
            

    def write_masked_tRNA_database(self, out_dir='tRNA_database_masked'):
        # Make out_dir folder:
        out_dir_abs = '{}/{}'.format(self.TM_dir_abs, out_dir)
        try:
            os.mkdir(out_dir_abs)
        except:
            shutil.rmtree(out_dir_abs)

        tRNA_database_masked = dict()
        for species in self.tr_muts:
            # Make out_dir species folder:
            out_dir_sp_abs = '{}/{}'.format(out_dir_abs, species)
            try:
                os.mkdir(out_dir_sp_abs)
            except:
                shutil.rmtree(out_dir_sp_abs)
            
            # Write the masked sequences as fasta:
            out_fnam = '{}-tRNAs.fa'.format(species)
            out_fnam_abs = '{}/{}'.format(out_dir_sp_abs, out_fnam)
            tRNA_database_masked[species] = out_fnam_abs
            with open(out_fnam_abs, 'w') as fh:
                for anno in self.tr_muts[species]:
                    if self.tr_muts[species][anno]['seq_masked_exp'] == '':
                        masked_seq = self.tr_muts[species][anno]['seq']
                    else:
                        masked_seq = self.tr_muts[species][anno]['seq_masked_exp']
                    fh.write('>{}\n{}\n'.format(anno, masked_seq))
            
            # Make BLAST DB on fasta:
            self.__makeblastdb(out_dir_sp_abs, out_fnam)

        return(tRNA_database_masked)
    
    def __makeblastdb(self, file_dir, fnam):
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



