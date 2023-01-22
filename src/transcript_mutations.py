import sys, os, shutil, bz2, copy, contextlib, gc
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
    '''
    def __init__(self, dir_dict, sample_df, tRNA_database, \
                 pull_default=False, common_seqs=None, ignore_common_count=False, \
                 overwrite_dir=False):
        # Input:
        self.sample_df, self.tRNA_database = sample_df, tRNA_database
        self.dir_dict = dir_dict
        self.common_seqs_fnam = common_seqs
        self.common_seqs_dict = dict() # map common sequences to index
        self.char_str = 'ACGTUN-'
        self.char_list = [c for c in self.char_str]
        self.char_dict = {c: i for i, c in enumerate(self.char_str)}
        self.tr_muts_masked = None
        
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
        self.__make_dir(overwrite_dir)

    def __make_dir(self, overwrite):
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

    def find_muts(self, unique_anno=True, match_score=1, mismatch_score=-1, \
                  open_gap_score=-2, extend_gap_score=-1, n_jobs=4, verbose=True, \
                  sample_list=None, fix_end=True):
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
        # Fill out the transcript mutations per sample:
        for unam_res in results:
            unam, res = unam_res
            self.tr_muts[unam] = res
        # Fix ends if requested:
        if fix_end:
            self.__fix_end()

    def __collect_transcript_muts(self, index, row):
        if self.verbose:
            print('  {}'.format(row['sample_name_unique']), end='')

        species = row['species']
        # Dictionary to store mutation info for each transcript:
        tr_muts_sp = copy.deepcopy(self.tr_muts_tmp)
        
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
            sample_stats = pd.read_csv(stats_fh, keep_default_na=False, low_memory=False)
        ID2anno = {rid: tan.split('@') for rid, tan, _3c in zip(sample_stats['readID'].values, sample_stats['tRNA_annotation'].values, sample_stats['3p_cover'].values) if _3c}
        sample_stats = None
        del sample_stats
        gc.collect()

        for seq in dedup_seq_count:
            # Check if common sequence, then readID has changed:
            if not self.common_seqs_fnam is None and seq in self.common_seqs_dict:
                readID = str(self.common_seqs_dict[seq])
            else:
                readID = dedup_seq_count[seq]['id']

            # Get list of annotations:
            if readID in ID2anno:
                anno_list = ID2anno[readID]
            else:
                # Skip unaligned reads:
                continue
            # Skip if multiple annotations found but unique requested:
            if self.unique_anno and len(anno_list) > 1:
                continue

            for anno in anno_list:
                # Generate alignments:
                target = tr_muts_sp[species][anno]['seq']
                alignments = aligner.align(target, seq)
                # If multiple alignments with the same score these should be weighted
                # so one read contributes with one observation:
                weight = 1.0 / len(alignments) * dedup_seq_count[seq]['count']
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
                    # weight variable, hence the dtype needs more than int8:
                    count_mat = np.zeros((len(target), len(self.char_list)), dtype=np.uint32)
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

        # Convert the count matrix to dataframe
        # calculate mutation/gap frequencies and return:
        for anno in tr_muts_sp[species]:
            tr_muts_sp[species][anno]['PSCM'] = pd.DataFrame(tr_muts_sp[species][anno]['PSCM'], columns=self.char_list)
            if tr_muts_sp[species][anno]['PSCM'].max().max() > 0:
                # Mutation frequencies:
                freq_mut = self.__calc_mut_freq(tr_muts_sp, anno, species, gap_only=False, min_count_show=0)
                tr_muts_sp[species][anno]['mut_freq'] = freq_mut
                # Gap frequencies:
                freq_gap = self.__calc_mut_freq(tr_muts_sp, anno, species, gap_only=True, min_count_show=0)
                tr_muts_sp[species][anno]['gap_freq'] = freq_gap

        return((row['sample_name_unique'], tr_muts_sp))
 
    def __combine_tr_muts(self, sample_list, freq_avg_weighted=True):
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
                else:
                    tr_muts_combi[species][anno]['PSCM'] += sp_muts[species][anno]['PSCM'].copy()
                    # Take the cumulative average for the frequencies:
                    tr_muts_combi[species][anno]['mut_freq'] = (sp_muts[species][anno]['mut_freq'] + avg_count*tr_muts_combi[species][anno]['mut_freq']) / (avg_count+1)
                    tr_muts_combi[species][anno]['gap_freq'] = (sp_muts[species][anno]['gap_freq'] + avg_count*tr_muts_combi[species][anno]['gap_freq']) / (avg_count+1)
            avg_count += 1

            # Get the frequency average weighted by total observations:
            if freq_avg_weighted:
                for anno in sp_muts[species]:
                    # Skip if no observations:
                    if sp_muts[species][anno]['PSCM'].max().max() == 0:
                        continue
                    tr_muts_combi[species][anno]['mut_freq'] = self.__calc_mut_freq(tr_muts_combi, anno, species, gap_only=False, min_count_show=0)
                    tr_muts_combi[species][anno]['gap_freq'] = self.__calc_mut_freq(tr_muts_combi, anno, species, gap_only=True, min_count_show=0)

        if len(sample_list_cp) > 0:
            print('Following samples could not be found and therefore not combined: {}'.format(str(sample_list_cp)))
        return(tr_muts_combi)

    def __fix_end(self):
        for unam in self.tr_muts:
            for species in self.tr_muts[unam]:
                for anno in self.tr_muts[unam][species]:
                    seq_len = self.tr_muts[unam][species][anno]['seq_len']
                    if self.tr_muts[unam][species][anno]['PSCM'].max().max() > 0:
                        end_obs = self.tr_muts[unam][species][anno]['PSCM']['C'].values[-2]
                        end_ar = np.zeros(len(self.char_list))
                        A_idx = self.char_list.index('A')
                        end_ar[A_idx] = end_obs
                        self.tr_muts[unam][species][anno]['PSCM'].loc[seq_len-1, :] = end_ar

    def plot_transcript_logo(self, topN=30, species='human', plot_name='tr-mut_logos', \
                             sample_list=None):
        # Get the mutations combined for the requested samples:
        if sample_list is None:
            sample_list = list(self.tr_muts.keys())
        tr_muts_combi = self.__combine_tr_muts(sample_list)

        # Sort according to observations:
        anno_sorted = self.__sort_anno(tr_muts_combi, species)
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
                            sort_rows=True, sample_list=None):
        # Get the mutations combined for the requested samples:
        if sample_list is None:
            sample_list = list(self.tr_muts.keys())
        tr_muts_combi = self.__combine_tr_muts(sample_list)

        # Sort according to observations:
        anno_sorted = self.__sort_anno(tr_muts_combi, species, mito=mito)
        if topN > len(anno_sorted):
            topN = len(anno_sorted)
        
        # Find the observations and insert them into a matrix:
        obs_mat = np.zeros((topN, self.longest_tRNA))
        anno_topN = [anno_sorted[i][0] for i in range(topN)]
        for i in range(topN):
            anno = anno_topN[i]
            if mito and 'mito' not in anno:
                continue
            # Count all observations for a given positions:
            counts_all = tr_muts_combi[species][anno]['PSCM'].sum(1)
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




    def plot_transcript_mut_compare(self, species='human', \
                                    plot_name='tr-mut_matrix_comp', \
                                    no_plot_return=False, \
                                    mito=False, gap_only=False, \
                                    min_count_show=100, \
                                    sample_pairs=None, sample_unique_pairs=None, \
                                    tr_compare_inp=None, \
                                    anno_substring_compare=None,\
                                    sample_list_exl=None, bc_list_exl=None,
                                    freq_avg_weighted=True, \
                                    topN=10, topN_select='max_diff'):

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
                mask1 = mask_exl & (self.sample_df['sample_name'] == sp1)
                mask2 = mask_exl & (self.sample_df['sample_name'] == sp2)
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
                tr_muts_combi_s1 = self.__combine_tr_muts(sp1l, freq_avg_weighted=freq_avg_weighted)
                tr_muts_combi_s2 = self.__combine_tr_muts(sp2l, freq_avg_weighted=freq_avg_weighted)

                # Find the transcripts to compare:
                if not anno_substring_compare is None:
                    tr_compare_list = [anno for anno in tr_muts_combi_s1[species] if anno_substring_compare in anno]
                elif tr_compare_inp is None:
                    # Sort annotations according to largest distance between samples:
                    topN_anno = self.__sort_freq_diff(tr_muts_combi_s1, tr_muts_combi_s2, species, min_count_show, gap_only, mito, topN_select)
                    tr_compare_list = topN_anno[:topN]

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
                    freq_mut_s1, min_count_mask_s1 = self.__get_mut_freq_filted(tr_muts_combi_s1, species, anno, min_count_show, gap_only)
                    freq_mut_s2, min_count_mask_s2 = self.__get_mut_freq_filted(tr_muts_combi_s2, species, anno, min_count_show, gap_only)
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

    def __sort_freq_diff(self, tr_muts_combi_s1, tr_muts_combi_s2, species, min_count_show, gap_only, mito, topN_select):
        topN_anno = list()
        for anno in tr_muts_combi_s1[species]:
            if not anno in tr_muts_combi_s2[species]:
                continue
            # If mito is specified, skip non-mito annotations:
            if mito and 'mito' not in anno:
                continue
            freq_mut_s1, min_count_mask_s1 = self.__get_mut_freq_filted(tr_muts_combi_s1, species, anno, min_count_show, gap_only)
            freq_mut_s2, min_count_mask_s2 = self.__get_mut_freq_filted(tr_muts_combi_s2, species, anno, min_count_show, gap_only)
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

    def __get_mut_freq_filted(self, tr_muts_combi, species, anno, min_count_show, gap_only):
        if gap_only:
            freq_mut = tr_muts_combi[species][anno]['gap_freq']
        else:
            freq_mut = tr_muts_combi[species][anno]['mut_freq']
        # Enforce a minimum number of observations:
        counts_all = tr_muts_combi[species][anno]['PSCM'].sum(1)
        min_count_mask = counts_all >= min_count_show
        if min_count_mask.sum() == 0:
            return(None, None)
        freq_mut[~min_count_mask] = 0
        return(freq_mut, min_count_mask)

    def plot_transcript_mut(self, topN=50, species='human', plot_name='tr-mut_matrix', \
                            png_dpi=False, no_plot_return=False, mito=False, gap_only=False, \
                            sort_rows=True, min_count_show=10, sample_list=None,
                            freq_avg_weighted=True):
        # Get the mutations combined for the requested samples:
        if sample_list is None:
            sample_list = list(self.tr_muts.keys())
        tr_muts_combi = self.__combine_tr_muts(sample_list, freq_avg_weighted=freq_avg_weighted)

        # Sort according to observations:
        anno_sorted = self.__sort_anno(tr_muts_combi, species, mito=mito)
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

    def __sort_anno(self, tr_muts_combi, species, mito=False):
        # Sort according to observations:
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

    def __calc_mut_freq(self, tr_muts_combi, anno, species, gap_only, min_count_show):
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

    def mask_tRNA_database(self, min_mut_freq=0.1, min_pos_count=10, min_tr_count=30, \
                           frac_max_score=0.90, match_score=1, mismatch_score=-1, \
                           open_gap_score=-2, extend_gap_score=-2, sample_list=None,
                           freq_avg_weighted=True):
        # Get the mutations combined for the requested samples:
        if sample_list is None:
            sample_list = list(self.tr_muts.keys())
        tr_muts_combi = self.__combine_tr_muts(sample_list, freq_avg_weighted=freq_avg_weighted)

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
        aligner.mode = 'local'
        aligner.match_score = match_score
        aligner.mismatch_score = mismatch_score
        aligner.open_gap_score = open_gap_score
        aligner.extend_gap_score = extend_gap_score

        for species in tr_muts_combi:
            for anno1 in tr_muts_combi[species]:
                target = tr_muts_combi[species][anno1]['seq']
                ac1 = anno1.split('-')[2] # anticodon
                tr_muts_combi[species][anno1]['seq_masked_exp'] = ''

                # Add masked characters from highly similar sequences:
                for anno2 in tr_muts_combi[species]:
                    # Skip sequences with different anticodon:
                    if ac1 != anno2.split('-')[2]:
                        continue
                    # Skip sequences with no masked characters to contribute:
                    elif not 'seq_masked' in tr_muts_combi[species][anno2]:
                        continue

                    query = tr_muts_combi[species][anno2]['seq']
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
                    q_masked = tr_muts_combi[species][anno2]['seq_masked']
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
                    if tr_muts_combi[species][anno1]['seq_masked_exp'] == '':
                        # If "N" is found in the masked query then pick it out
                        # otherwise stick to the target sequence:
                        seq_masked_exp = ''.join([qc if qc == 'N' else tc for tc, qc in zip(target, q_masked_trans)])
                        tr_muts_combi[species][anno1]['seq_masked_exp'] = seq_masked_exp
                    else:
                        # Expand the masked character:
                        t_masked = tr_muts_combi[species][anno1]['seq_masked_exp']
                        seq_masked_exp = ''.join([qc if qc == 'N' else tc for tc, qc in zip(t_masked, q_masked_trans)])
                        tr_muts_combi[species][anno1]['seq_masked_exp'] = seq_masked_exp

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



