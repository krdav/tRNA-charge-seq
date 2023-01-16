import sys, os, shutil, bz2, warnings, copy, contextlib
import Bio.Data.CodonTable
import pandas as pd
import numpy as np
from mpire import WorkerPool
from Bio import Align

### Plotting imports ###
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as mcolors
import matplotlib as mpl
from matplotlib.patches import StepPatch
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import logomaker as lm
palette = list(mcolors.TABLEAU_COLORS.keys())
sns.set_theme(style="ticks", palette="muted")
sns.set_context("talk")



def freq2ratio(freq):
    return(freq / (1 - freq))

def ratio2freq(ratio):
    return(ratio / (ratio + 1))


class TRNA_plot:
    '''
    This class is used to generate some standard plots
    from the tRNAseq data such as charge per codon,
    read per million (RPM) etc.
    '''
    def __init__(self, dir_dict, sample_df, pull_default=False):
        # Input:
        self.sample_df = sample_df
        self.dir_dict = dir_dict

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

        # Load aggregated CSV data:
        stats_fnam = '{}/ALL_stats_aggregate.csv'.format(self.stats_dir_abs)
        self.all_stats = pd.read_csv(stats_fnam, keep_default_na=False)
        
        # Create single codon filter, to filter out sequences that map to
        # tRNA sequences with different codon/anticodon:
        single_codon = list()
        for anno_str, anticodon in zip(self.all_stats['tRNA_annotation'].values, self.all_stats['anticodon'].values):
            sc = True
            anno_list = anno_str.split('@')
            for anno in anno_list:
                if not anno.split('-')[2] == anticodon:
                    sc = False
            single_codon.append(sc)
        self.all_stats['single_codon'] = single_codon
        self.all_stats['mito_codon'] = ['mito_tRNA' in anno for anno in self.all_stats['tRNA_annotation'].values]
        self.all_stats['Ecoli_ctr'] = ['Escherichia_coli' in anno for anno in self.all_stats['tRNA_annotation'].values]
        # Translate codon into single letter amino acid:
        ctab_cyto = copy.deepcopy(Bio.Data.CodonTable.generic_by_id[1].forward_table)
        # Selenocysteine has to be added:
        ctab_cyto['UGA'] = 'SeC'
        ctab_mito = copy.deepcopy(Bio.Data.CodonTable.generic_by_id[2].forward_table)
        ctab_bac = copy.deepcopy(Bio.Data.CodonTable.generic_by_id[11].forward_table)
        aa_letters = list()
        for ecoli, mito, codon in zip(self.all_stats['Ecoli_ctr'], self.all_stats['mito_codon'], self.all_stats['codon'].values):
            if mito:
                aa_letters.append(ctab_mito[codon])
            elif ecoli:
                aa_letters.append(ctab_bac[codon])
            else:
                aa_letters.append(ctab_cyto[codon])
        self.all_stats['AA_letter'] = aa_letters
                

    def make_dir(self, overwrite=True):
        # Create folder for files:
        self.plotting_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['plotting_dir'])
        try:
            os.mkdir(self.plotting_dir_abs)
        except:
            if overwrite:
                shutil.rmtree(self.plotting_dir_abs)
                os.mkdir(self.plotting_dir_abs)
            else:
                print('Folder exists and overwrite set to false... Doing nothing.')
        return(self.plotting_dir_abs)

    
    
    def plot_charge(self, compartment='cyto', plot_type='codon', plot_name='charge_plot_cyto_codon', sample_list=None):
        # Plotting the charge of all aa/codons/transcripts
        if plot_type == 'aa':
            pass
        elif plot_type == 'codon':
            pass
        elif plot_type == 'transcript':
            pass
        else:
            raise Exception('Unknown plot type specified: {}\nValid strings are either either "aa", "codon" or "transcript".'.format(plot_type))
        
        if sample_list is None:
            sample_list = [row['sample_name_unique'] for _, row in self.sample_df.iterrows()]

            
    def plot_charge_corr(self, sample_pair_list):
        # Plot the charge of sample vs sample
        pass
    
    
    def plot_RPM_corr(self, sample_pair_list):
        # Plot the RPM of sample vs sample
        pass
    
    
    
    def plot_abundance(self, compartment='cyto', plot_type='aa', plot_name='abs_plot_cyto_aa'):
        # Plotting the abundance (in RPM) of all aa/codons/transcripts
        if plot_type == 'aa':
            pass
        elif plot_type == 'codon':
            pass
        elif plot_type == 'transcript':
            pass
        else:
            raise Exception('Unknown plot type specified: {}\nValid strings are either either "aa", "codon" or "transcript".'.format(plot_type))
    
    
    
    
    def plot_mutations(self):
        # misincorporation rate matrix like Behrens figure 6
        pass
    
    
    
    def plot_coverage_matrix(self):
        # coverage as a heat plot matrix instead of stacks of bars
        pass
    
    
    
    
    
    
    def plot_coverage(self, compartment='cyto', plot_type='needle', aa_norm=False, y_norm=False, plot_name='cov_plot_cyto_needle', sample_list=None, verbose=True):
        # Check input:
        if plot_type != 'behrens' and plot_type != 'needle':
            raise Exception('Unknown plot type specified: {}\nValid strings are either either "behrens" or "needle".'.format(plot_type))
        if sample_list is None:
            sample_list = [row['sample_name_unique'] for _, row in self.sample_df.iterrows()]
        sample_list = set(sample_list)

        # Columns used for data aggregation:
        aa_cols = ['sample_name_unique', 'tRNA_annotation_len', 'align_5p_idx', 'align_3p_idx', 'AA_letter', 'count']
        # Require rows to be covered at the 3', have single codon annotation,
        # not be Ecoli controls and have a single letter amino acid code (this excludes selenocysteine):
        type_mask = (self.all_stats['3p_cover'] == True) & (self.all_stats['single_codon']) & (~self.all_stats['Ecoli_ctr']) & (self.all_stats['AA_letter'].apply(len) == 1)
        # Choose rows from requested compartment:
        if compartment == 'mito':
            type_mask &= (self.all_stats['mito_codon'])
        elif compartment == 'cyto':
            type_mask &= (~self.all_stats['mito_codon'])
        else:
            raise Exception('Unknown compartment specified: {}\nValid strings are either either "mito" or "cyto"'.format(compartment))

        if verbose:
            print('\nNow plotting sample:', end='')
        # Use a 20 color colormap:
        cmap = mpl.colormaps['tab20']
        # Print each plot to the same PDF file:
        cov_fnam = '{}/{}.pdf'.format(self.plotting_dir_abs, plot_name)
        with PdfPages(cov_fnam) as pp:
            # Loop through and generate plots for each sample:
            for _, row in self.sample_df.iterrows():
                if not row['sample_name_unique'] in sample_list:
                    continue
                print('  {}'.format(row['sample_name_unique']), end='')
                try: # Catch odd errors, like empty dataframe etc.
                    # Final rows selected:
                    sample_mask = (self.all_stats['sample_name_unique'] == row['sample_name_unique']) & type_mask

                    # Generate dataframe with coverage information.
                    # The coverage can be calculated from the count and align_5p_idx.
                    cov_df = self.all_stats.loc[sample_mask, aa_cols].copy()
                    cov_df = cov_df.groupby(aa_cols[:-1], as_index=False).agg({"count": "sum"}).reset_index(drop=True)

                    # tRNAs differ in length, so to plot coverage on the same x-axis,
                    # we map the coverage of each tRNA to the longest tRNA in the set.
                    # Mapping is done using nearest percentile indexing
                    # from the queried tRNA to the longest in the set:
                    max_len = cov_df['tRNA_annotation_len'].max()
                    min_len = cov_df['tRNA_annotation_len'].min()
                    len_map_len = dict() # prepare mapping for all tRNAs
                    for len_i in range(min_len, max_len+1):
                        len_map_len[len_i] = np.percentile(np.arange(max_len), np.linspace(0, 100, len_i), method='nearest')

                    # Get the amino acid one letter names sorted and their indexing:
                    aa_ordered_list = sorted(set(cov_df['AA_letter']), reverse=True)
                    aa_order = {aa: i for i, aa in enumerate(aa_ordered_list)}

                    # Add the coverage counts to a matrix:
                    cov_count = np.zeros((len(aa_order), max_len))
                    # Each row has an amino acid letter, align_5p_idx and a count.
                    # Insert these into the matrix:
                    for _, cov_row in cov_df.iterrows():
                        aa_idx = aa_order[cov_row['AA_letter']]
                        _5p_idx = cov_row['align_5p_idx'] - 1 # shift alignment index to 0 indexing
                        anno_len = cov_row['tRNA_annotation_len']
                        _5p_idx_trans = len_map_len[anno_len][_5p_idx]
                        # Amino acids have multiple codons, each with multiple
                        # transcripts that can vary in length, thus add instead of assign:
                        cov_count[aa_idx, _5p_idx_trans] += cov_row['count']

                    # Each count in the matrix repressents a number of reads
                    # and their length mapped to the tRNA. Now this is converted
                    # to coverage by summation from left to right i.e. 5p to 3p:
                    for i in range(cov_count.shape[0]):
                        for j in range(cov_count.shape[1]):
                            if j > 0:
                                cov_count[i, j] += cov_count[i, j-1]

                    # If "aa_norm" is requested the read count for each
                    # amino acid is weighted so all amino acids start of with
                    # the same coverage at the 3p and normalized so it sums to 1:
                    if aa_norm:
                        # Weigh each amino acid:
                        cov_count = cov_count.T / cov_count[:, -1]
                        # Normalize so coverage at 3p is 1:
                        cov_count = cov_count.T * (1/cov_count.shape[1])
                    # If "y_norm" is requested the y-axis is normalized to 1
                    elif y_norm:
                        # Normalize so coverage at 3p is summing to 100:
                        cov_count = 100 * cov_count / sum(cov_count[:, -1])

                    # "cov_count" is summed on each row from left to right,
                    # but for plotting we also need to do summation on each
                    # column from top to bottom:
                    cov_count_sum = copy.deepcopy(cov_count)
                    for i in range(cov_count_sum.shape[0]):
                        if i > 0:
                            cov_count_sum[i] += cov_count_sum[i-1]
                    # Add an additional row of zeroes:
                    cov_count_sum = np.vstack((np.zeros(cov_count_sum.shape[1]), cov_count_sum))
                    title = 'Sample: {}, coverage from {} obs'\
                    .format(row['sample_name_unique'], cov_df['count'].sum())                    
                    if plot_type == 'behrens': # Behrens type coverage plot
                        fig, axes = plt.subplots(1, 2, figsize=(13, 7), \
                                                 sharey=False, gridspec_kw={'width_ratios': [1, 7]})

                        # Plot small 10 nt. window from the 5p end:
                        for i in range(1, len(cov_count_sum)):
                            axes[0].stairs(cov_count_sum[i, 0:10], baseline=cov_count_sum[i-1, 0:10], \
                                           fill=True, alpha=0.7, color=cmap(i-1))
                        # Add a black line to outline the top:
                        axes[0].stairs(cov_count_sum[i, 0:10], baseline=None, fill=False, \
                                       edgecolor='black', alpha=1, linewidth=1)

                        # Plot the whole window, from 5p to 3p:
                        for i in range(1, len(cov_count_sum)):
                            axes[1].stairs(cov_count_sum[i], baseline=cov_count_sum[i-1], fill=True, \
                                           alpha=0.7, label=aa_ordered_list[i-1], color=cmap(i-1))
                        axes[1].stairs(cov_count_sum[i], baseline=None, fill=False, \
                                       edgecolor='black', alpha=1, linewidth=1)

                        # Invert the legend labels so they are the same order
                        # as the stacks appear on the plot:
                        handles, labels = axes[1].get_legend_handles_labels()
                        axes[1].legend(handles[::-1], labels[::-1], title='Amino acid', \
                                       bbox_to_anchor=(1.01,1), borderaxespad=0, ncol=2);

                        # Add title and axis labels:
                        if aa_norm:
                            axes[0].set_ylabel("Normalized coverage (amino acids equally weighed at 3')")
                        elif y_norm:
                            axes[0].set_ylabel("Normalized coverage (%)")
                        else:
                            axes[0].set_ylabel('Read count')
                        axes[1].set_xlabel("5' to 3' index (mapped to longest tRNA)");
                        axes[1].set_title(title)

                    elif plot_type == 'needle': # Needle-like type coverage plot
                        # Eeach "needle" follows the midpoint of
                        # the summed count for the 3p (last) position:
                        last_col = cov_count_sum[:, -1]
                        last_col_mid = np.zeros(last_col.shape[0]-1)
                        for i in range(1, last_col.shape[0]):
                            delta = last_col[i] - last_col[i-1]
                            last_col_mid[i-1] = last_col[i-1] + delta/2

                        # Generate top and bottom values at each position
                        # from 5p to 3p for each "needle":
                        cov_funnel_top = np.zeros(cov_count.shape)
                        cov_funnel_bot = np.zeros(cov_count.shape)
                        for i in range(0, cov_count.shape[1]):
                            cov = cov_count[:, i]
                            cov_funnel_top[:, i] = last_col_mid + cov/2
                            cov_funnel_bot[:, i] = last_col_mid - cov/2

                        fig, ax = plt.subplots(1, 1, figsize=(10, 7))
                        for i in range(0, len(cov_funnel_top)):
                            # Plot each "needle":
                            ax.stairs(cov_funnel_top[i, :], baseline=cov_funnel_bot[i, :], fill=True, \
                                      alpha=0.7, label=aa_ordered_list[i], color=cmap(i))

                            # Plot an outline around the needle:
                            funnel_mask = cov_funnel_top[i] > cov_funnel_bot[i]
                            edges = np.where(funnel_mask)[0]
                            edges = np.hstack((edges, edges[-1]+1))
                            ax.stairs(cov_funnel_top[i, funnel_mask], edges, \
                                      baseline=cov_funnel_bot[i, funnel_mask], \
                                      edgecolor='black', alpha=1, linewidth=1)
                        # Invert the legend labels so they are the same order
                        # as the stacks appear on the plot:
                        handles, labels = ax.get_legend_handles_labels()
                        ax.legend(handles[::-1], labels[::-1], title='Amino acid', \
                                  bbox_to_anchor=(1.01,1), borderaxespad=0, ncol=2);

                        # Add title and axis labels:
                        if aa_norm:
                            ax.set_ylabel("Normalized coverage (amino acids equally weighed at 3')")
                        else:
                            ax.set_ylabel('Read count')
                        ax.set_xlabel("5' to 3' index (mapped to longest tRNA)");
                        ax.set_title(title)

                    # Write to PDF file:
                    fig.tight_layout()
                    pp.savefig(fig, bbox_inches='tight')
                    plt.close(fig)
                except Exception as err:
                    pass
                    # print(err)

    def plot_non_temp(self, end, plot_name, seq_len_percentile=99, seq_len=None, _3p_cover=False):
        if end == '3p':
            col_sele = '3p_non-temp'
        elif end == '5p':
            col_sele = '5p_non-temp'
        else:
            raise Exception('end input has to be either 3p and 5p, not: {}'.format(end))

        if seq_len is None:
            # Limit the length of extracted 5p non-template bases to "seq_len_percentile":
            seq_len = np.percentile([len(s) for s in self.all_stats[col_sele]], seq_len_percentile, method='nearest')
        # Collect data for plotting for each sample:
        plot_data = list()
        for _, row in self.sample_df.iterrows():
            title = 'Sample: {}, logo from {} obs'.format(row['sample_name_unique'], row['N_mapped'])
            mask = np.array([len(s) <= seq_len for s in self.all_stats[col_sele].values])
            mask = mask & (self.all_stats['sample_name_unique'] == row['sample_name_unique'])
            if _3p_cover:
                mask = mask & (self.all_stats['3p_cover'] == True)
            if mask.sum() > 1:
                # The dash "-" sign is not counted in "counts_mat":
                if end == '3p':
                    seq_list = [s + '-'*(seq_len - len(s)) for s in self.all_stats.loc[mask, '3p_non-temp'].values if len(s) <= seq_len]
                else:
                    seq_list = ['-'*(seq_len - len(s)) + s for s in self.all_stats.loc[mask, col_sele].values if len(s) <= seq_len]
            # Will throw an exception if all gaps:
            try:
                with warnings.catch_warnings(): # ignoring a warning about array assignment
                    warnings.simplefilter(action='ignore', category=FutureWarning)
                    counts_mat = lm.alignment_to_matrix(seq_list, self.all_stats.loc[mask, 'count'].values)
            except:
                counts_mat = False
            plot_data.append([counts_mat, title])
        # Make logo:
        self.__run_logomaker(plot_data, plot_name)
        # Also make it as log10 transformed:
        self.__run_logomaker(plot_data, plot_name, log10=True)


    def plot_UMI_logo(self, plot_name='UMI_logo', n_jobs=4):
        # Check files exists before starting:
        for _, row in self.sample_df.iterrows():
            stats_fnam = '{}/{}_stats.csv.bz2'.format(self.stats_dir_abs, row['sample_name_unique'])
            assert(os.path.exists(stats_fnam))

        # Read each stats CSV file and generate its count matrix:
        data = list(self.sample_df.iterrows())
        with WorkerPool(n_jobs=n_jobs) as pool:
            results = pool.map(self.__collect_UMI_mat, data)
        self.__run_logomaker(results, plot_name)

    def __collect_UMI_mat(self, index, row):
        stats_fnam = '{}/{}_stats.csv.bz2'.format(self.stats_dir_abs, row['sample_name_unique'])
        with bz2.open(stats_fnam, 'rt', encoding="utf-8") as stats_fh:
            sample_stats = pd.read_csv(stats_fh, keep_default_na=False)
        try:
            UMI_mask = sample_stats['5p_UMI'] != ''
            with warnings.catch_warnings(): # ignoring a warning about array assignment
                warnings.simplefilter(action='ignore', category=FutureWarning)
                counts_mat = lm.alignment_to_matrix(sample_stats[UMI_mask]['5p_UMI'].values)
            title = 'Sample: {}, logo from {} obs (excluding common seqs), {} % of obs vs exp'.format(row['sample_name_unique'], UMI_mask.sum(), round(row['percent_UMI_obs-vs-exp']))
        except:
            counts_mat = False
            title = ''

        return([counts_mat, title])

    def __run_logomaker(self, plot_data, plot_name, log10=False):
        if log10:
            plot_name = plot_name + '_log10'
        # Print each to the same PDF file:
        logo_fnam = '{}/{}.pdf'.format(self.plotting_dir_abs, plot_name)
        with PdfPages(logo_fnam) as pp:
            for dat in plot_data:
                counts_mat, title = dat
                if counts_mat is False:
                    continue
                if log10:
                    counts_mat = np.log10(counts_mat+1)
                # logomaker prints a warning when encountering "N"
                # characters, we don't want that:
                with contextlib.redirect_stdout(None):
                    logo_plot = lm.Logo(counts_mat, color_scheme='classic');
                logo_plot.ax.set_title(title, fontsize=15)
                logo_plot.ax.set_xlabel("5' to 3'")
                logo_plot.ax.set_ylabel("Count");
                logo_plot.fig.tight_layout()
                pp.savefig(logo_plot.fig, bbox_inches='tight')
                plt.close(logo_plot.fig)


