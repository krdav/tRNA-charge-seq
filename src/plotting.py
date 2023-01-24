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
    def __init__(self, dir_dict, sample_df=None, stats_fnam=None, \
                 pull_default=False, overwrite_dir=False):
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
        self.stats_agg_cols = ['sample_name_unique', 'sample_name', 'replicate', 'barcode', \
                               'tRNA_annotation', 'tRNA_annotation_len', 'unique_annotation', \
                               '5p_cover', 'align_3p_nt', 'codon', 'anticodon', 'amino_acid', \
                               'count']
        self.stats_agg_cols_type = [str, str, int, str, str, int, bool, \
                                    bool, str, str, str, str, int]
        self.stats_agg_cols_td = {nam:tp for nam, tp in zip(self.stats_agg_cols, self.stats_agg_cols_type)}

        # Input:
        self.dir_dict = dir_dict
        self.charge_df = None
        self.charge_filt = dict()

        self.stats_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['stats_dir'])
        self.align_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['align_dir'])
        self.UMI_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['UMI_dir'])
        # Attempt to load sample_df from default path:
        if pull_default:
            sample_df_path = '{}/{}'.format(self.align_dir_abs, 'sample_stats.xlsx')
            if os.path.isfile(sample_df_path):
                sample_df = pd.read_excel(sample_df_path, index_col=0)
            else:
                raise Exception('Default "sample_df" could not be found: {}'.format(sample_df_path))
        elif sample_df is None:
            raise Exception('A "sample_df" must be specified or "pulled from default".')
        self.sample_df = sample_df

        # Load aggregated CSV data:
        if stats_fnam is None:
            stats_fnam = '{}/ALL_stats_aggregate.csv'.format(self.stats_dir_abs)
        self.all_stats = pd.read_csv(stats_fnam, keep_default_na=False, dtype=self.stats_agg_cols_td)
        # Get rid of 1/2 in mito Leu amino acid name:
        self.all_stats['amino_acid'] = [AA[:-1] if AA[-1]=='1' or AA[-1]=='2' else AA for AA in self.all_stats['amino_acid'].values]
        # Add amino acid - codon string:
        self.all_stats['AA_codon'] = [AA + '-' + codon for codon, AA in zip(self.all_stats['codon'].values, self.all_stats['amino_acid'].values)]
        # tRNA annotation short form:
        self.all_stats['tRNA_anno_short'] = ['-'.join(an.split('-')[1:]) for an in self.all_stats['tRNA_annotation'].values]

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
        # Create single amino acid filter, to filter out sequences that map to
        # tRNA sequences with different codon/anticodon:
        single_aa = list()
        for anno_str, amino_acid in zip(self.all_stats['tRNA_annotation'].values, self.all_stats['amino_acid'].values):
            sa = True
            anno_list = anno_str.split('@')
            for anno in anno_list:
                if not anno.split('-')[1] == amino_acid:
                    sa = False
            single_aa.append(sa)
        self.all_stats['single_aa'] = single_aa
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

        # Make output folder:
        self.__make_dir(overwrite_dir)

    def __make_dir(self, overwrite):
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

    def get_charge_df(self):
        # Filter and group to only extract rows
        # valid for charge determination:
        stats_agg_cols = ['sample_name_unique', 'sample_name', 'replicate', 'barcode', \
                          'tRNA_annotation', 'tRNA_anno_short', 'tRNA_annotation_len', \
                          'unique_annotation', 'codon', 'anticodon', 'AA_codon', 'amino_acid', \
                          'single_codon', 'single_aa', 'mito_codon', 'Ecoli_ctr', 'AA_letter', \
                          'align_3p_nt', 'count']
        row_mask = ((self.all_stats['align_3p_nt'] == 'A') | (self.all_stats['align_3p_nt'] == 'C'))
        all_stats_filt = self.all_stats.loc[row_mask, stats_agg_cols]
        all_stats_filt = all_stats_filt.groupby(stats_agg_cols[:-1], as_index=False).agg({"count": "sum"}).reset_index(drop=True)

        # Count A/C endings:
        ac_dict = dict()
        for snu, tan, _3nt, count in zip(all_stats_filt['sample_name_unique'], all_stats_filt['tRNA_annotation'], all_stats_filt['align_3p_nt'], all_stats_filt['count']):
            key = (snu, tan)
            if key in ac_dict:
                ac_dict[key][_3nt] = count
            else:
                ac_dict[key] = dict()
                ac_dict[key][_3nt] = count

        # Add 0 count if count is missing:
        for key in ac_dict.keys():
            if 'A' not in ac_dict[key]:
                ac_dict[key]['A'] = 0
            elif 'C' not in ac_dict[key]:
                ac_dict[key]['C'] = 0

        # Make charge dataframe:
        row_mask = np.zeros(len(all_stats_filt), dtype=bool)
        bag = set() # avoid duplicates
        a_count = list()
        c_count = list()
        for row_idx, key in enumerate(zip(all_stats_filt['sample_name_unique'], all_stats_filt['tRNA_annotation'])):    
            # Take the first encounter of a row:
            if key in bag:
                continue
            else:
                bag.add(key)
            row_mask[row_idx] = True

            # Get and add the A/C count:
            a_count.append(ac_dict[key]['A'])
            c_count.append(ac_dict[key]['C'])

        # Extract the dataframe and calculate the charge:
        charge_df = all_stats_filt.loc[row_mask, stats_agg_cols[:-2]].copy()
        charge_df['A_count'] = a_count
        charge_df['C_count'] = c_count
        charge_df['count'] = charge_df['A_count']+charge_df['C_count']
        charge_df['charge'] = 100*charge_df['A_count'] / charge_df['count']

        # Calculate read per million (RPM).
        # Notice, the Ecoli control is not counted in the total:
        df_count = charge_df[~charge_df['Ecoli_ctr']].groupby(['sample_name_unique'], as_index=False).agg({"count": "sum"}).reset_index(drop=True)
        # Add the sample total count to the rows:
        charge_df = charge_df.merge(df_count, on='sample_name_unique', suffixes=('', '_sample_tot'))
        # Calculated the RPM and get rid of the total count:
        charge_df['RPM'] = charge_df['count'] / (charge_df['count_sample_tot'] / 1e6)
        charge_df = charge_df.drop(columns=['count_sample_tot'])

        # Filter and group by same amino acid:
        aa_mask = charge_df['single_aa']
        charge_df_aa = charge_df[aa_mask].groupby(['sample_name_unique', 'sample_name', 'replicate', 'barcode', 'amino_acid', 'AA_letter', 'mito_codon', 'Ecoli_ctr'], as_index=False).agg({"count": "sum", "A_count": "sum", "C_count": "sum", "RPM": "sum"}).reset_index(drop=True)
        charge_df_aa['charge'] = 100 * charge_df_aa['A_count'] / charge_df_aa['count']
        self.charge_filt['aa'] = charge_df_aa

        # Filter and group by same codon:
        cd_mask = charge_df['single_codon']
        charge_df_cd = charge_df[cd_mask].groupby(['sample_name_unique', 'sample_name', 'replicate', 'barcode', 'codon', 'anticodon', 'AA_codon', 'amino_acid', 'AA_letter', 'mito_codon', 'Ecoli_ctr'], as_index=False).agg({"count": "sum", "A_count": "sum", "C_count": "sum", "RPM": "sum"}).reset_index(drop=True)
        charge_df_cd['charge'] = 100 * charge_df_cd['A_count'] / charge_df_cd['count']
        self.charge_filt['codon'] = charge_df_cd

        # Filter and group by same transcript:
        tr_mask = charge_df['unique_annotation']
        charge_df_tr = charge_df[tr_mask].groupby(['sample_name_unique', 'sample_name', 'replicate', 'barcode', 'tRNA_annotation', 'tRNA_anno_short', 'tRNA_annotation_len', 'codon', 'anticodon', 'AA_codon', 'amino_acid', 'AA_letter', 'mito_codon', 'Ecoli_ctr'], as_index=False).agg({"count": "sum", "A_count": "sum", "C_count": "sum", "RPM": "sum"}).reset_index(drop=True)
        charge_df_tr['charge'] = 100 * charge_df_tr['A_count'] / charge_df_tr['count']
        self.charge_filt['tr'] = charge_df_tr
        self.charge_df = charge_df

    def write_charge_df(self, fnam='charge_df.csv'):
        fnam_abs = '{}/{}.pdf'.format(self.plotting_dir_abs, fnam)
        if self.charge_df is None:
            self.get_charge_df()
        self.charge_df.to_csv(fnam_abs, header=True, index=False)

    def plot_Ecoli_ctr(self, plot_name='ecoli-ctr_charge_plot', sample_list=None,
                       charge_plot=False, min_obs=1, \
                       sample_list_exl=None, bc_list_exl=None):

        if charge_plot:
            y_axis = 'charge'
        else:
            y_axis = 'RPM'

        if sample_list is None:
            sample_list = list(set(self.sample_df['sample_name_unique']))

        # Sample rows selected:
        charge_df_type = self.charge_filt['tr'].copy()
        sample_mask = np.array([False]*len(charge_df_type))
        for unam in sample_list:
            sample_mask |= (charge_df_type['sample_name_unique'] == unam)

        # Exclude unique samples or barcodes:
        if not sample_list_exl is None:
            for unam in sample_list_exl:
                sample_mask &= (charge_df_type['sample_name_unique'] != unam)
        if not bc_list_exl is None:
            for bc in bc_list_exl:
                sample_mask &= (charge_df_type['barcode'] != bc)

        # Only take rows with Ecoli controls:
        sample_mask &= charge_df_type['Ecoli_ctr'] & (charge_df_type['count'] >= min_obs)
        charge_sample = charge_df_type[sample_mask].copy()
        
        # Different settings for different plot types:
        Ngrp = len(set(charge_sample['sample_name_unique']))
        x_list = sorted(set(charge_sample['sample_name_unique'].values), key=str.casefold)

        # Make the plot:
        fig, ax = plt.subplots(1, 1, figsize=(0.6*Ngrp, 6))
        g1 = sns.barplot(ax=ax, x='sample_name_unique', y=y_axis, order=x_list, hue='tRNA_anno_short', \
                         data=charge_sample, capsize=.1, errwidth=2, \
                         edgecolor='black', linewidth=2, alpha=0.85)

        # Set axis/title text:
        g1.set_title('E.coli control samples')
        g1.set_xticklabels(g1.get_xticklabels(), rotation=90)
        g1.grid(True, axis='y')
        g1.set_xlabel('');
        
        if charge_plot is True:
            g1.set_ylabel('Charge (%)');
            g1.set(ylim=[0, 100])
        else:
            g1.set_ylabel('Abundance (RPM)');

        # Write to PDF file:
        plot_fnam = '{}/{}.pdf'.format(self.plotting_dir_abs, plot_name)
        fig.tight_layout()
        fig.savefig(plot_fnam, bbox_inches='tight')
        plt.close(fig)

    def plot_abundance(self, plot_type='aa', group=False, charge_plot=False, \
        plot_name='abundance_plot_aa', min_obs=100, sample_list=None, verbose=True, \
        sample_list_exl=None, bc_list_exl=None):
        if self.charge_df is None:
            self.get_charge_df()

        if plot_type == 'aa':
            charge_df_type = self.charge_filt['aa'].copy()
        elif plot_type == 'codon':
            charge_df_type = self.charge_filt['codon'].copy()
        elif plot_type == 'transcript':
            charge_df_type = self.charge_filt['tr'].copy()
        else:
            raise Exception('Unknown plot type specified: {}\nValid strings are either either "aa", "codon" or "transcript".'.format(plot_type))
        
        if charge_plot:
            y_axis = 'charge'
        else:
            y_axis = 'RPM'

        if sample_list is None:
            sample_list = list(self.sample_df['sample_name'].drop_duplicates())

        if verbose:
            print('\nNow plotting sample/group:', end='')

        # Use a 40 color colormap:
        cmap_b = mpl.colormaps['tab20']
        cmap_c = mpl.colormaps['tab20']
        cmap_b.colors = cmap_b.colors + cmap_c.colors
        cmap_b.N = 40
        cmap_b._i_bad = 42
        cmap_b._i_over = 40
        cmap_b._i_under = 0

        snam_df = self.sample_df.loc[:, ['sample_name', 'plot_group', 'hue_name', 'hue_value', 'hue_order']].drop_duplicates()
        charge_df_type = charge_df_type.merge(snam_df, on='sample_name')
        if group:
            plot_iter = snam_df['plot_group'].drop_duplicates()
        else:
            plot_iter = snam_df['sample_name'].drop_duplicates()

        # Print each plot to the same PDF file:
        charge_fnam = '{}/{}.pdf'.format(self.plotting_dir_abs, plot_name)
        with PdfPages(charge_fnam) as pp:
            # Loop through and generate plots for each sample:
            for sg_name in plot_iter:
                if group is False and not sg_name in sample_list:
                    continue
                print('  {}'.format(sg_name), end='')

                # Sample rows selected:
                if group:
                    sample_mask = (charge_df_type['plot_group'] == sg_name) & (charge_df_type['count'] >= min_obs)
                else:
                    sample_mask = (charge_df_type['sample_name'] == sg_name) & (charge_df_type['count'] >= min_obs)

                # Exclude unique samples or barcodes:
                if not sample_list_exl is None:
                    for unam in sample_list_exl:
                        sample_mask &= (charge_df_type['sample_name_unique'] != unam)
                if not bc_list_exl is None:
                    for bc in bc_list_exl:
                        sample_mask &= (charge_df_type['barcode'] != bc)

                charge_sample = charge_df_type[sample_mask].copy()
                # Plot separate for mito/cyto:
                mask_cyto = (~charge_sample['Ecoli_ctr']) & (~charge_sample['mito_codon'])
                mask_mito = (~charge_sample['Ecoli_ctr']) & (charge_sample['mito_codon'])

                # Different settings for different plot types:
                Ngrp = len(set(charge_sample['sample_name']))
                hue_order = [t[0] for t in sorted(dict(zip(charge_sample['hue_value'].values, charge_sample['hue_order'].values)).items(), key=lambda x: x[1])]
                if plot_type == 'aa':
                    x_axis = 'amino_acid'
                    x_list_cyto = sorted(set(charge_sample['amino_acid'][mask_cyto].values), key=str.casefold)
                    x_list_mito = sorted(set(charge_sample['amino_acid'][mask_mito].values), key=str.casefold)
                    AA_set = set(x_list_cyto+x_list_mito)
                    AAi = {aa:i for i, aa in enumerate(sorted(AA_set))}
                    colors_cyto = {c: cmap_b(AAi[c]) for c in x_list_cyto}
                    colors_mito = {c: cmap_b(AAi[c]) for c in x_list_mito}
                    fig = plt.figure(figsize=(8*Ngrp, 9))
                    gs = fig.add_gridspec(2, 4)
                    ax1 = fig.add_subplot(gs[0, :])
                    ax2 = fig.add_subplot(gs[1, :])
                elif plot_type == 'codon':
                    x_axis = 'AA_codon'
                    x_list_cyto = sorted(set(charge_sample['AA_codon'][mask_cyto].values), key=str.casefold)
                    x_list_mito = sorted(set(charge_sample['AA_codon'][mask_mito].values), key=str.casefold)
                    AA_set = {c.split('-')[0] for c in x_list_cyto+x_list_mito}
                    AAi = {aa:i for i, aa in enumerate(sorted(AA_set))}
                    colors_cyto = {c: cmap_b(AAi[c.split('-')[0]]) for c in x_list_cyto}
                    colors_mito = {c: cmap_b(AAi[c.split('-')[0]]) for c in x_list_mito}
                    fig = plt.figure(figsize=(16*Ngrp, 9))
                    gs = fig.add_gridspec(2, 4)
                    ax1 = fig.add_subplot(gs[0, :])
                    ax2 = fig.add_subplot(gs[1, 1:3])
                elif plot_type == 'transcript':
                    x_axis = 'tRNA_anno_short'
                    mask_cyto = (~charge_sample['Ecoli_ctr']) & (~charge_sample['mito_codon'])
                    mask_mito = (~charge_sample['Ecoli_ctr']) & (charge_sample['mito_codon'])
                    x_list_cyto = sorted(set(charge_sample['tRNA_anno_short'][mask_cyto].values), key=str.casefold)
                    x_list_mito = sorted(set(charge_sample['tRNA_anno_short'][mask_mito].values), key=str.casefold)
                    colors_cyto = None
                    colors_mito = None
                    fig = plt.figure(figsize=(45*Ngrp, 15))
                    gs = fig.add_gridspec(2, 10)
                    ax1 = fig.add_subplot(gs[0, :])
                    ax2 = fig.add_subplot(gs[1, 4:6])

                # Plot either as grouped or single bars:
                if group:
                    # Cyto tRNAs
                    g1 = sns.barplot(ax=ax1, x=x_axis, y=y_axis, hue='hue_value', hue_order=hue_order, order=x_list_cyto, data=charge_sample[mask_cyto], capsize=.025, errwidth=2, edgecolor='black', linewidth=2, alpha=0.8)
                    # Mito tRNAs
                    g2 = sns.barplot(ax=ax2, x=x_axis, y=y_axis, hue='hue_value', hue_order=hue_order, order=x_list_mito, data=charge_sample[mask_mito], capsize=.025, errwidth=2, edgecolor='black', linewidth=2, alpha=0.8)
                    # Legend:
                    g2.legend_.remove()
                    old_legend = g1.legend_
                    handles = old_legend.legendHandles
                    labels = hue_order
                    title = charge_sample['hue_name'].values[0]
                    g1.legend(handles, labels, title=title, bbox_to_anchor=(1.01,1), borderaxespad=0)
                else:
                    # Cyto tRNAs
                    g1 = sns.barplot(ax=ax1, x=x_axis, y=y_axis, order=x_list_cyto, data=charge_sample[mask_cyto], capsize=.1, errwidth=2, edgecolor='black', linewidth=2, alpha=0.85, palette=colors_cyto)
                    # Mito tRNAs
                    g2 = sns.barplot(ax=ax2, x=x_axis, y=y_axis, order=x_list_mito, data=charge_sample[mask_mito], capsize=.1, errwidth=2, edgecolor='black', linewidth=2, alpha=0.85, palette=colors_mito)

                # Set axis/title text:
                g1.set_title('Cytoplasmic tRNA for sample/group {}'.format(sg_name))
                g1.set_xticklabels(g1.get_xticklabels(), rotation=90)
                g1.grid(True, axis='y')
                g1.set_xlabel('');
                g2.set_title('Mitochondrial tRNA for sample/group {}'.format(sg_name))
                g2.set_xticklabels(g2.get_xticklabels(), rotation=90)
                g2.grid(True, axis='y')
                g2.set_xlabel('');

                if charge_plot is True:
                    g1.set_ylabel('Charge (%)');
                    g2.set_ylabel('Charge (%)');
                    g1.set(ylim=[0, 100])
                    g2.set(ylim=[0, 100])
                else:
                    g1.set_ylabel('Abundance (RPM)');
                    g2.set_ylabel('Abundance (RPM)');

                # Write to PDF file:
                fig.tight_layout()
                pp.savefig(fig, bbox_inches='tight')
                plt.close(fig)

    def plot_abundance_corr(self, sample_pairs=None, sample_unique_pairs=None, \
        plot_type='aa', charge_plot=False, log=False, plot_name='abundance_corr_plot_aa', \
        min_obs=100, sample_list=None, verbose=True, sample_list_exl=None, \
        bc_list_exl=None):
        if self.charge_df is None:
            self.get_charge_df()

        if plot_type == 'aa':
            charge_df_type = self.charge_filt['aa'].copy()
            merge_on = ['amino_acid', 'mito_codon', 'Ecoli_ctr']
        elif plot_type == 'codon':
            charge_df_type = self.charge_filt['codon'].copy()
            merge_on = ['AA_codon', 'mito_codon', 'Ecoli_ctr']
        elif plot_type == 'transcript':
            charge_df_type = self.charge_filt['tr'].copy()
            merge_on = ['tRNA_annotation', 'mito_codon', 'Ecoli_ctr']
        else:
            raise Exception('Unknown plot type specified: {}\nValid strings are either either "aa", "codon" or "transcript".'.format(plot_type))

        # Enforce minimum observations,
        # and not Ecoli control:
        min_obs_mask = (charge_df_type['count'] > min_obs) & ~charge_df_type['Ecoli_ctr']
        charge_df_type = charge_df_type[min_obs_mask]

        if not sample_pairs is None:
            assert(len(sample_pairs[0]) == len(sample_pairs[1]))
            pair_list = sample_pairs

            # Exclude samples/barcodes:
            mask_exl = np.array([True]*len(charge_df_type))
            if not sample_list_exl is None:
                for snam in sample_list_exl:
                    mask_exl &= (charge_df_type['sample_name_unique'] != snam)
            if not bc_list_exl is None:
                for bc in bc_list_exl:
                    mask_exl &= (charge_df_type['barcode'] != bc)
            charge_df_type = charge_df_type[mask_exl]
            
            # For sample names, merge replicates by taking the mean:
            colnames = charge_df_type.columns
            gr_cols = colnames.drop(['sample_name_unique', 'replicate', 'barcode', 'count', 'A_count', 'C_count', 'RPM', 'charge'])
            charge_df_type = charge_df_type.groupby(list(gr_cols), as_index=False).agg({'RPM': 'mean', 'charge': 'mean'}).reset_index(drop=True)
            charge_df_type_stdev = charge_df_type.groupby(list(gr_cols), as_index=False).agg({'RPM': 'std', 'charge': 'std'}).reset_index(drop=True)
        elif not sample_unique_pairs is None:
            assert(len(sample_unique_pairs[0]) == len(sample_unique_pairs[1]))
            pair_list = sample_unique_pairs
        else:
            raise Exception('Neither sample_pairs nor sample_unique_pairs lists were provided.')

        if charge_plot:
            metric = 'charge'
        else:
            metric = 'RPM'

        if sample_list is None:
            sample_list = list(self.sample_df['sample_name'].drop_duplicates())

        if verbose:
            print('\nNow plotting sample pairs:', end='')


        # Print each plot to the same PDF file:
        corr_fnam = '{}/{}.pdf'.format(self.plotting_dir_abs, plot_name)
        with PdfPages(corr_fnam) as pp:
            # Plot each pair:
            for p1, p2 in zip(*pair_list):
                # Either as unique samples are with replicates:
                if not sample_pairs is None:
                    mask1 = charge_df_type['sample_name'] == p1
                    mask2 = charge_df_type['sample_name'] == p2
                else:
                    mask1 = charge_df_type['sample_name_unique'] == p1
                    mask2 = charge_df_type['sample_name_unique'] == p2

                if mask1.sum() == 0 or mask2.sum() == 0:
                    continue
                print('  ({} - {})'.format(p1, p2), end='')

                # If charge is plotted both samples
                # have to be present, if RPM is plotted
                # it can be zero in one sample and >0 in another:
                if not charge_plot and min_obs == 0:
                    merge_how = 'outer'
                else:
                    merge_how = 'inner'

                # Merge the two:
                m12 = pd.merge(
                    charge_df_type[mask1],
                    charge_df_type[mask2],
                    how=merge_how,
                    on=merge_on,
                    left_index=False,
                    right_index=False,
                    sort=True,
                    suffixes=("_1", "_2"),
                    copy=True,
                    indicator=False,
                    validate=None,
                )
                m12 = m12.fillna(0)

                # Handle if log transform is requested:
                if log:
                    mask0 = (m12[metric+'_1'] > 0) & (m12[metric+'_2'] > 0)
                    m12 = m12[mask0]
                    m12[metric+'_1'] = np.log10(m12[metric+'_1'])
                    m12[metric+'_2'] = np.log10(m12[metric+'_2'])

                # Plot the scatterplot:
                fig, ax = plt.subplots(1, 1, figsize=(7, 6))
                g1 = sns.regplot(ax=ax, data=m12, x=metric+'_1', y=metric+'_2')
                if not charge_plot:
                    xy_min = [min(g1.get_xlim()), min(g1.get_ylim())]
                    xy_max = [max(g1.get_xlim()), max(g1.get_ylim())]
                    ax.plot([max(xy_min), min(xy_max)], [max(xy_min), min(xy_max)], c='r')

                # Set axis/title text:
                g1.set_title('Correlation between {} and {}'.format(p1, p2))
                # g1.grid(True, axis='y')

                if charge_plot is True:
                    g1.set_xlabel('Charge (%), {}'.format(p1));
                    g1.set_ylabel('Charge (%), {}'.format(p2));
                elif log is True:
                    g1.set_xlabel('Abundance (log10 RPM), {}'.format(p1));
                    g1.set_ylabel('Abundance (log10 RPM), {}'.format(p2));
                else:
                    g1.set_xlabel('Abundance (RPM), {}'.format(p1));
                    g1.set_ylabel('Abundance (RPM), {}'.format(p2));

                # Write to PDF file:
                fig.tight_layout()
                pp.savefig(fig, bbox_inches='tight')
                plt.close(fig)

    def plot_non_temp(self, end='5p', plot_name='_5p-non-template_logo', n_jobs=4, seq_len_percentile=99, seq_len=None, _3p_cover=False, verbose=True):
        if end == '3p':
            col_sele = '3p_non-temp'
        elif end == '5p':
            col_sele = '5p_non-temp'
        else:
            raise Exception('end input has to be either 3p and 5p, not: {}'.format(end))
        self.verbose = verbose
        if verbose:
            print('\nNow collecting data for sample:', end='')

        # Check data exists:
        for _, row in self.sample_df.iterrows():
            stats_fnam = '{}/{}_stats.csv.bz2'.format(self.stats_dir_abs, row['sample_name_unique'])
            assert(os.path.exists(stats_fnam))

        # Read each stats CSV file and generate its count matrix:
        data = list()
        for index, row in self.sample_df.copy().iterrows():
            row['end'] = end
            row['col_sele'] = col_sele
            row['seq_len'] = seq_len
            row['seq_len_percentile'] = seq_len_percentile
            row['_3p_cover'] = _3p_cover
            data.append((index, row))
        # Collect data for plotting for each sample:
        with WorkerPool(n_jobs=n_jobs) as pool:
            results = pool.map(self.__collect_no_temp_mat, data)

        # Make logo:
        self.__run_logomaker(results, plot_name)
        # Also make it as log10 transformed:
        self.__run_logomaker(results, plot_name, log10=True)

    def __collect_no_temp_mat(self, index, row):
        if self.verbose:
            print('  {}'.format(row['sample_name_unique']), end='')
        stats_fnam = '{}/{}_stats.csv.bz2'.format(self.stats_dir_abs, row['sample_name_unique'])
        with bz2.open(stats_fnam, 'rt', encoding="utf-8") as stats_fh:
            sample_stats = pd.read_csv(stats_fh, keep_default_na=False, dtype=self.stats_csv_header_td)

        col_sele = row['col_sele']
        if row['seq_len'] is None:
            # Limit the length of extracted 5p non-template bases to "seq_len_percentile":
            seq_len = np.percentile([len(s) for s in sample_stats[col_sele]], row['seq_len_percentile'], method='nearest')
        else:
            seq_len = row['seq_len']

        title = 'Sample: {}, logo from {} obs'.format(row['sample_name_unique'], row['N_mapped'])
        mask = np.array([len(s) <= seq_len for s in sample_stats[col_sele].values])
        if row['_3p_cover']:
            mask = mask & (sample_stats['3p_cover'] == True)
        if mask.sum() > 1:
            # The dash "-" sign is not counted in "counts_mat":
            if row['end'] == '3p':
                seq_list = [s + '-'*(seq_len - len(s)) for s in sample_stats.loc[mask, col_sele].values if len(s) <= seq_len]
            else:
                seq_list = ['-'*(seq_len - len(s)) + s for s in sample_stats.loc[mask, col_sele].values if len(s) <= seq_len]
        # Will throw an exception if all gaps:
        try:
            with warnings.catch_warnings(): # ignoring a warning about array assignment
                warnings.simplefilter(action='ignore', category=FutureWarning)
                counts_mat = lm.alignment_to_matrix(seq_list, sample_stats.loc[mask, 'count'].values)
        except:
            counts_mat = False
        return([counts_mat, title])

    def plot_UMI_logo(self, plot_name='UMI_logo', n_jobs=4, verbose=True):
        self.verbose = verbose
        if verbose:
            print('\nNow collecting data for sample:', end='')
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
        if self.verbose:
            print('  {}'.format(row['sample_name_unique']), end='')
        stats_fnam = '{}/{}_stats.csv.bz2'.format(self.stats_dir_abs, row['sample_name_unique'])
        with bz2.open(stats_fnam, 'rt', encoding="utf-8") as stats_fh:
            sample_stats = pd.read_csv(stats_fh, keep_default_na=False, dtype=self.stats_csv_header_td)
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
        if self.verbose:
            print('\nNow plotting logo plot.', end='')
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

    def plot_coverage(self, compartment='cyto', plot_type='needle', \
        aa_norm=False, y_norm=False, plot_name='cov_plot_cyto_needle', \
        sample_list=None, verbose=True, n_jobs=4):
        # Check input:
        if compartment != 'mito' and compartment != 'cyto':
            raise Exception('Unknown compartment specified: {}\nValid strings are either either "mito" or "cyto"'.format(compartment))
        if plot_type != 'behrens' and plot_type != 'needle':
            raise Exception('Unknown plot type specified: {}\nValid strings are either either "behrens" or "needle".'.format(plot_type))
        if sample_list is None:
            sample_list = [row['sample_name_unique'] for _, row in self.sample_df.iterrows()]
        sample_list = set(sample_list)

        # Columns used for data aggregation:
        self.aa_cov_cols = ['sample_name_unique', 'tRNA_annotation_len', 'align_5p_idx', 'align_3p_idx', 'AA_letter', 'count']
        # Check data exists:
        for _, row in self.sample_df.iterrows():
            stats_fnam = '{}/{}_stats.csv.bz2'.format(self.stats_dir_abs, row['sample_name_unique'])
            assert(os.path.exists(stats_fnam))
        self.verbose = verbose
        if verbose:
            print('\nNow collecting data for sample:', end='')

        # Read each stats CSV file and generate its coverage matrix:
        data = list()
        for index, row in self.sample_df.copy().iterrows():
            if not row['sample_name_unique'] in sample_list:
                continue
            row['compartment'] = compartment
            row['aa_norm'] = aa_norm
            row['y_norm'] = y_norm
            data.append((index, row))
        # Collect data for plotting for each sample:
        with WorkerPool(n_jobs=n_jobs) as pool:
            results = pool.map(self.__collect_coverage_data, data)
        if verbose:
            print('\nNow plotting sample:', end='')

        # Use a 20 color colormap:
        cmap = mpl.colormaps['tab20']
        # Print each plot to the same PDF file:
        cov_fnam = '{}/{}.pdf'.format(self.plotting_dir_abs, plot_name)
        with PdfPages(cov_fnam) as pp:
            # Loop through and generate plots for each sample:
            for res in results:
                cov_count_sum, aa_ordered_list, title_info = res
                if verbose:
                    print('  {}'.format(title_info[0]), end='')

                try: # Catch odd errors, like empty dataframe etc. 
                    title = 'Sample: {}, coverage from {} obs'\
                    .format(title_info[0], title_info[1])
                    
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

    def __collect_coverage_data(self, index, row):
        print('  {}'.format(row['sample_name_unique']), end='')
        # Read data and add necessary columns:
        stats_fnam = '{}/{}_stats.csv.bz2'.format(self.stats_dir_abs, row['sample_name_unique'])
        with bz2.open(stats_fnam, 'rt', encoding="utf-8") as stats_fh:
            sample_stats = pd.read_csv(stats_fh, keep_default_na=False, dtype=self.stats_csv_header_td)
        # Create single codon filter, to filter out sequences that map to
        # tRNA sequences with different codon/anticodon:
        single_codon = list()
        for anno_str, anticodon in zip(sample_stats['tRNA_annotation'].values, sample_stats['anticodon'].values):
            sc = True
            anno_list = anno_str.split('@')
            for anno in anno_list:
                if not anno.split('-')[2] == anticodon:
                    sc = False
            single_codon.append(sc)
        sample_stats['single_codon'] = single_codon
        sample_stats['Ecoli_ctr'] = ['Escherichia_coli' in anno for anno in sample_stats['tRNA_annotation'].values]
        sample_stats['mito_codon'] = ['mito_tRNA' in anno for anno in sample_stats['tRNA_annotation'].values]
        # Translate codon into single letter amino acid:
        ctab_cyto = copy.deepcopy(Bio.Data.CodonTable.generic_by_id[1].forward_table)
        # Selenocysteine has to be added:
        ctab_cyto['UGA'] = 'SeC'
        ctab_mito = copy.deepcopy(Bio.Data.CodonTable.generic_by_id[2].forward_table)
        ctab_bac = copy.deepcopy(Bio.Data.CodonTable.generic_by_id[11].forward_table)
        aa_letters = list()
        for ecoli, mito, codon in zip(sample_stats['Ecoli_ctr'], sample_stats['mito_codon'], sample_stats['codon'].values):
            if mito:
                aa_letters.append(ctab_mito[codon])
            elif ecoli:
                aa_letters.append(ctab_bac[codon])
            else:
                aa_letters.append(ctab_cyto[codon])
        sample_stats['AA_letter'] = aa_letters

        # Choose rows from requested compartment:
        type_mask = (sample_stats['3p_cover'] == True) & (sample_stats['single_codon']) & (~sample_stats['Ecoli_ctr']) & (sample_stats['AA_letter'].apply(len) == 1)
        if row['compartment'] == 'mito':
            type_mask &= (sample_stats['mito_codon'])
        else:
            type_mask &= (~sample_stats['mito_codon'])

        # Generate dataframe with coverage information.
        # The coverage can be calculated from the count and align_5p_idx.
        cov_df = sample_stats.loc[type_mask, self.aa_cov_cols].copy()
        cov_df = cov_df.groupby(self.aa_cov_cols[:-1], as_index=False).agg({"count": "sum"}).reset_index(drop=True)

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
        for aa_let, _5pi, alen, count in zip(cov_df['AA_letter'], cov_df['align_5p_idx'], cov_df['tRNA_annotation_len'], cov_df['count']):
            aa_idx = aa_order[aa_let]
            _5p_idx = _5pi - 1 # shift alignment index to 0 indexing
            _5p_idx_trans = len_map_len[alen][_5p_idx]
            # Amino acids have multiple codons, each with multiple
            # transcripts that can vary in length, thus add instead of assign:
            cov_count[aa_idx, _5p_idx_trans] += count

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
        if row['aa_norm']:
            # Weigh each amino acid:
            cov_count = cov_count.T / cov_count[:, -1]
            # Normalize so coverage at 3p is 1:
            cov_count = cov_count.T * (1/cov_count.shape[1])
        # If "y_norm" is requested the y-axis is normalized to 1
        elif row['y_norm']:
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
        title_info = (row['sample_name_unique'], cov_df['count'].sum())
        return([cov_count_sum, aa_ordered_list, title_info])




