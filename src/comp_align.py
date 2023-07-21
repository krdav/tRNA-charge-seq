import numpy as np
import json
import pandas as pd
import bz2


def compare_SWalign(sample_list, data_dir, \
                    aln_dir1='SWalign', aln_dir2='SWalign_masked', \
                    verbose=True):
    if verbose:
        print('Now processing sample:', end='')
    # Count read annotations:
    anno2anno = dict()
    for snu in sample_list:
        if verbose:
            print(' {} '.format(snu), end='')
        try:
            SWres_fnam1 = '{}/{}/{}_SWalign.json.bz2'.format(data_dir, aln_dir1, snu)
            SWres_fnam2 = '{}/{}/{}_SWalign.json.bz2'.format(data_dir, aln_dir2, snu)
            with bz2.open(SWres_fnam1, 'rt') as fh_in1:
                json_aln1 = json.load(fh_in1)
            with bz2.open(SWres_fnam2, 'rt') as fh_in2:
                json_aln2 = json.load(fh_in2)
        except Exception as e:
            print('Skipping: {}'.format(snu))
            print(e)
            continue

        for id1 in json_aln1:
            try:
                anno1 = json_aln1[id1]['name']
            except:
                continue
            try:
                anno2 = json_aln2[id1]['name']
            except:
                continue

            try:
                anno2anno[anno1][anno2] += 1
            except:
                try:
                    anno2anno[anno1][anno2] = 1
                except:
                    anno2anno[anno1] = dict()
                    anno2anno[anno1][anno2] = 1

    # Collect data for dataframe:
    ch_anno = list()
    for anno1 in anno2anno:
        for anno2 in anno2anno[anno1]:
            # Extract for annotation 1:
            anno1_inf = [[], []]
            for sa1 in anno1.split('@'):
                anno1_inf[0].append('-'.join(sa1.split('_')[-1].split('-')[1:]))
                anno1_inf[1].append('-'.join(sa1.split('_')[-1].split('-')[1:3]))
            # Extract for annotation 2:
            anno2_inf = [[], []]
            for sa2 in anno2.split('@'):
                anno2_inf[0].append('-'.join(sa2.split('_')[-1].split('-')[1:]))
                anno2_inf[1].append('-'.join(sa2.split('_')[-1].split('-')[1:3]))

            ch_anno.append(('@'.join(anno1_inf[0]), \
                            '@'.join(anno2_inf[0]), \
                            '@'.join(sorted(set(anno1_inf[1]))), \
                            '@'.join(sorted(set(anno2_inf[1]))), \
                            anno2anno[anno1][anno2]))

    # Generate dataframe:
    anno_df = pd.DataFrame(ch_anno, columns=['anno1', 'anno2', 'codon1', 'codon2', 'count'])

    # Find the total count for each annotation:
    a1c = anno_df.groupby('anno1').agg({'anno1': 'max', \
                                        'count': 'sum'}).reset_index(drop=True)
    a2c = anno_df.groupby('anno2').agg({'anno2': 'max', \
                                        'count': 'sum'}).reset_index(drop=True)

    c1c = anno_df.groupby('codon1').agg({'codon1': 'max', \
                                         'count': 'sum'}).reset_index(drop=True)
    c2c = anno_df.groupby('codon2').agg({'codon2': 'max', \
                                         'count': 'sum'}).reset_index(drop=True)

    anno_df = anno_df.merge(a1c, on='anno1', suffixes=('', '_a1_tot'), copy=False)
    anno_df = anno_df.merge(a2c, on='anno2', suffixes=('', '_a2_tot'), copy=False)

    anno_df = anno_df.merge(c1c, on='codon1', suffixes=('', '_c1_tot'), copy=False)
    anno_df = anno_df.merge(c2c, on='codon2', suffixes=('', '_c2_tot'), copy=False)

    # Isolate the codon annotations:
    codon_df = anno_df.groupby(['codon1', 'codon2']).agg({'codon1': 'max', \
                                                          'codon2': 'max', \
                                                          'count_c1_tot': 'max', \
                                                          'count_c2_tot': 'max', \
                                                          'count': 'sum'}).reset_index(drop=True)

    return(anno_df, codon_df)








