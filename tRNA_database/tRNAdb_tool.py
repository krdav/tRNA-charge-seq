#! /usr/bin/env python3

import sys
from Bio import Seq, SeqIO

mito_fnam = sys.argv[1]
cyto_fnam = sys.argv[2]
try:
    tRNAscan_out = sys.argv[3]
    tRNAscan_rename = sys.argv[4]
except:
    tRNAscan_out = None
    tRNAscan_rename = None

eColiLys = '''>Escherichia_coli_str_K_12_substr_MG1655_tRNA-eColiLys-TTT-1-1
GGGTCGTTAGCTCAGTTGGTAGAGCAGTTGACTTTTAATCAATTGGTCGCAGGTTCGAATCCTGCACGACCCACCA'''

eColiThr = '''>Escherichia_coli_str_K_12_substr_MG1655_tRNA-eColiThr-CGT-1-1
GCCGATATAGCTCAGTTGGTAGAGCAGCGCATTCGTAATGCGAAGGTCGTAGGTTCGACTCCTATTATCGGCACCA'''



name_dict = dict()
name_dict_rv = dict()
if not tRNAscan_rename is None:
    with open(tRNAscan_rename, 'r') as fh_in:
        fh_in.readline()
        for l in fh_in:
            ch, db = l.split()
            name_dict[db] = ch
            name_dict_rv[ch] = db


# Code snippet from Behrens et al. mim-tRNAseq #
intron_dict = dict()
intron_count = 0
if not tRNAscan_out is None:
    with open(tRNAscan_out, 'r') as fh_in:
        for l in fh_in:
            if not l.startswith("chr"):
                continue
            tRNA_ID = '{}.trna{}'.format(*l.split()[0:2])
            if not tRNA_ID in name_dict_rv:
                print(tRNA_ID)
            assert(tRNA_ID in name_dict_rv)
            tRNA_start = int(l.split()[2])
            intron_start = int(l.split()[6])
            intron_stop = int(l.split()[7])
            # if inton boundaries are not 0, i.e. there is an intron then add to dict
            if (intron_start > 0) & (intron_stop > 0):
                if tRNA_start > intron_start: # tRNA is on reverse strand
                    intron_count += 1
                    intron_start = tRNA_start - intron_start
                    intron_stop = tRNA_start - intron_stop + 1
                else: # tRNA is on forward strand
                    intron_count += 1
                    intron_start -= tRNA_start
                    intron_stop -= tRNA_start
                    intron_stop += 1

                intron_dict[tRNA_ID] = (intron_start, intron_stop)



print(eColiLys)
print(eColiThr)

db_id_set = set()
seq_set = set()
for record in SeqIO.parse(mito_fnam, "fasta"):
    assert(record.id not in db_id_set)
    assert(record.seq not in seq_set)
    db_id_set.add(record.id)
    seq = str(record.seq).upper()
    seq_set.add(seq)
    ID, sp, numb, AA, AC = record.id.split('|')
    '''mtdbD00000547|Homo_sapiens|9606|Thr|TGT'''
    '''>Homo_sapiens_mito_tRNA-Leu-TAA-1-1'''
    print('>{}_mito_tRNA-{}-{}'.format(sp, AA, AC))
    # Histidine tRNA post translational modification.
    # See Heinemann et al. 2012.
    if AA == 'His':
        print('G{}CCA'.format(seq))
    else:
        print('{}CCA'.format(record.seq))


intron_found = 0
db_id_set = set()
seq_set = set()
for record in SeqIO.parse(cyto_fnam, "fasta"):
    '''>Homo_sapiens_tRNA-Ala-AGC-1-1 (tRNAscan-SE ID: chr6.trna112) Ala (AGC) 72 bp mature sequence Sc: 84.9 chr6:28763741-28763812 (-)'''
    '''>Homo_sapiens_tRNA-Ala-AGC-1-1'''
    ID = record.id.split(' ')[0]
    assert(ID not in db_id_set)
    seq = str(record.seq).upper().replace('U', 'T')
    if seq in seq_set or '-Und-NNN-' in ID:
        continue
    
    db_id_set.add(ID)
    seq_set.add(seq)
    print('>{}'.format(ID))
    
    # Cut out intro if any:
    tr_anno = ID.split('_')[2]
    if tr_anno in name_dict and name_dict[tr_anno] in intron_dict:
        intron_start, intron_stop = intron_dict[name_dict[tr_anno]]
        seq = seq[:intron_start] + seq[intron_stop:]
        intron_found += 1
    if 'His' in ID:
        print('G{}CCA'.format(str(seq)))
    else:
        print('{}CCA'.format(str(seq)))


