#! /usr/bin/env python3

import sys
from Bio import Seq, SeqIO

mito_fnam = sys.argv[1]
cyto_fnam = sys.argv[2]
# out_fnam = sys.argv[3]

eColiLys = '''>Escherichia_coli_str_K_12_substr_MG1655_tRNA-eColiLys-TTT-1-1
GGGTCGTTAGCTCAGTTGGTAGAGCAGTTGACTTTTAATCAATTGGTCGCAGGTTCGAATCCTGCACGACCCACCA'''

eColiThr = '''>Escherichia_coli_str_K_12_substr_MG1655_tRNA-eColiThr-CGT-1-1
GCCGATATAGCTCAGTTGGTAGAGCAGCGCATTCGTAATGCGAAGGTCGTAGGTTCGACTCCTATTATCGGCACCA'''


print(eColiLys)
print(eColiThr)

db_id_set = set()
seq_set = set()
for record in SeqIO.parse(mito_fnam, "fasta"):
    assert(record.id not in db_id_set)
    assert(record.seq not in seq_set)
    db_id_set.add(record.id)
    seq_set.add(record.seq)
    ID, sp, numb, AA, AC = record.id.split('|')
    '''mtdbD00000547|Homo_sapiens|9606|Thr|TGT'''
    '''>Homo_sapiens_mito_tRNA-Leu-TAA-1-1'''
    print('>{}_mito_tRNA-{}-{}'.format(sp, AA, AC))
    print('{}CCA'.format(str(record.seq)))


db_id_set = set()
seq_set = set()
for record in SeqIO.parse(cyto_fnam, "fasta"):
    '''>Homo_sapiens_tRNA-Ala-AGC-1-1 (tRNAscan-SE ID: chr6.trna112) Ala (AGC) 72 bp mature sequence Sc: 84.9 chr6:28763741-28763812 (-)'''
    '''>Homo_sapiens_tRNA-Ala-AGC-1-1'''
    ID = record.id.split(' ')[0]
    assert(ID not in db_id_set)
    if record.seq in seq_set:
        continue
    db_id_set.add(ID)
    seq_set.add(record.seq)
    seq = str(record.seq).replace('U', 'T')
    print('>{}'.format(ID))
    print('{}CCA'.format(str(seq)))






