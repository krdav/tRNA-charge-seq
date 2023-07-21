#! /usr/bin/env python3

import sys
from Bio import Seq, SeqIO

cyto_fnam = sys.argv[1]

db_id_set = set()
seq_set = set()
for record in SeqIO.parse(cyto_fnam, "fasta"):
    '''>Homo_sapiens_tRNA-Ala-AGC-1-1 (tRNAscan-SE ID: chr6.trna112) Ala (AGC) 72 bp mature sequence Sc: 84.9 chr6:28763741-28763812 (-)'''
    '''>Homo_sapiens_tRNA-Ala-AGC-1-1'''
    ID = record.id.split(' ')[0]
    ID = ID.replace('Escherichia_coli_str_K-12_substr', 'Escherichia_coli_str_K12_substr')
    assert(ID not in db_id_set)
    seq = str(record.seq).upper().replace('U', 'T')
    if seq in seq_set or '-Und-NNN-' in ID:
        continue
    db_id_set.add(ID)
    seq_set.add(seq)
    print('>{}'.format(ID))
    if 'His' in ID:
        print('G{}'.format(str(seq)))
    else:
        print('{}'.format(str(seq)))




