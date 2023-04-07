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
    assert(ID not in db_id_set)
    if record.seq in seq_set:
        continue
    db_id_set.add(ID)
    seq_set.add(record.seq)
    seq = str(record.seq).replace('U', 'T')
    print('>{}'.format(ID))
    print('{}'.format(str(seq)))




