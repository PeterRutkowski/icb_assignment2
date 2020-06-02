from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast import NCBIXML
import pandas as pd
import numpy as np

# translate DNA code to amino acid code
def translate(seq):
    seq = str(seq)
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    rest = len(seq) % 3
    seq = seq[:-rest]
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
    return Seq(protein)

# load ecoli genes from fasta file
db = []
for record in SeqIO.parse('data/genes_e_coli_new.fa', 'fasta'):
    db.append(record)

# translate ecoli genes to amino acids
for i in range(len(db)):
    db[i].seq = translate(db[i].seq)

# save ecoli genes in amino acid format
SeqIO.write(db, 'data/db.fa', 'fasta')

# create BLAST database
cline = NcbimakeblastdbCommandline(dbtype="prot", input_file='data/db.fa')
stdout, stderr = cline()

# execute BLASTp algorithm
cline = NcbiblastpCommandline(query='data/protein_fragments.fa', db='data/db.fa',
                              evalue=0.0001, outfmt=5, out='data/blastp.xml')
cline
stdout, stderr = cline()

# parse BLASTp result
result_handle = open('data/blastp.xml')
blast_records = NCBIXML.parse(result_handle)
blast_records = list(blast_records)

# save protein fragment and ecoli labels for convenience
protein_ids = []
for record in SeqIO.parse('data/protein_fragments.fa', 'fasta'):
    protein_ids.append(record.id)
ecoli_ids = []
for record in SeqIO.parse('data/db.fa', 'fasta'):
    ecoli_ids.append(record.id)

# calculate best matches based on blast records
results = []
ind = 0
for blast_record in blast_records:
    alignment = blast_record.alignments[0]
    index = ''
    for i in range(14,len(alignment.title)):
        if alignment.title[i] != ' ':
            index = index + alignment.title[i]
        else:
            break
    index = int(index)
    results.append([protein_ids[ind],ecoli_ids[index],alignment.hsps[0].expect])
    ind += 1

# save best matches as a csv file
df = pd.DataFrame(results, columns=['protein_id','ecoli_id','e_value'])
df.to_csv('results/matches.csv', index=False)

# in preperation for problem 2
# load labels of ecoli genes form both groups A and B
group_ids_A = np.asarray(pd.read_csv('results/matches.csv')['ecoli_id'][:63])
group_ids_B = np.asarray(pd.read_csv('results/matches.csv')['ecoli_id'][63:])

# choose corresponding promoters to saved group labels
proms_A, proms_B = [], []
for record in SeqIO.parse('data/proms_e_coli_fixed.fa', 'fasta'):
    if record.id in group_ids_A:
        proms_A.append(record)
    if record.id in group_ids_B:
        proms_B.append(record)

# save promoters corresponding to groups A and B
SeqIO.write(proms_A, 'data/proms_A.fa', 'fasta')
SeqIO.write(proms_B, 'data/proms_B.fa', 'fasta')