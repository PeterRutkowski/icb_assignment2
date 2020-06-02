import numpy as np
from Bio import SeqIO
from scipy.stats import binom_test
import pandas as pd

def convert_to_int(text):
    for i in range(len(text)):
        if text[0] == ' ':
            text = text[1:]
    return int(text)
def parse_line(line):
    return [convert_to_int(line[1:6]), convert_to_int(line[8:13]),
            convert_to_int(line[15:20]), convert_to_int(line[22:27])]

def log_odds_matrices(group):
    load_matrix = False
    line_counter = 0
    log_odds_matrix = []
    with open('data/meme_%s.txt'%(group), 'r') as file:
            for line in file:
                if load_matrix:
                    line_counter += 1
                    log_odds_line.append(parse_line(line))

                    if line_counter == 15:
                        line_counter = 0
                        load_matrix = False
                        log_odds_matrix.append(log_odds_line)
                if line.startswith('log-odds matrix'):
                    load_matrix = True
                    log_odds_line = []
    return log_odds_matrix

log_odds_A = log_odds_matrices('A')
log_odds_B = log_odds_matrices('B')

proms_A = []
for record in SeqIO.parse('data/proms_A.fa', 'fasta'):
    proms_A.append(record.seq)
proms_B = []
for record in SeqIO.parse('data/proms_B.fa', 'fasta'):
    proms_B.append(record.seq)

def log_odds_score(sequence, log_odds_matrix):
    score = 0
    for i in range(len(sequence)):
        if sequence[i] == 'A':
            score += log_odds_matrix[i][0]
        elif sequence[i] == 'C':
            score += log_odds_matrix[i][1]
        elif sequence[i] == 'G':
            score += log_odds_matrix[i][2]
        else:
            score += log_odds_matrix[i][3]
    return score

def log_odds_hits(sequence, log_odds_matrix):
    sequence = str(sequence)
    hits = 0
    for i in range(0,(len(sequence)-np.shape(log_odds_matrix)[0])):
        if log_odds_score(sequence[i:(i+np.shape(log_odds_matrix)[0])], log_odds_matrix) > 0:
            hits += 1
    return hits


motifs_enrichment = []

for i in range(len(log_odds_A)):
    sum_A = 0
    sum_B = 0
    for j in range(len(proms_A)):
        sum_A += log_odds_hits(proms_A[j], log_odds_A[i])
    for j in range(len(proms_B)):
        sum_B += log_odds_hits(proms_B[j], log_odds_A[i])

    motifs_enrichment.append(['A%d' % (i), sum_A, sum_B,
                              binom_test(sum_A, n=86 * len(proms_A), p=sum_A / (86 * len(proms_A))),
                              binom_test(sum_B, n=86 * len(proms_B), p=sum_A / (86 * len(proms_A)))])
for i in range(len(log_odds_B)):
    sum_A = 0
    sum_B = 0
    for j in range(len(proms_A)):
        sum_A += log_odds_hits(proms_A[j], log_odds_B[i])
    for j in range(len(proms_B)):
        sum_B += log_odds_hits(proms_B[j], log_odds_B[i])

    motifs_enrichment.append(['B%d' % (i), sum_A, sum_B,
                              binom_test(sum_A, n=86 * len(proms_A), p=sum_B / (86 * len(proms_B))),
                              binom_test(sum_B, n=86 * len(proms_B), p=sum_B / (86 * len(proms_B)))])

df = pd.DataFrame(motifs_enrichment, columns=['motif_id','hits_A','hits_B','evalue_A','evalue_B'])
df.to_csv('results/enrichments.csv', index=False)