# icb_assignment2
MIM UW 2019/20 Introduction to computational biology Assignment 2

1. Identify the closest related protein-coding sequences in the E. coli genome for each of the input sequences given in the fasta file. This should be done by downloading the E. coli protein-coding sequences, creating a BLAST database and performing a local BLAST search against this database. (Here’s a tutorial on running blast locally.) The result of this step should be a comma-separated list containing: input sequence id, best matching E. coli gene id and the associated e-value.

2. For each of the identified E. coli genes, you should find the associated promoter DNA sequence in the given file. Please note that the input file contained sequences from groups A and B,  we speculate (based on the empirical evidence), that these groups of genes should have different regulatory mechanisms. Identify 10 sequence motifs present in the promoters associated with group A and 10 motifs associated with group B, independently of each other. The result of this step should be two sets of 10 different motif possition-specific matrices in a .pfm format. You can do this with a consensus-like method you implement yourself, or you can use the MEME-suite and process it’s output.

3. Given the two sets of motifs, we would like to select only motifs that are specific to group A or group B. For this, we need to compute the occurrences of the motifs from .pfm in both groups. For this, we can scan the promoter sequences of genes from group A and B, and obtain a number of sequences in these promoters that have a log-odds score higher than 0. Now, given the combined lengths of the promoters in each group and total number of “hits” in each group, we can use the binomial test to identify these motifs that are significantly enriched in one of the groups than in the other. The result of this step should be a list of all 20 motifs with associated, number of hits in the promoters from group A, and promoters from group B and the associated p-values for enrichment in group A and group B respectively.

All of these results and your implementation, should be described in a short report.

Altogether, you should submit your solution in 4 weeks (until June 9th).
