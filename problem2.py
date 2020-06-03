from Bio import motifs
from Bio.Seq import Seq

# ASSUMPTION:
# we already have MEME reports generated by MEME Suite (http://meme-suite.org)
# they are saved as data/meme_A.txt and data/meme_B.txt
# they correpond to group A and group B promoters respectively

# MEME parser
# returns a list of biopython motifs of given group
# group is a string
# either 'A' or 'B'
def MEME_parser(group):
    group_motifs = []
    build_motif = False
    motif_counter = 0
    add_index = 0
    with open('data/meme_%s.txt'%(group), 'r') as file:
        for line in file:
            if line.startswith('MOTIF'):
                motif_counter += 1
                if motif_counter == 10:
                    add_index = 1
                build_motif = True
                sites = line[51+add_index:53+add_index]
                if sites[1] == ' ':
                    sites = sites[0]
                sites = int(sites)
                counter = 0
                motif = []
            else:
                if build_motif:

                    counter += 1
                    if counter > 32 and counter < 32 + sites + 1:
                        motif.append(Seq(line[59:74]))

                    if counter == 32 + sites:
                        build_motif = False
                        group_motifs.append(motifs.create(motif))
    return group_motifs

# save index-th motif from given group to a pfm file
def save_pfm(motif, group, index):
    with open('results/motifs/%s%s.pfm'%(group,index), 'w') as f:
        f.writelines(motif.format('pfm'))

motifs_A = MEME_parser('A')
index = 0
for motif in motifs_A:
    save_pfm(motif, 'A', index)
    index += 1

motifs_B = MEME_parser('B')
index = 0
for motif in motifs_B:
    save_pfm(motif, 'B', index)
    index += 1