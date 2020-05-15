from Bio.Blast.Applications import NcbiblastxCommandline
#help(NcbiblastxCommandline)

blastx_cline = NcbiblastxCommandline(query='genes_e_coli_new.fa', db='protein_fragments.fa',
                              evalue=0.001, outfmt=5, out='result_blast.xml')
blastx_cline
print(blastx_cline)
stdout, stderr = blastx_cline()
#b_cline = NcbiblastxCommandline(cmd='blastx', out='result_blast.xml', outfmt=5, query='genes_e_coli_new.fa',
#db='protein_fragments.fa', evalue=0.001)
#print(b_cline)

#stdout, stderr = blastx_cline()