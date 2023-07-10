#!/usr/bin/env python3
from collections import Counter
import helpers as hp


"""
Based on current annotations, calculate how many times each amino acid identity and codon sequence occur in the 
T. kodakarensis TS559 reference genome. 
"""

if __name__ == "__main__":

    TS559_genome_path = './supporting_data/TS559_reference_genome.fasta'

    TS559_genome_sequence = hp.load_genome(TS559_genome_path)
    annotation_path = './supporting_data/TS559_annotations.gtf'

    with open(annotation_path, "r") as f:
        all_aa_id = []
        all_aa_codons = []
        for row in f:
            chrom, method, element, start, stop, dot, strand, zero, descriptors = row.split('\t')
            if element == 'mRNA' and int(start) != 0 and int(stop) != 0 :
                desc = hp.get_descript(descriptors)
                gene_seq = TS559_genome_sequence[int(start)-1:int(stop)]
                if len(gene_seq)%3 != 0:
                    continue
                if strand == '-':
                    gene_seq = hp.rev_comp(gene_seq, type = 'DNA')
                aa_seq = hp.get_aaseq(gene_seq)
                aa_codons = hp.get_aacodon(gene_seq)
                for i in aa_seq:
                    all_aa_id.append(i)
                for i in aa_codons:
                    all_aa_codons.append(i)


    # print codon usage
    for k,v in Counter(all_aa_codons).items():
        print(k,v)


    # print amino acid usage
    for k,v in Counter(all_aa_id).items():
        print(k,v)

