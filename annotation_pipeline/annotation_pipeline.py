import pandas as pd
import helpers as hp

if __name__ == "__main__":
    #in/out: line separated file of TS559 positions
    #index_path = './supporting_data/unit_test_positions.lsv'
    index_path = './supporting_data/m5C_positions.lsv'

    #Name if output file, default: annotation_ready.tsv
    outfile_name = 'annotation_ready.tsv'

    #Supporting files to help with annotations (auxillary files = 7)
    TS559_genome_path = './supporting_data/TS559_reference_genome.fasta'
    KOD1_genome_path = './supporting_data/KOD1_reference_genome.fasta'
    annotation_path = './supporting_data/TS559_annotations.gtf'
    UTR5_path = './supporting_data/UTR5_TS559.gtf'
    UTR3_path = './supporting_data/UTR3_TS559.gtf'
    TSS_TS559_positions_path ='./supporting_data/TSSs_TS559.tsv'
    RNA_folds_path = "./supporting_data/Tko_allgenes_85C_fold.tsv"

    final_dict = {"TS559_position": [],
                  "KOD1_position": [],
                  "strand" : [],
                  "logos_sequence_41bp" : [],
                  "element_name" : [],
                  "element_type" : [],
                  "element_description": [],
                  "element_length" : [],
                  "position_in_transcript" : [],
                  "percent_position_in_transcript" : [],
                  "positional_enrichment":[],
                  "codon_position" : [],
                  "amino_acid_sequence": [],
                  "amino_acid_ID" : [],
                  "local_41bp_predicted_fold": [],
                  "m5C_position_fold":[],
                  "MFE": [],
                  "associated_TSS_id": [],
                  "TSS_direction": [],
                  "TSS_description": [],
                  "total_annotations": [],
                  "alternate_annotations": []
                  }

    #Load reference genomes into memory
    TS559_genome_sequence = hp.load_genome(TS559_genome_path)
    KOD1_genome_sequence = hp.load_genome(KOD1_genome_path)

    #Initiate dictionary for gene lengths
    gene_lens = {}
    gene_start = {}
    gene_stop = {}

    #Create dataframe with annotation files
    annotation_cols = ['chromosome', 'method', 'element', 'start', 'stop', 'dot', 'strand','zero', 'descriptors']
    annotation_df = pd.read_csv(annotation_path, sep='\t', header = None, names=annotation_cols)

    UTR5_df = pd.read_csv(UTR5_path, sep='\t', header = None, names=annotation_cols)
    UTR3_df = pd.read_csv(UTR3_path, sep='\t', header = None, names=annotation_cols)
    TSS_df = pd.read_csv(TSS_TS559_positions_path, sep='\t', header = None,
                         names=["promoter_id","promoter_start","promoter_downstream",
                                "promoter_direction","promoter_description"])
    rna_folds_df = pd.read_csv(RNA_folds_path, sep='\t', header = None,
                         names=["gene_id","MFE","RNA_predicted_fold"])

    #populate gene_lens dictionary
    for index, row in annotation_df.iterrows():
        element_type = row['element']
        descriptions = hp.get_descript(row['descriptors'])
        element_name = descriptions['gene_id']
        element_length = abs(row['stop'] - row['start'])
        gene_lens[element_name] = element_length
        gene_start[element_name] = row['start']
        gene_stop[element_name] = row['stop']

    with open(index_path, 'r') as index:
        for position in index:
            position = int(position.rstrip())  #make position 0-based indexing
            reference_nucelotide = TS559_genome_sequence[position]
            reference_strand = hp.ref_strand(reference_nucelotide)

            """kmer41 includes 40bp of adjacent sequence plus the target cytosine"""
            kmer41_start = position-20
            if kmer41_start < 1:
                kmer41_start = 1

            kmer41_stop = position + 21
            if kmer41_stop > len(TS559_genome_sequence)-1:
                kmer41_stop = len(TS559_genome_sequence)+1

            kmer41 = TS559_genome_sequence[kmer41_start:kmer41_stop].upper().replace("T", "U")

            if reference_strand == '-':
                kmer41 = hp.rev_comp(kmer41, type='RNA')

            #get KOD1 position
            query = TS559_genome_sequence[position:position+40]
            KOD1_position = hp.get_KOD1_position(query, KOD1_genome_sequence)
            if ';' in KOD1_position:
                query = TS559_genome_sequence[position:position + 60]
                KOD1_position = hp.get_KOD1_position(query, KOD1_genome_sequence)

            final_dict['TS559_position'].append(position)
            final_dict['KOD1_position'].append(KOD1_position)
            final_dict['strand'].append(reference_strand)
            final_dict['logos_sequence_41bp'].append(kmer41)

            element_name = "orphan"
            element_type = "ncRNA"
            element_description = "orphan"
            element_length = "."
            position_in_transcript = "."
            percent_position_in_transcript = '.'
            corrected_position_in_transcript = '.'
            codon_position = '.'
            amino_acid_sequence = '.'
            amino_acid_ID = '.'
            annotations = 0
            alternate_annotations = []

            for index, row in annotation_df.iterrows():
                if row['start'] <= position <= row['stop'] and reference_strand == row['strand']:
                    annotations +=1
                    descriptions = hp.get_descript(row['descriptors'])
                    element_length = abs(row['stop'] - row['start'] + 1)
                    element_name = descriptions['gene_id']
                    element_type = row['element']
                    element_description = row['element']
                    alternate_annotations.append(descriptions['Name'])

                    if 'product' in descriptions.keys():
                        element_description = descriptions['product']

                    position_in_transcript = position - row['start'] + 1


                    if reference_strand == '-':
                        position_in_transcript = row['stop'] - position + 1

                    percent_position_in_transcript = round((position_in_transcript/element_length)*100,1)

                    #corrected position (enrichment position) is converted to a number between 1-2.
                    corrected_position_in_transcript = round(1 + (position_in_transcript/element_length),2)

                    codon_position = hp.codon_pos(position_in_transcript)
                    if reference_nucelotide == 'G':
                        amino_acid_sequence = hp.codon_seq(codon_position, hp.rev_comp(seq=TS559_genome_sequence[position-2:position+3], type='RNA'))
                    else:
                        amino_acid_sequence = hp.codon_seq(codon_position, TS559_genome_sequence[position - 2:position + 3])

                    amino_acid_ID = hp.aa_id(amino_acid_sequence)

                    if element_type != "mRNA":
                        codon_position = '.'
                        amino_acid_sequence = '.'
                        amino_acid_ID = '.'

            if element_description == 'orphan':
                for index, row in UTR5_df.iterrows():
                    if row['start'] <= position <= row['stop'] and reference_strand == row['strand']:
                        annotations += 1
                        descriptions = hp.get_descript(row['descriptors'])
                        cor_gene = descriptions['gene_id']
                        cor_gene_len = gene_lens[cor_gene]
                        element_length = abs(row['stop'] - row['start'])

                        if row['strand'] == '+':
                            position_in_transcript = -1*(gene_start[cor_gene] - position)
                            percent_position_in_transcript = round((position_in_transcript / gene_lens[cor_gene]) * 100, 2)
                            # corrected position (enrichment position) is converted to a number between 0-1
                            corrected_position_in_transcript = round(((element_length+position_in_transcript) / element_length),2)

                        elif row['strand'] == '-':
                            position_in_transcript = -1*(position - gene_stop[cor_gene])
                            percent_position_in_transcript = round((position_in_transcript / gene_lens[cor_gene]) * 100,2)
                            # corrected position (enrichment position) is converted to a number between 0-1
                            corrected_position_in_transcript = round(((element_length+position_in_transcript) / element_length),2)

                        element_name = descriptions['Name']
                        element_type = row['element']
                        if element_type == 'UTR' or element_type == 'likely_UTR':
                            element_type = "5'UTR"
                        element_description =  "5-" + row['element']
                        alternate_annotations.append(descriptions['Name'])

            if element_description == 'orphan':
                for index, row in UTR3_df.iterrows():
                    if row['start'] <= position <= row['stop'] and reference_strand == row['strand']:
                        annotations += 1
                        descriptions = hp.get_descript(row['descriptors'])

                        cor_gene = descriptions['gene_id']
                        cor_gene_len = gene_lens[cor_gene]
                        element_length = abs(row['stop'] - row['start'])

                        if row['strand'] == '+':
                            position_in_transcript = position - gene_start[cor_gene]
                            percent_position_in_transcript = round((position_in_transcript / gene_lens[cor_gene]) * 100,2)
                            # corrected position (enrichment position) is converted to a number between 2-3
                            corrected_position_in_transcript = 2 + round((position - row['start']) / element_length,2)

                        elif row['strand'] == '-':
                            position_in_transcript = gene_stop[cor_gene] - position
                            percent_position_in_transcript = round((position_in_transcript / gene_lens[cor_gene]) * 100,2)
                            # corrected position (enrichment position) is converted to a number between 2-3
                            corrected_position_in_transcript = 2 + round((row['stop']-position) / element_length,2)

                        element_name = descriptions['Name']
                        element_type = row['element']
                        if element_type == 'UTR' or element_type == 'likely_UTR':
                            element_type = "3'UTR"
                        element_description = "3-" + row['element']
                        alternate_annotations.append(descriptions['Name'])

            if element_description == 'orphan':
                for index, row in annotation_df.iterrows():
                    if row['start'] <= position <= row['stop'] and reference_strand != row['strand']:
                        annotations +=1
                        descriptions = hp.get_descript(row['descriptors'])
                        element_length = abs(row['stop'] - row['start'])
                        element_name = "antisense-" + descriptions['gene_id']
                        element_type = "antisense-" + row['element']
                        element_description = "antisense-" + row['element']
                        alternate_annotations.append(descriptions['gene_id'])
                        if 'product' in descriptions.keys():
                            element_description = descriptions['product']
                        break

            #place holders
            m5C_position_fold = 'NA'
            MFE = '.'
            local_41bp_fold = '.'

            for index, row in rna_folds_df.iterrows():
                if row['gene_id'] == element_name and 0 < percent_position_in_transcript <= 100:
                    RNA_predicted_fold = row['RNA_predicted_fold']
                    MFE = row['MFE']

                    m5C_site_fold = RNA_predicted_fold[position_in_transcript-1]
                    if m5C_site_fold == '|':
                        m5C_position_fold = "base_paired"
                    elif m5C_site_fold == '.':
                        m5C_position_fold = "single_stranded"
                    else:
                        m5C_position_fold = "NA"

                    position_in_transcript_lower = position_in_transcript - 21
                    position_in_transcript_upper = position_in_transcript + 20

                    if position_in_transcript-20 < 0:
                        position_in_transcript_lower = 0
                        out_of_bounds_5 = "-"*((position_in_transcript-21)*-1)
                        local_41bp_fold = out_of_bounds_5 + \
                                          RNA_predicted_fold[position_in_transcript_lower:position_in_transcript_upper]

                    elif position_in_transcript + 20 > element_length:
                        position_in_transcript_upper = element_length + 1
                        out_of_bounds_3 = "-" * (20 - (element_length - position_in_transcript))
                        local_41bp_fold = RNA_predicted_fold[position_in_transcript_lower:position_in_transcript_upper] + \
                                          out_of_bounds_3
                    else:
                        local_41bp_fold = RNA_predicted_fold[position_in_transcript_lower:position_in_transcript_upper]

                    break

            promoter_id = '.'
            promoter_direction = '.'
            promoter_description = '.'

            for index, row in TSS_df.iterrows():
                if row['promoter_start'] <= position <= row['promoter_downstream'] and reference_strand == row['promoter_direction']:
                    promoter_id = row['promoter_id']
                    promoter_direction = row['promoter_direction']
                    promoter_description = row['promoter_description']
                    break

            if len(alternate_annotations) <= 1:
                alternate_annotations = "NA"
            else:
                alternate_annotations = hp.ls_to_csv(alternate_annotations)

            final_dict['element_name'].append(element_name)
            final_dict['element_description'].append(element_description)
            final_dict['element_type'].append(element_type)
            final_dict['element_length'].append(element_length)
            final_dict['position_in_transcript'].append(position_in_transcript)
            final_dict['percent_position_in_transcript'].append(percent_position_in_transcript)
            final_dict['positional_enrichment'].append(corrected_position_in_transcript)
            final_dict['codon_position'].append(codon_position)
            final_dict['amino_acid_sequence'].append(amino_acid_sequence)
            final_dict['amino_acid_ID'].append(amino_acid_ID)
            final_dict['local_41bp_predicted_fold'].append(local_41bp_fold)
            final_dict['m5C_position_fold'].append(m5C_position_fold)
            final_dict['MFE'].append(MFE)
            final_dict['total_annotations'].append(annotations)
            final_dict['alternate_annotations'].append(alternate_annotations)
            final_dict['associated_TSS_id'].append(promoter_id)
            final_dict['TSS_direction'].append(promoter_direction)
            final_dict['TSS_description'].append(promoter_description)


    df = pd.DataFrame(final_dict)
    print(df)


    df.to_csv(outfile_name, sep = '\t', index=False)