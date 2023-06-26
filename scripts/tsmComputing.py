""" :dna: returns the similarity matrix tsm+ or tsm scores for all pairs of homologous transcripts.
> Usage:
======
    python3 tsm-computing.py [-talg] [-gtot] [-tsm]

> Reference:
======
    https://github.com/UdeS-CoBIUS/TranscriptOrthology
"""
__authors__ = ("Wend Yam Donald Davy Ouedraogo")
__contact__ = ("wend.yam.donald.davy.usherbrooke.ca")
__copyright__ = "CoBIUS lab at Université de Sherbrooke, QC, CANADA"
__date__ = "2023-06-26"
__version__= "2.0.6"

import pandas as pd
import argparse
import time
import numpy as np

def build_arg_parser():
    '''Parsing function'''
    parser = argparse.ArgumentParser(description="parsor program parameter")
    parser.add_argument('-talg', '--tralignment', default=None)
    parser.add_argument('-gtot', '--genetotranscripts', default=None)
    parser.add_argument('-tsm', '--tsmvalue', default=False)
    parser.add_argument('-outf', '--outputfolder', default='.')
    return parser

def saving_files(file, outputfolder, name, ind):
    file.to_csv('{}/{}.csv'.format(outputfolder,name), header=True, index=ind, sep=';')
    return True

def get_matrix(transcripts_msa_path, gtot_path, tsm_condition, output_folder_path):
    '''Returns the matrix score'''
    try:
        print('\n\n++++++++++++++++Starting ....')
        start = time.time()
        df_transcripts, df_gtot, df_blocks_transcripts, df_blocks_genes = convert_files_to_dataframes(transcripts_msa_path, gtot_path, output_folder_path)
        print('+++++++ All data were retrieved & the representation of subtranscribed sequences of genes into blocks are available.')
    except:
        raise('Something wrong with the inputs! Please check it out or contact us. Thank you.')
    
    try:
        print('+++++ Computing matrix ...\t in progress')
        tsm_matrix = compute_matrix(df_transcripts, df_gtot, df_blocks_transcripts, df_blocks_genes, tsm_condition)
        #saving matrix
        saving_files(tsm_matrix, output_folder_path, 'matrix', True)
        end = time.time()
        print('+++++ Computing matrix ...\t status: Finished without errors in {} seconds'.format(str(end-start)))
    except:
        raise('The algorithm exits with errors!! Please contact us for help. Thank you.')
    return tsm_matrix, df_blocks_transcripts, df_blocks_genes

def convert_files_to_dataframes(transcripts_msa_path,gtot_path, output_folder_path):
    """Format inputs of users into Pandas DataFrames"""

    df_transcripts = pd.DataFrame(columns=['id_transcript','alg'])
    df_gtot = pd.DataFrame(columns=['id_transcript','id_gene'])

    file_transcripts_msa = open(transcripts_msa_path, 'r')
    transcripts_msa_rows = [str(_).split('\n')[0] for _ in file_transcripts_msa.readlines()]
    file_transcripts_msa.close()
    for i_number, row in enumerate(transcripts_msa_rows):
        if row.startswith('>'):
            df_transcripts.loc[row] = [row.split('>')[-1], transcripts_msa_rows[i_number+1]]
        
    file_gtot = open(gtot_path, 'r')
    gtot_rows = [str(_).split('\n')[0] for _ in file_gtot.readlines()]
    file_gtot.close()
    for row_gtot in gtot_rows:
        if row_gtot.startswith('>'):
            df_gtot.loc[row_gtot] = [row_gtot.split(':')[0].split('>')[-1], row_gtot.split(':')[-1]]
    

    df_blocks_transcripts = retrieve_blocks(transcripts_msa_path)

    df_blocks_genes = inferred_genes_blocks(df_blocks_transcripts, df_gtot)
    
    # save files
    saving_files(df_blocks_transcripts, output_folder_path, 'blocks_transcripts', False)
    saving_files(df_blocks_genes, output_folder_path, 'blocks_genes', False)

    return df_transcripts, df_gtot, df_blocks_transcripts, df_blocks_genes

def retrieve_blocks(transcripts_msa_path):
    """Retrieve blocks from transcripts and genes MSA"""

    list_of_ids_transcripts, blocks_data = convert_alignment_to_binary_sequence(transcripts_msa_path)
    transcripts_blocks = search_blocks_in_msa(list_of_ids_transcripts, blocks_data)
    return transcripts_blocks

def convert_alignment_to_binary_sequence(transcripts_msa_path):
    '''Return the binary sequences of an alignment'''
    file = open(transcripts_msa_path, 'r')
    lines = file.readlines()
    file.close()
    list_of_ids_transcripts = []
    blocks_data = []
    for (index,line) in enumerate(lines):
        if line.startswith('>'):
            j_index = index+1
            blocks_sequence = []
            nucleotides_in_sequence = lines[j_index].split('\n')[0]
            for nucleotide in nucleotides_in_sequence:
                if nucleotide == '-':
                    blocks_sequence.append(0)
                else:
                    blocks_sequence.append(1)
            blocks_data.append(blocks_sequence)
            list_of_ids_transcripts.append(line.split('>')[1].split('\n')[0])               
    
    return (list_of_ids_transcripts, blocks_data)

def search_blocks_in_msa(list_of_ids_transcripts, blocks_data):
    """Search blocks"""
    #mappings
    length_column = len(blocks_data[0])
    positions = []
    prec_value_block = ''
    gaps_zero = ''.join(['0' for i in range(len(blocks_data))])
    for i_column in range(length_column):
        current_value = []
        for i_row, row in enumerate(blocks_data):
            row = blocks_data[i_row]
            current_value.append(str(row[i_column]))
        current_block = ('').join(current_value)
        if current_block != gaps_zero:
            if current_block != prec_value_block:
                positions.append(i_column)
                prec_value_block = current_block
            else:
                prec_value_block = current_block
    positions_blocks = {}
    number = 0
    interval_stop_last_block = length_column
    for nt in range(len(blocks_data[0])-1,0,-1):
        last_nts = []
        for alg in blocks_data:
            last_nts.append(alg[nt])
        if 1 in last_nts:
            interval_stop_last_block = nt
            break

    for i in range(len(positions)-1):
        interval_start = positions[i]
        interval_stop = positions[i+1]
        interval = (interval_start, interval_stop)
        positions_blocks[i] = interval
        number = i
    interval_start_last_block = positions[-1]
    positions_blocks[number+1] = (interval_start_last_block, interval_stop_last_block)

    # transcripts
    sequences = {}
    for i in range(len(blocks_data)):
        transcript = list_of_ids_transcripts[i]
        list_sequence = [str(_) for _ in blocks_data[i]]
        sequence = ''.join(list_sequence)
        sequence_blocks = []
        for block in positions_blocks.keys():
            tuple_position = positions_blocks[block]
            start = tuple_position[0]
            stop = tuple_position[1]
            sequence_block = sequence[start:stop]
            boolCheck = False
            for nt in sequence_block:
                if nt == '1':
                    boolCheck = True
            if boolCheck:
                sequence_blocks.append(block)
        sequences[transcript] = sequence_blocks
    df_results = pd.DataFrame(columns=['id_transcript','blocks'])
    i_data = []
    j_data = []
    for transcript in sequences.keys():
        value_data = sequences[transcript]
        i_data.append(transcript)
        j_data.append('|'.join([str(_) for _ in value_data]))
    df_results['id_transcript'] = i_data
    df_results['blocks'] = j_data

    return df_results

def inferred_genes_blocks(transcripts_blocks, df_gtot):
    all_genes_blocks = []
    all_genes_id = []
    genes = list(df_gtot.id_gene.values)
    for gene in genes:
        transcripts = list(df_gtot[df_gtot['id_gene']==gene].id_transcript.values)
        list_of_nblocks = list(transcripts_blocks[transcripts_blocks['id_transcript'].isin(transcripts)].blocks.values)
        set_of_nblocks = [set(str(_).split('|')) for _ in list_of_nblocks]
        union_set_blocks = set()
        if len(set_of_nblocks) > 1:
            for i in range(len(set_of_nblocks)):
                union_set_blocks = union_set_blocks.union(set_of_nblocks[i])
        else:
            union_set_blocks = set_of_nblocks[0]
        list_union_set_blocks = list(union_set_blocks)
        list_union_set_blocks_sorted = [str(_) for _ in sorted([int(_) for _ in list_union_set_blocks])]
        gene_blocks = '|'.join(list_union_set_blocks_sorted)

        all_genes_blocks.append(gene_blocks)
        all_genes_id.append(gene)

    df = pd.DataFrame(data={'id_gene': all_genes_id, 'blocks': all_genes_blocks})
    return df

def compute_matrix(df_transcripts, df_gtot, df_blocks_transcripts, df_blocks_genes, tsm_condition):
    """ 
        ### Compute the tsm+ | tsm  score
        #### inputs :
        ##### arg1 : Pandas DataFrame = {columns=['id_transcript', 'alg']} containing ids of transcripts and their aligned sequences of nucleotides
        ##### arg2 : Pandas DataFrame = {columns=['id_gene', 'alg']} containing ids of genes and their aligned sequences of nucleotides
        ##### arg3 : Pandas DataFrame = {columns=['id_transcript', 'id_gene']} containing ids of transcripts and their corresponding genes (mappings set)
        ##### arg4 : Pandas DataFrame = {columns=['id_transcript', 'blocks']} containing ids of transcripts and their representation into blocks
        ##### arg5 : Pandas DataFrame = {columns=['id_gene', 'blocks']} containing ids of genes and their representation into blocks
        ##### arg6 : Bool : True | False. if True compute the tsm otherwise tsm+.
        --- 
        Returns the matrix similarity score
    """
    transcripts = list(df_transcripts['id_transcript'].values)
    matrix = np.zeros((len(transcripts), len(transcripts)))
    tsm_condition = int(tsm_condition)
    used_indice = {}

    if tsm_condition == 1:
        """ tsm+unitary"""
        for i, transcript_1 in enumerate(transcripts):
            for j, transcript_2 in enumerate(transcripts):
                if transcript_1 == transcript_2:
                    matrix[i,j] = 1.0
                else:
                    indice1 = transcript_1+'&'+transcript_2
                    indice2 = transcript_2+'&'+transcript_1
                    if indice1 not in list(used_indice.keys()) and indice2 not in list(used_indice.keys()):
                        # Compute the nucHOM score (unitary)
                        score_nuc_hom = round(compute_nuchom_score(transcript_1, transcript_2, df_transcripts),3)
                        matrix[i, j] = score_nuc_hom
                        used_indice[indice1] = score_nuc_hom
                        used_indice[indice2] = score_nuc_hom
                    else :
                        matrix[i, j] = used_indice[indice1]
    elif tsm_condition == 2:
        """ tsm+length"""
        for i, transcript_1 in enumerate(transcripts):
            for j, transcript_2 in enumerate(transcripts):
                if transcript_1 == transcript_2:
                    matrix[i,j] = 1.0
                else:
                    indice1 = transcript_1+'&'+transcript_2
                    indice2 = transcript_2+'&'+transcript_1
                    if indice1 not in list(used_indice.keys()) and indice2 not in list(used_indice.keys()):
                        # Compute the nucHOM score (unitary)
                        score_deg_hom = round(compute_deghom_score(transcript_1, transcript_2, df_blocks_transcripts),3)
                        matrix[i, j] = score_deg_hom
                        used_indice[indice1] = score_deg_hom
                        used_indice[indice2] = score_deg_hom
                    else :
                        matrix[i, j] = used_indice[indice1]
    elif tsm_condition == 3:
        """ tsm+mean"""
        for i, transcript_1 in enumerate(transcripts):
            for j, transcript_2 in enumerate(transcripts):
                if transcript_1 == transcript_2:
                    matrix[i,j] = 1.0
                else:
                    indice1 = transcript_1+'&'+transcript_2
                    indice2 = transcript_2+'&'+transcript_1
                    if indice1 not in list(used_indice.keys()) and indice2 not in list(used_indice.keys()):
                        # Compute the nucHOM score (unitary)
                        score_nuc_hom = round(compute_nuchom_score(transcript_1, transcript_2, df_transcripts),3)
                        # Compute the DegHOM score (length)
                        score_deg_hom = round(compute_deghom_score(transcript_1, transcript_2, df_blocks_transcripts),3)
                        # Compute tsm (mean)
                        score = round(((score_nuc_hom+score_deg_hom)/2),3)
                        matrix[i, j] = score
                        used_indice[indice1] = score
                        used_indice[indice2] = score
                    else :
                        matrix[i, j] = used_indice[indice1]
    elif tsm_condition == 4:
        """ tsm++unitary"""
        gene_model = get_gene_model(df_gtot, df_transcripts)
        for i, transcript_1 in enumerate(transcripts):
            for j, transcript_2 in enumerate(transcripts):
                if transcript_1 == transcript_2:
                    matrix[i,j] = 1.0
                else:
                    indice1 = transcript_1+'&'+transcript_2
                    indice2 = transcript_2+'&'+transcript_1
                    if indice1 not in list(used_indice.keys()) and indice2 not in list(used_indice.keys()):
                        gene_tr_1 = df_gtot[df_gtot['id_transcript'] == transcript_1].id_gene.values[0]
                        gene_tr_2 = df_gtot[df_gtot['id_transcript'] == transcript_2].id_gene.values[0]
                        score_nuc_hom_v2 = round(compute_nuchom_score_tsmplus(transcript_1, transcript_2, df_transcripts, gene_tr_1, gene_tr_2, gene_model),3)
                        matrix[i, j] = score_nuc_hom_v2
                        used_indice[indice1] = score_nuc_hom_v2
                        used_indice[indice2] = score_nuc_hom_v2
                    else :
                        matrix[i, j] = used_indice[indice1]
    elif tsm_condition == 5:
        """ tsm++length"""
        for i, transcript_1 in enumerate(transcripts):
            for j, transcript_2 in enumerate(transcripts):
                if transcript_1 == transcript_2:
                    matrix[i,j] = 1.0
                else:
                    indice1 = transcript_1+'&'+transcript_2
                    indice2 = transcript_2+'&'+transcript_1
                    if indice1 not in list(used_indice.keys()) and indice2 not in list(used_indice.keys()):
                        gene_tr_1 = df_gtot[df_gtot['id_transcript'] == transcript_1].id_gene.values[0]
                        gene_tr_2 = df_gtot[df_gtot['id_transcript'] == transcript_2].id_gene.values[0]
                        score_deg_hom_v2 = round(compute_deghom_score_tsmplus(transcript_1, transcript_2, df_blocks_transcripts, df_blocks_genes, gene_tr_1, gene_tr_2),3)
                        matrix[i, j] = score_deg_hom_v2
                        used_indice[indice1] = score_deg_hom_v2
                        used_indice[indice2] = score_deg_hom_v2
                    else :
                        matrix[i, j] = used_indice[indice1]
    elif tsm_condition == 6:
        """ tsm++mean"""
        gene_model = get_gene_model(df_gtot, df_transcripts)
        for i, transcript_1 in enumerate(transcripts):
            for j, transcript_2 in enumerate(transcripts):
                if transcript_1 == transcript_2:
                    matrix[i,j] = 1.0
                else:
                    indice1 = transcript_1+'&'+transcript_2
                    indice2 = transcript_2+'&'+transcript_1
                    if indice1 not in list(used_indice.keys()) and indice2 not in list(used_indice.keys()):
                        gene_tr_1 = df_gtot[df_gtot['id_transcript'] == transcript_1].id_gene.values[0]
                        gene_tr_2 = df_gtot[df_gtot['id_transcript'] == transcript_2].id_gene.values[0]
                        score_nuc_hom_v2 = round(compute_nuchom_score_tsmplus(transcript_1, transcript_2, df_transcripts, gene_tr_1, gene_tr_2, gene_model),3)
                        score_deg_hom_v2 = round(compute_deghom_score_tsmplus(transcript_1, transcript_2, df_blocks_transcripts, df_blocks_genes, gene_tr_1, gene_tr_2),3)
                        score_tsmplus = round(((score_nuc_hom_v2+score_deg_hom_v2)/2),3)
                        matrix[i, j] = score_tsmplus
                        used_indice[indice1] = score_tsmplus
                        used_indice[indice2] = score_tsmplus
                    else :
                        matrix[i, j] = used_indice[indice1]
    else:
        raise Exception('Error! the value of tsm parameter does not match!')
    df_matrix = pd.DataFrame(matrix, index=transcripts, columns=transcripts)
    return df_matrix

def get_gene_model(df_gtot, df_transcripts):
    list_all_distinct_genes = list(pd.Categorical(df_gtot['id_gene'].values))
    list_save_id_gene = []
    list_save_sequences = []
    for distinct_gene in list_all_distinct_genes:
        list_all_distinct_transcripts = list(pd.Categorical(df_gtot[df_gtot['id_gene']==distinct_gene].id_transcript.values))
        list_all_alg_distinct_transcripts = list(df_transcripts[df_transcripts['id_transcript'].isin(list_all_distinct_transcripts)].alg.values)
        length_seq = len(list_all_alg_distinct_transcripts[0])
        list_save_id_gene.append(distinct_gene)
        list_save_sequence = []
        for position in range(length_seq):
            is_conserved = True
            for alg_distinct_transcript in list_all_alg_distinct_transcripts:
                if alg_distinct_transcript[position] == '-':
                    list_save_sequence.append('-')
                    is_conserved = False
                    break
                if is_conserved:
                    list_save_sequence.append('c')
        list_save_sequences.append(''.join(list_save_sequence))
    df = pd.DataFrame(columns=['id_gene','gene_model'])
    df['id_gene'] = list_save_id_gene
    df['gene_model'] = list_save_sequences
    return df

def compute_deghom_score_tsmplus(tr_id, tr_ref, transcripts_blocks,data_block_alg_gene, g_id, g_ref ):
    """
        Compute the degHom score (do not account for blocks not appearing at the gene level)
    """
    block_id = transcripts_blocks[transcripts_blocks['id_transcript']==tr_id].blocks.values[0].split('|')
    block_ref = transcripts_blocks[transcripts_blocks['id_transcript']==tr_ref].blocks.values[0].split('|')
    block_g_id = data_block_alg_gene[data_block_alg_gene['id_gene']==g_id].blocks.values[0].split('|')
    block_g_ref = data_block_alg_gene[data_block_alg_gene['id_gene']==g_ref].blocks.values[0].split('|')

    intersect_set = list(set(block_id).intersection(set(block_ref)))
    union_set = list(set(block_id).union(set(block_ref)))
    denum_list = [_ for _ in union_set if _ in block_g_id and _ in block_g_ref]
    if len(denum_list) == 0:
        return 0
    score = len(intersect_set)/len(denum_list)   
    return score

def nucleotide_in_gene(gene_tr, position, gene_model):
    gene_model_seq = gene_model[gene_model['id_gene']==gene_tr].gene_model.values[0]
    if gene_model_seq[position] == '-':
        return False
    elif gene_model_seq[position] == 'c':
        return True
    else:
        raise ValueError('Could not find the conservation into the gene model')

def compute_nuchom_score_tsmplus(transcript_1, transcript_2, data_alg_seq, gene_tr_1, gene_tr_2, gene_model):
    """
        Compute nucHom score (do not account for nucleotides not appearing at the gene level)
    """
    alignment_tr1 = data_alg_seq[data_alg_seq['id_transcript']==transcript_1].alg.values[0]
    alignment_tr2 = data_alg_seq[data_alg_seq['id_transcript']==transcript_2].alg.values[0]
   
    length_alignment_tr1 = len(alignment_tr1)
    length_alignment_tr2 = len(alignment_tr2)
    if length_alignment_tr1 != length_alignment_tr2:
        return False
    else:
        matchs_mismatchs = 0
        indel = 0
        gap_gap = 0
        indel_denum = 0
        for position in range(length_alignment_tr1):
            nucleotide_tr1 = alignment_tr1[position]
            nucleotide_tr2 = alignment_tr2[position]
            if((nucleotide_tr1 != '-') and (nucleotide_tr2 != '-')):
                matchs_mismatchs += 1
            elif(((nucleotide_tr1 == '-') and (nucleotide_tr2 != '-'))) or ((nucleotide_tr1 != '-') and (nucleotide_tr2 == '-')):
                indel += 1
                if (nucleotide_tr1 != '-') and (nucleotide_tr2 == '-'):
                    nuc_g2 = nucleotide_in_gene(gene_tr_2, position, gene_model)
                    nuc_g1 = True
                else:
                    nuc_g1 = nucleotide_in_gene(gene_tr_1, position, gene_model)
                    nuc_g2 = True
                    
                if nuc_g1 and nuc_g2:
                    indel_denum += 1
            else:
                gap_gap += 1
        if (matchs_mismatchs + indel_denum) == 0:
            score = 0
            return score
        else:
            score = score = (matchs_mismatchs) / (matchs_mismatchs + indel_denum)
            return score

def compute_deghom_score(tr_id, tr_ref, transcripts_blocks):
    """
        Compute the degHom score (do account for blocks not appearing at the gene level)
    """
    block_id = transcripts_blocks[transcripts_blocks['id_transcript']==tr_id].blocks.values[0].split('|')
    block_ref = transcripts_blocks[transcripts_blocks['id_transcript']==tr_ref].blocks.values[0].split('|')
    blocks_in = 0
    blocks_out = 0
    for block in block_id:
        if block in block_ref:
            blocks_in += 1
        else:
            blocks_out += 1
    all_blocks = [_ for _ in block_id]
    for _ in block_ref:
        all_blocks.append(_)
    score = blocks_in / len(list(set(all_blocks)))             
    return score

def compute_nuchom_score(transcript_1, transcript_2, data_alg_seq):
    """
        Compute nucHom score (do account for blocks not appearing at the gene level)
    """
    alignment_tr1 = data_alg_seq[data_alg_seq['id_transcript']==transcript_1].alg.values[0]
    alignment_tr2 = data_alg_seq[data_alg_seq['id_transcript']==transcript_2].alg.values[0]
    length_alignment_tr1 = len(alignment_tr1)
    length_alignment_tr2 = len(alignment_tr2)
    if length_alignment_tr1 != length_alignment_tr2:
        return False
    else:
        matchs_mismatchs = 0
        indel = 0
        gap_gap = 0
        for position in range(length_alignment_tr1):
            nucleotide_tr1 = alignment_tr1[position]
            nucleotide_tr2 = alignment_tr2[position]
            if((nucleotide_tr1 != '-') and (nucleotide_tr2 != '-')):
                matchs_mismatchs += 1
            elif(((nucleotide_tr1 == '-') and (nucleotide_tr2 != '-'))) or ((nucleotide_tr1 != '-') and (nucleotide_tr2 == '-')):
                indel += 1
            else:
                gap_gap += 1
        score = (matchs_mismatchs) / (matchs_mismatchs + indel)
    return score

if __name__ == '__main__':
    # retrieve inputs given by user
    args = build_arg_parser().parse_args()
    transcripts_msa_path = args.tralignment
    gtot_path = args.genetotranscripts
    tsm_conditions_path = int(args.tsmuncorrected)
    output_folder_path = args.outputfolder

    # compute the algorithm
    get_matrix(transcripts_msa_path, gtot_path, tsm_conditions_path, output_folder_path)
