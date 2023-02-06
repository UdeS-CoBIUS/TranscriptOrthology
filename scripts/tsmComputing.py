""" :dna: returns the similarity matrix tsm+ or tsm scores for all pairs of homologous transcripts.
> Usage:
======
    python3 tsm-computing.py [-talg] [-galg] [-gtot] [-tsm] [-outf]

> Reference:
======
    https://github.com/UdeS-CoBIUS/TranscriptOrthology
"""
__authors__ = ("Wend Yam Donald Davy Ouedraogo")
__contact__ = ("wend.yam.donald.davy.usherbrooke.ca")
__copyright__ = "CoBIUS lab at Université de Sherbrooke, QC, CANADA"
__date__ = "2022-12-19"
__version__= "1.0.0"

import sys
import pandas as pd
import numpy as np
import argparse
import time

def build_arg_parser():
    '''Parsing function'''
    parser = argparse.ArgumentParser(description="parsor program parameter")
    parser.add_argument('-talg', '--transcriptsalignment', default=None)
    parser.add_argument('-galg', '--genesalignment', default=None)
    parser.add_argument('-gtot', '--genetotranscripts', default=None)
    parser.add_argument('-tsm', '--tsmuncorrected', default=False)
    parser.add_argument('-outf', '--outputfolder', default='.')
    return parser

def saving_files(file, outputfolder, name, ind):
    file.to_csv('{}/{}.csv'.format(outputfolder,name), header=True, index=ind, sep=';')
    return True

def get_matrix(transcripts_msa_path, genes_msa_path, gtot_path, tsm_condition, output_folder_path):
    '''Returns the matrix score'''
    try:
        print('++++++++++++++++Starting ....')
        start = time.time()
        df_transcripts, df_genes, df_gtot, df_blocks_transcripts, df_blocks_genes = convert_files_to_dataframes(transcripts_msa_path, genes_msa_path, gtot_path, output_folder_path)
        print('+++++++ All data retrieved & the representation of subtranscribed sequences of genes into blocks are available.')
    except:
        raise('Something wrong with the inputs! Please check it out or contact us. Thank you.')
    
    try:
        print('------- Computing matrix ... in progress')
        tsm_matrix = compute_matrix(df_transcripts, df_genes, df_gtot, df_blocks_transcripts, df_blocks_genes, tsm_condition)
        #saving matrix
        saving_files(tsm_matrix, output_folder_path, 'matrix', True)
        end = time.time()
        print('------- Computing matrix ... status: Finished without errors in {} seconds'.format(str(end-start)))
    except:
        raise('The algorithm exits with errors!! Please contact us for help. Thank you.')
    
    return tsm_matrix

def convert_files_to_dataframes(transcripts_msa_path,genes_msa_path,gtot_path, output_folder_path):
    """Format inputs of users into Pandas DataFrames"""

    df_transcripts = pd.DataFrame(columns=['id_transcript','alg'])
    df_genes = pd.DataFrame(columns=['id_gene','alg'])
    df_gtot = pd.DataFrame(columns=['id_transcript','id_gene'])

    file_transcripts_msa = open(transcripts_msa_path, 'r')
    transcripts_msa_rows = [str(_).split('\n')[0] for _ in file_transcripts_msa.readlines()]
    file_transcripts_msa.close()
    for i_number, row in enumerate(transcripts_msa_rows):
        if row.startswith('>'):
            df_transcripts.loc[row] = [row.split('>')[-1], transcripts_msa_rows[i_number+1]]
    
    file_genes_msa = open(genes_msa_path, 'r')
    genes_msa_rows = [str(_).split('\n')[0] for _ in file_genes_msa.readlines()]
    file_genes_msa.close()
    for g_number, row_g in enumerate(genes_msa_rows):
        if row_g.startswith('>'):
            df_genes.loc[row_g] = [row_g.split('>')[-1], genes_msa_rows[g_number+1]]
    
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
    
    return df_transcripts, df_genes, df_gtot, df_blocks_transcripts, df_blocks_genes

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
    for i in range(len(positions)-1):
        interval_start = positions[i]
        interval_stop = positions[i+1]
        interval = (interval_start, interval_stop)
        positions_blocks[i] = interval
        number = i
    interval_start_last_block = positions[-1]
    interval_stop_last_block = length_column
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
        gene_blocks = '|'.join(list_union_set_blocks)

        all_genes_blocks.append(gene_blocks)
        all_genes_id.append(gene)

    df = pd.DataFrame(data={'id_gene': all_genes_id, 'blocks': all_genes_blocks})
    return df

def compute_matrix(df_transcripts, df_genes, df_gtot, df_blocks_transcripts, df_blocks_genes, tsm_condition):
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
    matrix = []
    if tsm_condition:
        for transcript_1 in transcripts:
            row = []
            for transcript_2 in transcripts:
                # Compute the nucHOM score (unitary)
                score_nuc_hom = round(compute_nuchom_score(transcript_1, transcript_2, df_transcripts),3)
                # Compute the DegHOM score (length)
                score_deg_hom = round(compute_deghom_score(transcript_1, transcript_2, df_blocks_transcripts),3)
                # Compute tsm (mean)
                score = round(((score_nuc_hom+score_deg_hom)/2),3)
                row.append(score)
            matrix.append(row)
    else:
        for transcript_1 in transcripts:
            row = []
            for transcript_2 in transcripts:
                gene_tr_1 = df_gtot[df_gtot['id_transcript'] == transcript_1].id_gene.values[0]
                gene_tr_2 = df_gtot[df_gtot['id_transcript'] == transcript_2].id_gene.values[0]
                score_nuc_hom_v2 = round(compute_nuchom_score_tsmplus(transcript_1, transcript_2, df_transcripts, df_genes, gene_tr_1, gene_tr_2),3)
                score_deg_hom_v2 = round(compute_deghom_score_tsmplus(transcript_1, transcript_2, df_blocks_transcripts, df_blocks_genes, gene_tr_1, gene_tr_2),3)
                score_tsmplus = round(((score_nuc_hom_v2+score_deg_hom_v2)/2),3)
                row.append(score_tsmplus)
            matrix.append(row)
    df_matrix = pd.DataFrame(matrix, index=transcripts, columns=transcripts)
    return df_matrix

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

def compute_nuchom_score_tsmplus(transcript_1, transcript_2, data_alg_seq, data_alg_gene_seq, gene_tr_1, gene_tr_2):
    """
        Compute nucHom score (do not account for blocks not appearing at the gene level)
    """
    alignment_tr1 = data_alg_seq[data_alg_seq['id_transcript']==transcript_1].alg.values[0]
    alignment_tr2 = data_alg_seq[data_alg_seq['id_transcript']==transcript_2].alg.values[0]

    alignment_g1 = data_alg_gene_seq[data_alg_gene_seq['id_gene']==gene_tr_1].alg.values[0]
    alignment_g2 = data_alg_gene_seq[data_alg_gene_seq['id_gene']==gene_tr_2].alg.values[0]
   
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
                nuc_g1 = alignment_g1[position]
                nuc_g2 = alignment_g2[position]
                if nuc_g1 != '-' and nuc_g2 != '-':
                    indel_denum += 1
            else:
                gap_gap += 1
        score = (matchs_mismatchs) / (matchs_mismatchs + indel_denum)
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
    transcripts_msa_path = args.transcriptsalignment
    genes_msa_path = args.genesalignment
    gtot_path = args.genetotranscripts
    tsm_conditions_path = args.tsmuncorrected
    output_folder_path = args.outputfolder

    # compute the algorithm
    get_matrix(transcripts_msa_path, genes_msa_path, gtot_path, tsm_conditions_path, output_folder_path)