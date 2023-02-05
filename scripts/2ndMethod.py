import pandas as pd
import numpy as np
import os
import concurrent.futures

def tsm_measure_v2():
    return True

def compute_deghom_score_v2(tr_id, tr_ref, transcripts_blocks,data_block_alg_gene, g_id, g_ref ):
    """
        Compute the degHom score
    """
    block_id = transcripts_blocks[transcripts_blocks['id_transcript']==tr_id].block_sequence.values[0].split('.')
    block_ref = transcripts_blocks[transcripts_blocks['id_transcript']==tr_ref].block_sequence.values[0].split('.')
    block_g_id = data_block_alg_gene[data_block_alg_gene['id_gene']==g_id].block_sequence.values[0].split('.')
    block_g_ref = data_block_alg_gene[data_block_alg_gene['id_gene']==g_ref].block_sequence.values[0].split('.')

    intersect_set = list(set(block_id).intersection(set(block_ref)))
    union_set = list(set(block_id).union(set(block_ref)))
    denum_list = [_ for _ in union_set if _ in block_g_id and _ in block_g_ref]
    print('--------------------------------')
    print(tr_id, tr_ref)
    print(intersect_set)
    print(union_set)
    print(denum_list)
    if len(denum_list) == 0:
        return 0
    score = len(intersect_set)/len(denum_list)  
    print(score)   
    return score

def compute_nuchom_score_v2(transcript_1, transcript_2, data_alg_seq, data_alg_gene_seq, gene_tr_1, gene_tr_2):
    """
        Compute nucHom score
    """
    alignment_tr1 = data_alg_seq[data_alg_seq['id_transcript']==transcript_1].alg.values[0]
    alignment_tr2 = data_alg_seq[data_alg_seq['id_transcript']==transcript_2].alg.values[0]

    alignment_g1 = data_alg_gene_seq[data_alg_gene_seq['id_gene']==gene_tr_1].alg_seq.values[0].replace('!','-')
    alignment_g2 = data_alg_gene_seq[data_alg_gene_seq['id_gene']==gene_tr_2].alg_seq.values[0].replace('!','-')

    #print(alignment_g1)
   
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
                #print(length_alignment_tr1)
                #print(length_alignment_tr2)
                #print(len(alignment_g1))
                #print(len(alignment_g2))
                nuc_g2 = alignment_g2[position]
                if nuc_g1 != '-' and nuc_g2 != '-':
                    indel_denum += 1
            else:
                gap_gap += 1
        score = (matchs_mismatchs) / (matchs_mismatchs + indel_denum)
        #print(score)
    return score

def compute_matrix(data_alg_seq, data_alg_block, data_transcripts_to_gene, data_alg_gene_seq, data_block_alg_gene):
    """ 
        ### Compute the TDB transcript similarity score
        ##### inputs :
        ##### argument 1 : Pandas DataFrame = {columns=['id_transcript', 'alg']} containing ids of transcripts and their aligned sequences of nucleotides
        ##### argument 2 : Pandas DataFrame = {columns=['id_transcript', 'block_sequence']} containing ids of transcripts and their aligned sequences of blocks
        ##### argument 3 : Pandas DataFrame = {columns=['id_transcript', 'id_gene']} containing ids of transcripts and their corresponding genes
    """
    score_nuc= []
    score_deg=[]
    score_moy = []
    score_nucv2= []
    score_degv2=[]
    score_moyv2 = []
    transcripts = list(data_alg_seq['id_transcript'].values)
    matrix = []
    for transcript_1 in transcripts:
        row = []
        for transcript_2 in transcripts:
            gene_tr_1 = data_transcripts_to_gene[data_transcripts_to_gene['id_transcript'] == transcript_1].id_gene.values[0]
            gene_tr_2 = data_transcripts_to_gene[data_transcripts_to_gene['id_transcript'] == transcript_2].id_gene.values[0]

            '''if gene_tr_1 == gene_tr_2:
                #score = np.nan
                score_deg_hom = 0
                score_nuc_hom_v2 = 0
                score_nuc_hom = 0
                score_deg_hom_v2 =0
                scorev2 = np.nan
                score = np.nan
            else:'''
            # Compute the nucHOM score
            score_nuc_hom = round(compute_nuchom_score(transcript_1, transcript_2, data_alg_seq),3)
            #score_nuc_hom_v2 = round(compute_nuchom_score_v2(transcript_1, transcript_2, data_alg_seq, data_alg_gene_seq, gene_tr_1, gene_tr_2),3)
            # Compute the DegHOM score
            score_deg_hom = round(compute_deghom_score(transcript_1, transcript_2, data_alg_block),3)
            #score_deg_hom_v2 = round(compute_deghom_score_v2(transcript_1, transcript_2, data_alg_block, data_block_alg_gene, gene_tr_1, gene_tr_2),3)
            # Compute the TDB Transcript Similarity Measure 
            score = round(((score_nuc_hom+score_deg_hom)/2),3)
            #scorev2 = round(((score_nuc_hom_v2+score_deg_hom_v2)/2),3)

            row.append(score)
            #score_nuc.append(score_nuc_hom)
            #score_deg.append(score_deg_hom)
            #score_moy.append(score)

            #score_nucv2.append(score_nuc_hom_v2)
            #score_degv2.append(score_deg_hom_v2)
            #score_moyv2.append(scorev2)
        matrix.append(row)
    df_matrix = pd.DataFrame(matrix, index=transcripts, columns=transcripts)

    #df_m1 = pd.DataFrame(data={'nucHom':score_nuc, 'degHom':score_deg, 'tsm':score_moy})
    #df_m2 = pd.DataFrame(data={'nucHom':score_nucv2, 'degHom':score_degv2, 'tsm':score_moyv2})
    return [df_matrix]

def compute_deghom_score(tr_id, tr_ref, transcripts_blocks):
    """
        Compute the degHom score
    """
    block_id = transcripts_blocks[transcripts_blocks['id_transcript']==tr_id].block_sequence.values[0].split('.')
    block_ref = transcripts_blocks[transcripts_blocks['id_transcript']==tr_ref].block_sequence.values[0].split('.')
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
    #print(all_blocks)
    #print(set(all_blocks))
    score = blocks_in / len(list(set(all_blocks)))             
    return score

def compute_nuchom_score(transcript_1, transcript_2, data_alg_seq):
    """
        Compute nucHom score
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

def formatting_data(file_algn_transcripts_path, file_tr_to_gene_path, file_algn_genes_path, gt):
    """
        Format inputs of users into Pandas DataFrames
    """
    transcript_name = []
    tr_alignment = []
    data_cds_alg = pd.DataFrame(
        columns=['id_transcript','alg']
    )
    file_open = open(file_algn_transcripts_path,'r')
    lines = file_open.readlines()
    for (i,line) in enumerate(lines):
        if line.startswith('>'):
            transcript_name.append(line.split('\n')[0].split('>')[-1])
            tr_alignment.append(lines[i+1].split('\n')[0])
    file_open.close()
    data_cds_alg['id_transcript'] = transcript_name
    data_cds_alg['alg'] = tr_alignment

    transcript_name_two = []
    gene_name = []
    data_transcripts_to_gene = pd.DataFrame(
        columns=['id_transcript','id_gene']
    )
    
    file_open_two = open(file_tr_to_gene_path,'r')
    lines_two = file_open_two.readlines()
    for (i,line) in enumerate(lines_two):
        line_G = lines_two[i]
        if line_G.startswith('>'):
            transcript_name_two.append(line_G.split(':')[0].split('>')[-1])
            gene_name.append(line_G.split(':')[1].split('\n')[0])
    file_open_two.close()
    data_transcripts_to_gene['id_transcript'] = transcript_name_two
    data_transcripts_to_gene['id_gene'] = gene_name
    data_block_alg = retrieve_aligned_blocks_sequences(file_algn_transcripts_path)
    data_block_alg_gene = retrieve_aligned_blocks_gene_sequences(data_block_alg, data_transcripts_to_gene)
    data_gene_alg = retrieve_aligned_gene_sequences(file_algn_genes_path)

    #Next lines to delete
    data_block_alg.to_csv('./Confusion_matrix/tsm/Macse/m/{}/{}_blocks_alg.csv'.format(gt,gt), header=True, index=False, sep=';')
    data_cds_alg.to_csv('./Confusion_matrix/tsm/Macse/m/{}/{}_cds_alg.csv'.format(gt,gt), header=True, index=False, sep=';')
    return [data_cds_alg, data_block_alg, data_transcripts_to_gene, data_gene_alg, data_block_alg_gene]

def retrieve_aligned_gene_sequences(file_algn_genes_path):
    #print(file_algn_genes_path)
    file = open(file_algn_genes_path, 'r')
    lines = [str(_).split('\n')[0] for _ in file.readlines()]
    #print(len(lines))
    all_genes = []
    all_sequences = []
    for i in range(len(lines)):
        #print(i)
        line = lines[i]
        #print(line[i+1])
        if line.startswith('>'):
            gene = line.split('>')[1]
            #print(line)
            j = i+1
            sequence = ''
            linesj = lines[j]
            while not linesj.startswith('>'):
                sequence += linesj
                j += 1
                if j < len(lines):
                    linesj = lines[j]
                else:
                    break
            #print(len(sequence))
            all_genes.append(gene)
            all_sequences.append(sequence)
            #print('++++++++++++++++++++++++++++++===')    
    df = pd.DataFrame(data={'id_gene': all_genes, 'alg_seq': all_sequences})
    #print(df)
    file.close()
    return df

def retrieve_aligned_blocks_gene_sequences(data_block_alg, data_transcripts_to_gene):
    """
        Returns the blocks alignement of genes
    """
    #print(data_transcripts_to_gene)
    all_genes_blocks = []
    all_genes_id = []
    genes = list(data_transcripts_to_gene.id_gene.values)
    #print(genes)
    #print(data_block_alg)
    for gene in genes:
        transcripts = list(data_transcripts_to_gene[data_transcripts_to_gene['id_gene']==gene].id_transcript.values)
        #print(transcripts)
        list_of_nblocks = list(data_block_alg[data_block_alg['id_transcript'].isin(transcripts)].block_sequence.values)
        set_of_nblocks = [set(str(_).split('.')) for _ in list_of_nblocks]
        #print(set_of_nblocks)
        union_set_blocks = set()
        if len(set_of_nblocks) > 1:
            for i in range(len(set_of_nblocks)):
                union_set_blocks = union_set_blocks.union(set_of_nblocks[i])
        else:
            union_set_blocks = set_of_nblocks[0]
        list_union_set_blocks = list(union_set_blocks)
        gene_blocks = '.'.join(list_union_set_blocks)

        all_genes_blocks.append(gene_blocks)
        all_genes_id.append(gene)

    df = pd.DataFrame(data={'id_gene': all_genes_id, 'block_sequence': all_genes_blocks})
    return df

def retrieve_aligned_blocks_sequences(file):
    """ 
        Retrieve aligned of sequences of transcripts given at inputs
    """
    aligned_blocks = get_aligned_blocks(file)
    blocks = find_blocks_in_alg_transcripts(aligned_blocks)
    transcripts_blocks = map_blocks(aligned_blocks,blocks)
    return transcripts_blocks

def get_aligned_blocks(input_fasta_file_path):
    fasta_file = open(input_fasta_file_path, 'r')
    lines = fasta_file.readlines()
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
    fasta_file.close()
    return (list_of_ids_transcripts, blocks_data)

def find_blocks_in_alg_transcripts(aligned_blocks):
    """
        Find blocks in a Multiple Sequences Alignment
    """
    blocks_data = aligned_blocks[1]
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
    return positions_blocks

def map_blocks(aligned_blocks, blocks):
    """
        Map blocks into transcripts
    """
    transcripts = aligned_blocks[0]
    blocks_data = aligned_blocks[1]
    sequences = {}
    for i in range(len(blocks_data)):
        transcript = transcripts[i]
        list_sequence = [str(_) for _ in blocks_data[i]]
        sequence = ''.join(list_sequence)
        sequence_blocks = []
        for block in blocks.keys():
            tuple_position = blocks[block]
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
    df_results = pd.DataFrame(columns=['id_transcript','block_sequence'])
    i_data = []
    j_data = []
    for transcript in sequences.keys():
        value_data = sequences[transcript]
        i_data.append(transcript)
        j_data.append('.'.join([str(_) for _ in value_data]))
    df_results['id_transcript'] = i_data
    df_results['block_sequence'] = j_data
    return df_results

def main_function_tt(gt, inputs):
    """Main function"""
    data_alg_seq = inputs[0]
    data_alg_block = inputs[1]
    data_transcripts_to_gene = inputs[2]
    data_alg_gene_seq = inputs[3]
    data_block_alg_gene = inputs[4]
    dataframes = compute_matrix(data_alg_seq, data_alg_block, data_transcripts_to_gene, data_alg_gene_seq, data_block_alg_gene)
    #results = find_rbhs(matrix_out)
    #matrix_out_res = matrix_out
    #homologs_res = results
    #block_res = data_alg_block.set_index(data_alg_block['id_transcript'])
    #print(matrix_out)
    #df1 = dataframes[1]
    #df2 = dataframes[2]
    matrix = dataframes[0]
    
    #df1.to_csv('./Confusion_matrix/tsm/Kalign/nucHom/{}/{}_m1.csv'.format(gt,gt), header=True, sep=';', index=False)
    #df2.to_csv('./Confusion_matrix/tsm/Kalign/nucHom/{}/{}_m2.csv'.format(gt,gt), header=True, sep=';', index=False)
    matrix.to_csv('./Confusion_matrix/tsm/Macse/m/{}/{}_matrix.csv'.format(gt,gt), header=True, sep=';')
    data_transcripts_to_gene.to_csv('./Confusion_matrix/tsm/Macse/m/{}/{}_gtot.csv'.format(gt,gt), header=True, index=False, sep=';')
    return True

def parallelize_function(gt):
    print(gt)
    path_file_algn_transcripts = './macse_253/alignment_transcripts_sequences/{}.alg'.format(gt+'_98')
    path_file_algn_genes = './macse_253/alignment_genes_sequences/{}'.format(gt)
    path_file_gtot = './macse_253/gtot/{}_2.fasta'.format(gt)
    cmd = './Confusion_matrix/tsm/Macse/m/{}'.format(gt)
    os.system('mkdir {}'.format(cmd))
    data_transcript = formatting_data(path_file_algn_transcripts, path_file_gtot, path_file_algn_genes, gt)
    print(data_transcript)

    # Compute the main program
    results = main_function_tt(gt, data_transcript)
    return True


if __name__ == '__main__':
    file = open('./macse_253/trees_macse.txt','r')
    gts = [str(_).split('\n')[0] for _ in file.readlines()]
    # print(gts)
    #for gt in gts:
        #'ENSGT00390000008738'
        #if gt not in ['ENSGT00390000012186','ENSGT00390000012913','ENSGT00530000063156']:
        #if gt == 'ENSGT00390000008738':
        #parallelize_function(gt)
    '''
    for gt in gts:
        path_file_algn_transcripts = './alignements_transcripts/{}.alg'.format(gt)
        path_file_algn_genes = './alignements_genes/{}.alg'.format(gt)
        path_file_gtot = './gtot/{}_2.fasta'.format(gt)
        cmd = './results/{}'.format(gt)
        os.system('mkdir {}'.format(cmd))
        data_transcript = formatting_data(path_file_algn_transcripts, path_file_gtot, path_file_algn_genes, gt)
        #print(data_transcript)

        # Compute the main program
        results = main_function_tt(gt, data_transcript)'''
    file.close()
    with concurrent.futures.ProcessPoolExecutor(max_workers=5) as executor:
        executor.map(parallelize_function, gts)
    