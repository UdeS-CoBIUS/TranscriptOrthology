""" Performs the generation of graph clustering and returns connected components in the graph.
> Usage:
======
    python3 t-clustering.py [-m] [-nhxt] [-gtot] [-lowb] [-outf]

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
import networkx as nx
import time
from ete3 import Tree
import matplotlib.pyplot as plt
import argparse

def build_arg_parser():
    '''Parsing function'''
    parser = argparse.ArgumentParser(description="parsor program parameter")
    parser.add_argument('-m', '--matrix', default=None)
    parser.add_argument('-gtot', '--genetotranscripts', default=None)
    parser.add_argument('-nhxt', '--nhxgenetree', default=None)
    parser.add_argument('-lowb', '--lowerbound', default=0.5)
    parser.add_argument('-outf', '--outputfolder', default='.')
    return parser

def get_orthology_graph(matrix, gtot_path, gt_path, lower_bound, output_folder):
    df_gtot = convert_inputs_data(gtot_path)
    # matrix = pd.read_csv(matrix_path, sep=';', index_col=0)
    try:
        start1 = time.time()
        print('+++++ Searching for recent-paralogs ... \tstatus: processing')
        recent_paralogs = search_recent_paralogs(matrix, df_gtot)
        end1 = time.time()
        print('+++++ Searching for recent-paralogs ... \tstatus: finished in {} seconds'.format(str(end1-start1)))
    except:
        raise('Failed to retrieve recent-paralogs ! Errors occured.')
    
    try:
        start2 = time.time()
        print('+++++ Searching for RBHs ... \tstatus: processing')
        graphG = rhb_clustering(matrix, df_gtot, recent_paralogs, lower_bound)
        end2 = time.time()
        print('+++++ Searching for RBHs ... \tstatus: finished in {} seconds'.format(str(end2-start2)))
    except:
        raise('Failed to retrieve RBHs ! Errors occured.')
    
    try:
        start3 = time.time()
        print('+++++ Construction of the orthology graph (Adding nodes ...) ... \tstatus: processing')
        clusters = get_conserved_clusters(graphG, recent_paralogs, df_gtot, output_folder)
        end3 = time.time()
        print('+++++ Construction of the orthology graph (Adding nodes ...) ... \tstatus: finished in {} seconds'.format(str(end3-start3)))
    except:
        raise('Failed to retrieve connected components ! Errors occured.')
    
    try:
        start4 = time.time()
        print('+++++ Searching for connected components ... \tstatus: processing')
        df, df_orthology = write_results(clusters, df_gtot, matrix, gt_path, output_folder)
        end4 = time.time()
        print('+++++ Searching for connected components ... \tstatus: finished in {} seconds'.format(str(end4-start4)))

    except:
        raise('Failed to save results ! Errors occured.')
    return clusters, df, df_orthology

def convert_inputs_data(gtot_path):
    df_gtot = pd.DataFrame(columns=['id_transcript','id_gene'])
    file_gtot = open(gtot_path, 'r')
    gtot_rows = [str(_).split('\n')[0] for _ in file_gtot.readlines()]
    file_gtot.close()
    for row_gtot in gtot_rows:
        if row_gtot.startswith('>'):
            df_gtot.loc[row_gtot] = [row_gtot.split(':')[0].split('>')[-1], row_gtot.split(':')[-1]]
    return df_gtot

def search_recent_paralogs(matrix, g):
    all_transcripts = list(g.id_transcript.values)
    dict_inParalogs = {}

    genes = list(pd.Categorical(g.id_gene).categories.values)
    for gene in genes:
        transcripts = list(g[g['id_gene']==gene].id_transcript.values)
        if len(transcripts) == 1:
            transcript = transcripts[0]
            dict_inParalogs[transcript] = []
        else:
            for transcript in transcripts:
                transcripts_without_tr = [tr for tr in transcripts if tr != transcript]
                tr_matrix = matrix[matrix.index==transcript][transcripts_without_tr]
                max_value_inParalogs = max(tr_matrix.values[0])
                if max_value_inParalogs != 0:
                    tr_inParalogs = list(tr_matrix.apply(lambda row: row[row == max_value_inParalogs].index, axis=1)[0].values)

                    all_transcripts_without_tr = [tr for tr in all_transcripts if tr != transcript]
                    all_tr_matrix = matrix[matrix.index==transcript][all_transcripts_without_tr]
                    max_value_outParalogs = max(all_tr_matrix.values[0])
                    tr_outParalogs = list(all_tr_matrix.apply(lambda row: row[row == max_value_outParalogs].index, axis=1)[0].values)
                    if len(tr_inParalogs) > 1 or len(tr_inParalogs)==1:
                        list_inParalogs = []
                        for tr_inParalog in tr_inParalogs:
                            if tr_inParalog in tr_outParalogs:
                                transcripts_without_tr_ref = [tr for tr in transcripts if tr != tr_inParalog]
                                tr_matrix_ref = matrix[matrix.index==tr_inParalog][transcripts_without_tr_ref]
                                max_value_inParalogs_ref = max(tr_matrix_ref.values[0])
                                if max_value_inParalogs_ref != 0:
                                    tr_inParalogs_ref = list(tr_matrix_ref.apply(lambda row: row[row == max_value_inParalogs_ref].index, axis=1)[0].values)
                                    all_transcripts_without_tr_ref = [tr for tr in all_transcripts if tr != tr_inParalog]
                                    all_tr_matrix_ref = matrix[matrix.index==tr_inParalog][all_transcripts_without_tr_ref]
                                    max_value_outParalogs_ref = max(all_tr_matrix_ref.values[0])
                                    tr_outParalogs_ref = list(all_tr_matrix_ref.apply(lambda row: row[row == max_value_outParalogs_ref].index, axis=1)[0].values)
                                    intersect_trs = list(set(tr_inParalogs_ref).intersection(set(tr_outParalogs_ref)))
                                    if transcript in intersect_trs:
                                        # before adding we need to check if 
                                        list_inParalogs.append(tr_inParalog)
                        dict_inParalogs[transcript] = list_inParalogs
                    else:
                        dict_inParalogs[transcript] = []
                else:
                    dict_inParalogs[transcript] = []

    return dict_inParalogs

def rhb_clustering(matrix, g, inparalogs, lower_bound):
    """Returns a similarity graph G containing recent paralogs and the recent paralogs details"""
    
    #warning if lower_bound is greater than 1
    if lower_bound > 1 or lower_bound < 0:
        raise Exception('Each transcript is considered as a cluster. Check the lower bound! Enter a lower bound less greater than 1')
    elif lower_bound < 0:
        raise Exception('A negative value has been entered. Check the lower bound! Enter a valid lower bound')
    else:
        pass

    #initialize the graph G
    edges_rbhs = []
    pair_id = []
    pair_ref = []
    score_pair = []
    G = nx.Graph()
    G.add_nodes_from(matrix.index)
    for tr_recent in list(inparalogs.keys()):
        tr_recent_paralogs = inparalogs[tr_recent]
        if len(tr_recent_paralogs) != 0:
            for tr_recent_paralog in tr_recent_paralogs:
                G.add_edge(tr_recent, tr_recent_paralog)
    #nx.draw(G, with_labels=True)
    #plt.show(G)
    

    genes = list(pd.Categorical(g.id_gene).categories.values)
    for gene in genes:
        transcripts = list(g[g['id_gene']==gene].id_transcript.values)
        genes_without_gene = [_ for _ in genes if _ != gene]
        for transcript in transcripts:
            for gene_without_gene in genes_without_gene:
                other_transcripts = list(g[g['id_gene']==gene_without_gene].id_transcript.values)
                tr_matrix = matrix[matrix.index==transcript][other_transcripts]
                max_value = max(tr_matrix.values[0])
                if max_value != 0 and max_value >= lower_bound:
                    tr_orthologs = list(tr_matrix.apply(lambda row: row[row == max_value].index, axis=1)[0].values)
                    # true_rbhs = []
                    for tr_ortholog in tr_orthologs:
                        tr_ortholog_matrix = matrix[matrix.index==tr_ortholog][transcripts]
                        max_value_tr_ortholog = max(tr_ortholog_matrix.values[0])
                        if max_value_tr_ortholog >= lower_bound:
                            tr_reciprocal_orthologs = list(tr_ortholog_matrix.apply(lambda row: row[row == max_value_tr_ortholog].index, axis=1)[0].values)
                            if transcript in tr_reciprocal_orthologs:                               
                                if set([tr_ortholog, transcript]) not in edges_rbhs:
                                    pair_id.append(transcript)
                                    pair_ref.append(tr_ortholog)
                                    score_pair.append(max_value_tr_ortholog)
                                    edges_rbhs.append(set([tr_ortholog, transcript]))
                                
    df_edges_rbhs = pd.DataFrame(data={'pair_id': pair_id, 'pair_ref': pair_ref, 'score': score_pair})
    return G, df_edges_rbhs

def get_conserved_clusters(graphData, inParalogs, g, output_folder):
    """Returns the conserved clusters -- inference of transcripts homologies"""
    
    is_saving = False
    if(len(list(g.id_transcript.values))) <= 20:
       is_saving = True
    
    #Retrieve data
    G = graphData[0]
    edges_data = graphData[1]
    
    # save figures if number of transcripts is lower than 20
    if is_saving:
        nx.draw(G, with_labels=True)
        plt.savefig('{}/start_orthology_graph.pdf'.format(output_folder))
        #plt.show(G)


    edges_data = edges_data.sort_values(by=['score'], ascending=False)

    for index_edge, row_edge in edges_data.iterrows():
        transcript_id = row_edge['pair_id']
        transcript_ref = row_edge['pair_ref']

        adj_trID_with_weight = G.adj[transcript_id]
        adj_trID = [ _ for _ in list(adj_trID_with_weight.keys()) if _ != transcript_ref]
        adj_trREF_with_weight = G.adj[transcript_ref]
        adj_trREF = [ _ for _ in list(adj_trREF_with_weight.keys()) if _ != transcript_id]
        
        genes_adj_trID = list(set([g[g['id_transcript']==tr].id_gene.values[0] for tr in adj_trID]))
        genes_adj_trREF = list(set([g[g['id_transcript']==tr].id_gene.values[0] for tr in adj_trREF]))

        intersect_genes = list(set(genes_adj_trID).intersection(set(genes_adj_trREF)))

        if len(intersect_genes) == 0:
            G.add_edge(transcript_id, transcript_ref)
        else:
            for intersect_gene in intersect_genes:
                spes_adj_trID = [ _ for _ in adj_trID if g[g['id_transcript']==_].id_gene.values[0]==intersect_gene]
                spes_adj_trREF = [ _ for _ in adj_trREF if g[g['id_transcript']==_].id_gene.values[0]==intersect_gene]

                
                list_inParalogs_adj_trID = []
                inParalogs_adj_trID = [inParalogs[_] for _ in spes_adj_trID]
                for spe_adj_trID in spes_adj_trID:
                    inParalogs_adj_trID = inParalogs[spe_adj_trID]
                    for inparalog_adj_trID in inParalogs_adj_trID:
                        list_inParalogs_adj_trID.append(inparalog_adj_trID)
                
                list_inParalogs_adj_trREF = []
                inParalogs_adj_trREF = [inParalogs[_] for _ in spes_adj_trREF]
                for spe_adj_trREF in spes_adj_trREF:
                    inParalogs_adj_trREF = inParalogs[spe_adj_trREF]
                    for inparalog_adj_trREF in inParalogs_adj_trREF:
                        list_inParalogs_adj_trREF.append(inparalog_adj_trREF)

                list_inParalogs_adj_trID.extend(spes_adj_trID)
                list_inParalogs_adj_trREF.extend(spes_adj_trREF)
                set_id = set(list_inParalogs_adj_trID)
                set_ref =set(list_inParalogs_adj_trREF)

                inter_set = set_id.intersection(set_ref)

                if len(inter_set) == len(set_id) and len(inter_set) == len(set_ref):
                #if len(inter_set) == max(len(set_id), len(set_ref)):
                    G.add_edge(transcript_id, transcript_ref)
    if is_saving:
        plt.clf()
        nx.draw(G, with_labels=True)
        #plt.show(G)
        plt.savefig('{}/end_orthology_graph.pdf'.format(output_folder))
    clusters_components = [list(G.subgraph(c).copy()) for c in nx.connected_components(G)]
    return clusters_components

def write_results(clusters, df_gtot, matrix, gt_path, output_folder):
    '''Write the results'''

    transcripts = list(df_gtot['id_transcript'].values)
    tree = Tree(gt_path, format=1)
    df_orthology = pd.DataFrame(columns=['id_transcript_1','id_transcript_2','homology'])
    tsm_scores = []
    conserved_transcripts = []
    id_transcripts = []
    already_done = []
    for transcript in transcripts:
        id_transcripts.append(transcript)
        for cluster in clusters:
            if transcript in cluster:
                conserved_transcripts.append('|'.join([str(_) for _ in cluster if _ != transcript ]))
                tsm_scores.append('|'.join([str(matrix.loc[transcript,_]) for _ in cluster if _ != transcript]))
                for homolog in cluster:
                    if homolog != transcript and str(homolog)+'<=>'+str(transcript) not in already_done:
                        g1 = df_gtot[df_gtot.id_transcript==homolog].id_gene.values[0]
                        g2 = df_gtot[df_gtot.id_transcript==transcript].id_gene.values[0]
                        if g1 != g2:
                            node = tree.get_common_ancestor(g1, g2)
                            lca_node = ''
                            if hasattr(node, 'D'):
                                lca_node = list(node.D)[0]
                            elif hasattr(node, 'DD'):
                                lca_node = list(node.DD)[0]
                            if lca_node == 'Y':
                                df_orthology.loc[str(homolog)+'<=>'+str(transcript)] = [homolog, transcript, 'para-orthologs']
                            else:
                                df_orthology.loc[str(homolog)+'<=>'+str(transcript)] = [homolog, transcript, 'ortho-orthologs']
                        else:
                            df_orthology.loc[str(homolog)+str(transcript)] = [homolog, transcript, 'recent-paralogs']
                        already_done.append(str(homolog)+'<=>'+str(transcript))
                        already_done.append(str(transcript)+'<=>'+str(homolog))

  
    df = pd.DataFrame(data={'id_transcript':id_transcripts, 'isoorthologs|recent-paralogs': conserved_transcripts, 'tsm_scores': tsm_scores })
    df.to_csv('{}/groupsOfOrthologs.csv'.format(output_folder), sep=';', header=True)
    df_orthology.to_csv('{}/relationsOrthology.csv'.format(output_folder), sep=';', header=True)
    return df, df_orthology

if __name__ == '__main__':
    # retrieve inputs given by user
    args = build_arg_parser().parse_args()
    matrix_path = args.matrix
    gtot_path = args.genetotranscripts
    gt_path = args.nhxgenetree
    lower_bound = float(args.lowerbound)
    output_folder_path = args.outputfolder

    #compute the algorithm
    get_orthology_graph(matrix_path, gtot_path, gt_path, lower_bound, output_folder_path)
