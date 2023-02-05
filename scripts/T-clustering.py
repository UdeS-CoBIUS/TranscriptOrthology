"""
author: Ouedraogo Wend Yam Donald Davy
email: wend.yam.donald.davy.ouedraogo@usherbrooke.ca
"""
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import warnings

def write_function(clusters, g, matrix, gt, path_output):
    # modifier la function plus tard
    #set_index
    col0 = matrix.columns[0]
    matrix = matrix.rename(columns={col0: 'transcripts'})
    matrix = matrix.set_index('transcripts')

    #initialize
    dict_clusters = {}
    transcripts = list(g['id_transcript'].values)
    print(transcripts)
    tsm_scores = []
    conserved_transcripts = []
    id_transcripts = []
    for transcript in transcripts:
        id_transcripts.append(transcript)
        for cluster in clusters:
            if transcript in cluster:
                conserved_transcripts.append(':'.join([str(_) for _ in cluster if _ != transcript ]))
                tsm_scores.append(':'.join([str(matrix.loc[transcript,_]) for _ in cluster if _ != transcript]))
    print(conserved_transcripts)
    print(tsm_scores)
    df = pd.DataFrame(data={'id_transcript':id_transcripts, 'orthoorthologues': conserved_transcripts, 'tsm_scores': tsm_scores })
    df.to_csv('{}/{}/{}_clustering.csv'.format(path_output, gt, gt), sep=';', header=True)
    return True

def get_conserved_clusters(G, matrix, inParalogs, g):
    """Returns the conserved clusters"""

    #set_index
    col0 = matrix.columns[0]
    matrix = matrix.rename(columns={col0: 'transcripts'})
    matrix = matrix.set_index('transcripts')


    nx.draw(G, with_labels=True)
    plt.show(G)

    components = [G.subgraph(c).copy() for c in nx.connected_components(G)]
    for component in components:
        if len(list(component)) > 2:
            # algorithme progressif for the clustering
            edges = list(component.edges())
            edges_sorted = sorted(d for d in component.edges())
            for edge_sorted in edges_sorted:
                # stay or it will be removed
                transcript_id = edge_sorted[0]
                transcript_ref = edge_sorted[1]
                adj_trID_with_weight = G.adj[transcript_id]
                adj_trID = [ _ for _ in list(adj_trID_with_weight.keys()) if _ != transcript_ref]
                adj_trREF_with_weight = G.adj[transcript_ref]
                adj_trREF = [ _ for _ in list(adj_trREF_with_weight.keys()) if _ != transcript_id]

                gene_id = g[g['id_transcript']==transcript_id].id_gene.values[0]
                gene_ref = g[g['id_transcript']==transcript_ref].id_gene.values[0]
                
                genes_adj_trID = list(set([g[g['id_transcript']==tr].id_gene.values[0] for tr in adj_trID]))
                genes_adj_trREF = list(set([g[g['id_transcript']==tr].id_gene.values[0] for tr in adj_trREF]))
                print(genes_adj_trREF, genes_adj_trID)
                intersect_genes = list(set(genes_adj_trID).intersection(set(genes_adj_trREF)))
                cut = False
                print(edge_sorted)
                if len(intersect_genes) != 0:
                    print('ok')
                    
                    for intersect_gene in intersect_genes:
                        if cut == False:
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
                            print(set_id)
                            print(set_ref)
                            if len(inter_set) == 0:
                            #if len(inter_set) == max(len(set_id), len(set_ref)):
                                print('oh NON!!')
                                #supprimer l'arrete de G
                                G.remove_edge(transcript_id, transcript_ref)
                                cut = True
                        else:
                            break
                elif gene_id in genes_adj_trREF:
                    spes_adj_trREF = [ _ for _ in adj_trREF if g[g['id_transcript']==_].id_gene.values[0]==genes_adj_trREF]
                    list_inParalogs_adj_trREF = []
                    inParalogs_adj_trREF = [inParalogs[_] for _ in spes_adj_trREF]
                    for spe_adj_trREF in spes_adj_trREF:
                        inParalogs_adj_trREF = inParalogs[spe_adj_trREF]
                        for inparalog_adj_trREF in inParalogs_adj_trREF:
                            list_inParalogs_adj_trREF.append(inparalog_adj_trREF)
                    list_inParalogs_adj_trREF.extend(spes_adj_trREF)
                    if transcript_id not in list_inParalogs_adj_trREF:
                        print('oh NON!!')
                        #supprimer l'arrete de G
                        G.remove_edge(transcript_id, transcript_ref)
                elif gene_ref in genes_adj_trID:
                    spes_adj_trID = [ _ for _ in adj_trID if g[g['id_transcript']==_].id_gene.values[0]==genes_adj_trID]
                    list_inParalogs_adj_trID = []
                    inParalogs_adj_trID = [inParalogs[_] for _ in spes_adj_trID]
                    for spe_adj_trID in spes_adj_trID:
                        inParalogs_adj_trID = inParalogs[spe_adj_trID]
                        for inparalog_adj_trID in inParalogs_adj_trID:
                            list_inParalogs_adj_trID.append(inparalog_adj_trID)
                    
                    list_inParalogs_adj_trID.extend(spes_adj_trID)
                    if transcript_ref not in list_inParalogs_adj_trID:
                        print('oh NON!!')
                        #supprimer l'arrete de G
                        G.remove_edge(transcript_id, transcript_ref)

                else:
                    pass
    nx.draw(G, with_labels=True)
    plt.show(G)
    clusters_components = [list(G.subgraph(c).copy()) for c in nx.connected_components(G)]
    return clusters_components

def rhb_clustering(matrix, g, lower_bound):
    """Returns a similarity graph G"""

    #warning if lower_bound is greater than 1
    if lower_bound > 1:
        warnings.warn('Each transcript is considered as a cluster. Check the lower bound! Enter a lower bound less greater than 1')
    
    #set_index
    col0 = matrix.columns[0]
    matrix = matrix.rename(columns={col0: 'transcripts'})
    matrix = matrix.set_index('transcripts')

    #initialize the graph G
    G = nx.Graph()
    G.add_nodes_from(matrix.index)
    

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
                                # true_rbhs.append(tr_ortholog)
                                # a link with tr_ortholog and transcript
                                G.add_edge(tr_ortholog, transcript, weight=max_value_tr_ortholog)
    return G


def get_conserved_clusters_v2(graphData, inParalogs, g):
    """Returns the conserved clusters -- inference of transcripts homologies"""

    #Retrieve data
    G = graphData[0]
    edges_data = graphData[1]

    #nx.draw(G, with_labels=True)
    #plt.show(G)

    #sort edges by descending
    #print(edges_data)
    edges_data = edges_data.sort_values(by=['score'], ascending=False)
    #print(edges_data)

    for index_edge, row_edge in edges_data.iterrows():
        transcript_id = row_edge['pair_id']
        transcript_ref = row_edge['pair_ref']
        print(transcript_id, transcript_ref)
        adj_trID_with_weight = G.adj[transcript_id]
        adj_trID = [ _ for _ in list(adj_trID_with_weight.keys()) if _ != transcript_ref]
        adj_trREF_with_weight = G.adj[transcript_ref]
        adj_trREF = [ _ for _ in list(adj_trREF_with_weight.keys()) if _ != transcript_id]
        print(adj_trID_with_weight)
        
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
                else:
                    print('oh NON!!')

    #nx.draw(G, with_labels=True)
    #plt.show(G)
    clusters_components = [list(G.subgraph(c).copy()) for c in nx.connected_components(G)]
    return clusters_components

def rhb_clustering_v2(matrix, g, inparalogs, lower_bound):
    """Returns a similarity graph G containing recent paralogs and the recent paralogs details"""

    #warning if lower_bound is greater than 1
    if lower_bound > 1:
        warnings.warn('Each transcript is considered as a cluster. Check the lower bound! Enter a lower bound less greater than 1')
    
    #set_index
    col0 = matrix.columns[0]
    matrix = matrix.rename(columns={col0: 'transcripts'})
    matrix = matrix.set_index('transcripts')

    #initialize the graph G
    edges_rbhs = []
    pair_id = []
    pair_ref = []
    score_pair = []
    G = nx.Graph()
    G.add_nodes_from(matrix.index)
    for tr_recent in list(inparalogs.keys()):
        print(tr_recent)
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
                                # true_rbhs.append(tr_ortholog)
                                # a link with tr_ortholog and transcript
                                #G.add_edge(tr_ortholog, transcript, weight=max_value_tr_ortholog)
                                
                                if set([tr_ortholog, transcript]) not in edges_rbhs:
                                    #print(transcript, tr_ortholog)
                                    pair_id.append(transcript)
                                    pair_ref.append(tr_ortholog)
                                    score_pair.append(max_value_tr_ortholog)
                                    edges_rbhs.append(set([tr_ortholog, transcript]))
                                
    df_edges_rbhs = pd.DataFrame(data={'pair_id': pair_id, 'pair_ref': pair_ref, 'score': score_pair})
    return [G, df_edges_rbhs]

def find_internal_paralogs(matrix, g):
    """Returns a map of internal paralogs"""
    #set_index
    col0 = matrix.columns[0]
    matrix = matrix.rename(columns={col0: 'transcripts'})
    matrix = matrix.set_index('transcripts')

    #initialiaze
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

 

if __name__ == '__main__':
    file = open('./macse_253/trees_macse.txt','r')
    gts = [str(_).split('\n')[0] for _ in file.readlines()]
    #path_output = './results'

    paths_conf_matrix = ['Confusion_matrix/tsm/Macse/m/','Confusion_matrix/tsm/Macse/degHom', 'Confusion_matrix/tsm/Macse/nucHom', 'Confusion_matrix/tsm_corrigé/Macse/m/','Confusion_matrix/tsm_corrigé/Macse/degHom', 'Confusion_matrix/tsm_corrigé/Macse/nucHom']
    #print(gts)
    for path_output in paths_conf_matrix:
        for gt in gts:
            #if gt == 'ENSGT00390000000583':
            file_matrix = '{}/{}/{}_matrix.csv'.format(path_output, gt, gt)
            file_gtot = '{}/{}/{}_gtot.csv'.format(path_output, gt, gt)
            matrix = pd.read_csv(file_matrix, sep=';')
            gtot = pd.read_csv(file_gtot, sep=';')

            inparalogs = find_internal_paralogs(matrix, gtot)
            print(inparalogs)
            lower_bound = 0.1
            graphG = rhb_clustering_v2(matrix, gtot, inparalogs, lower_bound)
            print(graphG[1])
            print(graphG[0])
            clusters = get_conserved_clusters_v2(graphG, inparalogs, gtot)
            write_function(clusters, gtot, matrix, gt, path_output)
        