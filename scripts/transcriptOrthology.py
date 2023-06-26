""" Inferring transcript orthology.
> Usage:
======
    python3 transcriptOrthology.py [-h] [-talg transcripts alignment] [-gtot gene to transcripts mappings] [-nhxt NHX gene tree] [-lowb lower bound] [-outf output folder]

> Reference:
======
    https://github.com/UdeS-CoBIUS/TranscriptOrthology
"""
__authors__ = ("Wend Yam Donald Davy Ouedraogo")
__contact__ = ("wend.yam.donald.davy.usherbrooke.ca")
__copyright__ = "CoBIUS lab at Université de Sherbrooke, QC, CANADA"
__date__ = "2023-06-26"
__version__= "2.0.1"


import argparse
from Tclustering import get_orthology_graph
from tsmComputing import get_matrix


def build_arg_parser():
    '''Parsing function'''
    parser = argparse.ArgumentParser(description="program parameters")
    parser.add_argument('-talg', '--tralignment', default=None, help='Multiple Sequences Alignment of transcripts in FASTA format')
    parser.add_argument('-gtot', '--genetotranscripts', default=None, help="mappings transcripts to corresponding genes")
    parser.add_argument('-nhxt', '--nhxgenetree', default=None, help='NHX gene tree')
    parser.add_argument('-lowb', '--lowerbound', default=0.7, help='a lower bound for the selection of transcripts RBHs')
    parser.add_argument('-tsm', '--tsmvalue', default=1, help='an integer(1|2|3|4|5|6) that refers to the transcript similarity measure')
    parser.add_argument('-outf', '--outputfolder', default='.', help='the output folder to store the results')
    return parser

def inferring_transcripts_isoorthology(transcripts_msa_path, gtot_path, gt_path, tsm_conditions, lower_bound, output_folder):
    """inferring transcript isoorthologies"""
    try:
        tsm_matrix, df_blocks_transcripts, df_blocks_genes = get_matrix(transcripts_msa_path, gtot_path, tsm_conditions, output_folder)
    except:
        raise('Impossible to retrieve the matrix ! Errors occured ...')
    
    try:
        clusters, df, df_orthology = get_orthology_graph(tsm_matrix, gtot_path, gt_path, lower_bound, output_folder)
    except:
        raise('Impossible de retrieve transcripts orthologs! Errors occured ... ')
    
    
    return tsm_matrix, df_blocks_transcripts, df_blocks_genes, clusters, df, df_orthology

if __name__ == '__main__':
    # retrieve inputs given by user
    args = build_arg_parser().parse_args()
    gtot_path = args.genetotranscripts
    gt_path = args.nhxgenetree
    lower_bound = float(args.lowerbound)
    transcripts_msa_path = args.tralignment
    tsm_conditions = args.tsmvalue
    output_folder = args.outputfolder

    #compute the main algorithm
    tsm_matrix, df_blocks_transcripts, df_blocks_genes, clusters, df, df_orthology = inferring_transcripts_isoorthology(transcripts_msa_path, gtot_path, gt_path, tsm_conditions, lower_bound, output_folder)
    
    # finish
    print('++++++++++++++++Finished \n\n Succesful. Data results can be found in {}'.format(output_folder))

    
    


