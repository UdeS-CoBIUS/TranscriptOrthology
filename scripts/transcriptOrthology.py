""" Inferring transcript orthology.
> Usage:
======
    python3 transcriptOrthology.py [-h] [-talg transcripts alignment] [-galg genes alignment] [-gtot gene to transcripts mappings] [-nhxt NHX gene tree] [-lowb lower bound] [-outf output folder]

> Reference:
======
    https://github.com/UdeS-CoBIUS/TranscriptOrthology
"""
__authors__ = ("Wend Yam Donald Davy Ouedraogo")
__contact__ = ("wend.yam.donald.davy.usherbrooke.ca")
__copyright__ = "CoBIUS lab at Université de Sherbrooke, QC, CANADA"
__date__ = "2022-12-19"
__version__= "1.0.0"


import argparse

from Tclustering import get_orthology_graph
from tsmComputing import get_matrix


def build_arg_parser():
    '''Parsing function'''
    parser = argparse.ArgumentParser(description="parsor program parameter")
    parser.add_argument('-talg', '--tralignment', default=None)
    parser.add_argument('-galg', '--genesalignment', default=None)
    parser.add_argument('-gtot', '--genetotranscripts', default=None)
    parser.add_argument('-nhxt', '--nhxgenetree', default=None)
    parser.add_argument('-lowb', '--lowerbound', default=0.5)
    parser.add_argument('-tsm', '--tsmuncorrected', default=False)
    parser.add_argument('-outf', '--outputfolder', default='.')
    return parser

def inferring_transcripts_isoortholoy(transcripts_msa_path, genes_msa_path,gtot_path, gt_path, tsm_conditions, output_folder_path, lower_bound):
    """inferring transcript isoorthologies"""
    try:
        get_matrix(transcripts_msa_path, genes_msa_path, gtot_path, tsm_conditions, output_folder_path)
        matrix_path = '{}/{}.csv'.format(output_folder_path, 'matrix')
    except:
        raise('Impossible to retrieve the matrix ! Errors occured ...')
    
    try:
        get_orthology_graph(matrix_path, gtot_path, gt_path, lower_bound, output_folder_path)
    except:
        raise('Impossible de retrieve transcripts orthologs! Errors occured ... ')
    
    print('++++++++++++++++Finished \n\n Succesful. Data results can be found in {}'.format(output_folder_path))
    return True

if __name__ == '__main__':
    # retrieve inputs given by user
    args = build_arg_parser().parse_args()
    gtot_path = args.genetotranscripts
    gt_path = args.nhxgenetree
    lower_bound = float(args.lowerbound)
    transcripts_msa_path = args.tralignment
    genes_msa_path = args.genesalignment
    tsm_conditions = args.tsmuncorrected
    output_folder_path = args.outputfolder

    #compute the main algorithm
    inferring_transcripts_isoortholoy(transcripts_msa_path, genes_msa_path,gtot_path, gt_path, tsm_conditions, output_folder_path, lower_bound)
    
    