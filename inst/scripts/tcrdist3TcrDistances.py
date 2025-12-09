import sys
import ast
import re
import os
import numpy as np
import pandas as pd
from tcrdist.repertoire import TCRrep
from rpy2.robjects import default_converter
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import numpy2ri
from rpy2.robjects.packages import importr

def getTcrDistances(csv_path, 
                    organism='human', 
                    chainsString="alpha,beta", 
                    db_file='combo_xcr_2024-03-05',
                    debug = True):
    #read in the csv file
    df = pd.read_csv(csv_path)
    
    #convert debug strings from R into python booleans
    if debug == 'TRUE' or debug == 'true' or debug == 'True' or debug == 'T' or debug == 't' or debug == '1':
        debug = True
    else:
        debug = False
    
    if debug:
        print("DataFrame shape:\n", df.shape)
        print("DataFrame columns:\n", df.columns)

    #regex the chainsString argument to get the chains
    chains = []
    if re.search(r'alpha', chainsString):
        chains.append('alpha')
    if re.search(r'beta', chainsString):
        chains.append('beta')
    if re.search(r'gamma', chainsString):
        chains.append('gamma')
    if re.search(r'delta', chainsString):
        chains.append('delta')

    print("Computing distances for chains: ", chains)
    tr = TCRrep(cell_df = df, 
                organism = organism, 
                chains = chains, 
                db_file = db_file)
    #TODO: there seems to be an internal swap at n > 10,000 to returning pairwise distances as "rw"s instead 
    #of "pw"s, which we may want to use? 
    '''
    When TCRrep.<clone_df> size 18779 > 10,000.
        TCRrep.compute_distances() may be called explicitly by a user
        with knowledge of system memory availability.
        However, it's HIGHLY unlikely that you want to compute such
        a large numpy array. INSTEAD, if you want all pairwise distance,
        you will likely want to set an appropriate number of cpus with TCRrep.cpus = x,
        and then generate a scipy.sparse csr matrix of distances with:
        TCRrep.compute_sparse_rect_distances(radius=50, chunk_size=100), leaving df and df2 arguments blank.
        When you do this the results will be stored as TCRrep.rw_beta instead of TCRrep.pw_beta.
        This function is highly useful for comparing a smaller number of sequences against a bulk set
        In such a case, you can specify df and df2 arguments to create a non-square matrix of distances.
        See https://tcrdist3.readthedocs.io/en/latest/sparsity.html?highlight=sparse for more info.
    '''
    #TODO part two: for now, we'll just bypass the warning and compute the pairwise distances in a dense array.
    #check that the distances are fully pairwise 'pw' and not 'rw's, error out otherwise
    sample_attr = None
    for chain in chains:
        attr_name = f'pw_{chain}'
        if hasattr(tr, attr_name):
            sample_attr = getattr(tr, attr_name)
            break
    
    if sample_attr is not None and sample_attr.shape[0] > 10000:
       tr.compute_distances()
    
    #build the expected return dict dynamically based on available chains and distance matrices
    result = {}
    for chain in chains:
        attr_name = f'pw_{chain}'
        if hasattr(tr, attr_name):
            result[attr_name] = getattr(tr, attr_name)
    
    #look-up table for chain to attribute name mapping for CDR3-only distances
    chain_mapping = {
        'alpha': 'pw_cdr3_a_aa',
        'beta': 'pw_cdr3_b_aa', 
        'gamma': 'pw_cdr3_g_aa',
        'delta': 'pw_cdr3_d_aa'
    }
    #add the CDR3-only amino acid distance matrices to the result dict
    for chain in chains:
        if chain in chain_mapping:
            attr_name = chain_mapping[chain]
            if hasattr(tr, attr_name):
                result[attr_name] = getattr(tr, attr_name)
    
    return result

def writeTcrDistances(csv_path, 
                      organism='human', 
                      chainsString="alpha,beta", 
                      db_file='combo_xcr_2024-03-05',
                      rds_output_path='./tcrdist3DistanceMatrices/', 
                      debug=True):
    """
    Compute TCR distance matrices using tcrdist3 and save them as RDS files for use in R.
    Matrices saved with row/colnames matching cell_df CloneId order.
    """
    os.makedirs(rds_output_path, exist_ok=True)
    
    # Read CloneIds to preserve as row/colnames (cell_df order)
    clone_ids = pd.read_csv(csv_path)['CloneId'].astype(str).tolist()
    
    distances = getTcrDistances(csv_path, organism, chainsString, db_file, debug)

    base = importr('base')
    from rpy2 import robjects
    with localconverter(default_converter + numpy2ri.converter):
        for matrix_name, matrix_data in distances.items():
            output_file = os.path.join(rds_output_path, f'{matrix_name}.rds')
            r_matrix = numpy2ri.py2rpy(matrix_data)
            r_names = robjects.StrVector(clone_ids)
            #set dimnames(matrix) <- list(cloneIds, cloneIds)
            r_matrix = robjects.r['`dimnames<-`'](r_matrix, robjects.ListVector({'row': r_names, 'col': r_names}))
            base.saveRDS(r_matrix, output_file)
            if debug:
                print(f"Saved {matrix_name} to {output_file}")
   
#this file serves as a template. RunTcrdist3.R will copy this function, add a string with arguments to the end, and then call the whole file.

