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
                    db_file='alphabeta_gammadelta_db.tsv', 
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
    if tr.pw_beta.shape[0] > 10000:
       TCRrep.compute_distances()
    return {
        'pw_alpha': tr.pw_alpha,
        'pw_beta': tr.pw_beta,
        'pw_cdr3_a_aa': tr.pw_cdr3_a_aa,
        'pw_cdr3_b_aa': tr.pw_cdr3_b_aa
    }

def writeTcrDistances(csv_path, 
                      organism='human', 
                      chainsString="alpha,beta", 
                      db_file='alphabeta_gammadelta_db.tsv',
                      rds_output_path='./tcrdist3DistanceMatrices/', 
                      debug=True):
    """
    Compute TCR distance matrices using tcrdist3 and save them as RDS files for use in R.
    """
    #instantiate output directory if it doesn't exist.
    os.makedirs(rds_output_path, exist_ok=True)

    #compute distances
    distances = getTcrDistances(csv_path, organism, chainsString, db_file, debug)

    base = importr('base')
    #context manager to handle numpy <-> R conversion
    with localconverter(default_converter + numpy2ri.converter):
        base.saveRDS(distances['pw_alpha'], os.path.join(rds_output_path, 'pw_alpha.rds'))
        base.saveRDS(distances['pw_beta'], os.path.join(rds_output_path, 'pw_beta.rds'))
        base.saveRDS(distances['pw_cdr3_a_aa'], os.path.join(rds_output_path, 'pw_cdr3_a_aa.rds')
        base.saveRDS(distances['pw_cdr3_b_aa'], os.path.join(rds_output_path, 'pw_cdr3_b_aa.rds'))
   
#this file serves as a template. tcrDist3Wrapper.R will copy this function, add a string with arguments to the end, and then call the whole file. 

