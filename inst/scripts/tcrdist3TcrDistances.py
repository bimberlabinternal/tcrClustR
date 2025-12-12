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
        print("debug: DataFrame shape:\n", df.shape)
        print("debug: DataFrame columns:\n", df.columns)

    #store original data for validation after TCRrep processing
    if 'CloneId' not in df.columns:
        raise ValueError("Input CSV must contain a 'CloneId' column for order validation")
    original_clone_ids = df['CloneId'].tolist()
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
    
    if debug:
        print(f"debug: Original CloneId order (first 5): {original_clone_ids[:5]}")
        print(f"debug: Input dataframe length: {len(df)}")

    print("info: Computing distances for chains: ", chains)
    tr = TCRrep(cell_df = df, 
                organism = organism, 
                chains = chains, 
                db_file = db_file)
    
    #cell_df is the input after TCRrep receives it - should be identical to our input
    cell_df = getattr(tr, 'cell_df', None)
    if cell_df is None:
        raise ValueError("TCRrep object missing cell_df attribute")

    #get clone_df which contains grouped/collapsed clones
    clone_df = getattr(tr, 'clone_df', None)
    if clone_df is None:
        raise ValueError("TCRrep object missing clone_df attribute - cannot validate clone order")

    #validate CloneId to guard against reordering across the three dataframes
    #we deduplicate in df upstream, so all three dataframes should have the same CloneId order
    #CloneId must match exactly between the input df, the TCRrep cell_df, and TCRrep clone_df
    #Note: In this pipeline, we expect 1:1 mapping (no collapsing), so clone_df must also match df
    
    # first, check df vs cell_df
    if not df['CloneId'].equals(cell_df['CloneId']):
        if df['CloneId'].tolist() != cell_df['CloneId'].tolist():
            raise ValueError(f"CloneId mismatch between input df and TCRrep cell_df! \nInput head: {df['CloneId'].head().tolist()}\nCell_df head: {cell_df['CloneId'].head().tolist()}")
        print("warning: CloneId Series found unequal (likely index mismatch) but values/order are identical between df and cell_df.")

    # second, check df vs clone_df
    if 'CloneId' not in clone_df.columns:
         raise ValueError(f"TCRrep clone_df missing CloneId column. Available columns: {clone_df.columns.tolist()}")

    if not df['CloneId'].equals(clone_df['CloneId']):
        if df['CloneId'].tolist() != clone_df['CloneId'].tolist():
             raise ValueError(f"CloneId mismatch between input df and TCRrep clone_df! \nInput head: {df['CloneId'].head().tolist()}\nClone_df head: {clone_df['CloneId'].head().tolist()}")
        print("warning: CloneId Series found unequal (likely index mismatch) but values/order are identical between df and clone_df.")
    
    if debug:
        print(f"debug: cell_df length: {len(cell_df)}")
        print(f"debug: clone_df length: {len(clone_df)}")
        print(f"debug: clone_df columns: {list(clone_df.columns)}")
    
    clone_df_ids = clone_df['CloneId'].tolist()

    if debug:
        print(f"debug: clone_df CloneId order (first 5): {clone_df_ids[:5]}")
    
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
    
    #look-up table for chain to attribute name mapping for CDR3-only distances
    chain_cdr3_mapping = {
        'alpha': 'pw_cdr3_a_aa',
        'beta': 'pw_cdr3_b_aa', 
        'gamma': 'pw_cdr3_g_aa',
        'delta': 'pw_cdr3_d_aa'
    }
    
    for chain in chains:
        #full-length distance matrix
        attr_name = f'pw_{chain}'
        if hasattr(tr, attr_name):
            mat = getattr(tr, attr_name)
            #secondary validation: validate that full length dimensions match clone_df
            if len(clone_df) != len(mat):
                raise ValueError(
                    f"Matrix {attr_name} row length did not match clone_df! "
                    f"clone_df: {len(clone_df)}, matrix rows: {len(mat)}"
                )
            if len(clone_df) != len(mat[0]):
                raise ValueError(
                    f"Matrix {attr_name} column length did not match clone_df! "
                    f"clone_df: {len(clone_df)}, matrix cols: {len(mat[0])}"
                )
            result[attr_name] = mat
            if debug:
                print(f"validation: {attr_name} dimensions validated, {len(mat)} x {len(mat[0])}")
        
        #CDR3-only distance matrix
        if chain in chain_cdr3_mapping:
            cdr3_attr_name = chain_cdr3_mapping[chain]
            if hasattr(tr, cdr3_attr_name):
                mat = getattr(tr, cdr3_attr_name)
                #secondary validation: validate that cdr3 dimensions match clone_df
                if len(clone_df) != len(mat):
                    raise ValueError(
                        f"Matrix {cdr3_attr_name} row length did not match clone_df! "
                        f"clone_df: {len(clone_df)}, matrix rows: {len(mat)}"
                    )
                if len(clone_df) != len(mat[0]):
                    raise ValueError(
                        f"Matrix {cdr3_attr_name} column length did not match clone_df! "
                        f"clone_df: {len(clone_df)}, matrix cols: {len(mat[0])}"
                    )
                result[cdr3_attr_name] = mat
                if debug:
                    print(f"validation: {cdr3_attr_name} dimensions validated, {len(mat)} x {len(mat[0])}")
    
    #return the clone_ids alongside matrices so dimnames can be set correctly
    result['clone_ids'] = clone_df_ids
    
    if debug:
        print(f"\nDebug Summary:")
        print(f"debug:   Input df length: {len(df)}")
        print(f"debug:   Output cell_df length: {len(cell_df)}")
        print(f"debug:   Output clone_df length: {len(clone_df)}")
        print(f"debug:   Matrices returned: {[k for k in result.keys() if k != 'clone_ids']}")
    
    return result

def writeTcrDistances(csv_path, 
                      organism='human', 
                      chainsString="alpha,beta", 
                      db_file='combo_xcr_2024-03-05',
                      rds_output_path='./tcrdist3DistanceMatrices/', 
                      debug=True):
    """
    Compute TCR distance matrices using tcrdist3 and save them as RDS files for use in R.
    
    The RDS files will have dimnames set to the CloneId values from clone_df,
    ensuring proper row/column identification when loaded in R.
    """
    #instantiate output directory if it doesn't exist.
    os.makedirs(rds_output_path, exist_ok=True)

    #compute distances (returns dict with matrices + 'clone_ids' key)
    distances = getTcrDistances(csv_path, organism, chainsString, db_file, debug)
    
    #extract clone_ids (these come from clone_df, already validated to match matrix dimensions)
    clone_ids = distances.pop('clone_ids')
    
    if debug:
        print(f"\ndebug: Saving {len(clone_ids)} clone IDs from clone_df")

    base = importr('base')
    import rpy2.robjects as ro
    
    #context manager to handle numpy to R conversion
    with localconverter(default_converter + numpy2ri.converter):
        #save all available distance matrices
        for matrix_name, matrix_data in distances.items():
            output_file = os.path.join(rds_output_path, f'{matrix_name}.rds')
            
            base.saveRDS(matrix_data, output_file)
            if debug:
                print(f"debug: Saved {matrix_name} to {output_file}")
        
        #save clone_ids as separate RDS file for R to set dimnames in R
        clone_ids_file = os.path.join(rds_output_path, 'clone_ids.rds')
        r_clone_ids = ro.StrVector([str(cid) for cid in clone_ids])
        base.saveRDS(r_clone_ids, clone_ids_file)
        if debug:
            print(f"debug: Saved clone_ids to {clone_ids_file}")
   
#this file serves as a template. RunTcrdist3.R will copy this function, add a string with arguments to the end, and then call the whole file.

