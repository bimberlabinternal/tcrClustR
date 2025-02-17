import conga 
import csv
import sys
import ast
import numpy as np
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri

temp_tcr_file = '/Users/mcelfreg/Library/CloudStorage/OneDrive-OregonHealth&ScienceUniversity/scRNASeq_OneDrive/tcrClustR/congaInput.csv'

#this function, unlike tcrdist3TcrDistances.py, only computes a single distance matrix
def getTcrDistances(tcrFile, 
  organism = 'human', 
  chain = "TRA", 
  rds_output_path='./congaDistanceMatrices/'):
  tcrdist = conga.tcrdist.tcr_distances.TcrDistCalculator(organism)

  with open(tcrFile) as file:
    content = file.readlines()

  #import the tcrs as a list of dictionaries
  tcrs = ast.literal_eval(str(content[1]))

  #compute distances
  D = np.array([tcrdist.single_chain_distance(t1,t2) for t1 in tcrs for t2 in tcrs]).reshape((len(tcrs),len(tcrs)))
  
  #convert and save distance matrix to an RDS file for the R-based clustering steps
  numpy2ri.activate()
  base = importr('base')
  base.saveRDS(D, rds_output_path + '/congaTcrDistances_' + chain + '.rds')


#this file serves as a template. congaWrapper.R will copy this function, add a string with arguments to the end, and then call the whole file. 

