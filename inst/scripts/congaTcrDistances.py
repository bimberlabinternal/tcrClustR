#import conga 
import csv
import sys
import ast
import numpy as np
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri

def getTcrDistances(tcrFile, organism = 'human', rdsOutputPath = 'TCR_distance_matrix.rds'):
  tcrdist = conga.tcrdist.tcr_distances.TcrDistCalculator(organism)

  with open('tcrs_for_conga.csv') as file:
    content = file.readlines()

  #import the tcrs as a list of dictionaries
  tcrs = ast.literal_eval(str(content[1]))
  #compute distances
  D = np.array([tcrdist.single_chain_distance(t1,t2) for t1 in tcrs for t2 in tcrs]).reshape((len(tcrs),len(tcrs)))
  
  #convert and save distance matrix to an RDS file for the R-based clustering steps
  numpy2ri.activate()
  base = importr('base')
  base.saveRDS(D, file="TCR_distance_matrix.rds")


