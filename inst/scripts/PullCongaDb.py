import conga


def PullCongaDb(organism = 'human', 
                outputFilePath = 'conga_db.csv'): 
  tcrdist = conga.tcrdist.tcr_distances.TcrDistCalculator(organism)
  conga_db = conga.tcrdist.all_genes.df 
  conga_db = conga_db.loc[conga_db['organism'] == organism,:]
  #write conga_db df to file 
  conga_db.to_csv(outputFilePath, index=False)
