import inspect
import importlib.util
import os
import pandas as pd
from tcrdist.repertoire import TCRrep

def PullTcrdist3Db(organism = 'human', 
                   outputFilePath = 'tcrdist3_gene_segments.txt'):
  """
  Pull the default database file from TCRdist3.
  """
  #find the tcrdist package
  print("Reading tcrdist3 database file...")
  spec = importlib.util.find_spec("tcrdist")
  if spec and spec.submodule_search_locations:
      module = importlib.import_module("tcrdist.repertoire")
      signature = inspect.signature(TCRrep.__init__)
      tcrdist_path = spec.submodule_search_locations[0]
      default_db_file = signature.parameters["db_file"].default
      db_file_path = os.path.join(tcrdist_path, "db", default_db_file)

      #get the db file 
      if os.path.exists(db_file_path):
          df = pd.read_csv(db_file_path, sep="\t")
          #subset to organism specific TCRs
          df = df.loc[df['organism'] == organism,:]

          #gene segments
          gene_segments = df["id"].dropna().unique().tolist()

          #write gene segments
          with open(outputFilePath, "w") as f:
              f.write("gene_segments\n" + "\n".join(gene_segments) + "\n\n")
          print(f"Gene segments written to {outputFilePath}")
      else:
          print("DB file not found:", db_file_path)
  else:
      print("tcrdist package not found.")



