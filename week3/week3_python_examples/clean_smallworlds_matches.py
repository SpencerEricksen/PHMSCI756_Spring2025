
import sys
import pandas as pd

# input/output file arguments
intsv = sys.argv[1]
outsmi = sys.argv[2]

outcsv = intsv.split('.')[0] + '.csv'

# read in tab-delimited data file downloaded from SmallWorld query
#df = pd.read_csv('smallworld-results-1680712496020_CID25181377.tsv', delim_whitespace=True )
df = pd.read_csv( intsv, delim_whitespace=True )

# clean up fields in dataframe
df = df.reset_index()
df.rename( columns={'index':'sw_smiles','alignment':'molid'}, inplace=True )
df.head(10).to_csv( outcsv, index=False )

# dump SMILES for 10 closest matches
df2 = df[['sw_smiles','molid']].head(10)
#df2.to_csv( 'CID25181377_sw_matches.smi', index=False, header=False, sep=" ")
df2.to_csv( outsmi, index=False, header=False, sep=" ")

