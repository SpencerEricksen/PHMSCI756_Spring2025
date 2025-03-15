
import sys
import pandas as pd

vendors_csv = sys.argv[1]	# 'CID2790671_analogs_vendors.csv'
smallworld_csv = sys.argv[2]	# 'smallworld-results-1680712735959_CID2790671.csv'
swissadme_csv = sys.argv[3]     # 'swissadme_CID2790671_sw_matches.csv'
outcsv = sys.argv[4]		# 'CID2790671_analogs_merged_data.csv'

df1 = pd.read_csv( vendors_csv )
df2 = pd.read_csv( smallworld_csv )
df3 = pd.read_csv( swissadme_csv )
df3.rename( columns={'Molecule':'molid'}, inplace=True )

# merge
df = df2.merge( df1, on='molid' )
df = df.merge( df3, on='molid')
# keep columns for report
df_final = df[['molid','sw_smiles','ecfp4','Lipinski #violations','PAINS #alerts','BBB permeant','vendors']]

# dump
df_final.to_csv( outcsv, index=False )

