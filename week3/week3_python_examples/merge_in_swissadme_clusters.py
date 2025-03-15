
import pandas as pd

# read in assay data, clustering data, and swissADME data
df1 = pd.read_csv('merged_retest_data.csv')
df2 = pd.read_csv('chemmine_clustering_top40.csv')
df3 = pd.read_csv('top40_swissadme.csv')

# change molids to "PUBCHEM_CID"
df2.rename( columns={'ids':'PUBCHEM_CID'}, inplace=True )
df3.rename( columns={'Molecule':'PUBCHEM_CID'}, inplace=True )

# add the clustering results to merged assay data
df4 = df1.merge( df2, on='PUBCHEM_CID', how='left' )
df5 = df4.merge( df3, on='PUBCHEM_CID', how='left' )

# keep just the important fields in merged data file
df6 = df5[['PUBCHEM_CID','PUBCHEM_EXT_DATASOURCE_SMILES', 
           'Potency_1558', 'Efficacy_1558', 'Fit_CurveClass_1558', 'Max_Response_1558', 'rank_1558', 
           'Potency_1559', 'Efficacy_1559', 'Fit_CurveClass_1559', 'Max_Response_1559', 'rank_1559', 
           'quencher', 'rank_1558_1559', 'CLSZ_0.4', 'CLID_0.4','Lipinski #violations','PAINS #alerts','BBB permeant']]

# output merged dataframe as CSV
df6.to_csv('top40_merged_retest_clust_adme.csv', index=False )

