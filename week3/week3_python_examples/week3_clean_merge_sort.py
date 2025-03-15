
import pandas as pd

#  AID_1558_datatable_week3_exercise.csv
#  AID_1559_datatable_week3_exercise.csv
#  AID_1694_datatable_week3_exercise.csv  (counter)

df1 = pd.read_csv('AID_1558_datatable_week3_exercise.csv', low_memory=False )
df2 = pd.read_csv('AID_1559_datatable_week3_exercise.csv', low_memory=False )
df3 = pd.read_csv('AID_1694_datatable_week3_exercise.csv', low_memory=False )

def clean_data( df, aid ):
    # remove extra header lines (rows 0:4)
    df = df.iloc[5:]
    # make relevant columns numeric "floats" as their datatypes
    df['PUBCHEM_CID'] = df['PUBCHEM_CID'].fillna(0).astype('int')
    df['Fit_CurveClass'] = df['Fit_CurveClass'].astype('float')
    df['Potency'] = df['Potency'].astype('float')
    df['Max_Response'] = df['Max_Response'].astype('float')
    # set activators with positive curve class to -10.0
    df.loc[ df.Fit_CurveClass > 0.0, 'Fit_CurveClass'] = -10.0
    # prioritize samples with multi-tiered sort
    df = df.sort_values( by=['Fit_CurveClass','Max_Response','Potency'], ascending=[False,True,True] )
    df = df[['PUBCHEM_CID', 'PUBCHEM_EXT_DATASOURCE_SMILES', 'Potency', 'Efficacy', 'Fit_HillSlope',
             'Fit_R2', 'Fit_InfiniteActivity', 'Fit_ZeroActivity', 'Fit_CurveClass', 'Max_Response']]
    # add rank score (best=0, worst=138)
    df['rank'] = df.reset_index().index
    # add suffix with AID number
    new_col_names = [ c if "PUBCHEM" in c else c+'_'+str(aid) for c in df.columns.tolist() ]
    df.columns = new_col_names
    df.drop_duplicates( subset='PUBCHEM_CID', keep='first', inplace=True )
    return df

# clean 1558
df_1558 = clean_data( df1, 1558 )

# clean 1559
df_1559 = clean_data( df2, 1559 )
df_1559.drop( columns=['PUBCHEM_EXT_DATASOURCE_SMILES'], inplace=True )

# clean_1694
df_1694 = clean_data( df3, 1694 )
# make new column to label molecules as fluorescence quenchers
df_1694['quencher'] = 1
# only CurveClass=4 will be considered non-quenchers (quencher=0)
df_1694.loc[ df_1694['Fit_CurveClass_1694'] == -10.0, 'quencher'] = 0
# from 1694 keep only the quencher label and PUBCHEM_CID for this counter-screen data
df_1694 = df_1694[['PUBCHEM_CID','quencher']]


# merge 1558 and 1559 data
df = df_1558.merge( df_1559, on='PUBCHEM_CID' )
# add counter screen data
df = df.merge( df_1694, on='PUBCHEM_CID' )


# make sum of ranks score (based on 1558 and 1559 results)
df['rank_1558_1559'] = df[['rank_1558','rank_1559']].sum(axis=1)
df = df.sort_values( by=['quencher','rank_1558_1559'], ascending=[True,True] )


# dump merged, sorted data set
df.to_csv('merged_retest_data.csv', index=False )
  
# dump top 40 cpd in SMILES format
df[['PUBCHEM_EXT_DATASOURCE_SMILES','PUBCHEM_CID']].head(40).to_csv( 'top40_retest_cpds_sorted.smi', sep=" ", header=None, index=False)

