
import sys
import os
import pandas as pd

# input download from smallworld
smallworld_csv = sys.argv[1]
# write output to CSV file
outcsv = open( sys.argv[2], 'w' )

# read in smallworld search results
df_sw = pd.read_csv( smallworld_csv, low_memory=False )

# get ZINC molids for top 10 matches from
analog_list = df_sw['molid'].head(10).tolist()

# write field headers to output CSV
outcsv.write('molid,vendors\n')

# for each ZINC molecule, get vendor data and write to output CSV
for zincid in analog_list:
    # use curl to download ZINC20 data table, save as 'temp_zinc'
    os.system( "curl -s https://zinc20.docking.org/substances/"+zincid+"/catitems/subsets/for-sale/table.html -o temp_zinc" )
    # read in temp_zinc and parse out the vendor info
    with open('temp_zinc', 'r') as fh:
        data = fh.readlines()
    # delete the temp_zinc file
    os.remove( 'temp_zinc' )
    # store vendor URLs in a list
    vendor_list = []
    for line in data:
        if 'vendor' in line:
            vendor_url = line.split('"')[1]
            #vendor_molid = line.split('"')[4].split('<')[0][1:]
            vendor_list.append( vendor_url )
    # write vendor URLs, use "|" as delimiter
    outcsv.write( zincid+","+"|".join( vendor_list )+"\n" )

outcsv.close()


