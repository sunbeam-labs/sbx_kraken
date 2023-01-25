biom_to_tsv_skbio(snakemake.input[0], snakemake.output[0])

# Function to compile a tsv report from multiple kraken2 reports



# Function to convert biom files to tsv files.
import skbio as skb

#Input: biom file
#Output: tsv file with the following columns: OTU ID, Kingdom, Phylum, Class, Order, Family, Genus, Species
def biom_to_tsv_skbio(biom_file, tsv_file):

    table = skb.Table.read(biom_file)

    table.to_tsv(tsv_file)


#import pandas as pd
#import numpy as np
#from biom import load_table, Table

def biom_to_tsv_biom(biom_file):
    table = load_table(biom_file)

    # Convert to pandas dataframe.  This is the easiest way to get the data into a tsv file.  You can also use the write command from biom library but it is not as flexible.  See https://biom-format.org/documentation/format_versions/biom-2.0.html for more info on how to use write command in biom library.

    df = table.to_dataframe()

    # Transpose dataframe so that samples are columns and features are rows (this is how qiime2 likes it).  The transpose function in pandas does not work with sparse matrices, so we need to convert to dense matrix first using .toarray() function and then transpose using .T function (see https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.T.html for more info).  

    df = df.toarray().T

    # Get sample IDs from table object and convert to list of strings (the sample IDs are stored as a numpy array of bytes, so we need to convert them back into strings).  We will use these sample IDs later when we save the dataframe as a tsv file (see below).  

    samples = list(map(str, table.ids('sample')))

    # Get feature IDs from table object and convert them into a list of strings (the feature IDs are stored as a numpy array of bytes, so we need to convert them back into strings).  We will use these feature IDs later when we save the dataframe as a tsv file (see below).  

    features = list(map(str, table.ids('observation')))

    # Create new column names for our dataframe by appending 'X' at beginning of each sample ID string in samples list above and adding '.1' at end of each string in features list above (this is how qiime2 likes it).  

    #newcol