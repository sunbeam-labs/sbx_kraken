biom_to_tsv_skbio(snakemake.input[0], snakemake.output[0])

# Create BIOM-format tables from Kraken output
# Input:
#    1) A file containing the output from the kraken program.
#    2) A file containing a list of taxon names and their taxonomic lineage.
# Output:
#    A BIOM-format table containing the abundance of each taxon in each sample.
import sys, re, os, argparse, json, numpy as np
from collections import defaultdict


def parse_kraken(kraken_file):
    """Parse kraken output file"""

    # initialize variables to store data from kraken file.
    sample_id = ""  # sample id for each line of kraken output file.
    taxon = ""  # taxon name for each line of kraken output file.

    # initialize dictionary to store data from kraken file.
    dic_kraken = {}

    with open(kraken_file) as f:
        for line in f:
            if (
                line[0] == "C"
            ):  # if first character is 'C', then this is a comment line and can be skipped.
                continue

            else:  # otherwise, this is a data line and should be parsed.
                row = re.split(
                    "\t|\n", line
                )  # split the line into its components using tab or newline as delimiter.

                sample_id = row[
                    1
                ]  # second column contains sample id (first column contains read id).

                taxon = row[
                    4
                ]  # fifth column contains taxon name (fourth column contains read length).

                if (
                    sample_id not in dic_kraken
                ):  # if this is the first time seeing this sample id...
                    dic_kraken[
                        sample_id
                    ] = {}  # create a new dictionary entry for this sample id...

                    dic_kraken[sample_id][
                        taxon
                    ] = 1  # and initialize count for this taxon to 1 (since we've seen it once).

                else:  # otherwise, we've already seen this sample id before...
                    if (
                        taxon not in dic_kraken[sample_id]
                    ):  # if we haven't seen this particular taxon before...
                        dic_kraken[sample_id][
                            taxon
                        ] = 1  # initialize count for that particular taxon to 1 (since we've seen it once).
                    else:  # otherwise, we've already seen both the sample id and the particular taxon before...
                        dic_kraken[sample_id][
                            taxon
                        ] += 1  ##increment count by one since we've just seen that particular combination again.

    return dic_kraken


def parseTaxonomy(taxonomyFile):
    """Parse NCBI Taxonomy File"""
    with open(taxonomyFile) as f:
        lines = f.readlines()
        lines = [x.strip() for x in lines]
        return lines


def makeTaxaDict(lines):
    """Make Dictionary of Taxa Names"""
    d = {}
    for i in range(0, len(lines), 3):
        key = lines[i].replace(">", "")
        value = lines[i + 2].replace("\t", "").replace("|", ".")[:-1].split(".")[1:]
        value = [x for x in value if x != ""]
        value = ";".join(value) + ";"
        d[key] = value
        return d


def makeTableFromKrakOutAndTaxDict(dKrakOut, dTaxDict):
    """Make BIOM-format Table From Kraken Output File and Taxonomy Dictionary"""
    matrix = []
    samples = []
    rows = []
    cols = []
    rowsums = {}
    colsums = {}
    ##Get all unique samples from Kraken output dictionary##
    samples = [x for x in list(set([y for y in [z for z in list((dKrakOut))]]))]
    ##Get all unique rows from Kraken output dictionary##
    rows = [x for x in list((set([y for y in [z for z in list((dKrakOut))]])))]
    ##Get all unique columns from Kraken output dictionary##
    cols = [
        x
        + ";"
        + dTaxDict[x]
        + ";"
        + "unclassified"
        + ";"
        + "unclassified;"
        + "unclassified;"
        + "unclassified;"
        + "unclassified;"
        + "unclassified;unknown;unknown;unknown;unknown;unknown;unknown;unknown;unknown;unknown;unknown;;unknown;;unknown;;unknown;;unknown;;unknown;;unknown;;unknown;;unknown;;unknown;;Unknown;;;;Unknown;;;;Unknown;;;;Unknown;;;;Unknown;;;;Unknown;;;;Unknown;;;;Unknown;;;;Unknown;;;;Unknown;;;;Unknown;;;;Unassigned;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Unassigned;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Unassigned;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;"
    ]
    for i in range(len(rows)):
        matrixrow = []
        rowsum = 0
        rowlabel = rows[i]
        print("Processing Sample " + str((i + 1)) + " out of " + str((len(rows))))
        ##For each row (each unique sample), get all unique columns##
        cols = [
            x
            + ";"
            + (dTaxDict[x])
            + ";unclassified ; unclassified ; unclassified ; unclassified ; unclassified ; unknown ; unknown ; unknown ; unknown ; unknown ; unknown ; unknown ; unknown ; unknown ;; unknown ;; unknown ;; unknown ;; unknown ;; unknown ;; unkn"
        ]


# Function to convert biom files to tsv files.
import skbio as skb


# Input: biom file
# Output: tsv file with the following columns: OTU ID, Kingdom, Phylum, Class, Order, Family, Genus, Species
def biom_to_tsv_skbio(biom_file, tsv_file):
    table = skb.Table.read(biom_file)

    table.to_tsv(tsv_file)


# import pandas as pd
# import numpy as np
# from biom import load_table, Table


def biom_to_tsv_biom(biom_file):
    table = load_table(biom_file)

    # Convert to pandas dataframe.  This is the easiest way to get the data into a tsv file.  You can also use the write command from biom library but it is not as flexible.  See https://biom-format.org/documentation/format_versions/biom-2.0.html for more info on how to use write command in biom library.

    df = table.to_dataframe()

    # Transpose dataframe so that samples are columns and features are rows (this is how qiime2 likes it).  The transpose function in pandas does not work with sparse matrices, so we need to convert to dense matrix first using .toarray() function and then transpose using .T function (see https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.T.html for more info).

    df = df.toarray().T

    # Get sample IDs from table object and convert to list of strings (the sample IDs are stored as a numpy array of bytes, so we need to convert them back into strings).  We will use these sample IDs later when we save the dataframe as a tsv file (see below).

    samples = list(map(str, table.ids("sample")))

    # Get feature IDs from table object and convert them into a list of strings (the feature IDs are stored as a numpy array of bytes, so we need to convert them back into strings).  We will use these feature IDs later when we save the dataframe as a tsv file (see below).

    features = list(map(str, table.ids("observation")))

    # Create new column names for our dataframe by appending 'X' at beginning of each sample ID string in samples list above and adding '.1' at end of each string in features list above (this is how qiime2 likes it).

    # newcol
