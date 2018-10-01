#! /usr/bin/env python
""" 
Script that takes a matrix of counts
where the columns are genes and the rows
are spot coordinates
        gene    gene    
XxY
XxY

And removes the columns of genes
whose type is not in the allowed types given (Ensembl annotation
gene type). For this, the script needs to be given an annotation
file in GFF format. 

@Author Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>
"""

import argparse
import sys
import os
import pandas as pd
from stpipeline.common.gff_reader import *
              
def main(counts_matrix, gene_types_keep, outfile, annotation, ensembl_ids, gene_type_attribute):

    if not os.path.isfile(counts_matrix) or not os.path.isfile(annotation):
        sys.stderr.write("Error, input file not present or invalid format\n")
        sys.exit(1)
     
    if not outfile:
        outfile = "filtered_{}".format(os.path.basename(counts_matrix))
    
    gene_types = dict()
    for line in gff_lines(annotation):
        gene_name = line["gene_id"] if ensembl_ids else line["gene_name"]
        try:
            gene_types[gene_name] = line[gene_type_attribute]
        except KeyError:
            sys.stderr.write("Error, the gene_type_attribute={} is not present for gene {}\n".format(gene_type_attribute,gene_name))
            sys.exit(1)
    assert(len(gene_types) > 0)

    # Read the data frame (genes as columns)
    counts_table = pd.read_table(counts_matrix, sep="\t", header=0, index_col=0)
    genes = counts_table.columns
    # Filter out genes that match any of the allowed types
    genes_drop = list()
    for gene in genes:
        try:
            if gene_types[gene] not in gene_types_keep:
                genes_drop.append(gene)
        except KeyError:
            sys.stdout.write("Warning, {} was not found in the annotation\n".format(gene))
    if len(genes_drop) > 0:
        counts_table.drop(genes_drop, axis=1, inplace=True)
    else:
        print "Not a single gene could be discarded..."
    # Write filtered table
    counts_table.to_csv(outfile, sep='\t')
               
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--counts-matrix", required=True,
                        help="Matrix with gene counts (genes as columns)")
    parser.add_argument("--outfile", help="Name of the output file")
    parser.add_argument("--gene-types-keep", required=True, nargs='+', type=str,
                        help="List of Ensembl gene types to keep (E.x protein_coding lincRNA")
    parser.add_argument("--gene-type-attribute", required=False, default="gene_type", type=str,
                        help='Name of the gene type attribute as specified in the Ensembl annotation file (default is "gene_type", "gene_biotype" is also a common value).')
    parser.add_argument("--annotation", help="The Ensembl annotation file", required=True, type=str)
    parser.add_argument("--ensembl-ids", action="store_true", 
                        default=False, help="Pass this parameter if the genes in the matrix" \
                        "are named with Ensembl Ids instead of gene names")
    args = parser.parse_args()
    main(args.counts_matrix, args.gene_types_keep, args.outfile, args.annotation, args.ensembl_ids, args.gene_type_attribute)
