"""
Module Name:    gene_pair_gen
Author:         Simon Hegele
Date:           2025-09-19
Version:        1.0
License:        GPL-3
Description:    Pairing of predicted and assembled gene features by overlapping
"""

from itertools       import chain
from typing          import Generator
from multiprocessing import Pool
from pandas          import DataFrame, Series

from .gffutils import overlapping_features, feature_length

def gene_pairs(gff_prediction: DataFrame,
               gff_assembly:   DataFrame) -> Generator[tuple[Series, DataFrame], None, None]:
    
    predicted_genes = gff_prediction.loc[gff_prediction["type"]=="gene"]
    assembled_genes = gff_assembly.loc[gff_assembly["type"]=="gene"]

    for _, predicted_gene in predicted_genes.iterrows():

        matches = assembled_genes.loc[assembled_genes["strand"]==predicted_gene["strand"]]
        matches = overlapping_features(matches, predicted_gene)

        yield predicted_gene, matches

def gene_pairs_wrapper(args: tuple[DataFrame, DataFrame]) -> list[tuple[Series, DataFrame]]:

    gff_prediction, gff_assembly = args

    return list(gene_pairs(gff_prediction, gff_assembly))

def gene_pairs_multiprocessing(gff_prediction: dict[str, DataFrame],
                               gff_assembly:   dict[str, DataFrame],
                               seqnames:       list[str],
                               threads:        int) -> list[tuple[Series,DataFrame]]:
    
    args     = zip([gff_prediction[s] for s in seqnames],
                   [gff_assembly[s] for s in seqnames])
    
    with Pool(min(len(seqnames), threads)) as pool:

        total_gene_pairs: list[tuple[Series, DataFrame]] = []
        
        for gene_pair in list(chain(pool.map(gene_pairs_wrapper, args))):

            total_gene_pairs += gene_pair
            
        key = lambda gp: feature_length(gp[0])

        return sorted(total_gene_pairs, key=key, reverse=True)