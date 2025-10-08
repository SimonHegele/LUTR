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

def gene_pairs(predicted_gff: DataFrame,
               assembled_gff:   DataFrame) -> Generator[tuple[Series, DataFrame], None, None]:
    
    predicted_genes = predicted_gff.loc[predicted_gff["type"]=="gene"]
    assembled_genes = assembled_gff.loc[assembled_gff["type"]=="gene"]

    for _, predicted_gene in predicted_genes.iterrows():

        matches = assembled_genes.loc[assembled_genes["strand"]==predicted_gene["strand"]]
        matches = overlapping_features(matches, predicted_gene)

        yield predicted_gene, matches
        
        assembled_genes.drop(matches.index)

def gene_pairs_wrapper(args: tuple[DataFrame, DataFrame]) -> list[tuple[Series, DataFrame]]:

    predicted_gff, assembled_gff = args

    return list(gene_pairs(predicted_gff, assembled_gff))

def gene_pairs_multiprocessing(predicted_gff: dict[str, DataFrame],
                               assembled_gff: dict[str, DataFrame],
                               seqnames:      list[str],
                               threads:       int) -> list[tuple[Series,DataFrame]]:
    
    args     = zip([predicted_gff[s] for s in seqnames],
                   [assembled_gff[s] for s in seqnames])
    
    with Pool(min(len(seqnames), threads)) as pool:

        total_gene_pairs: list[tuple[Series, DataFrame]] = []
        
        for gene_pair in list(chain(pool.map(gene_pairs_wrapper, args))):

            total_gene_pairs += gene_pair
            
        key = lambda gp: feature_length(gp[0])

        return sorted(total_gene_pairs, key=key, reverse=True)