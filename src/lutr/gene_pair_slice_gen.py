"""
Module Name:    gps_gen_lt
Author:         Simon Hegele
Date:           2025-09-19
Version:        1.0
License:        GPL-3
Description:    Generation of GFF-slices for pairs of predicted and assembled genes
                when a low number of threads is given
"""

from multiprocessing import Pool
from os              import path
from pandas          import DataFrame

from .gffutils       import *
from .gene_pair_gen  import gene_pairs

def gene_pair_slices(predicted_gff: DataFrame,
                     assembled_gff: DataFrame,
                     tmpdir:        str):
    
    predicted_outdir   = path.join(tmpdir, "genes_predicted")
    assembled_outdir   = path.join(tmpdir, "genes_assembled")
    predicted_children = get_map_parent2children(predicted_gff)
    assembled_children = get_map_parent2children(assembled_gff)

    for predicted_gene, assembled_genes in gene_pairs(predicted_gff, assembled_gff):
        
        predicted_outfile = path.join(predicted_outdir, "gp_"+ predicted_gene["ID"] +".gff")
        assembled_outfile = path.join(assembled_outdir, "ga_"+ predicted_gene["ID"] +".gff")
        
        write_gff(predicted_gff.iloc[get_subtree(predicted_gff,
                                                  predicted_gene.name,
                                                  predicted_children)],
                  predicted_outfile)
        
        open(assembled_outfile, "w").close()
        
        for _, assembled_gene in assembled_genes.iterrows():
            
            write_gff(assembled_gff.iloc[get_subtree(assembled_gff,
                                                     assembled_gene.name,
                                                     assembled_children)],
                  assembled_outfile,
                  mode="a")
            
def generate_gene_pair_slices_wrapper(args: tuple[DataFrame,DataFrame,str]):
    
    predicted_gff, assembled_gff, tmpdir = args
    
    gene_pair_slices(predicted_gff,assembled_gff,tmpdir)
        
def gene_pair_slices_multiprocessed(prediction: dict[str, DataFrame],
                           assembly:   dict[str, DataFrame],
                           threads:    int,
                           tmpdir:     str):
    
    tasks    = [(prediction[seqname], assembly[seqname], tmpdir)
                for seqname in prediction.keys()]
    
    with Pool(min(threads, len(prediction.keys()))) as pool:
        
        pool.map(generate_gene_pair_slices_wrapper, tasks)