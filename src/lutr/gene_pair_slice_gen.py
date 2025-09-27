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
from pandas          import concat, DataFrame

from .gffutils       import *
from .gene_pair_gen  import gene_pairs

def gene_pair_slices(gff_prediction: DataFrame,
                     gff_assembly:   DataFrame,
                     tmpdir:         str):
    
    out_prediction = path.join(tmpdir, "genes_predicted")
    out_assembly   = path.join(tmpdir, "genes_assembled")
    
    map_i2i_prediction  = get_map_id2index(gff_prediction)
    map_i2i_assembly    = get_map_id2index(gff_assembly)
    map_p2ch_prediction = get_map_parent2children(gff_prediction, map_i2i_prediction)
    map_p2ch_assembly   = get_map_parent2children(gff_assembly,   map_i2i_assembly)
    
    for predicted_gene, matches in gene_pairs(gff_prediction, gff_assembly):
        
        gene_id               = predicted_gene["ID"]
        gene_slice_prediction = gff_prediction.iloc[get_subtree(gff_prediction,
                                                                predicted_gene.name,
                                                                map_p2ch_prediction)]
        
        if len(matches) > 0:
            gene_slice_assembly = concat([gff_assembly.iloc[get_subtree(gff_assembly,
                                                                        i,
                                                                        map_p2ch_assembly)]
                                            for i, _ in matches.iterrows()])
        else:
            gene_slice_assembly = empty_gff()
        
        write_gff(gene_slice_prediction, path.join(out_prediction, f"gp_{gene_id}.gff"))
        write_gff(gene_slice_assembly,   path.join(out_assembly, f"ga_{gene_id}.gff"))
            
def generate_gene_pair_slices_wrapper(args: tuple[DataFrame,DataFrame,str,int]):
    
    gff_prediction, gff_assembly, tmpdir = args
    
    gene_pair_slices(gff_prediction,gff_assembly,tmpdir)
        
def gene_pair_slices_multiprocessed(prediction: dict[str, DataFrame],
                           assembly:   dict[str, DataFrame],
                           threads:    int,
                           tmpdir:     str):
    
    tasks    = [(prediction[seqname], assembly[seqname], tmpdir)
                for seqname in prediction.keys()]
    
    with Pool(min(threads, len(prediction.keys()))) as pool:
        
        pool.map(generate_gene_pair_slices_wrapper, tasks)