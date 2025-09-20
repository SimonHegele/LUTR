"""
Module Name:    gps_gen_lt
Author:         Simon Hegele
Date:           2025-09-19
Version:        1.0
License:        GPL-3
Description:    Generation of GFF-slices for pairs of predicted and assembled genes
                when a low number of threads is given
"""

from logging         import info
from multiprocessing import Pool
from os              import listdir, path
from pandas          import concat, DataFrame
from time            import time

from .gffutils import empty_gff, get_subtree, write_gff
from .gene_pair_gen import gene_pairs

def generate_gene_pair_slices(gff_prediction: DataFrame,
                              gff_assembly:   DataFrame,
                              tmpdir:         str,
                              verbosity:      int):
    
    out_prediction = path.join(tmpdir, "genes_predicted")
    out_assembly   = path.join(tmpdir, "genes_assembled")
    
    for predicted_gene, matches in gene_pairs(gff_prediction, gff_assembly):
        
        start = time()
        
        gene_id               = predicted_gene["ID"]
        gene_slice_prediction = get_subtree(gff_prediction, predicted_gene)
        
        if len(matches) > 0:
            gene_slice_assembly = concat([get_subtree(gff_assembly, assembled_gene)
                                            for _, assembled_gene in matches.iterrows()])
        else:
            gene_slice_assembly = empty_gff()
        
        gff_prediction.drop(gene_slice_prediction.index, inplace=True)
        gff_assembly.drop(gene_slice_assembly.index, inplace=True)
        
        write_gff(gene_slice_prediction, path.join(out_prediction, f"gp_{gene_id}.gff"))
        write_gff(gene_slice_assembly,   path.join(out_assembly, f"ga_{gene_id}.gff"))
        
        progress = len(listdir(out_prediction))
        
        if progress % 10**verbosity == 0:
        
            info(f"{progress:>6} {time()-start:>7.2f} seconds {gene_id}")
            
def generate_gene_pair_slices_wrapper(args: tuple[DataFrame,DataFrame,str,int]):
    
    gff_prediction, gff_assembly, tmpdir, verbosity = args
    
    generate_gene_pair_slices(gff_prediction,gff_assembly,tmpdir,verbosity)
        
def gene_pair_slices_mp_lt(prediction: dict[str, DataFrame],
                           assembly:   dict[str, DataFrame],
                           threads:    int,
                           verbosity:  int,
                           tmpdir:     str):
    
    info("Generating gene slices")
    
    tasks    = [(prediction[seqname], assembly[seqname], tmpdir, verbosity)
                for seqname in prediction.keys()]
    
    with Pool(min(threads, len(prediction.keys()))) as pool:
        
        pool.map(generate_gene_pair_slices_wrapper, tasks)