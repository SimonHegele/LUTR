"""
Module Name:    gps_gen_mt
Author:         Simon Hegele
Date:           2025-09-19
Version:        1.0
License:        GPL-3
Description:    Generation of GFF-slices for pairs of predicted and assembled genes
                when a high number of threads is given
"""
from logging         import info
from multiprocessing import Pool 
from pandas          import concat, DataFrame, Series
from os              import listdir, path
from time            import time

from .gffutils       import empty_gff, get_subtree, write_gff
from .gene_pair_gen  import gene_pairs_multiprocessing
from .lutr_logging   import log_wrapper

shared_prediction = dict({})
shared_assembly   = dict({})

def init_globals(prediction, assembly):

    global shared_prediction, shared_assembly
    shared_prediction = prediction
    shared_assembly   = assembly

def gene_pair_slices(gff_prediction:  DataFrame,
                     gff_assembly:    DataFrame,
                     predicted_gene:  Series,
                     assembled_genes: DataFrame,
                     tmpdir:          str):
    
    predicted_gene_slice = get_subtree(gff_prediction,
                                      predicted_gene)
    gene_id = predicted_gene["ID"]
    
    if len(assembled_genes) == 0:
        assembled_gene_slice = empty_gff()
    else:
        assembled_gene_slice = concat([get_subtree(gff_assembly, gene)
                                      for _, gene in assembled_genes.iterrows()])
    
    write_gff(predicted_gene_slice, path.join(path.join(tmpdir, "genes_predicted"), f"gp_{gene_id}.gff"))
    write_gff(assembled_gene_slice, path.join(path.join(tmpdir, "genes_assembled"), f"ga_{gene_id}.gff"))
    
    return gene_id

def gene_pair_slices_wrapper(args: tuple):
    
    predicted_gene, assembled_genes, verbosity, n_genes, tmpdir = args
    
    seqname    = predicted_gene["seqname"]
    prediction = shared_prediction[seqname]
    assembly   = shared_assembly[seqname]
    start      = time()
    progress   = len(listdir(path.join(tmpdir, "genes_predicted")))
    duration   = time()-start
    percent    = round(100*progress/n_genes,1)
    gene_id    = log_wrapper(gene_pair_slices,
                            (prediction,assembly,predicted_gene,assembled_genes,tmpdir))
    if progress % 10**verbosity == 0:
        info(f"{progress:>6} / {n_genes} (~{percent:6>}%), {duration:>7.2f} seconds {gene_id}")

def gene_pair_slices_mp_mt(prediction: dict[str, DataFrame],
                           assembly:   dict[str, DataFrame],
                           threads:    int,
                           verbosity:  int,
                           tmpdir:     str,
                           seqnames:   list[str]):

    info("Identifying gene pairs")
    gene_pairs = gene_pairs_multiprocessing(prediction,
                                            assembly,
                                            seqnames,
                                            threads)
    
    np_genes   = len(gene_pairs)
    na_genes   = sum([1 for gene_pair in gene_pairs if not gene_pair[1].empty])
    ratio      = na_genes / np_genes
    info(f"{np_genes} genes predicted, {na_genes} ({ratio:.2%}) overlapped by assembled genes")
    
    info("Generating gene slices")
    
    tasks = [(gene_pair[0], gene_pair[1], verbosity, len(gene_pairs), tmpdir)
             for gene_pair in gene_pairs]

    with Pool(threads, initializer=init_globals, initargs=(prediction, assembly)) as pool:
        
        results = [pool.apply_async(gene_pair_slices_wrapper, (t,)) for t in tasks]
        results = [r.get() for r in results]