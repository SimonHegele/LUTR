"""
Module Name:    setup_logging.py
Author:         Simon Hegele
Date:           2025-09-19
Version:        1.0
License:        GPL-3

Description:    Facilitates the matching of and creation of UTR-variants from predicted
                and assembled transcripts
"""

from logging                     import info, error, debug
from multiprocessing             import Queue, Process
from numpy                       import array, sum, max, min, ndarray
from os                          import listdir, path
from pandas                      import concat, DataFrame
from time                        import time

from .gffutils                   import *
from .transcript_matching        import tmatch
from .utr_variant_generation     import generate_utr_variant
from .lutr_logging               import log_wrapper

def select_utr_variant(utr_variants: list[DataFrame],
                       select:       str) -> list[DataFrame]:
    """Length based selection of UTR-variants

    Args:
        utr_variants (list[DataFrame]): UTR-variant GFF-slices
        select (str):                   "shortest", "longest" or "all"

    Raises:
        Exception: If selction could not be performed

    Returns:
        list[DataFrame]: Selection from the input list
    """
    
    if select == "all" or len(utr_variants) <= 1:
        return utr_variants
    
    exons   = [v.loc[v["type"]=="exon"] for v in utr_variants]
    lengths = [sum(e["end"]-e["start"]-1) for e in exons]
        
    match select:
        case "shortest":
            return [utr_variants[lengths.index(min(lengths))]]
        case "longest":
            return [utr_variants[lengths.index(max(lengths))]]
        
    raise Exception(f"Could not select {select} from {utr_variants}")

def sort_by_number_of_exons(slices: list[DataFrame]) -> list[DataFrame]:

    return sorted(slices,
                  key=lambda slice: len(slice.loc[slice["type"]=="exon"]),
                  reverse=True)

def get_match_matrix(p_trans:      list[Series],
                     a_trans:      list[Series],
                     p_tran_exons: list[DataFrame],
                     a_tran_exons: list[DataFrame],
                     mftm:         float,
                     mtbm:         int,
                     match_all:    bool,
                     match_middle) -> ndarray:
    
    debug("Creating transcript match matrix:")
    
    m = array([[tmatch(p_tran, a_tran, p_tran_exons[i], a_tran_exons[j], mftm, mtbm, match_middle)
                for j, a_tran in enumerate(a_trans)] for i, p_tran in enumerate(p_trans)])
    
    if not match_all:
        # Use assembled transcripts for UTR-extensions of the predicted transcript they
        # overlap with the most
        m = array([[m[i][j] if m[i][j] == max(m[:,j]) else 0
            for j, _ in enumerate(a_trans)] for i, _ in enumerate(p_trans)])
    
    debug("%s", m)
    
    return m
    
def potential_utr_extension(p_tran:       Series,
                            a_tran:       Series,
                            matches:      int,
                            p_tran_slice: DataFrame) -> int:
    
    if not matches:
        return 0
    if not ((a_tran["start"] < p_tran["start"]) or (p_tran["end"] < a_tran["end"])):
        return 0
    if not len(p_tran_slice.loc[p_tran_slice["type"]=="CDS"]):
        return 0
    return 1
    
def get_utr_matrix(p_trans:        list[Series], 
                   a_trans:        list[Series],
                   p_trans_slices: list[DataFrame],
                   match_matrix:   ndarray) -> list[list[int]]:
    """
    Args:
        p_trans (list[Series]):           Predicted transcript features
        a_trans (list[Series]):           Assembled transcript features
        p_trans_slices (list[DataFrame]): Predicted transcript and associated features
        match_matrix (ndarray):           If predicted transcript i and assembled
                                          transcript j are matching:
                                              match_matrix[i][j] = number of shared bases
                                          else:
                                              match_matrix[i][j] = 0
 
    Returns:
        list[list[int]]: 
    """
    
    m = [[potential_utr_extension(p_tran, a_tran, match_matrix[i][j], p_trans_slices[i])
          for j, a_tran in enumerate(a_trans)] for i, p_tran in enumerate(p_trans)]
    
    m = [[sum(m[i][:j+1]) if m[i][j] == 1 else 0
          for j, _ in enumerate(a_trans)] for i, _ in enumerate(p_trans)]
    
    return m
   
def processgenepair(tmpdir:     str,
                    p_path:     str,
                    a_path:     str,
                    mftm:       float,
                    mtbm:       int,
                    remove:     bool,
                    select:       str,
                    match_all:    bool,
                    match_middle: bool) -> str:
    
    def transcript_slices(gene_slice: DataFrame,
                          transcripts: list[Series]) -> Generator:
        
        id2index     = get_map_id2index(gene_slice)
        parent2child = get_map_parent2children(gene_slice, id2index)
    
        for transcript in transcripts:
            
            indices          = get_subtree(gene_slice, transcript.name, parent2child)
            transcript_slice = gene_slice.iloc[indices].reset_index(drop=True)
            
            yield transcript_slice
       
    # Loading Data
    p_gene_slice = load_gff(p_path)
    a_gene_slice = load_gff(a_path)
    gene         = p_gene_slice.loc[p_gene_slice["type"]=="gene"]
    gene_id      = gene.iloc[0]["ID"]
    
    # Extracting transcript data
    p_trans        = [tran for _, tran in get_trans(p_gene_slice).iterrows()]
    a_trans        = [tran for _, tran in get_trans(a_gene_slice).iterrows()]
    p_trans_slices = list(transcript_slices(p_gene_slice, p_trans))
    a_trans_slices = list(transcript_slices(a_gene_slice, a_trans))
    p_trans_exons  = [slice.loc[slice["type"]=="exon"] for slice in p_trans_slices]
    a_trans_exons  = [slice.loc[slice["type"]=="exon"] for slice in a_trans_slices]
    
    # Creating the match matrix
    match_matrix = get_match_matrix(p_trans, a_trans, p_trans_exons, a_trans_exons, mftm, mtbm, match_all, match_middle)
    utr_matrix   = get_utr_matrix(p_trans, a_trans, p_trans_slices, match_matrix)
    
    # Creating UTR-variants and selecting transcrips for output
    lutr_trans_slices = []
    for i, p_tran_slice in enumerate(p_trans_slices):
        
        # Here, we decide what to do with unmatched transcripts (remove or keep)
        if not any(match_matrix[i]):
            if remove:
                continue
                
        # Here, we decide what to do with matched transcripts
        # (keep or create and select from UTR-variants):
        if not any(utr_matrix[i]):
            lutr_trans_slices.append(p_tran_slice)
        else:
            utr_variants = [generate_utr_variant(p_tran_slice,
                                                 a_tran_slice,
                                                 p_trans_exons[i],
                                                 a_trans_exons[j],
                                                 p_trans[i],
                                                 a_trans[j],
                                                 utr_matrix[i][j],
                                                 gene_id)
                            for j, a_tran_slice in enumerate(a_trans_slices)
                            if utr_matrix[i][j] >= 1]
            
            lutr_trans_slices += select_utr_variant(utr_variants, select)
            
    # Outputting results to tmpdir
    outfile = path.join(path.join(tmpdir, "genes_lutr"), f"lutr_{gene_id}.gff")
    
    if len(lutr_trans_slices) == 0:
        lutr_gene_slice = empty_gff()
    else:
        lutr_gene_slice = concat([slice[gff_columns] for slice in lutr_trans_slices], axis=0)
        gene.iloc[0, 3] = lutr_gene_slice["start"].min()
        gene.iloc[0, 4] = lutr_gene_slice["end"].max()
        lutr_gene_slice = concat([gene, lutr_gene_slice], )
    
    write_gff(lutr_gene_slice, outfile)
        
    return gene_id

def processgenepair_wrapper(args: tuple[str, str, str, float, int, int, int, bool, str, bool, bool]):

    start = time()

    tmpdir, p_path, a_path, mftm, mtbm, verbosity, n_genes, remove, select, match_all, match_middle = args

    gene_id  = log_wrapper(processgenepair,
                           (tmpdir,p_path,a_path,mftm,mtbm,remove,select,match_all,match_middle))

    progress = len(listdir(path.join(tmpdir, "genes_lutr")))
    duration = time()-start
    percent  = round(100*progress/n_genes,1)
    
    if progress % 10**verbosity == 0:
        info(f"{progress:>6} / {n_genes} (~{percent:6>}%), {duration:>6.2f} seconds {gene_id}")

def processgenepairs_multiprocessing(tmpdir: str,
                                     mftm: float,
                                     mtbm: int,
                                     threads: int,
                                     v: int,
                                     remove: bool,
                                     select: str,
                                     match_all: bool):

    p_dir = path.join(tmpdir, "genes_predicted")
    a_dir = path.join(tmpdir, "genes_assembled")

    p_paths = sorted([path.join(p_dir, f) for f in listdir(p_dir)])
    a_paths = sorted([path.join(a_dir, f) for f in listdir(a_dir)])

    if not len(p_paths) == len(a_paths):
        error(f"Incomplete gene pairs in tmpdir! Predicted: {len(p_paths)}, Assembled: {len(a_paths)}")
        exit(1)

    # Sort by estimated size (to process longer jobs first)
    gene_pairs = sorted(zip(p_paths, a_paths),
                        key=lambda gp: max([path.getsize(gp[0]), path.getsize(gp[1])]),
                        reverse=True)

    tasks = [(tmpdir, p, a, mftm, mtbm, v, len(gene_pairs), remove, select, match_all, False)
             for p, a in gene_pairs]

    # Create a multiprocessing Queue and enqueue tasks in order
    task_queue = Queue()
    for task in tasks:
        task_queue.put(task)

    # Add stop signals (one per worker)
    for _ in range(threads):
        task_queue.put(None)

    # Define worker function
    def worker():
        while True:
            task = task_queue.get()
            if task is None:
                break
            processgenepair_wrapper(task)

    # Launch workers
    workers = [Process(target=worker) for _ in range(threads)]
    for w in workers:
        w.start()

    # Wait for workers to complete
    for w in workers:

        w.join()
