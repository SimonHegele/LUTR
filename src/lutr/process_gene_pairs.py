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
from numpy                       import array, sum, max, min, ndarray, inf
from os                          import listdir, path
from pandas                      import concat, DataFrame
from time                        import time

from .gffutils                   import *
from .transcript_matching        import tmatch
from .utr_variant_generation     import generate_utr_variant
from .lutr_logging               import log_wrapper

def processgenepair(tmpdir:       str,
                    p_path:       str,
                    a_path:       str,
                    mftm:         float,
                    mtbm:         int,
                    remove:       bool,
                    select:       str,
                    match_all:    bool,
                    match_middle: bool) -> str:
    
    def extract_transcript_data(gene_slice: DataFrame) -> tuple[list[Series],
                                                                list[DataFrame],
                                                                list[DataFrame]]:
    
        def transcript_slices(transcripts: list[Series]) -> Generator:
        
            parent2child = get_map_parent2children(gene_slice)

            for transcript in transcripts:
                
                indices          = get_subtree(gene_slice, transcript.name, parent2child)
                transcript_slice = gene_slice.iloc[indices].reset_index(drop=True)
                
                yield transcript_slice
                
        trans  = [t for _, t in get_trans(gene_slice).iterrows()]
        slices = list(transcript_slices(trans))
        exons  = [s.loc[s["type"]=="exon"] for s in slices]
        
        return trans, slices, exons
                
    def get_match_matrix() -> ndarray:
    
        m = array([[tmatch(p_tran,
                           a_tran,
                           p_tran_exons[i],
                           a_tran_exons[j],
                           mftm, mtbm,
                           match_middle)
                    for j, a_tran in enumerate(a_trans)]
                   for i, p_tran in enumerate(p_trans)])
        
        if not match_all:

            m = array([[m[i][j] if m[i][j] == max(m[:,j]) else 0
                for j, _ in enumerate(a_trans)] for i, _ in enumerate(p_trans)])
        
        debug("%s", m)
        
        return m
    
    def get_utr_matrix() -> list[list[int]]:
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
        
        def potential_utr_extension(p_tran:        Series,
                                    a_tran:        Series,
                                    matched_bases: int,
                                    p_tran_slice:  DataFrame) -> int:
    
            if not matched_bases:
                return 0
            if not ((a_tran["start"] < p_tran["start"]) or (p_tran["end"] < a_tran["end"])):
                return 0
            if not len(p_tran_slice.loc[p_tran_slice["type"]=="CDS"]):
                return 0
            return 1
        
        m = [[potential_utr_extension(p_tran, a_tran, match_matrix[i][j], p_tran_slices[i])
            for j, a_tran in enumerate(a_trans)] for i, p_tran in enumerate(p_trans)]
        
        m = [[sum(m[i][:j+1]) if m[i][j] == 1 else 0
            for j, _ in enumerate(a_trans)] for i, _ in enumerate(p_trans)]
        
        return m

    def select_utr_variant(utr_variants: list[DataFrame]) -> list[DataFrame]:
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
    
    # 1. Preparing data
    
    p_gene_slice = load_gff(p_path)
    a_gene_slice = load_gff(a_path)
    
    gene         = p_gene_slice.loc[p_gene_slice["type"]=="gene"]
    gene_id      = gene.iloc[0]["ID"]
    outfile      = path.join(path.join(tmpdir, "genes_lutr"), f"lutr_{gene_id}.gff")
    
    p_trans, p_tran_slices, p_tran_exons = extract_transcript_data(p_gene_slice)
    a_trans, a_tran_slices, a_tran_exons = extract_transcript_data(a_gene_slice)
    
    # 2. Matching
    
    match_matrix = get_match_matrix()
    utr_matrix   = get_utr_matrix()
    
    # 3. Processing of predicted transcripts
    # (UTR-variant creation and transcript selection)
    
    l_tran_slices = []
    
    for i, p_tran_slice in enumerate(p_tran_slices):
        
        # Here, we decide what to do with unmatched transcripts
        # (remove or keep)
        if remove and not any(match_matrix[i]):
            continue
                
        # Here, we decide what to do with matched transcripts
        # (keep or create and select from UTR-variants):
        if not any(utr_matrix[i]):
            l_tran_slices.append(p_tran_slice)
        else:
            utr_variants = [generate_utr_variant(p_tran_slice,
                                                 a_tran_slice,
                                                 p_tran_exons[i],
                                                 a_tran_exons[j],
                                                 p_trans[i],
                                                 a_trans[j],
                                                 utr_matrix[i][j],
                                                 gene_id)
                            for j, a_tran_slice in enumerate(a_tran_slices)
                            if utr_matrix[i][j] >= 1]
            
            l_tran_slices += select_utr_variant(utr_variants)
            
    # 4. Outputting
    if len(l_tran_slices) > 0:

        lutr_gene_slice = concat([slice[gff_columns] for slice in l_tran_slices], axis=0)
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

    # Sort by estimated size (to process longer jobs first)
    gene_pairs = sorted(zip(p_paths, a_paths),
                        key=lambda gp: path.getsize(gp[0])*path.getsize(gp[1]),
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
