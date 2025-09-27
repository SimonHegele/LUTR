"""
Module Name:    main.py
Author:         Simon Hegele
Date:           2025-09-19
Version:        1.0
License:        GPL-3
"""

from logging                     import info
from multiprocessing             import Pool
from pandas                      import concat, DataFrame
from os                          import path, listdir, remove

from .argument_parser            import LUTRArgparser
from .lutr_logging               import logging_setup, log_wrapper
from .gffutils                   import load_gff, seqname_split_wrapper, write_gff
from .gene_pair_slice_gen        import gene_pair_slices_multiprocessed
from .process_gene_pairs         import processgenepairs_multiprocessing

def delete_all_files(folder):
    for filename in listdir(folder):
        filepath = path.join(folder, filename)
        if path.isfile(filepath):
            remove(filepath)
            
def load_data(path_prediction: str, path_assembly: str) -> tuple[DataFrame, DataFrame]:
    
    info("Loading GFF-file (Prediction) ...")
    prediction = log_wrapper(load_gff, path_prediction)
    info("Loading GFF-file (Assembly) ...")
    assembly   = log_wrapper(load_gff, path_assembly)
    
    return prediction, assembly

def main():

    args    = LUTRArgparser().parse_args()
    tmpdir  = path.join(args.outdir,"tmpdir")
    
    logging_setup(args)
    
    ######################################################################################
    # Step 1: Generating GFF-slices for pairs of overlapping genes                       #
    ######################################################################################
    
    if not path.isfile(path.join(args.outdir, "step_1")):
        
        info("Step 1:")
        
        prediction, assembly = load_data(args.prediction, args.assembly)
            
        info("Generating chromosome slices ...")
        seqnames   = sorted(set(list(prediction["seqname"].unique()) +
                                list(assembly["seqname"].unique())))
        with Pool(min(2,args.threads)) as pool:
            prediction, assembly = pool.map(seqname_split_wrapper,
                                            [(prediction, seqnames),(assembly,seqnames)])
            
        info("Generating gene slices ...")
        gene_pair_slices_multiprocessed(prediction,
                                        assembly,
                                        args.threads,
                                        tmpdir)

        open(path.join(args.outdir, "step_1"),"x").close()
        
        del prediction
        del assembly
    
    else:
        
        info("Skipping step 1")
    
    ######################################################################################
    # Step 2: Gene-pair wise transcript matching and UTR-extension                       #
    ######################################################################################
        
    info("Step 2:")
    if any(listdir(path.join(tmpdir, "genes_lutr"))):
        info("Delete files from previous run")
        delete_all_files(path.join(tmpdir, "genes_lutr"))
        
    info("Processing gene pairs ...")
    processgenepairs_multiprocessing(tmpdir,
                                     args.mftm,
                                     args.mtbm,
                                     args.threads,
                                     args.verbosity,
                                     args.remove_unsupported,
                                     args.select,
                                     args.match_all)
    
    info("Concatenating and sorting ...")
    
    tmpdir = path.join(tmpdir, "genes_lutr")
    
    with Pool(args.threads) as pool:
        lutr_gff = pool.map(load_gff,
                            [path.join(tmpdir,file) for file in listdir(tmpdir)])
    
    key      = lambda x: (x.iloc[0]["seqname"], x.iloc[0]["start"])
    lutr_gff = sorted(lutr_gff, key=key)
    lutr_gff = concat(lutr_gff)
    
    write_gff(lutr_gff, path.join(args.outdir, "lutr.gff"))
    
    info("############################################")
    info("#    Simon says: Thanks for using LUTR!    #")

    info("############################################")




