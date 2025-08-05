from logging                     import info
from multiprocessing             import Pool
from pandas                      import concat
from os                          import path, listdir, remove

from .argument_parser            import LUTRArgparser
from .setup_logging              import logging_setup
from .gffutils                   import load_gff, seqname_split_wrapper, write_gff
from .gene_pair_slice_gen_lt     import gene_pair_slices_mp_lt
from .gene_pair_slice_gen_mt     import gene_pair_slices_mp_mt
from .process_gene_pairs         import processgenepairs_multiprocessing

def delete_all_files(folder):
    for filename in listdir(folder):
        filepath = path.join(folder, filename)
        if path.isfile(filepath):
            remove(filepath)

def main():

    args    = LUTRArgparser().parse_args()
    tmpdir  = path.join(args.outdir,"tmpdir")
    
    logging_setup(args)
    
    ######################################################################################
    # Step 1: Generating GFF-slices for pairs of overlapping genes                       #
    ######################################################################################
    
    if not path.isfile(path.join(args.outdir, "step_1")):
        
        info("Step 1:")
        
        info("Loading GFF-files ...")
        with Pool(2) as pool:
            prediction, assembly = pool.map(load_gff, [args.prediction, args.assembly])
            
        info("Generating chromosome slices ...")
        seqnames   = sorted(set(list(prediction["seqname"].unique()) +
                                list(assembly["seqname"].unique())))
        with Pool(min(2,args.threads)) as pool:
            prediction, assembly = pool.map(seqname_split_wrapper,
                                            [(prediction, seqnames),(assembly,seqnames)])
            
        # Depending on the number of threads we will be using a different strategy for
        # parallelizing the gene pair slice generation
        
        if len(seqnames)/2 > args.threads:
            # If given a low number of threads we will be processing each chromosome
            # seperately without shared memory across all processes
            gene_pair_slices_mp_lt(prediction, assembly, args.threads, args.verbosity, tmpdir)
            
        else:
            # If given a high number of threads we will be processing each gene separately
            # with shared memory. Takes more time per gene but is overall faster as we 
            # can use all processes fully
            gene_pair_slices_mp_mt(prediction, assembly, args.threads, args.verbosity, tmpdir, seqnames)

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
    
    name = f"lutr_{args.select}_mtbm{args.mtbm}_mftm{args.mftm:.2f}.gff"
    write_gff(lutr_gff, path.join(args.outdir, name))
    
    info("############################################")
    info("#    Simon says: Thanks for using LUTR!    #")
    info("############################################")