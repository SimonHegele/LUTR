from argparse        import ArgumentParser
from multiprocessing import Pool
from os              import path, mkdir
from pandas          import concat, DataFrame
from subprocess      import run

from lutr.gffutils import *

class MyArgumentParser(ArgumentParser):

    prog        =   "noUTR"

    description =   """
                    Generating UTR-clipped annotations of coding genes
                    Caution: Does not create exons, use agat_convert_gxf2gxf.pl!
                    """
    
    def __init__(self) -> None:

        super().__init__(prog=self.prog, description=self.description)

        self.add_argument("GFF_in")
        self.add_argument("GFF_out")
        self.add_argument("-t","--threads",
                          type=int,
                          default=8)
        
def process(args):

    gff: DataFrame = args[0]
    seqname: str   = args[1]
    gff            = gff.loc[gff["type"].isin(["gene","mRNA","CDS", "exon"])]
    
    # Process genes
    for i, gene_split in enumerate(type_split(gff, "gene")):
        
        if len(gene_split.loc[gene_split["type"]=="mRNA"]) == 0:
            continue
        
        # Process gene feature
        gene                             = gene_split.loc[gene_split["type"]=="gene"]
        gene.loc[gene.index[0], "start"] = gene_split.loc[gene_split["type"]=="CDS"]["start"].min()
        gene.loc[gene.index[0], "end"]   = gene_split.loc[gene_split["type"]=="CDS"]["end"].max()
        
        write_gff(gene, path.join("coding", f"{seqname}.gff"), mode="a")
        
        # process mRNAs
        for j, mRNA_split in enumerate(type_split(gene_split, "mRNA")):
            
            mRNA       = mRNA_split.loc[mRNA_split["type"]=="mRNA"]
            cds        = mRNA_split.loc[mRNA_split["type"]=="CDS"]
            mRNA_start = cds["start"].min()
            mRNA_end   = cds["end"].max()
            exons      = mRNA_split.loc[mRNA_split["type"]=="exon"]
            exons      = exons.loc[(exons["end"]   > mRNA_start) & 
                                   (exons["start"] < mRNA_end)]
            
            mRNA_split = concat([mRNA, cds, exons])
            
            mRNA_split["start"] = [max(mRNA_start, start) for start in mRNA_split["start"]]
            mRNA_split["end"]   = [min(mRNA_end,   end)   for end   in mRNA_split["end"]]
            
            write_gff(mRNA_split, path.join("coding", f"{seqname}.gff"), mode="a")
            
        print(f"{seqname:>10}: {i:6}")

def main():

    args = MyArgumentParser().parse_args()
    print("Loading ...")
    gff  = load_gff(args.GFF_in)
    gff  = seqname_split(gff)
    seqs = sorted(gff.keys())

    mkdir("coding")

    with Pool(args.threads) as pool:

        pool.map(process, zip([gff[s] for s in seqs], seqs))
        
    run(f"cat coding/* > {args.GFF_out}",
        shell=True,
        check=True)

if __name__ == "__main__":
    main()
