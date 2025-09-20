from argparse   import ArgumentParser, RawTextHelpFormatter
from logging    import error
from subprocess import run

class LUTRprepArgparser(ArgumentParser):

    prog        =  "assembly"

    description =   """
                    -------------------------------------------------------------------
                    Preparing assembly annotation for LUTR
                    
                    1. Alignment:
                        a) Short reads: STAR
                        b) Long reads:  Minimap2
                    2. Alignment conversion and sorting:
                        Samtools
                    3. Assembly:
                        StringTie
                    4. Conversion from GTF to GFF
                        AGAT
                   -------------------------------------------------------------------
                    """
                       
    def __init__(self) -> None:

        super().__init__(prog=self.prog,
                         description=self.description,
                         formatter_class=RawTextHelpFormatter)
        
        self.add_argument("GFF")
        
        grp1 = self.add_argument_group(title="Input")
        grp1.add_argument("-sr","--short_reads",
                          help="Illumina reads ( .fasta / .fastq)",
                          nargs="+",
                          metavar="")
        grp1.add_argument("-sb","--short_bam",
                          help="Illumina read alignments, sorted ( .bam )",
                          metavar="")
        grp1.add_argument("-lr","--long_reads",
                          help="Long reads ( .fasta / .fastq)",
                          metavar="")
        grp1.add_argument("-lb","--long_bam",
                          help="Illumina read alignments, sorted ( .bam )",
                          metavar="")
        grp1.add_argument("-g","--genome",
                          help="Reference genome ( .fasta )",
                          metavar="")
        grp1.add_argument("-si","--star_idx",
                          help="Preconstructed STAR index of your genome for short read mapping",
                          metavar="")
        grp1.add_argument("-mi","--mm2_idx",
                          help="Preconstructed Minimap2 index of your genome for long read mapping",
                          metavar="")
        
        grp3 = self.add_argument_group("StringTie")
        grp3.add_argument("-j","--junction_coverage",
                          help="Minimum coverage to accept a splice junction [default: 2]",
                          type=int,
                          metavar="",
                          default=2)
        grp3.add_argument("-c","--bp_coverage",
                          help="Minimum base pair coverage to accept a transcript [default: 2]",
                          type=int,
                          metavar="",
                          default=2)
        
        grp3 = self.add_argument_group("Performance")
        grp3.add_argument("-t", "--threads",
                          help="Number of parallel processes to use [Default:8]",
                          type=int,
                          metavar="",
                          default=8)
        
        grp4 = self.add_argument_group("Others")
        grp4.add_argument("-o", "--outdir",
                          help="Path to output directory [default: \".\"]",
                          metavar="")
        
    def check_input(self):
        
        if (not self.args.short_reads) and (not self.args.short_bam):
            error("Options --short_reads and --short_bam are mutually exclusive")
            exit(1)
        if (not self.args.long_reads ) and (not self.args.long_bam ):
            error("Options --long_reads and --long_bam are mutually exclusive")
            exit(1)
        if self.args.short_reads or self.args.long_reads:
            if self.args.short_reads and not (self.args.genome or self.args.star_index):
                error("Genome or STAR index needed for short reads")
            if self.args.long_reads and not (self.args.genome or self.args.minimap2_index):
                error("Genome or Minimap2 index needed for long reads")
                
    def parse_args(self):

        self.args = super().parse_args()

        self.check_input()  
        
        return self.args
    
def cmd(command: str):
    
    run(command,
        shell=True,
        check=True)
    
def assembly_type(short: str | None, long: str | None) -> str:
    
    if not short is None:
        if not long is None:
            return "hybrid"
        else:
            return "short"
    else:
        return "long"
 
def main():
    
    args    = LUTRprepArgparser().parse_args()
    gff     = args.GFF
    outdir  = args.outdir
    threads = args.threads
    
    # Short read alignment 
    if args.short_reads:
        
        # Mapping
        if args.star_idx is None:
            star_index = outdir + "/star_index"   
            cmd(f"STAR --runThreadN {threads} --genomeDir {star_index} --runMode genomeGenerate --genomeFastaFiles {args.genome}")
        else:
            star_index = args.star_idx
        input  = " ".join(args.short_reads)
        output = outdir + "/star"
        if input.endswith(".gz"):
            cmd(f"STAR --runThreadN {threads} --genomeDir {star_index} --outFileNamePrefix {output} --outSAMtype BAM Unsorted --readFilesIn {input} --readFilesCommand zcat")
        else:
            cmd(f"STAR --runThreadN {threads} --genomeDir {star_index} --outFileNamePrefix {output} --outSAMtype BAM Unsorted --readFilesIn {input}")
        if not args.star_idx is None:
            cmd(f"rm -r {args.star_idx}")

        # SAM -> BAM conversion
        input  = outdir + "/starAligned.out.bam"
        output = outdir + "/short.sorted.bam"
        cmd(f"samtools sort -@ {threads} -o {output} {input}")
        cmd(f"rm {input}")
        
        short = output
        
    elif args.short_bam:
        
        long = args.short_bam
        
    # Long reads
    if args.long_reads:
        
        # Mapping
        input  = args.long_reads
        output = outdir + "/long.sam"
        cmd(f"minimap2 -ax splice -o {output} -t {threads} {args.genome} {input}")

        # SAM -> BAM conversion
        input  = outdir + "/long.sam"
        output = outdir + "/long.bam"
        cmd(f"samtools view -Sb -o {output} {input}")
        cmd(f"rm {input}")
        
        # BAM sorting
        input  = outdir + "/long.bam"
        output = outdir + "/long.sorted.bam"
        cmd(f"samtools sort -@ {threads} -o {output} {input}")
        cmd(f"rm {input}")
        
        long = output
    
    elif args.long_bam:
        
        long = args.long_bam
           
    # Assembly
    output = outdir + "/stringtie.gtf"
    match assembly_type(short, long):
        case "short":
            cmd(f"stringtie -G {gff} -c {args.bp_coverage} -j {args.junction_coverage} -o {output} -p {threads} {short}")
            cmd(f"rm {short}")
        case "long":
            cmd(f"stringtie -G {gff} -c {args.bp_coverage} -j {args.junction_coverage} -o {output} -p {threads} -L {long}")
            cmd(f"rm {long}")
        case "hybrid":
            cmd(f"stringtie -G {gff} -c {args.bp_coverage} -j {args.junction_coverage}  -o {output} -p {threads} --mix {short} {long}")
            cmd(f"rm {short}")
            cmd(f"rm {long}")
            
    # GTF -> GFF conversion
    input  = outdir + "/stringtie.gtf"
    output = outdir + "/assembly.gff"
    cmd(f"agat_convert_sp_gxf2gxf.pl --gtf {input} -o {output}")
    cmd(f"rm {input}")