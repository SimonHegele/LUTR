"""
Module Name:    argumentparser.py
Description:    Provides class LUTRArgparser(ArgumentParser)
                - Arguments for LUTR added on initialization
                - Extended parse_args() to 
                    a) Check input 
                    b) Create the output directory
                    c) Write a file with parameter selection to the output directory
Author:         Simon Hegele
Date:           2025-09-19
Version:        1.0
License:        GPL-3
"""

from argparse import ArgumentParser, RawTextHelpFormatter
from logging  import error
from os       import mkdir, path

class LUTRArgparser(ArgumentParser):

    prog        =  "LUTR"

    description =   """
                    ----------------------------------------------------------------------------------------------------
                      ██░   ██░   ██░   ██░   ██░   ██░   ██░   ██░   ██░   ██░   ██░   ██░   ██░   ██░   ██░   ██░   █░
                     █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ 
                      ██░   ██░   ██░   ██░   ██░                                          ██░   ██░   ██░   ██░   ██░
                     █░ █░ █░ █░ █░ █░ █░ █░ █░ █░   █░      █░     █░███████░  ██░       █░ █░ █░ █░ █░ █░ █░ █░ █░ █░
                    █░   ██░   ██░   ██░   ██░   █░  |█░    █|█░   █░   |█░   |█░ █░     █░   ██░   ██░   ██░   ██░   █░  
                     █░ █░ █░ █░ █░ █░ █░ █░ █░ █░   |█░     |█░   █░   |█░   |█░  █░     █░ █░ █░ █░ █░ █░ █░ █░ █░ █░   
                      ██░   ██░   ██░   ██░   ██░    |█░     |█░   █░   |█░   |█░ █░       ██░   ██░   ██░   ██░   ██░  
                      ██░   ██░   ██░   ██░   ██░    |█░     |█░   █░   |█░   |█░█░        ██░   ██░   ██░   ██░   ██░  
                     █░ █░ █░ █░ █░ █░ █░ █░ █░ █░   |█░     |█░   █░   |█░   |█░ █░      █░ █░ █░ █░ █░ █░ █░ █░ █░ █░
                    █░   ██░   ██░   ██░   ██░   █░  |█░   █░|█░   █░   |█░   |█░  █░█░  █░   ██░   ██░   ██░   ██░   █░ 
                     █░ █░ █░ █░ █░ █░ █░ █░ █░ █░   █ ████░   ████░    |█░  |█░    █░    █░ █░ █░ █░ █░ █░ █░ █░ █░ █░
                      ██░   ██░   ██░   ██░   ██░                                          ██░   ██░   ██░   ██░   ██░
                     █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ █░ 
                      ██░   ██░   ██░   ██░   ██░   ██░   ██░   ██░   ██░   ██░   ██░   ██░   ██░   ██░   ██░   ██░   █░  
                    ----------------------------------------------------------------------------------------------------
                    
                    UTR-extensions for transcripts in functional annotations from orthology-based gene prediction tools
                    using matching transcripts in structural annotations from reference-based transcriptome assembly.
                    
                    Additional parameter explanation:
                    *1 By default, LUTR only allows transcript matches in which all other exons between the first and
                       last matched exon must also be matched exactly. This is based on the assumption that UTRs can
                       affect splicing within the coding regions of a transcript, e.g. by changing its secondary
                       structure.
                    *2 By default, LUTR assignes assembled transcripts only to the matching predicted transcripts they
                       share the most bases with. Assigning them to all matching predicted transcripts can increase the
                       number of annotated UTRs but bears the risk of an increased rate of false positives.
                    ----------------------------------------------------------------------------------------------------
                    """
            
    def __init__(self) -> None:

        super().__init__(prog=self.prog,
                         description=self.description,
                         formatter_class=lambda prog: RawTextHelpFormatter(prog, max_help_position=35, width=100))

        # Input files
        self.add_argument("prediction",
                          help="Annotation from gene prediction (GFF)")
        self.add_argument("assembly",
                          help="Annotation from transcriptome assembly (GFF)")
        self.add_argument("outdir",
                          help="Output directory (Must not exist already)")
        
        grp1 = self.add_argument_group(title="Transcript matching")
        grp1.add_argument("-mf","--mftm",
                          help="Minimum fraction of the predicted transcript to be matched [default: 0.9]",
                          metavar="",
                          type=float,
                          default=0.9,)
        grp1.add_argument("-mb","--mtbm",
                          help="Minimum bases of the predicted transcript to be matched [default: 10,000]",
                          metavar="",
                          type=int,
                          default=10000,)
        grp1.add_argument("-nm","--no_match_middle",
                          action="store_true",
                          help="Allows missing / additional middle exons in the transcrpt matching *1")
        grp1.add_argument("-ma","--match_all",
                          action="store_true",
                          help="Assign assembled transcripts to all matching predicted transcripts *2")
        
        grp2 = self.add_argument_group(title="Transcript selection")
        grp2.add_argument("-s","--select",
                          help="How to select from multiple UTR-variants [choices: shortest, longest, all] [default: all]",
                          choices=["shortest", "longest", "all"],
                          default="all",
                          metavar="")
        grp1.add_argument("-mupl","--max_utr_piece_len",
                          metavar="",
                          help="Limits the length of UTR-pieces allowed in UTR-variants # TODO")
        grp2.add_argument("-r","--remove_unsupported",
                          action="store_true",
                          help="Remove unmatched transcripts")
        
        grp3 = self.add_argument_group("Performance")
        grp3.add_argument("-t", "--threads",
                          help="Number of parallel processes to use [Default:8]",
                          type=int,
                          metavar="",
                          default=8)
        
        grp4 = self.add_argument_group(title="Others")
        grp4.add_argument("-v","--verbosity",
                          help="Report progress every 10^v genes [default: 0]",
                          type=int,
                          default=0,
                          metavar="")
        grp4.add_argument("-l","--log_level",
                          help="[default: info]",
                          default="info",
                          metavar="")
               
    def write_parameter_file(self):

        with open(path.join(self.args.outdir, "lutr.param"), "w") as param_file:
            
            param_file.write(f"prediction           {self.args.prediction:<15}\n")
            param_file.write(f"assembly             {self.args.assembly:<15}\n")
            param_file.write(f"--mftm               {self.args.mftm:<15}\n")
            param_file.write(f"--mtbm               {self.args.mtbm:<15}\n")
            param_file.write(f"--select             {self.args.select:<15}\n")
            param_file.write(f"--threads            {self.args.threads:<15}\n")
            
            if self.args.remove_unsupported:
                param_file.write(f"--remove_unsupported")

    def check_input(self):
        
        tmpdir = path.join(self.args.outdir,"tmpdir")

        if not path.isfile(self.args.prediction):
            error(f"{self.args.prediction} is not a file!")
            exit(1)
        if not path.isfile(self.args.assembly):
            error(f"{self.args.assembly} is not a file!")
            exit(1)
        if not path.isdir(self.args.outdir):
            mkdir(self.args.outdir)
        if not path.isdir(tmpdir):
            mkdir(path.join(self.args.outdir,"tmpdir"))
            mkdir(path.join(tmpdir,"genes_predicted"))
            mkdir(path.join(tmpdir,"genes_assembled"))
            mkdir(path.join(tmpdir,"genes_lutr"))
        if not 0 < self.args.mftm and self.args.mftm <= 1:
            error("Value for mftm must be from in (0,1]!")
        if not 0 < self.args.mtbm:
            error("Value for mtbm must be positive!")

    def parse_args(self):

        self.args = super().parse_args()

        self.check_input()
        self.write_parameter_file()  
        
        return self.args
