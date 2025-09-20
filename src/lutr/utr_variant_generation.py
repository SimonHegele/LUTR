"""
Module Name:    setup_logging.py
Author:         Simon Hegele
Date:           2025-09-19
Version:        1.0
License:        GPL-3

Description:    Provides function generate_utr_variant for UT-variant creation 
"""

from pandas  import concat, DataFrame, Series

from logging import debug

try:
    from .gffutils import attributes_dict, attributes_str, gff_columns
except:
    from gffutils import attributes_dict, attributes_str, gff_columns
    
def get_UTRvar_attrs(p_tran:   Series,
                     a_tran:   Series,
                     gene_id:  str,
                     tran_id:  str) -> dict[str, str]:
    
    attrs = attributes_dict(p_tran)
    
    attrs["gene_id"]       = gene_id
    attrs["transcript_id"] = tran_id
    attrs["ID"]            = tran_id 
    attrs["predictor"]     = p_tran["source"]
    attrs["assembler"]     = a_tran["source"]
    attrs["predicted"]     = str(p_tran["start"]) + "-" + str(p_tran["end"])
    attrs["assembled"]     = str(a_tran["start"]) + "-" + str(a_tran["end"])
    
    return attrs

def get_UTRvar(p_tran: Series, a_tran: Series, tran_id: str, gene_id: str) -> DataFrame:

        return DataFrame({"seqname":    [p_tran["seqname"]],
                          "source":     ["LUTR"],
                          "type":       [p_tran["type"]],
                          "start":      [min(a_tran["start"],p_tran["start"])],
                          "end":        [max(a_tran["end"],p_tran["end"])],
                          "score":      ["."],
                          "strand":     [p_tran["strand"]],
                          "frame":      ["."],
                          "attributes": [attributes_str(get_UTRvar_attrs(p_tran, a_tran, gene_id, tran_id))]},
                         columns=gff_columns)

def is_assembled(a_tran: Series, exon_start: int, exon_end: int) -> bool:
    """
    Args:
        a_tran (Series):  GFF-row of the assembled transcript used for UTR-variant generation
        exon_start (int): Start position of an exon from the UTR-variant
        exon_end (int):   End position of an exon from the UTR-variant

    Returns:
        bool: True  if the exon of the UTR-variant was (at least partially) assembled,
                    i.e. it overlaps the assembled transcript
              False else
    """
    
    if (a_tran["start"]<=exon_end) and (a_tran["end"]>=exon_start):
        return True
    if (exon_start<=a_tran["end"]) and (exon_end>=a_tran["start"]):
        return True
    return False
            
def is_predicted(p_tran: Series, exon_start: int, exon_end: int) -> bool:
    """
    Args:
        p_tran (Series):  GFF-row of the assembled transcript used for UTR-variant generation
        exon_start (int): Start position of an exon from the UTR-variant
        exon_end (int):   End position of an exon from the UTR-variant

    Returns:
        bool: True  if the exon of the UTR-variant was (at least partially) predicted,
                    i.e. it overlaps the predicted transcript
              False else
    """
    if (p_tran["start"]<=exon_end) and (p_tran["end"]>=exon_start):
        return True
    if (exon_start<=p_tran["end"]) and (exon_end>=p_tran["start"]):
        return True
    return False

def exon_attrs(exon_start: int,
               exon_end: int,
               i: int,
               gene_id: str,
               tran_id: str,
               p_tran: Series,
               a_tran: Series) -> str:
            
    attrs = dict()
    
    attrs["ID"]             = f"exon-{tran_id}-{i}"
    attrs["gene_id"]        = gene_id
    attrs["Parent"]         = tran_id
    attrs["transcript_id"]  = tran_id
    
    if is_assembled(a_tran, exon_start, exon_end):
        assembly_start     = max(exon_start, a_tran["start"])
        assembly_end       = min(exon_end,   a_tran["end"])
        attrs["assembler"] = a_tran["source"]
        attrs["assembled"] = f"{assembly_start}-{assembly_end}"
        
    if is_predicted(a_tran, exon_start, exon_end):
        prediction_start     = max(exon_start, p_tran["start"])
        prediction_end       = min(exon_end,   p_tran["end"])
        attrs["predictor"]   = p_tran["source"]
        attrs["predicted"]   = f"{prediction_start}-{prediction_end}"
        
    return attributes_str(attrs)
        
def get_exons(p_tran: Series, 
              a_tran: Series,
              p_tran_exons: DataFrame,
              a_tran_exons: DataFrame,
              gene_id: str,
              tran_id: str):
        
        # Case 1: Assembled transcript extends predicted transcripts to both sides
        if a_tran["start"] < p_tran["start"] and p_tran["end"] < a_tran["end"]:
            debug("Assembled transcript extends predicted transcripts to both sides")
            exon_starts = list(a_tran_exons["start"])
            exon_ends   = list(a_tran_exons["end"])
        # Case 2: Assembled transcript extends predicted transcripts only to the left
        elif a_tran["start"] < p_tran["start"] and a_tran["end"] <= p_tran["end"]:
            debug("Assembled transcript extends predicted transcripts only to the left")
            exon_starts  = list(a_tran_exons.loc[a_tran_exons["start"] < p_tran["start"]]["start"]) 
            exon_starts += list(p_tran_exons.loc[p_tran_exons["start"] >= p_tran["start"]]["start"])
            exon_ends  = list(a_tran_exons.loc[a_tran_exons["start"] < p_tran["start"]]["end"]) 
            exon_ends += list(p_tran_exons.loc[p_tran_exons["start"] >= p_tran["start"]]["end"])
        # Case 3: Assembled transcript extends predicted transcripts only to the right
        elif p_tran["start"] <= a_tran["start"] and p_tran["end"] < a_tran["end"]:
            debug("Assembled transcript extends predicted transcripts only to the right")
            exon_starts  = list(a_tran_exons.loc[a_tran_exons["end"] > p_tran["end"]]["start"]) 
            exon_starts += list(p_tran_exons.loc[p_tran_exons["end"] <= p_tran["end"]]["start"])
            exon_ends    = list(a_tran_exons.loc[a_tran_exons["end"] > p_tran["end"]]["end"]) 
            exon_ends   += list(p_tran_exons.loc[p_tran_exons["end"] <= p_tran["end"]]["end"])
        else:
            raise Exception(f"UTR-variant cannot be created from {list(p_tran)} and {list(a_tran)}!")
        
        exon_starts = sorted(exon_starts)
        exon_ends   = sorted(exon_ends)

        return DataFrame({"seqname":    [p_tran["seqname"]  for _ in exon_starts],
                          "source":     ["LUTR"             for _ in exon_starts],
                          "type":       ["exon"             for _ in exon_starts],
                          "start":      exon_starts,
                          "end":        exon_ends,
                          "score":      ["."                for _ in exon_starts],
                          "strand":     [p_tran["strand"]   for _ in exon_starts],
                          "frame":      ["."                for _ in exon_starts],
                          "attributes": [exon_attrs(exon_start, exon_end, i, gene_id, tran_id, p_tran, a_tran)
                                         for i, (exon_start, exon_end) in enumerate(zip(exon_starts, exon_ends))]},
                         columns=gff_columns)
        
def get_non_exons(p_tran_slice: DataFrame,
                  p_tran:       Series,
                  gene_id:      str,
                  tran_id:      str) -> DataFrame:
        
    non_exons = p_tran_slice.loc[(p_tran_slice["type"] != "exon")]
    non_exons = non_exons[~(non_exons == p_tran).all(axis=1)].reset_index(drop=True)
    
    for i, non_exon in non_exons.iterrows():
        
        attrs = attributes_dict(non_exon)
        type  = non_exon["type"]
        
        attrs["ID"]             = f"{type}-{tran_id}-{i}"
        attrs["gene_id"]        = gene_id
        attrs["Parent"]         = tran_id
        attrs["transcript_id"]  = tran_id
        
        non_exons.iloc[i, 8] = attributes_str(attrs)

    return non_exons

def generate_utr_variant(p_tran_slice: DataFrame,
                         a_tran_slice: DataFrame,
                         p_tran_exons: DataFrame,
                         a_tran_exons: DataFrame,
                         p_tran:       Series,
                         a_tran:       Series,
                         index:        int,
                         gene_id:      str) -> DataFrame:
    """
    Creating an UTR-variant for a predicted transcript by transfering extending exons from
    a matching assembled transcript.

    Args:
        p_tran_slice (DataFrame): GFF-slice, features of a predicted transcript
        a_tran_slice (DataFrame): GFF-slice, features of a matching assembled transcript
        p_tran_exons (DataFrame): 
        a_tran_exons (DataFrame): 
        p_tran (Series):          GFF-row of the predicted transcript
        a_tran (Series):          GFF-row of the assembled transcript
        index (int):              Index of the UTR-variant 0,1,...
        gene_id (str):            ID of the gene of the predicted transcript

    Returns:
        DataFrame: GFF-slice, features of the UTR-variant
    """
    
    debug("Creating UTR-Variant from:")
    debug(f"Predicted: {p_tran_slice}")
    debug(f"Assembled: {a_tran_slice}")
    debug(f"Predicted exons:\n{p_tran_exons}")
    debug(f"Assembled exons:\n{a_tran_exons}")
    
    tran_id = attributes_dict(p_tran)["ID"] + f"_UTRvar_{index}"
    
    transcript = get_UTRvar(p_tran, a_tran, tran_id, gene_id)
    exons      = get_exons(p_tran, a_tran, p_tran_exons, a_tran_exons, gene_id, tran_id)
    non_exons  = get_non_exons(p_tran_slice, p_tran, gene_id, tran_id)
  
    df = concat([transcript, exons, non_exons])

    return df