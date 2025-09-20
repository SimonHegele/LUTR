"""
Module Name:    transcript_matching.py
Author:         Simon Hegele
Date:           2025-09-19
Version:        1.0
License:        GPL-3

Description:    Provides functions:
                    tmatch() to determine if and for how many bases predicted and
                    assembled transcript match.
"""

from logging  import warning
from pandas   import DataFrame, Series

from .gffutils import overlapping_features, features_overlap, overlap_length

def sufficently_covered(covered: int,
                        p_exons: DataFrame,
                        mftm:    float,
                        mtbm:    int) -> bool:
    """
    Args:
        covered (int):       The number of bases of a predicted transcript that are
                             covered by a matching assembled transcript.
        p_exons (DataFrame): GFF-slice with the exons of the predicted transcript to
                             calculate its length.
        mftm (float):        Minimum fraction the predicted transcript to be covered.
                             Applied only to 
        mtbm (int):          Minimum  bases the predicted transcript to be covered.

    Returns:
        bool: True if the predicted transcript is sufficiently covered by the assembled one
    """

    p_tran_len   = sum([e["end"]-e["start"] for _, e in p_exons.iterrows()])
    min_coverage = min(p_tran_len, max((mftm * p_tran_len), mtbm))
    
    return covered >= min_coverage

def match_single_exon(p_exons: DataFrame, a_exons: DataFrame) -> int:
    
    if len(p_exons) == 1:
        single_exon    = p_exons.iloc[0]
        matching_exons = overlapping_features(a_exons, single_exon).reset_index(drop=True)
    elif len(a_exons) == 1:
        single_exon    = a_exons.iloc[0]
        matching_exons = overlapping_features(p_exons, single_exon).reset_index(drop=True)
    else:
        raise Exception("Illegal function call, neither transcript is single-exon!")
    
    if not len(matching_exons) == 1:
        return 0
    
    overlap_start = max(single_exon["start"], matching_exons.iloc[0]["start"])
    overlap_end   = min(single_exon["end"],   matching_exons.iloc[0]["end"])
    
    return overlap_end -overlap_start

def match_multi_exon(p_exons:           DataFrame,
                     a_exons:           DataFrame, 
                     no_match_middle:   bool) -> int:
    
    # Check 1: Transcript start and ends
    #          a) For the first exon of one of the two transcripts there is an exon with
    #             the same end in the other transcript
    #          b) For the last exon of one of the two transcripts there is an exon with
    #             the same start in the other transcript
    #          c) An exon whose end was matched in a) starts after the corresponding exon
    #             in the other transcript must not have preceeding exons
    #          d) An exon whose end was matched in b) ends before the corresponding exon
    #             in the other transcript must not have succeeding exons
    
    first_exon_end_to_match  = int(max(p_exons.iloc[0]["end"],    a_exons.iloc[0]["end"]))
    last_exon_start_to_match = int(min(p_exons.iloc[-1]["start"], a_exons.iloc[-1]["start"]))
    p_start                  = p_exons.loc[p_exons["end"]   == first_exon_end_to_match]
    a_start                  = a_exons.loc[a_exons["end"]   == first_exon_end_to_match]
    p_end                    = p_exons.loc[p_exons["start"] == last_exon_start_to_match]
    a_end                    = a_exons.loc[a_exons["start"] == last_exon_start_to_match]

    # a)
    if (len(p_start) != 1) or (len(a_start) != 1):
        return 0
    
    # b)
    if (len(p_end) != 1) or (len(a_end) != 1):
        return 0
    
    # c)
    if p_start.iloc[0]["start"] > a_start.iloc[0]["start"]:
        # p_start must be the first exon in p_exons
        if p_exons.index.get_loc(p_start.index[0]) != 0:
            return 0
    elif a_start.iloc[0]["start"] > p_start.iloc[0]["start"]:
        # a_start must be the first exon in a_exons
        if a_exons.index.get_loc(a_start.index[0]) != 0:
            return 0
    
    # d) 
    if p_end.iloc[0]["end"] < a_end.iloc[0]["end"]:
        # p_end must be the last exon in p_exons
        if p_exons.index.get_loc(p_end.index[0]) != len(p_exons) - 1:
            return 0
    elif a_end.iloc[0]["end"] < p_end.iloc[0]["end"]:
        # a_end must be the last exon in a_exons
        if a_exons.index.get_loc(a_end.index[0]) != len(a_exons) - 1:
            return 0

    # Check 2: Check middle exons and compute number of shared bases
    middle_p_exons = p_exons.loc[(p_exons["start"] > first_exon_end_to_match) &
                                 (p_exons["end"]   < last_exon_start_to_match)].reset_index(drop=True)
    middle_a_exons = a_exons.loc[(a_exons["start"] > first_exon_end_to_match) &
                                 (a_exons["end"]   < last_exon_start_to_match)].reset_index(drop=True)
    
    first_matched_start = max(p_exons.iloc[0]["start"],  a_exons.iloc[0]["start"])
    last_matched_end    = min(p_exons.iloc[-1]["end"],   a_exons.iloc[-1]["end"])
    
    first_exon_covered   = first_exon_end_to_match - first_matched_start
    last_exon_covered    = last_matched_end - last_exon_start_to_match
    
    if no_match_middle:
        
        middle_exons_covered = 0
        
        i = 0
        j = 0
        
        while i < len(middle_p_exons) and j < len(middle_a_exons):
            
            middle_exon_covered = overlap_length(middle_p_exons.iloc[i],middle_a_exons.iloc[j])
            
            if middle_exons_covered > 0:
                
                middle_exons_covered += middle_exon_covered
                
            if middle_p_exons.iloc[i]["start"] <= middle_a_exons.iloc[j]["start"]:
                i += 1
            if middle_p_exons.iloc[i]["start"] >= middle_a_exons.iloc[j]["start"]:
                j += 1
            
    else:
    
        if not middle_p_exons['start'].equals(middle_a_exons['start']):
            return 0
        if not middle_p_exons['end'].equals(middle_a_exons['end']):
            return 0
    
        middle_exons_covered = sum([exon["end"]-exon["start"]
                                    for _, exon in middle_p_exons.iterrows()])
        
    covered              = first_exon_covered + last_exon_covered + middle_exons_covered
    
    return covered

def tmatch(p_tran:       Series,
           a_tran:       Series,
           p_exons:      DataFrame,
           a_exons:      DataFrame,
           mftm:         float,
           mtbm:         int,
           match_middle: bool) -> int:
    
    # Pre-check 1: Both transcripts have exons
    if len(p_exons) == 0:
        id = p_tran["ID"]
        warning(f"Predicted transcript {id} has no exons (Treated as single exon transcript)")
        p_exons = DataFrame([p_tran]).reset_index(drop=True)
    if len(a_exons) == 0:
        id = "ID"
        warning(f"Assembled transcript {id} has no exons (Treated as single exon transcript)")
        a_exons = DataFrame([a_tran]).reset_index(drop=True)

    # Single exon transcripts
    if len(p_exons) == 1 or len(a_exons) == 1:
        covered = match_single_exon(p_exons,a_exons)
    # Multi-exon transcripts
    else:
        covered = match_multi_exon(p_exons, a_exons, match_middle)
        
    if sufficently_covered(covered, p_exons, mftm, mtbm):
        return covered
    else:
        return 0
    
##########################################################################################
# END OF FILE                                                                            #
##########################################################################################
    
def transcript_match_matrix(trans_predicted: list[Series],
                            trans_assembled: list[Series],
                            exons_predicted: list[DataFrame],
                            exons_assembled: list[DataFrame],
                            min_exons_match: int,
                            min_fract_match: int,
                            min_bases_match: int,
                            best_match_only: int):
    
    def transcript_match(trans_predicted: list[Series],
                         trans_assembled: list[Series],
                         exons_predicted: list[DataFrame],
                         exons_assembled: list[DataFrame],
                         min_exons_match: int,
                         min_fract_match: int,
                         min_bases_match: int):
        
        pass
    
    match_matrix = [[transcript_match(trans_predicted,
                                      trans_assembled,
                                      exons_predicted,
                                      exons_assembled,
                                      min_exons_match,
                                      min_fract_match,
                                      min_bases_match,)
                    for j, a_tran in enumerate(trans_assembled)]
                    for i, p_tran in enumerate(trans_predicted)]