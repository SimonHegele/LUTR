"""
Module Name:    gffutils.py
Author:         Simon Hegele
Date:           2025-09-19
Version:        1.0
License:        GPL-3
Description:    Functions and custom Exceptions for pandas.Dataframe represented GFF-files
                (alphabetically sorted)
"""

from collections     import defaultdict
from itertools       import chain
from pandas          import DataFrame, Series, read_csv
from typing          import Callable, Generator, Iterable, Hashable


gff_columns = ["seqname",
               "source",
               "type",
               "start",
               "end",
               "score",
               "strand",
               "frame",
               "attributes"]

class AttributesError(Exception):
    
    def __init__(self, feature) -> None:
        
        error_msg = f"Invalid format of attributes for feature: {feature_string(feature)}"
        
        super().__init__(error_msg)
        
class NoIDError(Exception):
    
    def __init__(self, feature) -> None:
        
        error_msg = f"No ID for feature: {feature_string(feature)}"
        
        super().__init__(error_msg)

class MultipleParentsError(Exception):
    
    def __init__(self,
                 gff:     DataFrame,
                 feature: Series,
                 parents: list[Hashable]):
        
        error_msg  = f"Multiple Parents:\n"
        error_msg += f"Feature:  {feature_string(feature)}"
        
        for i, parent_index in enumerate(parents):
            error_msg += f"Parent i: {feature_string(gff.iloc[parent_index])}"
            
        super().__init__(error_msg)

def attributes_dict(feature: Series) -> dict[str, str]:
    """
    Parsing the key=value pairs from a features attributes fields into a hashmap
    """
    
    try:
        attributes = {a.split("=")[0]: a.split("=")[1]
                      for a in feature["attributes"].split(";")} 
    except:
        raise AttributesError(feature)
    
    if "ID" in attributes.keys():
        return attributes
    else:
        raise NoIDError(feature)

def attributes_str(attributes: dict[str, str]) -> str:
    
    return ";".join([f"{key}={attributes[key]}" for key in attributes.keys()])

def check_strands(feature1: Series,
                  feature2: Series,
                  know_strand=False) -> bool:

    if know_strand:
        if not feature1["strand"]=="." or feature2["strand"]==".":
            if feature1["strand"]==feature2["strand"]:
                return True
    else:
        if feature1["strand"]=="." or feature2["strand"]==".":
            return True
        if feature1["strand"]==feature2["strand"]:
            return True
    return False

def empty_gff() -> DataFrame:
    """
    Returns:
        DataFrame: An empty GFF-file
    """
    return DataFrame({}, columns = gff_columns)

def get_map_id2index(gff: DataFrame) -> defaultdict[str, list[Hashable]]:
    
    map_id2index = defaultdict(list)
        
    for i, feature in gff.iterrows():
        
        map_id2index[feature["ID"]].append(i)
        
    return map_id2index

def get_map_parent2children(gff: DataFrame,
                            map_id2index: defaultdict[str, list[Hashable]]) -> defaultdict[Hashable, list[Hashable]]:
    
    map_parent2children = defaultdict(list)
    
    for i, feature in gff.iterrows():
        
        if feature["Parent"] is None:
            continue
        
        parents = map_id2index[feature["Parent"]]
        
        if len(parents) > 1:
            raise MultipleParentsError(gff, feature, parents)
        if len(parents) == 1:
            map_parent2children[parents[0]].append(i)
            
    return map_parent2children
            
def get_subtree(gff:                 DataFrame,
                index:               Hashable,
                map_parent2children: defaultdict[Hashable, list[Hashable]]) -> list:
    
    children = map_parent2children[index]
    
    if children == 0:
        return []
    
    return [index] + list(chain.from_iterable([get_subtree(gff, child, map_parent2children)
                                               for child in children]))

def get_trans(gff: DataFrame):
    
    return gff.loc[(gff["type"]=="transcript") |
                   (gff["type"].str.contains("RNA"))].reset_index(drop=True)

def feature_length(feature: Series) -> int:
    
    return int(feature["end"] - feature["start"])

def feature_includes(feature_1: Series,
                     feature_2: Series) -> bool:
    """
    Checks if the genomic location of feature_1 includes the genomic location of feature_2
    """

    if not feature_1["seqname"] == feature_2["seqname"]:
        return False
    if not feature_1["start"] <= feature_2["start"]:
        return False
    if not feature_1["end"] >= feature_2["end"]:
        return False
    return True

def features_overlap(feature_1: Series,
                     feature_2: Series) -> bool:
    """
    Checks if two features overlap
    """

    if feature_1["seqname"] == feature_2["seqname"]:
        if (feature_1["start"]<=feature_2["end"]) and (feature_1["end"]>=feature_2["start"]):
            return True
        if (feature_2["start"]<=feature_1["end"]) and (feature_2["end"]>=feature_1["start"]):
            return True
    return False

def feature_pairs(gff_1: DataFrame,
                  gff_2: DataFrame,
                  pairing: Callable[[DataFrame, Series],DataFrame],
                  type="") -> Generator:
    
    gff_1_features = gff_1.loc[gff_1["type"].str.contains(type, regex=False)]
    gff_2_features = gff_2.loc[gff_2["type"].str.contains(type, regex=False)]

    for i, gff_1_feature in gff_1_features.iterrows():

        for j, gff_2_feature in pairing(gff_2_features, gff_1_feature).iterrows():

            yield gff_1_feature, gff_2_feature
            
def feature_string(feature: Series) -> str:
    return "\t".join([str(c) for c in list(feature)])
    
def load_gff(file_path: str) -> DataFrame:
    """
    Loading a GFF-file from the file-system
    """
    
    gff           = read_csv(file_path, sep="\t", header=None, comment="#", names=gff_columns)
    
    # We will make extensive use of the ID and Parent attributes so we will add them as 
    # columns to avoid repeated parsing of the attributes string
    attrs = [attributes_dict(feature) for _, feature in gff.iterrows()]
    gff["ID"]     = [a["ID"] for a in attrs]
    gff["Parent"] = [a["Parent"] if "Parent" in a.keys() else None for a in attrs]
    
    return gff

def included_features(gff: DataFrame,
                      feature: Series,
                      type="") -> DataFrame:
    """
    Returns features of the specified type from the input GFF included by the input feature 
    """
    mask_seqname  = gff["seqname"] == feature["seqname"]
    mask_type     = gff["type"].str.contains(type, regex=False)
    prefiltered   = gff[mask_seqname & mask_type]
    mask_includes = prefiltered.apply(lambda f: feature_includes(feature, f), axis=1)

    return prefiltered[mask_includes]

def including_features(gff: DataFrame,
                      feature: Series,
                      type="") -> DataFrame:
    """
    Returns features of the specified type from the input GFF included by the input feature 
    """
    mask_seqname  = gff["seqname"] == feature["seqname"]
    mask_type     = gff["type"].str.contains(type, regex=False)
    prefiltered   = gff[mask_seqname & mask_type]
    mask_includes = prefiltered.apply(lambda f: feature_includes(f, feature), axis=1)

    return prefiltered[mask_includes]
    
def overlap_length(feature_1: Series,
                   feature_2: Series) -> int:
    """
    Checks if two features overlap
    """

    if feature_1["seqname"] == feature_2["seqname"]:
        start = min(feature_1["end"],feature_2["end"])
        end   = max(feature_1["start"],feature_2["start"])
        return max(0, end-start)
    
    return 0

def overlapping_features(gff: DataFrame,
                         feature: Series,
                         type="") -> DataFrame:
    """
    Returns features of the specified type from the input GFF overlapping the input feature 
    """
    
    mask_seqname = gff["seqname"] == feature["seqname"]
    mask_type    = gff["type"].str.contains(type, regex=False)
    prefiltered  = gff[mask_seqname & mask_type]
    mask_overlap = prefiltered.apply(lambda f: features_overlap(feature, f), axis=1)

    return prefiltered[mask_overlap]

def seqname_split(gff: DataFrame,
                  seqnames=None) -> dict[str, DataFrame]:

    if seqnames is None:
        seqnames = gff["seqname"].unique()
        
    splits = dict()
        
    for seqname in seqnames:
        
        splits[seqname] = gff.loc[gff["seqname"]==seqname].sort_values("start").reset_index(drop=True)
        gff.drop(splits[seqname].index)

    return splits
    
def seqname_split_wrapper(args: tuple[DataFrame, Iterable[str] | None]) -> dict[str, DataFrame]:
    
    return seqname_split(args[0], args[1])
    
def to_string(feature: Series):

    return "\t".join(str(v) for v in list(feature))

def type_split(gff: DataFrame, type: str) -> Generator[DataFrame, None, None]:
    
    features = gff.loc[gff["type"].str.contains(type, regex=False)]
    map_p2ch = get_map_parent2children(gff, get_map_id2index(gff))
    
    for i, feature in features.iterrows():

        yield gff.iloc[get_subtree(gff, i, map_p2ch)]

def write_gff(gff: DataFrame, file_path: str, mode='w'):

    gff.to_csv(file_path, sep="\t", index=False, header=None, mode=mode, columns=gff_columns)