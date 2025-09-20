from os import path

from src.lutr.transcript_matching import tmatch
from src.lutr.gffutils            import load_gff

def test_match():

    data_path = path.join(path.dirname(__file__), "test_data_transcriptmatching")
    
    predicted_exons   = load_gff(path.join(data_path,"predicted_transcript_exons.gff"))
    assembled_exons_1 = load_gff(path.join(data_path,"assembled_transcript_1_exons.gff"))
    assembled_exons_2 = load_gff(path.join(data_path,"assembled_transcript_2_exons.gff"))
    assembled_exons_3 = load_gff(path.join(data_path,"assembled_transcript_3_exons.gff"))
    assembled_exons_4 = load_gff(path.join(data_path,"assembled_transcript_4_exons.gff"))
    assembled_exons_5 = load_gff(path.join(data_path,"assembled_transcript_5_exons.gff"))
    assembled_exons_6 = load_gff(path.join(data_path,"assembled_transcript_6_exons.gff"))
    assembled_exons_7 = load_gff(path.join(data_path,"assembled_transcript_7_exons.gff"))

    assert     tmatch(predicted_exons, assembled_exons_1, 0.8, 1000)
    assert not tmatch(predicted_exons, assembled_exons_2, 0.8, 1000)
    assert not tmatch(predicted_exons, assembled_exons_3, 0.8, 1000)
    assert     tmatch(predicted_exons, assembled_exons_4, 0.8, 1000)
    assert not tmatch(predicted_exons, assembled_exons_4, 0.95, 1000)
    assert     tmatch(predicted_exons, assembled_exons_5, 0.7, 1000)
    assert     tmatch(predicted_exons, assembled_exons_6, 0.5, 1000)
    assert not tmatch(predicted_exons, assembled_exons_7, 0.95, 1000)
