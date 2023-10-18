from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import substitution_matrices


def biopython_align(seq1, seq2, matrix="BLOSUM62", gap_open=-10, gap_extend=-0.5):
    matrix = substitution_matrices.load("BLOSUM62")
    alignments = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)
    top_alignment = alignments[0]
    aligned_seq1, aligned_seq2, score, begin, end = top_alignment
    return aligned_seq1, aligned_seq2, score, begin, end