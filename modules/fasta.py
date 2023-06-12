from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt, molecular_weight as mw, GC, seq3

def get_fasta_data(query):
    protein3l = ""

    fasta = SeqIO.read(query.name, "fasta")

    seq = fasta.seq
    
    mrna = seq.transcribe()
    protein = [seq3(i) for i in seq.translate().split("*")]
    for i in protein:
        protein3l = protein3l + i + "*"

    gc = round(GC(seq), 2)
    at = round(100 - gc, 2)

    melt_temp = round(mt.Tm_Wallace(seq), 2)
    weight = round(mw(seq), 2)

    #fasta_seq, fasta_mrna, fasta_protein_3l, fasta_gc_percent, fasta_at_percent, fasta_mt, fasta_weight
    return seq, mrna, protein3l, gc, at, melt_temp, weight