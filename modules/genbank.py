import gradio as gr

from Bio import Entrez, SeqIO

def getDesc(query):
    Entrez.email = "your_email@example.com"

    handle = Entrez.esearch(db="nucleotide", term=query, retmax=1)
    record = Entrez.read(handle)
    handle.close()

    sequence_ids = record["IdList"]

    handle = Entrez.efetch(db="nucleotide", id=sequence_ids, rettype="fasta", retmode="text")
    sequences = SeqIO.parse(handle, "fasta")

    try:
        for seq_record in sequences:
            seq_desc = seq_record.description.split(seq_record.id)[1]
        handle.close()
        return seq_desc
    except:
        seq_desc = "Not Found"
        handle.close()
        return seq_desc

def getDNA(query):
    Entrez.email = "your_email@example.com"

    handle = Entrez.esearch(db="nucleotide", term=query, retmax=1)
    record = Entrez.read(handle)
    handle.close()

    sequence_ids = record["IdList"]

    handle = Entrez.efetch(db="nucleotide", id=sequence_ids, rettype="fasta", retmode="text")
    sequences = SeqIO.parse(handle, "fasta")

    try:
        for seq_record in sequences:
            seq_id = seq_record.id
            seq_desc = seq_record.description.split(seq_record.id)[1]
            seq_len = len(seq_record.seq)
            seq = seq_record.seq
        handle.close()
        return seq_id, seq_desc, seq_len, seq
    except:
        seq_id = "Not Found"
        seq_desc = "Not Found"
        seq_len = "Not Found"
        seq = "Not Found"
        handle.close()
        return seq_id, seq_desc, seq_len, seq

def get_all_seq(query):
    Entrez.email = "your_email@example.com"

    handle = Entrez.esearch(db="nucleotide", term=query, retmax=10)
    record = Entrez.read(handle)
    handle.close()

    sequence_ids = record["IdList"]
    seq_list = []

    handle = Entrez.efetch(db="nucleotide", id=sequence_ids, rettype="fasta", retmode="text")
    sequences = SeqIO.parse(handle, "fasta")

    try:
        for seq_record in sequences:
            seq_list.append(seq_record.id)
        handle.close()
    except:
        seq_list = ["Not Found!"]
        handle.close()

    return gr.Dropdown.update(choices=seq_list)

def get_seq_details(seq_id):
    return gr.Tabs.update(selected=2), gr.Textbox.update(value=seq_id)