import gradio as gr
import os
import wget
import shutil

from Bio import Entrez
from Bio import SeqIO
from pypdb import *

resources = [
    "https://github.com/williamgilpin/pypdb", 
    "https://www.ncbi.nlm.nih.gov/nuccore/",
    "https://www.rcsb.org/",
    "https://huggingface.co/blog/spaces_3dmoljs",
    "https://3dmol.org/doc/global.html#createViewer",
    "https://gradio.app/docs",
]

resources_list = ""
for rs in resources:
    resources_list = resources_list + f"{rs}\n"

## PDB DATABASE
def read_mol(molpath):
    with open(molpath, "r") as fp:
        lines = fp.readlines()
    mol = ""
    for l in lines:
        mol += l

    return mol

def molecule(input_pdb):
    mol = read_mol(input_pdb)
    mol = mol.replace("'", "")

    x = (
        """<!DOCTYPE html>
        <html>
        <head>    
    <meta http-equiv="content-type" content="text/html; charset=UTF-8" />
    <style>
    body{
        font-family:sans-serif
    }
    .mol-container {
    width: 100%;
    height: 600px;
    position: relative;
    }
    .mol-container select{
        background-image:None;
    }
    </style>
     <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.3/jquery.min.js" integrity="sha512-STof4xm1wgkfm7heWqFJVn58Hm3EtS31XFaagaa8VMReCXAkQnJZ+jEy8PCC/iT18dFy95WcExNHFTqLyp72eQ==" crossorigin="anonymous" referrerpolicy="no-referrer"></script>
    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    </head>
    <body>  
    <div id="container" class="mol-container"></div>
  
            <script>
               let pdb = `"""
        + mol
        + """`  
      
             $(document).ready(function () {
                let element = $("#container");
                let config = { backgroundColor: "white" };
                let viewer = $3Dmol.createViewer(element, config);
                viewer.addModel(pdb, "pdb");
                viewer.setStyle({ cartoon: { color:"spectrum" } });
                viewer.zoomTo();
                viewer.render();
                viewer.zoom(0.8, 2000);
              })
        </script>
        </body></html>"""
    )

    return f"""<iframe style="width: 100%; height: 600px" name="result" sandbox="allow-modals allow-forms 
    allow-scripts allow-same-origin allow-popups 
    allow-top-navigation-by-user-activation allow-downloads" allowfullscreen="" 
    allowpaymentrequest="" frameborder="0" srcdoc='{x}'></iframe>"""

def get_pdb(pdb_code):
    path = f"./pdbs/{pdb_code}.pdb"
    if os.path.isfile(path) is False:
        wget.download(f"https://files.rcsb.org/download/{pdb_code}.pdb")
        shutil.move(f"{pdb_code}.pdb", f"./pdbs/")
    return f"./pdbs/{pdb_code}.pdb"

def update(query):
    err = "Success!"
    try:
        pdb_path = get_pdb(query)
        return molecule(pdb_path), err
    except:
        err = "File not found or server error!"
        return None, err

def get_all_pdb(query):
    choice = ['Not Found!']
    try:
        found_pdbs = Query(query).search()
        text = ""
        for pdb in found_pdbs[:10]:
            info = get_info(pdb)
            title = info['citation'][0]['title']

            text = text + f"Title: {title}\nPDB ID: {pdb}\n\n"
        return text, gr.Dropdown.update(choices=found_pdbs[:10])
    except:
        text = "Not Found"
    return text, gr.Dropdown.update(choices=choice)

## GENBANK DATABASE
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

    handle = Entrez.efetch(db="nucleotide", id=sequence_ids, rettype="fasta", retmode="text")
    sequences = SeqIO.parse(handle, "fasta")

    text = ""
    try:
        for seq_record in sequences:
            text = text + f'Sequence ID: {seq_record.id}\nSequence Description: {seq_record.description.split(seq_record.id)[1]}\n\n'
    except:
        text = "Not Found"

    handle.close()

    return text

def change_tab(value):
    return gr.Tabs.update(selected=2), gr.Textbox.update(value=value)

with gr.Blocks() as ui:
    gr.Markdown("""
        # Bioinformatics Search Engine
        ### Dibuat oleh Kelompok 12
    """)
    with gr.Tabs() as tabs:
        with gr.TabItem(label="Genbank Search"):
            with gr.Tab(label="Genbank ID Search"):
                with gr.Row():
                    query = gr.Textbox(label="Search query", placeholder="enter your search query here, ex: Zika or H1N1")
                with gr.Row():
                    gr.Markdown("""
                        ## Search result
                    """)
                with gr.Row():
                    id_results = gr.Textbox(label="ID and description results")
                with gr.Row():        
                    genbank_btn = gr.Button("Search for IDs")
                    genbank_btn.click(fn=get_all_seq, inputs=[query], outputs=id_results)   
            with gr.Tab(label="Genbank Detail Search"):
                with gr.Row():
                    query = gr.Textbox(label="Search query", placeholder="enter Genbank id here, ex: OQ180929.1")
                with gr.Row():
                    gr.Markdown("""
                        ## Search result
                    """)
                with gr.Row():
                    seq_id = gr.Textbox(label="Sequence ID")
                    seq_len = gr.Textbox(label="Sequence length")
                with gr.Row():
                    seq_desc = gr.Textbox(label="Sequence description")
                with gr.Row():
                    seq = gr.Textbox(label="Sequence").style(show_copy_button=True)
                # with gr.Row():
                #     openedai = gr.Textbox(label="AI analysis")    
                with gr.Row():        
                    genbank_btn = gr.Button("Find sequence")
                    # find_button.click(fn=getDNA, inputs=[query], outputs=[seq_id, seq_desc, seq_len, seq, openedai])
                    genbank_btn.click(fn=getDNA, inputs=query, outputs=[seq_id, seq_desc, seq_len, seq])
        with gr.TabItem(label="PDB Search"):
            with gr.Tabs() as pdb_s:
                with gr.TabItem(label="PDB ID Search", id=1):
                    with gr.Row():
                        with gr.Column():
                            query = gr.Textbox(label="Search query", placeholder="enter your search query here, ex: Zika or H1N1")
                    with gr.Row():
                        pdb_btn = gr.Button("Search for IDs")
                    with gr.Row():
                        gr.Markdown("""
                            ## Search result
                        """)
                    with gr.Row():
                        pdb_results = gr.Textbox(label="PDB results")
                    with gr.Row():
                        pdb_dropdown = gr.Dropdown(label="PDB results selection", interactive=True)
                    with gr.Row():
                        model_button = gr.Button("View 3D model")
                    pdb_btn.click(fn=get_all_pdb, inputs=query, outputs=[pdb_results, pdb_dropdown])
                with gr.TabItem(label="PDB 3D View", id=2):
                    with gr.Row():
                        with gr.Column():
                            query_d = gr.Textbox(label="Search query", placeholder="enter your PDB id here, ex: 5PTY or 7R98")
                    with gr.Row():
                        gr.Markdown("""
                            ## Search result
                        """)
                    with gr.Row():
                        mol = gr.HTML()
                    with gr.Row():
                        err_out = gr.Textbox(show_label=False).style(container=False)
                    with gr.Row():
                        pdb_btn = gr.Button("View structure")
                        pdb_btn.click(fn=update, inputs=query_d, outputs=[mol, err_out])    
                model_button.click(change_tab, pdb_dropdown, [pdb_s, query_d])
    with gr.Row():
        # resources = gr.TextArea(value=resources_list)
        resources = gr.Markdown("""
### Resources and guides used
- [PyPDB](https://github.com/williamgilpin/pypdb)
- [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/)
- [RSCB](https://www.rcsb.org/)
- [Visualize proteins on Hugging Face Spaces](https://huggingface.co/blog/spaces_3dmoljs)
- [3dmol.js](https://3dmol.org/doc/global.html)
- [Biopython](https://biopython.org/docs)
- [Gradio WebUI](https://gradio.app/docs)
        """)

ui.launch(server_name="0.0.0.0", server_port=7870)