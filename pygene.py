import gradio as gr

import modules.genbank as genbank
import modules.pdb as pdb
import modules.fasta as fasta

with gr.Blocks(title="pygene") as ui:
    gr.Markdown("""
        # Bioinformatics Search Engine
    """)
    with gr.Tabs() as tabs:
        with gr.TabItem(label="Genbank Search"):
            with gr.Tabs() as gb_s:
                with gr.TabItem(label="Genbank ID Search", id=1):
                    with gr.Row():
                        query = gr.Textbox(label="Search query", placeholder="enter your search query here, ex: Zika or H1N1")
                    with gr.Row():        
                        genbank_btn = gr.Button("Search for IDs")
                    with gr.Row():
                        gr.Markdown("""
                            ## Search result
                        """)
                    with gr.Row():
                        gb_dropdown = gr.Dropdown(label="IDs List", interactive=True)
                    with gr.Row():
                        gb_detail = gr.Textbox(label="ID Description")
                    with gr.Row():        
                        detail_btn = gr.Button("Find selected ID details")
                    genbank_btn.click(fn=genbank.get_all_seq, inputs=[query], outputs=gb_dropdown)
                with gr.TabItem(label="Genbank Detail Search", id=2):
                    with gr.Row():
                        query = gr.Textbox(label="Search query", placeholder="enter Genbank id here, ex: OQ180929.1")
                    with gr.Row():        
                        genbank_btn = gr.Button("Find sequence")
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
                    genbank_btn.click(fn=genbank.getDNA, inputs=query, outputs=[seq_id, seq_desc, seq_len, seq])
                gb_dropdown.select(fn=genbank.getDesc, inputs=gb_dropdown, outputs=gb_detail)
                detail_btn.click(fn=genbank.get_seq_details, inputs=gb_dropdown, outputs=[gb_s, query])
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
                    pdb_btn.click(fn=pdb.get_all_pdb, inputs=query, outputs=[pdb_results, pdb_dropdown])
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
                    pdb_btn.click(fn=pdb.update, inputs=query_d, outputs=[mol, err_out])    
                model_button.click(pdb.change_pdb3d, pdb_dropdown, [pdb_s, query_d])
        with gr.TabItem(label="File Analysis"):
            with gr.Tabs() as file_ly:
                with gr.TabItem(label="FASTA Analysis", id=1):
                    with gr.Row():
                        with gr.Column():
                            with gr.Row():
                                fasta_file = gr.File(file_types=[".fasta"])
                            with gr.Row():
                                fasta_button = gr.Button("Analyze").style(full_width=True)
                        with gr.Column():
                            fasta_gc_percent = gr.Textbox(label="GC Content Percentage")
                            fasta_at_percent = gr.Textbox(label="AT Content Percentage")
                            fasta_mt = gr.Textbox(label="Melting Temperature (Wallace Method)")
                            fasta_weight = gr.Textbox(label="Molecular Weight")
                    with gr.Row():
                        fasta_seq = gr.Textbox(label="Sequence").style(show_copy_button=True)
                    with gr.Row():
                        fasta_mrna = gr.Textbox(label="Transcribed mRNA").style(show_copy_button=True)
                        with gr.Column():
                            fasta_protein_3l = gr.Textbox(label="Translated Protein (3-letter code)").style(show_copy_button=True)
                    fasta_button.click(fn=fasta.get_fasta_data, inputs=fasta_file, outputs=[fasta_seq, fasta_mrna, fasta_protein_3l, fasta_gc_percent, fasta_at_percent, fasta_mt, fasta_weight])
    with gr.Row():
        resources = gr.Markdown("""
#### [Github page](https://github.com/Simplegram/pygene)
### Resources and guides used
- [PyPDB](https://github.com/williamgilpin/pypdb)
- [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/)
- [RSCB](https://www.rcsb.org/)
- [3dmol.js](https://3dmol.org/doc/global.html)
- [Biopython](https://biopython.org/docs)
- [Gradio WebUI](https://gradio.app/docs)
- [Visualize proteins on Hugging Face Spaces - Simon Duerr](https://huggingface.co/blog/spaces_3dmoljs)
        """)

ui.launch(server_name="0.0.0.0", server_port=7870)