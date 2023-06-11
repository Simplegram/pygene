import wget
import shutil
import os
import gradio as gr

from pypdb import *

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

def change_pdb3d(pdb_id):
    return gr.Tabs.update(selected=2), gr.Textbox.update(value=pdb_id)