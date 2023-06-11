# pygene
A simple Gradio WebUI to find details about a gene and its 3D visualization using Python

Tested using Windows 11. I'm not a Python pro and there will be a lot of bugs in the program. A bug report and a fix is appreaciated!

## Requirements
- Conda (optional, recommended)
- Gradio
- Biopython
- PyPDB
- wget

## Setup Conda Environment
```
conda create --name pygene
```

## Clone Repo
```
git clone https://github.com/Simplegram/pygene.git
```

## Install Requirements
```
pip install -r requirements
```

## Run
Run with
```py
python pygene.py
```

All PDB files are saved to `pdbs` folder. Create one in the root directory.

## Resources and guides used
- [PyPDB](https://github.com/williamgilpin/pypdb)
- [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/)
- [RSCB](https://www.rcsb.org/)
- [3dmol.js](https://3dmol.org/doc/global.html)
- [Biopython](https://biopython.org/docs)
- [Gradio WebUI](https://gradio.app/docs)
- [Visualize proteins on Hugging Face Spaces - Simon Duerr](https://huggingface.co/blog/spaces_3dmoljs)
### Many thanks to Simon Duerr for his implementation of py3dmol and protein visualization with HTML
