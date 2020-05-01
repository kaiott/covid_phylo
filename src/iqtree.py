
import os
from src.config import MEDIA_DIR, TREE_DIR, FASTA_DIR


def tree_creator(selectname):
    print('Executing tree inference')
    filename = selectname
    subfolder = selectname.split('.')[0]
    os.system(f'iqtree -s {TREE_DIR / subfolder / filename}')
    print('Tree inference completed')


def align_selector(origname, destname, n_genomes):
    # Part to take the alignments with lowest number of gaps
    file = open(FASTA_DIR / origname, 'r')

    # Read the information
    data = file.read()
    file.close()
    data = data.split('>')[1::]  # The first element is an empty string

    # Take the best n models
    gaps = list(enumerate([model.count('-') for model in data]))
    data = '\n'.join(['>' + data[element[0]] for element in sorted(gaps, key=lambda x: x[1])[0:n_genomes]])

    # Write the selected data into another file
    sel_dir = TREE_DIR / destname.split('.')[0]
    sel_dir.mkdir(exist_ok=True)
    file = open(sel_dir / destname, 'w')
    file.write(data)
    file.close()

