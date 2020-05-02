
import os
import subprocess

from config import MEDIA_DIR, TREE_DIR, FASTA_DIR


def tree_creator(selectname):
    """
    DESCRIPTION:
    A function create the tree inference and store the results in a subfolder within covid_phylo/tree/
    :param selectname: [string] name of the file to be put in the tree folder. The same as the name of the subfolder.
    :return: None
    """
    print('Executing tree inference')
    filename = selectname
    subfolder = selectname.split('.')[0]
    process = subprocess.run(['iqtree', '-s', f'{TREE_DIR / subfolder / filename}', '-bnni', '-nt', 'AUTO'])
    print('Tree inference completed with exit code %d' % process.returncode)


def align_selector(origname, destname, n_genomes):
    """
    DESCRIPTION:
    Function to select the n alignments with the lowest number of gaps.
    :param origname: [string] name of the file with the complete list of alignments in the fasta folder.
    :param destname: [string] name of the file to be put in the tree folder. The same as the name of the subfolder.
    :param n_genomes: [integer] number of alignments to be taken.
    :return: None. It writes the selected aignments in the destname folder.
    """
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

