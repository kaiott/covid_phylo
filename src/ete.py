import os
from config import MEDIA_DIR
from ete3 import Tree, TreeStyle


def tree_viewer(treefolder):
    """
    DESCRIPTION:
    A function to create the representation of the tree with ete. Given the folder with the files in the tree folder it
    creates the tree representation in the media folder.
    :param treefolder: [pathlib] route to the subfolder in the tree folder which contains the data of the inference.
    :return: None. It writes the image in the media folder.
    """
    filename = os.path.basename(os.path.normpath(treefolder)) + '.txt.treefile'
    imagefile = os.path.basename(os.path.normpath(treefolder)) + '.png'
    treeroute = treefolder / filename
    file = open(treeroute, 'r')
    tree = file.read()
    file.close()
    tree = Tree(tree)
    circular_style = TreeStyle()
    circular_style.mode = "c"
    circular_style.scale = 20
    tree.render(str(MEDIA_DIR / imagefile), w=1024, units='mm', tree_style=circular_style)
