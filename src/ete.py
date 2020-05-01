
import os
from config import MEDIA_DIR
from ete3 import Tree, TreeStyle


def newick_extractor(fileroute):
    pass

def tree_viewer(treefolder):
    filename = os.path.basename(os.path.normpath(treefolder)) + '.txt.treefile'
    treeroute = treefolder / filename
    file = open(treeroute, 'r')
    tree = file.read()
    file.close()
    tree = Tree(tree)
    circular_style = TreeStyle()
    circular_style.mode = "c"
    circular_style.scale = 20
    tree.render(str(MEDIA_DIR / 'mytree.png'), w=1024, units='mm', tree_style=circular_style)
