from Bio import Phylo
import numpy as np
from io import StringIO
import ete3
import tempfile

# RelTime-like implementation from https://github.com/adamhockenberry/dca-weighting/tree/master/Code/supporting_functions.py
def RelTime_AJH(my_tree, modifier=1.):
    my_tree.add_features(rate=modifier)
    if len(my_tree.get_children()) == 0:
        return
    else:
        l_child = my_tree.children[0]
        r_child = my_tree.children[1]
        l_side_dists = [i.dist*len(i.get_leaves()) for i in l_child.traverse()]
        r_side_dists = [i.dist*len(i.get_leaves()) for i in r_child.traverse()]
        

        total_modifier = np.sum(l_side_dists+r_side_dists)/len(my_tree.get_leaves())
        l_modifier = np.sum(l_side_dists)/len(l_child.get_leaves())/total_modifier * modifier
        r_modifier = np.sum(r_side_dists)/len(r_child.get_leaves())/total_modifier * modifier

        my_tree = RelTime_AJH(l_child, l_modifier)
        my_tree = RelTime_AJH(r_child, r_modifier)
    return

def getRelTimeTree(tree_string):
    
    # Read in tree from R as string
    return(RelTime_AJH(Phylo.read(StringIO(tree_string),format='newick',rooted=True)).write())
