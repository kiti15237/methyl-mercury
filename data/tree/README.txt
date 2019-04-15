Tree downloaded from Silva living tree site. Due to some weird formatting issues, can only be read into python. 

import dendropy
tree = dendropy.Tree.get(path = "C:/Users/ctata/Documents/Lab/methyl-mercury/data/tree/LTPs132_SSU_tree2.newick.txt", schema="newick")
#dists = tree.phylogenetic_distance_matrix()

tree.write(path = "C:/Users/ctata/Documents/Lab/methyl-mercury/data/tree/LTPs132_SSU_tree.tre", schema = "newick")