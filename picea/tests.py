from .tree import Tree

def test_empty_init():
    tree = Tree()

def test_parsing():
    newick = '(((a,b),(c,d)),e);'
    tree = Tree.from_newick(newick)

def test_input_output():
    newick = '(((a,b),(c,d)),e);'
    tree = Tree.from_newick(newick)
    print(newick,tree.to_newick(branch_lengths = False))
    assert newick == tree.to_newick(branch_lengths = False)