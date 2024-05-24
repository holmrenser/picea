from unittest import TestCase

from picea import Tree


class TreeTests(TestCase):
    def setUp(self):
        # Very simple tree
        self.newick1 = '(((a,b),(c,d)),e);'
        # Tree with internal node names and branch lengths
        self.newick2 = '((p_patens:1,a_trichopoda:1)N1:1,((z_mays:1,o_sativa:1)N3:1,\
            (a_thaliana:1,(m_truncatula:1,g_max:1)N5:1)N4:1)N2:1)N0:0;'

    def test_empty_init(self):
        Tree()

    def test_parsing_from_string(self):
        Tree.from_newick(string=self.newick1)
        Tree.from_newick(self.newick2)

    def test_newick_writing(self):
        tree1 = Tree.from_newick(string=self.newick1)
        newick1= tree1.to_newick()
        self.assertEqual(newick1, self.newick1)

    def test_input_output(self):
        tree = Tree.from_newick(string=self.newick1)
        self.assertEqual(self.newick1, tree.to_newick(branch_lengths=False))

    def test_root(self):
        tree = Tree.from_newick(string=self.newick1)
        deep_node = tree.loc["a"]
        self.assertEqual(deep_node.root, tree)

    def test_number_of_elements(self):
        tree1 = Tree.from_newick(string=self.newick1)
        self.assertEqual(9, len(tree1.nodes))
        self.assertEqual(5, len(tree1.leaves))
        tree2 = Tree.from_newick(string=self.newick2)
        self.assertEqual(13, len(tree2.nodes))
        self.assertEqual(7, len(tree2.leaves))

    """
    def test_quoted_fasttree_newick(self):
        print(__file__)
        print('hi')
        Tree.from_newick(
            filename='./tests/data/fasttree.quoted_labels.newick'
        )
    """
