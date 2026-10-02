import math
import warnings
from unittest import TestCase

import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt  # noqa: E402

from picea import Tree, calculate_tree_layout, treeplot  # noqa: E402
from picea.tree import TreeStyle  # noqa: E402


class TreeTests(TestCase):
    def setUp(self):
        # Very simple tree
        self.newick1 = "(((a,b),(c,d)),e);"
        # Tree with internal node names and branch lengths
        self.newick2 = "((p_patens:1,a_trichopoda:1)N1:1,((z_mays:1,o_sativa:1)N3:1,\
            (a_thaliana:1,(m_truncatula:1,g_max:1)N5:1)N4:1)N2:1)N0:0;"

    def test_empty_init(self):
        Tree()

    def test_parsing_from_string(self):
        Tree.from_newick(string=self.newick1)
        Tree.from_newick(self.newick2)

    def test_newick_writing(self):
        tree1 = Tree.from_newick(string=self.newick1)
        newick1 = tree1.to_newick()
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


class TreeBugfixTests(TestCase):
    def setUp(self):
        self.newick = "((a:1,b:2)ab:1,c:3)root;"

    def tearDown(self):
        plt.close("all")

    def test_root_without_branch_length_does_not_warn(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            tree = Tree.from_newick(self.newick)
        self.assertEqual(tree.cumulative_length, 0.0)

    def test_zero_length_branch_cumulative_length(self):
        tree = Tree.from_newick("((a:1,(b:1)x:0)y:2,c:1);")
        self.assertEqual(tree.loc["b"].cumulative_length, 3.0)
        self.assertEqual(tree.loc["x"].cumulative_length, 2.0)

    def test_subtree_to_newick(self):
        tree = Tree.from_newick(self.newick)
        self.assertEqual(tree.loc["ab"].to_newick(branch_lengths=True), "(a:1.0,b:2.0)ab;")

    def test_rename_leaves_copy(self):
        tree = Tree.from_newick(self.newick)
        renamed = tree.rename_leaves(str.upper, inplace=False)
        self.assertEqual(sorted(leaf.name for leaf in renamed.leaves), ["A", "B", "C"])
        self.assertEqual(sorted(leaf.name for leaf in tree.leaves), ["a", "b", "c"])
        self.assertIsNone(tree.rename_leaves(str.upper))

    def test_layout_uses_distance_to_root(self):
        tree = Tree.from_newick(self.newick)
        layout = calculate_tree_layout(tree)
        for node in tree.nodes:
            self.assertEqual(layout[node.ID].x, node.cumulative_length)
        rtl = calculate_tree_layout(tree, ltr=False)
        self.assertEqual(rtl[tree.loc["b"].ID].x, -3.0)

    def test_cladogram_layout_aligns_leaves(self):
        for tree, branchlengths in ((Tree.from_newick(self.newick), False), (Tree.from_newick("((a,b),c);"), True)):
            layout = calculate_tree_layout(tree, branchlengths=branchlengths)
            self.assertEqual({layout[leaf.ID].x for leaf in tree.leaves}, {2.0})

    def test_radial_layout(self):
        tree = Tree.from_newick(self.newick)
        layout = calculate_tree_layout(tree, style="radial")
        for leaf in tree.leaves:
            x, y = layout[leaf.ID]
            self.assertAlmostEqual(math.hypot(x, y), leaf.cumulative_length)
        angles = sorted(math.atan2(layout[leaf.ID].y, layout[leaf.ID].x) % (2 * math.pi) for leaf in tree.leaves)
        self.assertAlmostEqual(angles[1] - angles[0], 2 * math.pi / 3)

    def test_invalid_style(self):
        with self.assertRaises(ValueError):
            calculate_tree_layout(Tree.from_newick(self.newick), style="rectangular")

    def test_treeplot_styles(self):
        for newick in (self.newick, "(((a,b),(c,d)),e);"):
            tree = Tree.from_newick(newick)
            for style in ("square", "triangular", "radial", TreeStyle.square):
                for ltr in (True, False):
                    treeplot(tree, style=style, ltr=ltr)

    def test_treeplot_options(self):
        tree = Tree.from_newick(self.newick)
        linestyles = []

        def branch_linestyle(branch):
            linestyles.append(branch)
            return {"color": "red"}

        ax = treeplot(
            tree,
            leaf_labels=False,
            leaf_marker=None,
            leaf_marker_fill=lambda leaf: "red",
            branch_linestyle=branch_linestyle,
        )
        self.assertEqual(len(linestyles), len(tree.links))
        self.assertEqual([text.get_text() for text in ax.texts], ["root", "ab"])
        ax, layout = treeplot(tree, node_labels=False, return_layout=True)
        self.assertEqual(sorted(text.get_text() for text in ax.texts), ["a", "b", "c"])
        self.assertEqual(layout[tree.loc["c"].ID].x, 3.0)
