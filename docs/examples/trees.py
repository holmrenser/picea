# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Trees
#
# A {class}`~picea.Tree` is a recursive node: every node has a name, a branch length and a list of children. Trees
# are read from and written to Newick, and can be plotted with {func}`~picea.treeplot`.

# %%
from matplotlib import pyplot as plt

from picea import Tree, treeplot

# %% [markdown]
# ## Parsing and navigating

# %%
tree = Tree.from_newick("((a:1,b:2)ab:1,(c:1,d:1)cd:2)root:0;")
tree.to_newick(branch_lengths=True)

# %%
len(tree.nodes), [leaf.name for leaf in tree.leaves]

# %% [markdown]
# Nodes can be looked up by name with `loc`, and every node knows its parent and root.

# %%
a = tree.loc["a"]
a.parent.name, a.root.name, a.cumulative_length

# %% [markdown]
# Traverse the tree depth first (pre- or post-order) or breadth first.

# %%
[node.name for node in tree.depth_first(post_order=False)]

# %%
[node.name for node in tree.breadth_first()]

# %% [markdown]
# ## Plotting
#
# This gene tree of hydroxycinnamoyl transferase (HCT) homologs has branch lengths and support values as internal
# node names.

# %%
hct = Tree.from_newick(filename="data/tree.newick")
len(hct.leaves)

# %%
hct.rename_leaves(lambda name: name.removesuffix(".1"))

fig, ax = plt.subplots(figsize=(8, 9))
treeplot(hct, style="square", ax=ax);

# %% [markdown]
# The radial style puts the root in the center.

# %%
fig, ax = plt.subplots(figsize=(9, 9))
treeplot(hct, style="radial", node_labels=False, ax=ax);

# %% [markdown]
# Leaf markers can be styled per leaf with a function, for example to color leaves by species. Without
# branch lengths the tree is drawn as a cladogram.

# %%
from matplotlib.lines import Line2D

species = {
    "AT": ("Arabidopsis thaliana", "tab:blue"),
    "Potri": ("Populus trichocarpa", "tab:orange"),
    "Eucgr": ("Eucalyptus grandis", "tab:green"),
    "Glyma": ("Glycine max", "tab:red"),
    "Medtr": ("Medicago truncatula", "tab:purple"),
    "Fvesca": ("Fragaria vesca", "tab:brown"),
    "PanWU": ("Parasponia andersonii", "tab:pink"),
    "TorRG": ("Trema orientalis", "tab:olive"),
}


def leaf_color(leaf: Tree) -> str:
    return next(color for prefix, (_, color) in species.items() if leaf.name.startswith(prefix))


fig, ax = plt.subplots(figsize=(8, 9))
treeplot(hct, style="square", branchlengths=False, node_labels=False, leaf_marker_fill=leaf_color, ax=ax)
ax.legend(
    handles=[Line2D([], [], marker="o", linestyle="", color=color, label=name) for name, color in species.values()],
    loc="upper left",
    bbox_to_anchor=(1.01, 1),
);

# %% [markdown]
# ## From hierarchical clustering
#
# Trees can also be created from a fitted scikit-learn
# [`AgglomerativeClustering`](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.AgglomerativeClustering.html)
# model.

# %%
import numpy as np
from sklearn.cluster import AgglomerativeClustering

X = np.array([[1, 2], [1, 4], [1, 0], [4, 2], [4, 4], [4, 0]])
clustering = AgglomerativeClustering().fit(X)
cluster_tree = Tree.from_sklearn(clustering)
cluster_tree.to_newick()

# %%
fig, ax = plt.subplots(figsize=(4, 3))
treeplot(cluster_tree, style="square", branchlengths=False, ax=ax);
