from picea.tree import Tree, treeplot
from typing import \
    Iterable, Callable, List, Optional, Generator, Dict, Union, Tuple, Set
import itertools
import numpy as np
from dataclasses import dataclass, field, InitVar, asdict, replace
from caretta import neighbor_joining as nj


@dataclass
class Clade:
    internode: Tree
    clade: Set[Tuple[str]]
    depth: int

    @classmethod
    def from_internode(cls, internode: Tree, banned: [None, Set[str]] = None):
        assert internode.children is not None
        if banned is None:
            banned = set()
        clade: Set[Tuple[str, ...]] = set()
        for children in internode.children:
            clade.add(tuple(sorted(list({x.name for x in children.leaves if x.name not in banned}))))
        return cls(internode, clade, internode.depth)

    def compare_to(self, other: 'Clade'):
        if self.clade == other.clade:
            return True
        else:
            return False


@dataclass
class Clades:
    tree: Tree
    clades: List[Clade]
    bootstrap_scores: Dict[str, int]

    @classmethod
    def from_tree(cls, tree: Tree, banned: [None, Set[str]] = None):
        tree: Tree = tree.root
        queue: List[Tree] = [tree]
        clades: List[Clade] = list()
        while queue:
            node = queue.pop(0)
            if node.children is not None:
                clades.append(Clade.from_internode(node, banned))
            else:
                continue
            queue += node.children
        return cls(tree.root, clades, {x.name: 0 for x in tree.breadth_first() if x.children is not None})

    def redo_with_banned(self, banned: Set[str]):
        new_clades: Clades = Clades.from_tree(self.tree, banned)
        new_clades.bootstrap_scores = self.bootstrap_scores
        return new_clades

    def compare_to_other(self, other: 'Clades'):
        for this_clade in self.clades:
            for other_clade in other.clades:
                if this_clade.clade == other_clade.clade:
                    self.bootstrap_scores[this_clade.internode.name] += 1
                else:
                    continue


def prepare_bootstrap_trees(distance_matrix: np.ndarray,
                            names: [None, List[str]] = None,
                            iteration: int = 10) -> Tuple[Tree, List[Tree]]:
    if names is None:
        names = [str(x) for x in range(distance_matrix.shape[0])]
    tree: Tree = build_nj_tree_from_distance_matrix(distance_matrix, names)
    names = np.array(names)
    other_trees: List[Tree] = list()
    for i in range(iteration):
        selected_ids = np.unique(np.random.choice(np.arange(names.shape[0]), size=names.shape[0], replace=True))
        other_trees.append(build_nj_tree_from_distance_matrix(distance_matrix[selected_ids, :][:, selected_ids],
                                                              list(names[selected_ids])).root)
    return tree.root, other_trees


def bootstrap(tree: Tree, bootstrap_trees: List[Tree]) -> Clades:
    tree_clades: Clades = Clades.from_tree(tree)
    leaf_names: Set[str] = {x.name for x in tree.leaves}
    other_tree_clades: Iterable[Clades] = (Clades.from_tree(x) for x in bootstrap_trees)
    for other_clades in other_tree_clades:
        banned_names = leaf_names - {x.name for x in other_clades.tree.leaves}
        tree_clades.redo_with_banned(banned_names).compare_to_other(other_clades)
    return tree_clades


def build_nj_tree_from_distance_matrix(distance_matrix: np.ndarray, names: List[str]) -> Tree:
    list_of_edges: List[Tuple[int, int]] = [x[::-1] for x in nj.neighbor_joining(distance_matrix)[0]]
    names = names[::-1]
    return Tree.from_edge_list(list_of_edges, {i: names[i] for i in range(distance_matrix.shape[0])})


def main():
    iteration=100
    from sklearn.metrics import pairwise_distances
    import matplotlib.pyplot as plt
    some_arrays = np.random.random_sample(size=(10, 5))
    print(build_nj_tree_from_distance_matrix(pairwise_distances(some_arrays, metric="euclidean"),
                                             [str(x) for x in range(some_arrays.shape[0])]))
    tree, other_trees = prepare_bootstrap_trees(pairwise_distances(some_arrays), iteration=iteration)

    tree_clade = bootstrap(tree, other_trees)
    for clade in tree_clade.clades:
        clade.internode.name = str(int(tree_clade.bootstrap_scores[clade.internode.name]/iteration*100))
    print(tree_clade.bootstrap_scores)
    treeplot(tree_clade.tree, internode_names=True)
    plt.show()


if __name__ == "__main__":
    main()
