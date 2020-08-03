import re
from typing import \
    Iterable, Callable, List, Optional, Generator, Dict, Union, Tuple, Set
import json
import itertools
import numpy as np

from dataclasses import dataclass, field, InitVar, asdict, replace
from collections import defaultdict
from matplotlib import pyplot as plt


TreeDict = Dict[str, Union[str, int, float, List[Optional['TreeDict']]]]


def index_equals(
    node: 'Tree',
    index: int
) -> bool:
    """[summary]

    Args:
        node (Tree): [description]
        index (int): [description]

    Returns:
        bool: [description]
    """
    return node.ID == index


def name_equals(
    node: 'Tree',
    name: str
) -> bool:
    """[summary]

    Args:
        node (Tree): [description]
        name (str): [description]

    Returns:
        bool: [description]
    """
    return node.name == name


@dataclass
class Tree:
    ID: int = 0
    depth: int = 0
    name: str = ""
    length: float = 0.
    cumulative_length: float = 0.
    children: [None, List['Tree']] = None
    parent: InitVar[[None, 'Tree']] = None

    def __post_init__(self, parent: [None, 'Tree']):
        self.parent = parent
        if not len(self.name):
            self.name: str = str(self.ID)

    @property
    def loc(self) -> 'Tree':
        """Name based index

        Example:
            >>> from picea import Tree
            >>> newick = '(((a,b),(c,d)),e);'
            >>> tree = Tree.from_newick(newick)
            >>> tree.loc['a']
            Tree(name='a', length=0.0, children=[])

        Returns:
            Tree: tree node matching name

        Raises:
            IndexError
        """
        return TreeIndex(iterator=self.depth_first, eq_func=name_equals)

    @property
    def iloc(self) -> 'Tree':
        """Index based index

        Example:
            >>> from picea import Tree
            >>> newick = '(((a,b),(c,d)),e);'
            >>> tree = Tree.from_newick(newick)
            >>> tree.iloc[2]
            Tree(name='', length=0.0, children=[Tree(name='a', length=0.0, \
children=[]), Tree(name='b', length=0.0, children=[])])

        Returns:
            Tree: tree node matching index
        """
        return TreeIndex(iterator=self.depth_first, eq_func=index_equals)

    @property
    def root(self) -> 'Tree':
        """Root node of the (sub)tree

        Returns:
            Tree: Root node
        """
        root = self
        while root.parent:
            root = root.parent
        return root

    @property
    def nodes(self) -> List['Tree']:
        """A list of all tree nodes in breadth-first order

        Returns:
            list: A list of all tree nodes
        """
        return list(self.breadth_first())

    @property
    def leaves(self) -> Iterable['Tree']:
        """A list of leaf nodes only
dict_factory function
        Returns:
            list: A list of leaf nodes only
        """
        return (n for n in self.breadth_first() if not n.children)

    @property
    def links(self) -> List[Tuple['Tree', 'Tree']]:
        """A list of all (parent, child) combinations

        Returns:
            list: All (parent,child) combinations
        """
        _links = []
        for node in self.nodes:
            if node.children:
                for child in node.children:
                    _links.append((node, child))
        return _links

    def find_leaf_by_name(self, name: str) -> [None, 'Tree']:
        for leaf in self.leaves:
            if leaf.name == name:
                return leaf
        return None

    @classmethod
    def from_newick(
        cls,
        string: Optional[str] = None,
        filename: Optional[str] = None
    ) -> 'Tree':
        """Parse a newick formatted string into a Tree object

        Arguments:
            newick_string (string): Newick formatted tree string

        Returns:
            Tree: Tree object
        """
        assert filename or string
        assert not (filename and string)
        if filename:
            with open(filename) as filehandle:
                string = filehandle.read()
        tokens = re.split(r'\s*(;|\(|\)|,|:)\s*', string)
        ID = 0
        tree = cls(ID=ID, depth=0, length=0.0, cumulative_length=0.0)
        ancestors = list()
        for i, token in enumerate(tokens):
            if token == '(':
                ID += 1
                subtree = cls(ID=ID, depth=0)
                tree.children = [subtree]
                ancestors.append(tree)
                tree = subtree
            elif token == ',':
                ID += 1
                subtree = cls(ID=ID, depth=0)
                ancestors[-1].add_as_child(subtree)
                tree = subtree
            elif token == ')':
                tree = ancestors.pop()
            else:
                previous_token = tokens[i - 1]
                if previous_token in ('(', ')', ','):
                    tree.name = token
                elif previous_token == ':':
                    tree.length = float(token)
        tree.depth = 0
        queue = [tree]
        while queue:
            node = queue.pop(0)
            for child in node.children:
                child.parent = node
                child.depth = node.depth + 1
                child.cumulative_length = node.cumulative_length \
                    + abs(child.length)
            queue += node.children

        return tree

    def to_newick(
        self,
        branch_lengths: bool = True
    ) -> str:
        """Make a Newick formatted string

        Args:
            branch_lengths (bool, optional): Whether to include branch lengths\
             in the Newick string. Defaults to True.

        Returns:
            String: Newick formatted tree string
        """
        if self.name:
            name = str(self.name)
        else:
            name = ''

        if self.children:
            subtree_string = ','.join([
                c.to_newick(branch_lengths=branch_lengths)
                for c in self.children
            ])
            newick = f'({subtree_string}){name}'
        else:
            newick = name

        if branch_lengths and self.ID != 0:
            length = self.length
            if length == 0:
                length = int(0)
            newick += f':{length}'

        if self == self.root:
            newick += ';'

        return newick

    @classmethod
    def from_sklearn(
        cls,
        clustering
    ) -> 'Tree':
        """Read a tree from sklearn agglomerative clustering

        Args:
            clustering (sklearn object): sklearn agglomerative clustering\
                 object.

        Returns:
            Tree: Tree object
        """
        nodes = clustering.children_
        n_leaves = nodes.shape[0] + 1
        tree = cls(ID=nodes.shape[0] * 2, depth=0)

        queue = [tree]
        while queue:
            node = queue.pop(0)
            if node.ID < n_leaves:
                node.name = str(node.ID)
                continue
            for child_ID in nodes[node.ID - n_leaves]:
                child = cls(ID=child_ID, depth=node.depth + 1)
                child.parent = node
                node.add_as_child(child)
            queue += node.children

        return tree

    def to_sklearn(self):
        # TODO
        raise NotImplementedError()

    @classmethod
    def from_edge_dict(cls, tree_dict: Dict[int, List[int]], names: [None, Dict[int, str]] = None) -> 'Tree':
        if names is None:
            names = dict()
        # Find the root node
        nodes_with_parents: Set[int] = set(itertools.chain.from_iterable(tree_dict.values()))
        parent_nodes: Set[int] = set(tree_dict.keys())
        nodes_without_parents: Set[int] = parent_nodes - nodes_with_parents
        assert len(nodes_without_parents) == 1, "Not a tree structure, multiple roots found."
        root_id: int = list(nodes_without_parents)[0]

        # start building the tree
        tree: 'Tree' = cls(ID=root_id, depth=0, name="root")
        tree.depth = 0
        queue: List['Tree'] = [tree]
        while queue:
            node = queue.pop(0)
            if node.ID in tree_dict:
                for child_ID in tree_dict[node.ID]:
                    if child_ID in names:
                        child = cls(ID=child_ID, name=names[child_ID], depth=node.depth + 1)
                    else:
                        child = cls(ID=child_ID, depth=node.depth + 1)
                    child.parent = node
                    node.add_as_child(child)
                queue += node.children
            else:
                continue

        return tree

    def add_as_child(self, node: 'Tree'):
        if not self.children:
            self.children = [node]
        else:
            self.children.append(node)

    @classmethod
    def from_edge_list(cls, edge_list: List[Tuple[int, int]], names: [None, Dict[int, str]] = None) -> 'Tree':
        edge_dict: Dict[int, List[int]] = dict()
        for n1, n2 in edge_list:
            if n1 not in edge_dict:
                edge_dict[n1] = [n2]
            else:
                edge_dict[n1].append(n2)
        return cls.from_edge_dict(edge_dict, names=names)

    @classmethod
    def from_json(cls):
        # TODO
        raise NotImplementedError()

    def to_json(self, indent: Optional[int] = None) -> str:
        return json.dumps(self.to_dict(), indent=indent)

    @classmethod
    def from_dict(cls, tree_dict):
        # TODO
        # raise NotImplementedError()
        tree = cls()
        return tree

    def to_dict(self) -> TreeDict:
        """[summary]

        Returns:
            TreeDict: [description]
        """
        return asdict(self)

    def breadth_first(self) -> Generator['Tree', None, None]:
        """Generator implementing breadth first search starting at root node
        """
        queue = [self]
        while queue:
            node = queue.pop(0)
            if node.children is not None:
                queue += node.children
            yield node

    def depth_first(
        self,
        post_order: bool = True
    ) -> Generator['Tree', None, None]:
        """Generator implementing depth first search in either post- or
        pre-order traversel

        Keyword Arguments:
            post_order (bool, optional): Depth first search in post-order
            traversal or not. Defaults to True
        """
        if not post_order:
            yield self
        if self.children is not None:
            for child in self.children:
                yield from child.depth_first(post_order=post_order)
        if post_order:
            yield self


class TreeIndex(object):
    def __init__(
        self,
        iterator: Iterable[Tree],
        eq_func: Callable[[int, str], bool]
    ):
        """[summary]

        Args:
            object ([type]): [description]
            iterator ([type]): [description]
            eq_func ([type]): [description]
        """
        self.iterator = iterator
        self.eq_func = eq_func

    def __getitem__(self, key):
        for element in self.iterator():
            if self.eq_func(element, key):
                return element
        raise IndexError(f'{key} is not valid index')


def unequal_separation(
    node_a: 'Tree',
    node_b: 'Tree',
    sep_1: float = 1.0,
    sep_2: float = 2.0
) -> float:
    """[summary]

    Args:
        node_a (Tree): [descripti        print(node_a.parent)on]
        node_b (Tree): [description]
        sep_1 (float, optional): [description]. Defaults to 1.0.
        sep_2 (float, optional): [description]. Defaults to 2.0.

    Returns:
        float: [description]
    """
    if node_a.parent == node_b.parent:
        return sep_1
    return sep_2


def equal_separation(
    node_a: 'Tree',
    node_b: 'Tree',
    separation: float = 1.0
) -> float:
    """[summary]

    Args:
        node_a (Tree): [description]
        node_b (Tree): [description]
        separation (float, optional): [description]. Defaults to 1.0.

    Returns:
        float: [description]
    """
    return separation


@dataclass
class TwoDCoordinate():
    x: float = 0.0
    y: float = 0.0


def treeplot(
    tree: Tree,
    separation: Callable = equal_separation,
    d_x: float = 1.0,
    d_y: float = 1.0,
    ltr: bool = True,
    ax: Optional[Callable] = None,
    internode_names: bool = False
):
    coords = defaultdict(TwoDCoordinate)
    previous_node = None
    y = 0
    for node in tree.depth_first(post_order=True):
        if node.children:
            coords[node.ID].y = sum([coords[c.ID].y for c in node.children]) \
                / len(node.children)
            if ltr:
                coords[node.ID].x = 1 + max([
                    coords[c.ID].x for c in node.children
                ])
            else:
                coords[node.ID].x = min([
                    coords[c.ID].x for c in node.children
                ]) - 1
        else:
            if previous_node:
                y += separation(node, previous_node)
                coords[node.ID].y = y
            else:
                coords[node.ID].y = 0
            coords[node.ID].x = 0
            previous_node = node

    for node in tree.depth_first(post_order=True):
        coords[node.ID].x = (coords[tree.ID].x - coords[node.ID].x) * d_x
        coords[node.ID].y = (coords[node.ID].y - coords[tree.ID].y) * d_y

    if not ax:
        fig, ax = plt.subplots(figsize=(10, 10))

    for node1, node2 in tree.links:
        ax.plot(
            (coords[node1.ID].x, coords[node2.ID].x),
            (coords[node1.ID].y, coords[node2.ID].y),
            c='k'
        )
        if internode_names:
            ax.text(coords[node1.ID].x,
                    coords[node1.ID].y,
                    node1.name,
                    fontsize=15,
                    # in_layout=True,
                    clip_on=True,
                    color="green")
    for leaf in tree.leaves:
        ax.text(
            coords[leaf.ID].x + .1,
            coords[leaf.ID].y - .1,
            leaf.name,
            fontsize=12,
            # in_layout=True,
            clip_on=True
        )
    ax.set_xticks(())
    ax.set_yticks(())
    fig.tight_layout()
    return ax


