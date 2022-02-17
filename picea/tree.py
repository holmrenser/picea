import re
from typing import (
    Iterable,
    Callable,
    List,
    Optional,
    Generator,
    Dict,
    Union,
    Tuple,
    Type,
    DefaultDict,
)
from enum import Enum
import json
from dataclasses import dataclass, field, InitVar, asdict
from collections import defaultdict
from matplotlib import pyplot as plt
from matplotlib.axes import SubplotBase
import numpy as np
from copy import deepcopy

TreeDict = Dict[str, Union[str, int, float, List[Optional["TreeDict"]]]]


@dataclass
class Tree:
    """Recursive Tree object"""

    name: Optional[str] = None
    length: float = 0.0
    children: Optional[List["Tree"]] = field(default_factory=list)

    ID: InitVar[Optional[int]] = None
    depth: InitVar[Optional[int]] = None
    parent: InitVar[Optional["Tree"]] = None
    cumulative_length: InitVar[float] = 0.0

    def __post_init__(self, ID, *args, **kwargs):
        """[summary]

        Args:
            ID ([type]): [description]
        """
        self.ID = ID

    @property
    def loc(self) -> "Tree":
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
        return TreeIndex(
            iterator=self.depth_first, eq_func=lambda node, name: node.name == name
        )

    @property
    def iloc(self) -> "Tree":
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
        return TreeIndex(
            iterator=self.depth_first, eq_func=lambda node, index: node.ID == index
        )

    @property
    def root(self) -> "Tree":
        """Root node of the (sub)tree

        Returns:
            Tree: Root node
        """
        root = self
        while root.parent:
            root = root.parent
        return root

    @property
    def nodes(self) -> List["Tree"]:
        """A list of all tree nodes in breadth-first order

        Returns:
            list: A list of all tree nodes
        """
        return list(self.breadth_first())

    @property
    def leaves(self) -> List["Tree"]:
        """A list of leaf nodes only

        Returns:
            list: A list of leaf nodes only
        """
        return [n for n in self.nodes if not n.children]

    @property
    def links(self) -> List[Tuple["Tree", "Tree"]]:
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

    @classmethod
    def from_newick(
        cls, string: Optional[str] = None, filename: Optional[str] = None
    ) -> "Tree":
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
        tokens = re.split(r"\s*(;|\(|\)|,|:)\s*", string)
        ID = 0
        tree = cls(ID=ID, length=0.0, cumulative_length=0.0)
        ancestors = list()
        for i, token in enumerate(tokens):
            if token == "(":
                ID += 1
                subtree = cls(ID=ID)
                tree.children = [subtree]
                ancestors.append(tree)
                tree = subtree
            elif token == ",":
                ID += 1
                subtree = cls(ID=ID)
                ancestors[-1].children.append(subtree)
                tree = subtree
            elif token == ")":
                tree = ancestors.pop()
            else:
                previous_token = tokens[i - 1]
                if previous_token in ("(", ")", ","):
                    tree.name = token
                elif previous_token == ":":
                    tree.length = float(token)
        tree.depth = 0
        queue = [tree]
        while queue:
            node = queue.pop(0)
            for child in node.children:
                child.parent = node
                child.depth = node.depth + 1
                child.cumulative_length = node.cumulative_length + abs(child.length)
            queue += node.children

        return tree

    def to_newick(self, branch_lengths: bool = True) -> str:
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
            name = ""

        if self.children:
            subtree_string = ",".join(
                [c.to_newick(branch_lengths=branch_lengths) for c in self.children]
            )
            newick = f"({subtree_string}){name}"
        else:
            newick = name

        if branch_lengths and self.ID != 0:
            length = self.length
            if length == 0:
                length = int(0)
            newick += f":{length}"

        if self == self.root:
            newick += ";"

        return newick

    @classmethod
    def from_sklearn(cls, clustering) -> "Tree":
        """Read a tree from sklearn agglomerative clustering

        Args:
            clustering (sklearn object): sklearn agglomerative clustering\
                 object.

        Returns:
            Tree: Tree object
        """
        nodes = clustering.children_
        n_leaves = nodes.shape[0] + 1
        tree = cls(ID=nodes.shape[0] * 2)

        queue = [tree]
        while queue:
            node = queue.pop(0)
            if node.ID < n_leaves:
                node.name = str(node.ID)
                continue
            for child_ID in nodes[node.ID - n_leaves]:
                child = cls(ID=child_ID)
                child.parent = node
                node.children.append(child)
            queue += node.children

        return tree

    def to_sklearn(self):
        # TODO
        raise NotImplementedError()

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

    def breadth_first(self) -> Generator["Tree", None, None]:
        """Generator implementing breadth first search starting at root node"""
        queue = [self]
        while queue:
            node = queue.pop(0)
            queue += node.children
            yield node

    def depth_first(self, post_order: bool = True) -> Generator["Tree", None, None]:
        """Generator implementing depth first search in either post- or
        pre-order traversel

        Keyword Arguments:
            post_order (bool, optional): Depth first search in post-order
            traversal or not. Defaults to True
        """
        if not post_order:
            yield self
        for child in self.children:
            yield from child.depth_first(post_order=post_order)
        if post_order:
            yield self

    def rename_leaves(
        self, rename_func: Callable, inplace: bool = True
    ) -> Optional["Tree"]:
        """[summary]"""
        tree = self if inplace else deepcopy(self)
        for leaf in tree.leaves:
            leaf.name = rename_func(leaf.name)


class TreeIndex(object):
    def __init__(self, iterator: Iterable[Tree], eq_func: Callable[[int, str], bool]):
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
        raise IndexError(f"{key} is not valid index")


def unequal_separation(
    node_a: "Tree", node_b: "Tree", sep_1: float = 1.0, sep_2: float = 2.0
) -> float:
    """[summary]

    Args:
        node_a (Tree): [description]
        node_b (Tree): [description]
        sep_1 (float, optional): [description]. Defaults to 1.0.
        sep_2 (float, optional): [description]. Defaults to 2.0.

    Returns:
        float: [description]
    """
    if node_a.parent == node_b.parent:
        return sep_1
    return sep_2


def equal_separation(node_a: "Tree", node_b: "Tree", separation: float = 1.0) -> float:
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
class TwoDCoordinate:
    x: float = 0.0
    y: float = 0.0

    def __iter__(self):
        yield from (self.x, self.y)

    def to_polar(self):
        return TwoDCoordinate(x=self.x * np.cos(self.y), y=self.x * np.sin(self.y))

    def to_cartesian(self):
        return TwoDCoordinate(
            x=np.sqrt(self.x**2 + self.y**2), y=np.arctan2(self.y, self.x)
        )


Ax = Type[SubplotBase]
TreeStyle = Enum("TreeStyle", ("square", "radial", "triangular"))
LayoutDict = DefaultDict[int, TwoDCoordinate]


def calculate_tree_layout(
    tree: Tree,
    style: TreeStyle = TreeStyle.square,
    ltr: bool = True,
    branchlengths: bool = True,
) -> LayoutDict:
    """[summary]

    Args:
        tree (Tree): [description]
        style (Treestyle, optional): [description]. Defaults to \
            Treestyle.square.
        ltr (bool, optional): [description]. Defaults to True.

    Returns:
        LayoutDict: [description]
    """
    layout = defaultdict(TwoDCoordinate)
    previous_node = None
    y = 0
    # separation = equal_separation
    for node in tree.depth_first(post_order=True):
        node_coords = layout[node.ID]
        if node.children:
            child_x_coords, child_y_coords = zip(*(layout[c.ID] for c in node.children))
            node_coords.y = sum(child_y_coords) / len(node.children)
            increment = node.length if branchlengths else 1.0
            if ltr:
                node_coords.x = increment + max(child_x_coords)
            else:
                node_coords.x = min(child_x_coords) - increment
        else:
            if previous_node:
                y += 1.0
                layout[node.ID].y = y
            else:
                layout[node.ID].y = 0
            layout[node.ID].x = node.length if branchlengths else 0.0
            previous_node = node

    for node in tree.depth_first(post_order=True):
        layout[node.ID].x = (layout[tree.ID].x - layout[node.ID].x) * 1.0
        layout[node.ID].y = (layout[node.ID].y - layout[tree.ID].y) * 1.0

    if style == "radial":
        for node_id in layout.keys():
            layout[node_id] = layout[node_id].to_polar()
    return layout


def treeplot(
    tree: Tree,
    style: TreeStyle = TreeStyle.square,
    branchlengths: bool = True,
    ltr: bool = True,
    node_labels: bool = True,
    leaf_labels: bool = True,
    leaf_marker: Union[str, Callable, None] = "o",
    leaf_marker_fill: Union[str, Callable[[Tree], "str"], None] = "white",
    leaf_marker_edge: Union[str, Callable[[Tree], "str"], None] = "black",
    branch_linestyle: Union[dict, Callable[[Tree], dict], None] = None,
    ax: Optional[Ax] = None,
    return_layout: bool = False,
) -> Union[Ax, Tuple[Ax, LayoutDict]]:
    """[summary]

    Args:
        tree (Tree): [description]
        style (TreeStyle, optional): [description]. Defaults to TreeStyle.\
            square.
        branchlengths (bool, optional): [description]. Defaults to True.
        ltr (bool, optional): [description]. Defaults to True.
        node_labels (bool, optional): [description]. Defaults to True.
        leaf_labels (bool, optional): [description]. Defaults to True.
        leaf_marker (Union[str, Callable, None], optional): [description].\
             Defaults to 'o'.
        leaf_marker_fill (Union[str, Callable[[Tree],, optional): \
            [description]. Defaults to 'white'.
        leaf_marker_edge (Union[str, Callable[[Tree],, optional): \
            [description]. Defaults to 'black'.
        branch_linestyle (Union[dict, Callable[[Tree], dict], None], \
            optional): [description]. Defaults to None.
        ax (Optional[Ax], optional): [description]. Defaults to None.
        return_layout (bool, optional): [description]. Defaults to False.

    Returns:
        Union[Ax, Tuple[Ax, LayoutDict]]: [description]
    """
    layout = calculate_tree_layout(
        tree=tree, style=style, ltr=ltr, branchlengths=branchlengths
    )

    if not ax:
        _, ax = plt.subplots(figsize=(6, 6))

    default_linestyle = dict(linewidth=1, color="black", zorder=1)
    if branch_linestyle is None:

        def linestyle_fun(_) -> dict:
            return default_linestyle

    elif isinstance(branch_linestyle, dict):

        def linestyle_fun(_) -> dict:
            return {**default_linestyle, **branch_linestyle}

    elif callable(branch_linestyle):

        def linestyle_fun(branch: Tuple[Tree, Tree]) -> dict:
            return {**default_linestyle, **branch_linestyle(leaf)}

    else:
        raise TypeError(f"{type(branch_linestyle)} is not valid branch_linestyle type")

    for node1, node2 in tree.links:
        node1_x, node1_y = node1_coords = layout[node1.ID]
        node2_x, node2_y = node2_coords = layout[node2.ID]
        if node_labels:
            ax.text(
                node1_x,
                node1_y,
                node1.name,
                fontsize=8,
                verticalalignment="center_baseline",
            )
        if style == "square":
            ax.plot(
                (node1_x, node1_x), (node1_y, node2_y), **linestyle_fun((node1, node2))
            )
            ax.plot(
                (node1_x, node2_x), (node2_y, node2_y), **linestyle_fun((node1, node2))
            )
        elif style == "radial":
            if node2.root == node1:
                ax.plot(
                    (node1_x, node2_x),
                    (node1_y, node2_y),
                    **linestyle_fun((node1, node2)),
                )
            else:
                corner = TwoDCoordinate(
                    x=node1_coords.to_cartesian().x, y=node2_coords.to_cartesian().y
                ).to_polar()

                ax.plot(
                    (node1_x, corner.x),
                    (node1_y, corner.y),
                    **linestyle_fun((node1, node2)),
                )
                ax.plot(
                    (corner.x, node2_x),
                    (corner.y, node2_y),
                    **linestyle_fun((node1, node2)),
                )
        else:
            ax.plot(
                (node1_x, node2_x), (node1_y, node2_y), **linestyle_fun((node1, node2))
            )

    xmin, xmax = ax.get_xlim()
    xspacer = 0.25  # 0.01 * (xmax - xmin)

    if isinstance(leaf_marker, str):

        def leaf_marker_fun(_):
            return leaf_marker

    elif callable(leaf_marker):
        leaf_marker_fun = leaf_marker
    elif leaf_marker is None:
        pass
    else:
        raise TypeError(f"{type(leaf_marker)} is not a valid leaf_marker type")

    if isinstance(leaf_marker_fill, str) or leaf_marker_fill is None:

        def leaf_marker_fill_fun(_):
            return leaf_marker_fill

    elif callable(leaf_marker_fill):
        leaf_marker_fill_fun = leaf_marker_fill
    else:
        raise TypeError(f"{type(leaf_marker)} is not a valid leaf_marker type")

    if isinstance(leaf_marker_edge, str) or leaf_marker_edge is None:

        def leaf_marker_edge_fun(_):
            return leaf_marker_edge

    elif callable(leaf_marker_edge):
        leaf_marker_edge_fun = leaf_marker_edge
    else:
        raise TypeError(f"{type(leaf_marker)} is not a valid leaf_marker type")

    for leaf in tree.leaves:
        x, y = leaf_coords = layout[leaf.ID]
        if leaf_marker:
            ax.scatter(
                x,
                y,
                c=leaf_marker_fill_fun(leaf),
                alpha=1,
                edgecolors=leaf_marker_edge_fun(leaf),
                marker=leaf_marker_fun(leaf),
                zorder=2,
            )

        # x, y = leaf_coords = layout[leaf.ID]
        if style == "radial":
            # pass
            polar_coords = leaf_coords.to_polar()
            polar_coords.x *= 1.1
            x, y = leaf_coords = polar_coords.to_cartesian()
        else:
            x = x + xspacer if ltr else x - xspacer
            y += 0.05
        horizontalalignment = "left" if ltr else "right"
        ax.text(
            x,
            y,
            leaf.name,
            fontsize=12,
            in_layout=True,
            clip_on=True,
            verticalalignment="center_baseline",
            horizontalalignment=horizontalalignment,
        )

    ax.set_xticks(())
    xmin, xmax = ax.get_xlim()
    if style != "radial":
        if ltr:
            ax.set_xlim((0.8 * xmin, 1.8 * xmax))
        else:
            ax.set_xlim((1.8 * xmin, 1.2 * xmax))

    ax.set_yticks(())

    plt.tight_layout()

    if return_layout:
        return (ax, layout)
    return ax
