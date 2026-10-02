import json
import re
from collections import defaultdict
from copy import deepcopy
from dataclasses import InitVar, asdict, dataclass, field
from enum import Enum
from typing import (
    Callable,
    DefaultDict,
    Dict,
    Generator,
    Iterable,
    List,
    Optional,
    Tuple,
    Type,
    Union,
)
from warnings import warn

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import SubplotBase

TreeDict = Dict[str, Union[str, int, float, List[Optional["TreeDict"]]]]


@dataclass
class Tree:
    """Tree (e.g. a phylogeny). Trees are recursive: every node is a :class:`Tree` with a list of child nodes, so
    every node is also the root of its own subtree.

    Trees are usually created with :meth:`from_newick` or :meth:`from_sklearn`. These also set the derived node
    attributes ``parent``, ``depth``, ``cumulative_length`` (distance to the root) and ``ID`` (see :attr:`iloc`).

    Examples:
        >>> tree = Tree.from_newick('((a,b)ab,c);')
        >>> [leaf.name for leaf in tree.leaves]
        ['c', 'a', 'b']
        >>> tree.loc['a'].parent.name
        'ab'

    Args:
        name (Optional[str]): Node name
        length (Optional[float]): Length of the branch to the parent node
        children (Optional[List[Tree]]): Child nodes
        ID (Optional[int]): Node ID
    """

    name: Optional[str] = None
    length: Optional[float] = None
    children: Optional[List["Tree"]] = field(default_factory=list)

    ID: InitVar[Optional[int]] = None
    depth: InitVar[Optional[int]] = None
    parent: InitVar[Optional["Tree"]] = None
    cumulative_length: InitVar[Optional[float]] = None

    def __post_init__(self, ID, *args, **kwargs):
        """Store the node ID. The other init-only fields are set by the parsers."""
        self.ID = ID

    @property
    def loc(self) -> "TreeIndex":
        """Name based index: ``tree.loc[name]`` returns the first node with that name (in post-order depth first
        traversal)

        Example:
            >>> from picea import Tree
            >>> newick = '(((a,b),(c,d)),e);'
            >>> tree = Tree.from_newick(newick)
            >>> tree.loc['a']
            Tree(name='a', length=None, children=[])

        Raises:
            IndexError: If no node has this name
        """
        return TreeIndex(iterator=self.depth_first, eq_func=lambda node, name: node.name == name)

    @property
    def iloc(self) -> "TreeIndex":
        """ID based index: ``tree.iloc[ID]`` returns the node with that ID. :meth:`from_newick` numbers nodes in
        pre-order, starting with 0 for the root.

        Example:
            >>> from picea import Tree
            >>> newick = '(((a,b),(c,d)),e);'
            >>> tree = Tree.from_newick(newick)
            >>> [leaf.name for leaf in tree.iloc[2].leaves]
            ['a', 'b']

        Raises:
            IndexError: If no node has this ID
        """
        return TreeIndex(iterator=self.depth_first, eq_func=lambda node, index: node.ID == index)

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
    def from_newick(cls, string: Optional[str] = None, filename: Optional[str] = None) -> "Tree":
        """Parse a Newick formatted file or string. Exactly one of ``string`` or ``filename`` must be given.

        Node names (including internal node names such as support values) and branch lengths are read when present.
        If some nodes have a branch length, nodes without one get length 0.0 with a warning.

        Examples:
            >>> tree = Tree.from_newick('((a:1,b:2)ab:1,c:3)root:0;')
            >>> tree.loc['b'].cumulative_length
            3.0

        Args:
            string (Optional[str]): Newick formatted string
            filename (Optional[str]): Newick filename

        Returns:
            Tree: Root node
        """
        assert filename or string
        assert not (filename and string)
        if filename:
            with open(filename) as filehandle:
                string = filehandle.read()
        tokens: list[str] = re.split(r"\s*(;|\(|\)|,|:)\s*", string)
        ID = 0
        tree = cls(ID=ID)
        ancestors: list[Tree] = list()
        found_branchlengths = False
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
                    found_branchlengths = True
                    tree.length = float(token)
                    tree.cumulative_length = 0.0
        tree.depth = 0
        queue: list[Tree] = [tree]
        while queue:
            node = queue.pop(0)
            if found_branchlengths:
                if node.length is None:
                    warn(
                        "Found branchlengths on some parts of the tree, but node "
                        f"{node.ID} has no branchlength specified, setting to "
                        "branchlength 0.0"
                    )
                    node.length = 0.0
                    node.cumulative_length = 0.0
            for child in node.children:
                child.parent = node
                child.depth = node.depth + 1
                if child.length:
                    child.cumulative_length = node.cumulative_length + abs(child.length)
            queue += node.children

        return tree

    def to_newick(self, branch_lengths: bool = False) -> str:
        """Newick formatted string of the (sub)tree

        Examples:
            >>> tree = Tree.from_newick('((a:1,b:2)ab:1,c:3)root:0;')
            >>> tree.to_newick(branch_lengths=True)
            '((a:1.0,b:2.0)ab:1.0,c:3.0)root;'

        Args:
            branch_lengths (bool, optional): Include branch lengths. Nodes without a branch length are written with
                length 0, with a warning. Defaults to False.

        Returns:
            str: Newick formatted string
        """
        if self.name:
            name = str(self.name)
        else:
            name = ""

        if self.children:
            subtree_string = ",".join([c.to_newick(branch_lengths=branch_lengths) for c in self.children])
            newick = f"({subtree_string}){name}"
        else:
            newick = name

        if branch_lengths and self.ID != 0:
            length = self.length
            if length is None:
                warn(
                    "Trying to write branch length for node that has no branch length \
                     set, defaulting to zero length branch."
                )
                length = 0
            if length == 0:
                length = int(0)
            newick += f":{length}"

        if self == self.root:
            newick += ";"

        return newick

    @classmethod
    def from_sklearn(cls, clustering) -> "Tree":
        """Create a tree from a fitted scikit-learn ``AgglomerativeClustering`` model. Leaves are named by sample
        index, and the tree has no branch lengths.

        Args:
            clustering (sklearn.cluster.AgglomerativeClustering): Fitted clustering model

        Returns:
            Tree: Root node
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
        """Not implemented yet"""
        raise NotImplementedError()

    @classmethod
    def from_json(cls):
        # TODO
        """Not implemented yet"""
        raise NotImplementedError()

    def to_json(self, indent: Optional[int] = None) -> str:
        """json formatted string of :meth:`to_dict`

        Args:
            indent (Optional[int]): Indentation, passed to :func:`json.dumps`

        Returns:
            str: json formatted string
        """
        return json.dumps(self.to_dict(), indent=indent)

    @classmethod
    def from_dict(cls, tree_dict):
        # TODO
        """Not implemented yet"""
        raise NotImplementedError()
        # tree = cls()
        # return tree

    def to_dict(self) -> TreeDict:
        """Nested dictionary with the ``name``, ``length`` and ``children`` of every node

        Returns:
            TreeDict: Tree dictionary
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

    def rename_leaves(self, rename_func: Callable, inplace: bool = True) -> Optional["Tree"]:
        """Rename all leaves by calling ``rename_func`` on every leaf name

        Args:
            rename_func (Callable[[str], str]): Function that takes a leaf name and returns a new name
            inplace (bool): Rename the leaves of this tree. If False, rename the leaves of a copy.
        """
        tree = self if inplace else deepcopy(self)
        for leaf in tree.leaves:
            leaf.name = rename_func(leaf.name)


class TreeIndex(object):
    """Index into a tree, see :attr:`Tree.loc` and :attr:`Tree.iloc`

    Args:
        iterator (Callable[[], Iterable[Tree]]): Function that returns an iterator over all nodes
        eq_func (Callable[[Tree, Any], bool]): Function that tests whether a node matches an index key
    """

    def __init__(self, iterator: Iterable[Tree], eq_func: Callable[[int, str], bool]):
        self.iterator = iterator
        self.eq_func = eq_func

    def __getitem__(self, key):
        for element in self.iterator():
            if self.eq_func(element, key):
                return element
        raise IndexError(f"{key} is not valid index")


def unequal_separation(node_a: "Tree", node_b: "Tree", sep_1: float = 1.0, sep_2: float = 2.0) -> float:
    """Separation between two neighbouring nodes: ``sep_1`` for siblings, ``sep_2`` otherwise

    Args:
        node_a (Tree): First node
        node_b (Tree): Second node
        sep_1 (float, optional): Separation between siblings. Defaults to 1.0.
        sep_2 (float, optional): Separation between non-siblings. Defaults to 2.0.

    Returns:
        float: Separation
    """
    if node_a.parent == node_b.parent:
        return sep_1
    return sep_2


def equal_separation(node_a: "Tree", node_b: "Tree", separation: float = 1.0) -> float:
    """Constant separation between two neighbouring nodes

    Args:
        node_a (Tree): First node
        node_b (Tree): Second node
        separation (float, optional): Separation. Defaults to 1.0.

    Returns:
        float: Separation
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
        return TwoDCoordinate(x=np.sqrt(self.x**2 + self.y**2), y=np.arctan2(self.y, self.x))


Ax = Type[SubplotBase]
TreeStyle = Enum("TreeStyle", ("square", "radial", "triangular"))
LayoutDict = DefaultDict[int, TwoDCoordinate]


def calculate_tree_layout(
    tree: Tree,
    style: TreeStyle = TreeStyle.square,
    ltr: bool = True,
    branchlengths: bool = True,
) -> LayoutDict:
    """Calculate 2D coordinates of all nodes, as used by :func:`treeplot`

    Leaves get consecutive y coordinates in depth first order, and internal nodes are centered on their children.
    The x coordinate is based on branch lengths, or on the number of levels below a node if ``branchlengths`` is
    False.

    Args:
        tree (Tree): Tree
        style (TreeStyle, optional): ``"square"``, ``"triangular"`` or ``"radial"``. Only ``"radial"`` changes
            the layout (to polar coordinates). Defaults to ``TreeStyle.square``.
        ltr (bool, optional): Left to right layout (root on the left). Defaults to True.
        branchlengths (bool, optional): Use branch lengths. Defaults to True.

    Returns:
        LayoutDict: Coordinates of every node, by node ID
    """
    layout = defaultdict(TwoDCoordinate)
    previous_node = None
    y = 0
    # separation = equal_separation
    for node in tree.depth_first(post_order=True):
        node_coords = layout[node.ID]
        if node.children:
            child_x_coords, child_y_coords = zip(*(layout[c.ID] for c in node.children), strict=True)
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
    """Plot a tree with matplotlib. See the :doc:`/examples/trees` example for usage.

    Args:
        tree (Tree): Tree to plot
        style (TreeStyle, optional): Branch style: ``"square"`` (right angles), ``"triangular"`` (straight lines
            from parent to child), or ``"radial"``. Defaults to ``TreeStyle.square``.
        branchlengths (bool, optional): Scale branches by their length. Use False to plot a cladogram, or to
            plot a tree without branch lengths. Defaults to True.
        ltr (bool, optional): Plot left to right (root on the left). Defaults to True.
        node_labels (bool, optional): Show the names of internal nodes, e.g. support values. Defaults to True.
        leaf_labels (bool, optional): Currently unused: leaf names are always shown. Defaults to True.
        leaf_marker (Union[str, Callable, None], optional): Matplotlib marker for leaves, a function that takes a
            leaf and returns a marker, or None for no markers. Defaults to ``"o"``.
        leaf_marker_fill (Union[str, Callable[[Tree], str], None], optional): Marker fill color, or a function
            that takes a leaf and returns a color. Defaults to ``"white"``.
        leaf_marker_edge (Union[str, Callable[[Tree], str], None], optional): Marker edge color, or a function
            that takes a leaf and returns a color. Defaults to ``"black"``.
        branch_linestyle (Union[dict, Callable[[Tree], dict], None], optional): Keyword arguments for
            :meth:`matplotlib.axes.Axes.plot` used to draw branches, or a function that takes a
            ``(parent, child)`` tuple and returns them. Defaults to None (thin black lines).
        ax (Optional[Ax], optional): Axes to plot on. A new figure is created when not given.
        return_layout (bool, optional): Also return the node coordinates (see :func:`calculate_tree_layout`).
            Defaults to False.

    Returns:
        Union[Ax, Tuple[Ax, LayoutDict]]: The axes, or an ``(axes, layout)`` tuple if ``return_layout`` is True
    """
    layout = calculate_tree_layout(tree=tree, style=style, ltr=ltr, branchlengths=branchlengths)

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
            ax.plot((node1_x, node1_x), (node1_y, node2_y), **linestyle_fun((node1, node2)))
            ax.plot((node1_x, node2_x), (node2_y, node2_y), **linestyle_fun((node1, node2)))
        elif style == "radial":
            if node2.root == node1:
                ax.plot(
                    (node1_x, node2_x),
                    (node1_y, node2_y),
                    **linestyle_fun((node1, node2)),
                )
            else:
                corner = TwoDCoordinate(x=node1_coords.to_cartesian().x, y=node2_coords.to_cartesian().y).to_polar()

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
            ax.plot((node1_x, node2_x), (node1_y, node2_y), **linestyle_fun((node1, node2)))

    xmin, xmax = ax.get_xlim()
    xspacer = 0.0  # 25  # 0.01 * (xmax - xmin)

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
