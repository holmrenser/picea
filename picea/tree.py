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
        If some nodes have a branch length, other nodes (except the root) without one get length 0.0 with a warning.

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
        if found_branchlengths:
            tree.cumulative_length = 0.0
        queue: list[Tree] = [tree]
        while queue:
            node = queue.pop(0)
            for child in node.children:
                child.parent = node
                child.depth = node.depth + 1
                if found_branchlengths:
                    if child.length is None:
                        warn(
                            "Found branchlengths on some parts of the tree, but node "
                            f"{child.ID} has no branchlength specified, setting to "
                            "branchlength 0.0"
                        )
                        child.length = 0.0
                    child.cumulative_length = node.cumulative_length + abs(child.length)
            queue += node.children

        return tree

    def to_newick(self, branch_lengths: bool = False) -> str:
        """Newick formatted string of the (sub)tree. The node this is called on is the root of the Newick tree, so its
        branch length is not written.

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
        return f"{self._to_newick(branch_lengths=branch_lengths, include_length=False)};"

    def _to_newick(self, branch_lengths: bool, include_length: bool = True) -> str:
        """Newick string of the (sub)tree without the closing semicolon"""
        name = str(self.name) if self.name else ""

        if self.children:
            subtree_string = ",".join(child._to_newick(branch_lengths=branch_lengths) for child in self.children)
            newick = f"({subtree_string}){name}"
        else:
            newick = name

        if branch_lengths and include_length:
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
            inplace (bool): Rename the leaves of this tree. If False, rename and return a copy.

        Returns:
            Optional[Tree]: The renamed copy if ``inplace`` is False, otherwise None
        """
        tree = self if inplace else deepcopy(self)
        for leaf in tree.leaves:
            leaf.name = rename_func(leaf.name)
        return None if inplace else tree


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
    if node_a.parent is node_b.parent:
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


Ax = Type[SubplotBase]
TreeStyle = Enum("TreeStyle", ("square", "radial", "triangular"))
LayoutDict = DefaultDict[int, TwoDCoordinate]


def _tree_style(style: Union[str, TreeStyle]) -> TreeStyle:
    """Convert a style name to a TreeStyle"""
    if isinstance(style, TreeStyle):
        return style
    try:
        return TreeStyle[style]
    except KeyError:
        raise ValueError(f"Unknown tree style {style!r}, must be one of {[s.name for s in TreeStyle]}") from None


def _rectangular_layout(tree: Tree, ltr: bool, branchlengths: bool) -> Tuple[LayoutDict, int]:
    """Rectangular layout (see :func:`calculate_tree_layout`) and the number of leaves"""
    layout: LayoutDict = defaultdict(TwoDCoordinate)

    n_leaves = 0
    for node in tree.depth_first(post_order=True):
        if node.children:
            layout[node.ID].y = sum(layout[child.ID].y for child in node.children) / len(node.children)
        else:
            layout[node.ID].y = float(n_leaves)
            n_leaves += 1

    if branchlengths and any(node.length is not None for node in tree.nodes if node is not tree):
        # distance to the root
        stack = [(tree, 0.0)]
        while stack:
            node, x = stack.pop()
            layout[node.ID].x = x
            stack.extend((child, x + abs(child.length or 0.0)) for child in node.children)
    else:
        # cladogram: levels below the root, with all leaves aligned
        heights = dict()
        for node in tree.depth_first(post_order=True):
            heights[node.ID] = 1 + max(heights[child.ID] for child in node.children) if node.children else 0
        for node in tree.nodes:
            layout[node.ID].x = float(heights[tree.ID] - heights[node.ID])

    if not ltr:
        for coordinate in layout.values():
            coordinate.x = -coordinate.x
    return layout, n_leaves


def _polar(coordinate: TwoDCoordinate, n_leaves: int) -> Tuple[float, float]:
    """Radius and angle of a rectangular layout coordinate in the radial layout"""
    return abs(coordinate.x), 2 * np.pi * coordinate.y / max(n_leaves, 1)


def calculate_tree_layout(
    tree: Tree,
    style: Union[str, TreeStyle] = TreeStyle.square,
    ltr: bool = True,
    branchlengths: bool = True,
) -> LayoutDict:
    """Calculate 2D coordinates of all nodes, as used by :func:`treeplot`

    Leaves get consecutive y coordinates (0, 1, 2, ...) in depth first order, and internal nodes are centered on
    their children. With branch lengths, x is the distance to the root. Without branch lengths, or for trees that
    have none, x is the number of levels below the root, with all leaves aligned (a cladogram). In the radial layout
    the root is in the center, the radius is the distance to the root, and leaves are spread evenly over the circle.

    Examples:
        >>> tree = Tree.from_newick('((a:1,b:2)ab:1,c:3)root;')
        >>> layout = calculate_tree_layout(tree)
        >>> [(node.name, layout[node.ID].x, layout[node.ID].y) for node in tree.leaves]
        [('c', 3.0, 2.0), ('a', 2.0, 0.0), ('b', 3.0, 1.0)]

    Args:
        tree (Tree): Tree
        style (Union[str, TreeStyle], optional): ``"square"``, ``"triangular"`` or ``"radial"``. Only ``"radial"``
            changes the node coordinates. Defaults to ``TreeStyle.square``.
        ltr (bool, optional): Left to right layout (root on the left). If False, x coordinates are negative.
            Ignored for the radial layout. Defaults to True.
        branchlengths (bool, optional): Use branch lengths. Defaults to True.

    Returns:
        LayoutDict: Coordinates of every node, by node ID

    Raises:
        ValueError: If ``style`` is not a valid style
    """
    style = _tree_style(style)
    layout, n_leaves = _rectangular_layout(tree, ltr=ltr, branchlengths=branchlengths)
    if style == TreeStyle.radial:
        for coordinate in layout.values():
            radius, angle = _polar(coordinate, n_leaves)
            coordinate.x, coordinate.y = radius * np.cos(angle), radius * np.sin(angle)
    return layout


def _as_function(value: Union[str, dict, Callable, None], name: str) -> Callable:
    """Wrap a constant plotting option in a function, so that constant and per-node options can be used the same way"""
    if callable(value):
        return value
    if value is None or isinstance(value, (str, dict)):
        return lambda _: value
    raise TypeError(f"{type(value)} is not a valid {name} type")


def treeplot(
    tree: Tree,
    style: Union[str, TreeStyle] = TreeStyle.square,
    branchlengths: bool = True,
    ltr: bool = True,
    node_labels: bool = True,
    leaf_labels: bool = True,
    leaf_marker: Union[str, Callable, None] = "o",
    leaf_marker_fill: Union[str, Callable[[Tree], "str"], None] = "white",
    leaf_marker_edge: Union[str, Callable[[Tree], "str"], None] = "black",
    branch_linestyle: Union[dict, Callable[[Tuple[Tree, Tree]], dict], None] = None,
    ax: Optional[Ax] = None,
    return_layout: bool = False,
) -> Union[Ax, Tuple[Ax, LayoutDict]]:
    """Plot a tree with matplotlib. See the :doc:`/examples/trees` example for usage.

    Args:
        tree (Tree): Tree to plot
        style (Union[str, TreeStyle], optional): Branch style: ``"square"`` (right angles), ``"triangular"``
            (straight lines from parent to child), or ``"radial"`` (root in the center). Defaults to
            ``TreeStyle.square``.
        branchlengths (bool, optional): Scale branches by their length. Use False to plot a cladogram. Trees without
            branch lengths are always plotted as a cladogram. Defaults to True.
        ltr (bool, optional): Plot left to right (root on the left). Ignored for the radial style. Defaults to True.
        node_labels (bool, optional): Show the names of internal nodes, e.g. support values. Defaults to True.
        leaf_labels (bool, optional): Show leaf names. Defaults to True.
        leaf_marker (Union[str, Callable, None], optional): Matplotlib marker for leaves, a function that takes a
            leaf and returns a marker, or None for no markers. Defaults to ``"o"``.
        leaf_marker_fill (Union[str, Callable[[Tree], str], None], optional): Marker fill color, or a function
            that takes a leaf and returns a color. Defaults to ``"white"``.
        leaf_marker_edge (Union[str, Callable[[Tree], str], None], optional): Marker edge color, or a function
            that takes a leaf and returns a color. Defaults to ``"black"``.
        branch_linestyle (Union[dict, Callable[[Tuple[Tree, Tree]], dict], None], optional): Keyword arguments for
            :meth:`matplotlib.axes.Axes.plot` used to draw branches, or a function that takes a
            ``(parent, child)`` tuple and returns them. Defaults to None (thin black lines).
        ax (Optional[Ax], optional): Axes to plot on. A new figure is created when not given.
        return_layout (bool, optional): Also return the node coordinates (see :func:`calculate_tree_layout`).
            Defaults to False.

    Returns:
        Union[Ax, Tuple[Ax, LayoutDict]]: The axes, or an ``(axes, layout)`` tuple if ``return_layout`` is True

    Raises:
        ValueError: If ``style`` is not a valid style
        TypeError: If a styling option has an invalid type
    """
    style = _tree_style(style)
    radial = style == TreeStyle.radial
    rectangular_layout, n_leaves = _rectangular_layout(tree, ltr=ltr, branchlengths=branchlengths)
    layout = calculate_tree_layout(tree=tree, style=style, ltr=ltr, branchlengths=branchlengths)

    default_linestyle = dict(linewidth=1, color="black", zorder=1)
    branch_linestyle_fun = _as_function(branch_linestyle, "branch_linestyle")
    leaf_marker_fun = _as_function(leaf_marker, "leaf_marker")
    leaf_marker_fill_fun = _as_function(leaf_marker_fill, "leaf_marker_fill")
    leaf_marker_edge_fun = _as_function(leaf_marker_edge, "leaf_marker_edge")

    if not ax:
        _, ax = plt.subplots(figsize=(6, 6))

    for parent, child in tree.links:
        linestyle = {**default_linestyle, **(branch_linestyle_fun((parent, child)) or dict())}
        parent_x, parent_y = layout[parent.ID]
        child_x, child_y = layout[child.ID]
        if style == TreeStyle.square:
            ax.plot((parent_x, parent_x, child_x), (parent_y, child_y, child_y), **linestyle)
        elif style == TreeStyle.triangular:
            ax.plot((parent_x, child_x), (parent_y, child_y), **linestyle)
        else:
            # arc at the radius of the parent, followed by a radial line to the child
            parent_radius, parent_angle = _polar(rectangular_layout[parent.ID], n_leaves)
            child_radius, child_angle = _polar(rectangular_layout[child.ID], n_leaves)
            angles = np.linspace(parent_angle, child_angle, max(2, int(abs(child_angle - parent_angle) / 0.02)))
            radii = np.append(np.full(angles.shape, parent_radius), child_radius)
            angles = np.append(angles, child_angle)
            ax.plot(radii * np.cos(angles), radii * np.sin(angles), **linestyle)

    if node_labels:
        for node in tree.nodes:
            if node.children and node.name:
                x, y = layout[node.ID]
                ax.annotate(
                    node.name,
                    (x, y),
                    xytext=(-2, 2) if ltr or radial else (2, 2),
                    textcoords="offset points",
                    fontsize=8,
                    horizontalalignment="right" if ltr or radial else "left",
                    verticalalignment="bottom",
                )

    for leaf in tree.leaves:
        x, y = layout[leaf.ID]
        if leaf_marker is not None:
            ax.scatter(
                x,
                y,
                c=leaf_marker_fill_fun(leaf),
                alpha=1,
                edgecolors=leaf_marker_edge_fun(leaf),
                marker=leaf_marker_fun(leaf),
                zorder=2,
            )
        if not leaf_labels:
            continue
        if radial:
            _, angle = _polar(rectangular_layout[leaf.ID], n_leaves)
            degrees = np.degrees(angle) % 360
            flip = 90 < degrees < 270
            text_options = dict(
                xytext=(6 * np.cos(angle), 6 * np.sin(angle)),
                rotation=degrees - 180 if flip else degrees,
                rotation_mode="anchor",
                horizontalalignment="right" if flip else "left",
            )
        else:
            text_options = dict(xytext=(6, 0) if ltr else (-6, 0), horizontalalignment="left" if ltr else "right")
        ax.annotate(
            leaf.name,
            (x, y),
            textcoords="offset points",
            fontsize=12,
            verticalalignment="center",
            **text_options,
        )

    ax.set_xticks(())
    ax.set_yticks(())
    xmin, xmax = ax.get_xlim()
    width = xmax - xmin
    if radial:
        # leave room for the leaf labels on all sides
        ax.set_aspect("equal")
        ymin, ymax = ax.get_ylim()
        ax.set_xlim((xmin - 0.4 * width, xmax + 0.4 * width))
        ax.set_ylim((ymin - 0.4 * (ymax - ymin), ymax + 0.4 * (ymax - ymin)))
    elif ltr:
        ax.set_xlim((xmin, xmax + 0.8 * width))
    else:
        ax.set_xlim((xmin - 0.8 * width, xmax))

    ax.figure.tight_layout()

    if return_layout:
        return (ax, layout)
    return ax
