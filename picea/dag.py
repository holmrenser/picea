"""Direct Acyclic Graphs"""

from collections import defaultdict
from copy import deepcopy
from typing import Callable, DefaultDict, Dict, Hashable, Iterable, List, Optional
from warnings import warn


class DAGElement:
    """
    Element in a DAG, intended to be used as BaseClass
    """

    def __init__(
        self,
        ID: str,
        container: "DirectedAcyclicGraph",
        children: List[str] = None,
        parents: List[str] = None,
    ):
        self._ID = ID
        self._original_ID = ID
        self._container = container
        if children is None:
            children = []
        self._children = children
        if parents is None:
            parents = []
        self._parents = parents

    def __repr__(self):
        classname = type(self).__name__
        return f"<{classname} ID={self.ID} at {hex(id(self))}>"

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __deepcopy__(self, memo):
        """
        Deep copy of a DAG element, container (DAG) is explicitly not copied
        User has to reassign copied DAG element to original DAG if desired
        """
        new_element = self.__class__(ID=self.ID, container=self._container.__class__())
        for key, value in vars(self).items():
            # _container is strictly a reference, do not deepcopy
            if key != "_container":
                new_element[key] = deepcopy(value)
        return new_element

    @property
    def ID(self):
        return self._ID

    @ID.setter
    def ID(self, value):
        old_ID = self._ID
        if not self._original_ID:
            self._original_ID = old_ID
        self._ID = value
        for child in self.children:
            if child._parents:
                child._parents = [value if p == old_ID else p for p in child._parents]
        if self._container:
            self._container[value] = self._container.pop(old_ID)

    @property
    def parents(self) -> "DirectedAcyclicGraph":
        graph = self._container.__class__()
        for element in self._traverse(direction="parents"):
            if element == self:
                continue
            graph[element.ID] = element
        return graph

    @property
    def children(self) -> "DirectedAcyclicGraph":
        graph = self._container.__class__()
        for element in self._traverse(direction="children"):
            if element == self:
                continue
            graph[element.ID] = element
        return graph

    def _traverse(
        self,
        direction: str = "children",
        visited: Optional[set] = None,
    ) -> Iterable["DAGElement"]:
        if visited is None:
            visited = set()
        if self not in visited:
            yield self
            visited.add(self)
        if not getattr(self, f"_{direction}"):
            return
        for next_ID in getattr(self, f"_{direction}"):
            next_element = self._container[next_ID]
            yield from next_element._traverse(direction=direction, visited=visited)


class DirectedAcyclicGraph:
    """Base DAG class"""

    def __init__(self) -> None:
        self._elements: Dict[Hashable, DAGElement] = dict()

    def __len__(self) -> int:
        return len(self._elements)

    def __getitem__(self, ID: Hashable) -> DAGElement:
        return self._elements[ID]

    def __setitem__(self, ID: Hashable, element: DAGElement) -> None:
        if ID in self._elements:
            warn(f"Turning duplicate ID {ID} into unique ID")
            element._original_ID = ID
            modifier = 0
            new_ID = ID
            while new_ID in self._elements:
                modifier += 1
                new_ID = f"{ID}_{modifier}"
            ID = new_ID
            element._ID = ID
        self._elements[ID] = element

    def __contains__(self, ID):
        return ID in self._elements

    def __iter__(self) -> Iterable[DAGElement]:
        try:
            yield from self._elements.values()
        except StopIteration:
            return

    @property
    def elements(self) -> List[DAGElement]:
        return list(self)

    def add(self, element: DAGElement) -> None:
        self[element.ID] = element

    def pop(self, ID: Hashable) -> DAGElement:
        return self._elements.pop(ID)

    def groupby(self, group_func: Callable[[DAGElement], Hashable]) -> DefaultDict[Hashable, "DirectedAcyclicGraph"]:
        """
        Group elements in the current collection by calling 'group_func' and using \
        the results as a dictionary key

        Args:
            group_func: Callable[[DagElement], Hashable] receives as input an element \
                of the current collection and should return something that can be used \
                as a dictionary key (i.e. it should be hashable)

        Returns:
            Dictionary with keys for groups and values individual DirectedAcyclicGraph \
                subclasses
        """
        grouped = defaultdict(self.__class__)
        for element in self:
            grouped[group_func(element)][element.ID] = element
            # element._container = self
        return grouped

    def filter(self, filter_func: Callable[[DAGElement], bool]) -> "DirectedAcyclicGraph":
        """
        Filter elements in the current collection based on the output of filter_func

        Args:
            filter_func: Callable[[DagElement], bool] receives as input an element \
                of the current collection and returns a boolean indicating whether the \
                element should be kept (i.e. an element with a 'True' return value will\
                be kept)

        Returns:
            DirectedAcyclicGraph subclass with unwanted elements removed
        """
        result = self.__class__()
        for element in self:
            if filter_func(element):
                result[element.ID] = element
        return result
