from typing import List
from itertools import groupby
from collections import defaultdict
from .dag import DirectedAcyclicGraph, DAGElement


class OntologyTerm(DAGElement):
    pass


class Ontology(DirectedAcyclicGraph):
    def __init__(self):
        super().__init__()
        self._header: List[str] = []

    @classmethod
    def from_obo(cls, filename: str = None, string: str = None) -> "Ontology":
        assert filename or string
        assert not (filename and string)
        ontology = cls()
        if filename:
            with open(filename) as filehandle:
                string = filehandle.read()

        obo_iter = (
            el
            for _, el in groupby(
                string.strip().split("\n"), lambda line: line[:1] == "["
            )
        )

        ontology._header = list(next(obo_iter))

        for element in obo_iter:
            element = next(element)
            if element != "[Term]":
                continue
            attributes = defaultdict(list)
            for attribute in next(obo_iter):
                if not attribute:
                    continue
                attr_key, attr_value = attribute.split(":", 1)
                attributes[attr_key].append(attr_value.strip())

            ID = attributes.pop("id")[0].strip()
            parents = [p.split("!")[0].strip() for p in attributes.get("is_a", "")]

            ontology[ID] = OntologyTerm(ID=ID, parents=parents, container=ontology)

        for ontology_term in ontology:
            for parent_id in ontology_term._parents:
                parent_term = ontology[parent_id]
                parent_term._children.append(ontology_term.ID)

        return ontology
