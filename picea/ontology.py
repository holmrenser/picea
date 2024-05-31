import warnings
from collections import defaultdict
from itertools import groupby
from typing import List

from .dag import DAGElement, DirectedAcyclicGraph


class OntologyTerm(DAGElement):
    _predefined_ontology_elements = ("name", "def")

    def __init__(self, ID, parents, container, **kwargs):
        super().__init__(ID=ID, parents=parents, container=container)

        for attr in self._predefined_ontology_elements:
            self[attr] = kwargs.get(attr, None)

        for key, value in kwargs.items():
            self[key] = value


class Ontology(DirectedAcyclicGraph):
    def __init__(self):
        super().__init__()
        self._header: List[str] = []

    def __getitem__(self, ID) -> OntologyTerm:
        term = self._elements[ID]
        if not term._children and not term._parents and term.__dict__.get("alt_id"):
            alt_id = term.__dict__.get("alt_id")[0]
            term = self._elements[alt_id]
            warnings.warn(f"Accessed GO term by alt ID {ID}, " f"returning main GO term with ID {alt_id}")
        return term

    @classmethod
    def from_obo(cls, filename: str = None, string: str = None, skip_obsolete=True) -> "Ontology":
        assert filename or string
        assert not (filename and string)
        ontology = cls()
        if filename:
            with open(filename) as filehandle:
                string = filehandle.read()

        obo_iter = (el for _, el in groupby(string.strip().split("\n"), lambda line: line[:1] == "["))

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
            if skip_obsolete and attributes.get("is_obsolete"):
                continue

            ID = attributes.pop("id")[0].strip()
            parents = [p.split("!")[0].strip() for p in attributes.get("is_a", "")]

            for relationship in attributes.get("relationship", []):
                relation_type, go_id = relationship.split("!")[0].strip().split(" ")
                if relation_type == "part_of":
                    parents.append(go_id)

            alt_ids = {*attributes.pop("alt_id", []), ID}
            for alt_id in alt_ids:
                ontology[alt_id] = OntologyTerm(
                    ID=alt_id,
                    parents=parents,
                    container=ontology,
                    alt_id=[id for id in alt_ids if id != alt_id],
                    **attributes,
                )

        for ontology_term in ontology:
            for parent_id in ontology_term._parents:
                parent_term = ontology[parent_id]
                parent_term._children.append(ontology_term.ID)

        return ontology
