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
# # Ontologies
#
# Biological ontologies such as the [Sequence Ontology](http://www.sequenceontology.org/) (SO) or the
# [Gene Ontology](https://geneontology.org/) are read from OBO files into an {class}`~picea.Ontology`: a directed
# acyclic graph of {class}`~picea.OntologyTerm` objects linked by `is_a` and `part_of` relationships.

# %%
from urllib.request import urlopen

from picea import Ontology

# %%
SO_URL = "https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/master/Ontology_Files/so.obo"
with urlopen(SO_URL) as response:
    so = Ontology.from_obo(string=response.read().decode())
len(so)

# %% [markdown]
# Terms are accessed by ID. OBO tags become attributes, with list values.

# %%
mrna = so["SO:0000234"]
mrna.name, mrna["def"]

# %% [markdown]
# `parents` and `children` are **transitive**: they contain all ancestors and all descendants of a term.

# %%
[(term.ID, term.name[0]) for term in mrna.parents]

# %%
len(mrna.children)

# %% [markdown]
# Both return a new {class}`~picea.Ontology`, so they can be filtered and grouped.

# %%
coding_rnas = mrna.children.filter(lambda term: "coding" in term.name[0])
sorted(term.name[0] for term in coding_rnas)
