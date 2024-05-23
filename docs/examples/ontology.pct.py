# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# This notebook shows how to work with biological ontologies such as the sequence ontology or the gene ontology.

# %%
import sys

import requests

sys.path.insert(0, '../../')
import picea

picea.__version__

# %%
obo_url = (
    'https://raw.githubusercontent.com/The-Sequence-Ontology/'
    'SO-Ontologies/master/Ontology_Files/so.obo'
)
r = requests.get(obo_url)
r

# %% tags=[]
r.text.split('\n')[:100]

# %%
so = picea.Ontology.from_obo(string=r.text)

# %%
ids = [el.ID for el in so['SO:0000866'].parents.elements]

# %%
'SO:0000866' in {el.ID for so_id in ids for el in so[so_id].children.elements}

# %%
len(so)

# %%
url = 'http://purl.obolibrary.org/obo/go.obo'
# url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
r = requests.get(url)
go = picea.Ontology.from_obo(string=r.text)
len(go.elements)

# %%
[(term.ID, term.name, len(term.parents)) for term in go['GO:0048316'].parents]

# %%
go['GO:0048316'].children

# %%
import networkx as nx

nx.__version__

# %%
graph = nx.DiGraph()
for term in [go['GO:0048316'], *go['GO:0048316'].children]:
    graph.add_node(term.ID, name=term.name)
    for child_ID in term._children:
        graph.add_edge(term.ID, child_ID)
layout = nx.planar_layout(graph)
nx.draw(graph, pos=layout, node_shape='s')

# %%
import sys

# !{sys.executable} -m pip install pygraphviz
nx.nx_agraph.to_agraph(graph)

# %%
[(term.ID, term.name) for term in go['GO:0048316'].children]

# %%
go['GO:0010431'].__dict__

# %%
go['GO:0048316'].__dict__

# %%
go['GO:0048316'].children._elements.keys()

# %% tags=[]
[(term.ID,term.name) for term in go if term.__dict__.get('alt_id') and term._parents]

# %%
[(term.ID,term.name) for term in go if not term.parents and term.children]

# %%
go['GO:0005554'].__dict__

# %%
go['GO:0003674'].__dict__
