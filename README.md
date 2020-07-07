_picea_
=======

Lightweight python library for working with trees and sequence collections

![tests](https://github.com/holmrenser/picea/workflows/tests/badge.svg)
![docs](https://github.com/holmrenser/picea/workflows/docs/badge.svg?branch=master)
[![Coverage Status](https://coveralls.io/repos/github/holmrenser/picea/badge.svg?branch=master)](https://coveralls.io/github/holmrenser/picea?branch=master)
[![PyPI version](https://badge.fury.io/py/picea.svg)](https://badge.fury.io/py/picea)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/picea.svg)](https://pypi.python.org/pypi/ansicolortags/)
[![PyPI status](https://img.shields.io/pypi/status/picea.svg)](https://pypi.python.org/pypi/ansicolortags/)



```
pip install picea
```

![example figure](https://github.com/holmrenser/picea/raw/master/docs/example1.png)

The above figure can be generated with the following code

```python
from picea import Tree
import matplotlib.pyplot as plt

newick = '(((a,b),(c,d)),e)'
tree = Tree.from_newick(newick)

fig, (ax1, ax2) = plt.subplots(ncols = 2, figsize = (10, 4))

#left-to-right layout with direct links
tree.layout(ltr = True)
for node1, node2 in tree.links:
    ax1.plot(
        (node1.x, node2.x),
        (node1.y, node2.y),
        c = 'k'
    )
for leaf in tree.leaves:
    ax1.text(
        leaf.x + .1, 
        leaf.y - .1, 
        leaf.name,
        fontsize = 18
    )

#right-to-left layout with square links
tree.layout(ltr = False)
for node1, node2 in tree.links:
    ax2.plot(
        (node1.x, node1.x),
        (node1.y, node2.y)
    )
    ax2.plot(
        (node1.x, node2.x),
        (node2.y, node2.y)
    )

#clean up plots
ax1.set_xlim((-.5, 3.5))
ax2.set_xlim((-3.5, .5))
for ax in (ax1,ax2):
    ax.set_xticks([],[])
    ax.set_yticks([],[])
```