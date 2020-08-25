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
from picea import Tree, treeplot
import matplotlib.pyplot as plt

newick = '(((a,b),(c,d)),e)'
tree = Tree.from_newick(newick)

fig, (ax1, ax2) = plt.subplots(ncols = 2, figsize = (10, 4))

#left-to-right layout with direct links
treeplot(tree, style='rectangular', ltr=True, ax=ax1)

#right-to-left layout with square links
treeplot(tree, style='square', ltr=False, ax=ax2)
```