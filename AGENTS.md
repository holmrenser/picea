# AGENTS.md

picea: small pure-Python library of bioinformatics data structures: sequences, annotations (GFF3/GTF), trees (Newick), ontologies (OBO). Runtime deps are only numpy and matplotlib.

## Layout
- `picea/sequence.py`: `Alphabet`, `Sequence`, readers, sequence collections, `SequenceAnnotation`/`SequenceInterval`
- `picea/dag.py`: `DirectedAcyclicGraph`/`DAGElement`, the base for annotations and ontologies
- `picea/ontology.py`: `Ontology`/`OntologyTerm`
- `picea/tree.py`: `Tree`, layout, matplotlib `treeplot`
- Public API = whatever `picea/__init__.py` re-exports.

## Design fundamentals
- **Symmetric I/O.** Each format is `Cls.from_<fmt>(filename=None, string=None)` (pass exactly one) plus `obj.to_<fmt>() -> str`. Whole files are read into memory; there is no streaming.
- **Storage-agnostic collections.** `AbstractSequenceCollection` subclasses implement only storage (`__get/set/delitem__`, `pop`, `headers`, `n_seqs`) and inherit I/O, iteration, `iloc`, renaming. `SequenceCollection` stores `{header: str}`; `MultipleSequenceAlignment` stores a gap-padded `uint8` numpy matrix plus `{header: row}`.
- **`Sequence` is a throwaway view.** Collections store raw strings and build a new `Sequence` on every access. Transforms return new objects. The alphabet (DNA/AminoAcid) is auto-detected by `guess_alphabet` unless given; derived sequences keep it.
- **No silent overwrite.** `__setitem__` on collections and DAGs never replaces: a duplicate header/ID becomes `X_1`, `X_2`, … and a warning is issued.
- **DAG by ID reference.** A DAG is `{ID: element}`; elements hold parent/child *IDs* plus a reference to their container. Files give child→parent links, and `_link_parents()` fills in children after load. `.children`/`.parents` are **transitive** (all descendants/ancestors). They return a new container of the same class, so `groupby`/`filter` chain.
- **GFF column 9 → instance attributes.** Keys are lowercased (except `ID`), and values are always lists (`iv.name == iv["name"]`). Keys that clash with the 8 fixed columns get a `_` prefix. Predefined keys are re-capitalized on write.
- **Tree = recursive dataclass node** (`name`, `length`, `children`). The parser sets `parent`, `depth`, `ID` and `cumulative_length`, which are not dataclass fields. `ID` is the parse-order index used by `iloc` and as the layout key. Equality is structural, so nodes are unhashable.

## Working here
- uv (`uv_build` backend, flat layout), Python ≥3.11, numpy 1.26 and 2.x: `uv sync`, `uv run pytest`, `uv run ruff check`. Black/ruff line length is 120.
- Tests are `tests/*_tests.py` (unittest, classes `*Tests`) **plus doctests** in `picea/` (`--doctest-modules`). Docstring examples are tests.
- Docs: Sphinx + MyST-NB + Furo in `docs/`. Build with `uv sync --group docs`, then `uv run sphinx-build -W --keep-going docs docs/_build/html`.
  - API pages (`docs/api/*.rst`) are autodoc of the public API. Docstrings are Google style and rendered as reStructuredText, so use `*emphasis*`, not `_emphasis_`.
  - Examples (`docs/examples/*.py`) are jupytext percent notebooks, executed on every build. A failing cell fails the build. Paths are relative to `docs/examples/`.
- `__version__` comes from installed metadata. A stale `picea.egg-info/` in the repo root shadows it, so delete that directory if you see one.
- Release: `uv version --bump <part>`, then a `vX.Y.Z` commit (see CONTRIBUTING.md).
