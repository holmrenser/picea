## Development setup

Requires [uv](https://docs.astral.sh/uv/).

```
uv sync
uv run ruff check
uv run pytest
```

## Documentation

Docs are built with [Sphinx](https://www.sphinx-doc.org/) and [MyST-NB](https://myst-nb.readthedocs.io/). API pages
are generated from docstrings (Google style). Example notebooks live in `docs/examples/` as
[jupytext](https://jupytext.readthedocs.io/) percent-format `.py` files, and are executed on every docs build.

```
uv sync --group docs
uv run sphinx-build -W --keep-going docs docs/_build/html
```

To edit a notebook interactively, open it in VS Code (`# %%` cells), or run `uv run jupyter lab` and open the `.py`
file with *Open With → Notebook*.

## Releasing to PyPI

```
uv lock --check
uv run coverage run
uv run coverage report
uv version --bump <major,minor,patch>
uv build
uv publish
```
