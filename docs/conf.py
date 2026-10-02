"""Sphinx configuration for the picea documentation"""

from importlib.metadata import version as get_version

project = "picea"
author = "Rens Holmer"
copyright = "Rens Holmer"
release = get_version("picea")
version = release

extensions = [
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
    "sphinx_copybutton",
]

# Example notebooks are stored as jupytext percent-format .py files, so conf.py must be excluded
exclude_patterns = ["_build", "conf.py", "**/.ipynb_checkpoints", ".DS_Store"]

# -- API docs ----------------------------------------------------------------
autodoc_default_options = {
    "members": True,
    "inherited-members": "set",  # Alphabet subclasses set, don't document set methods
    "member-order": "bysource",
}
autodoc_typehints = "description"
autodoc_typehints_description_target = "documented"
napoleon_google_docstring = True
napoleon_numpy_docstring = False
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
}

# -- Notebooks ---------------------------------------------------------------
myst_enable_extensions = ["colon_fence"]
nb_custom_formats = {".py": ["jupytext.reads", {"fmt": "py:percent"}]}
nb_execution_mode = "cache"
nb_execution_raise_on_error = True
nb_execution_timeout = 120

# -- HTML --------------------------------------------------------------------
html_theme = "furo"
html_title = f"picea {release}"
html_theme_options = {
    "source_repository": "https://github.com/holmrenser/picea",
    "source_branch": "master",
    "source_directory": "docs/",
}
