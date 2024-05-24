from importlib.metadata import version

__version__ = version(__name__)

from .tree import (  # noqa
    Tree,
    calculate_tree_layout,
    treeplot,
)
from .sequence import (  # noqa
    Alphabet,
    alphabets,
    Sequence,
    SequenceReader,
    BatchSequenceReader,
    AbstractSequenceCollection,
    SequenceCollection,
    MultipleSequenceAlignment,
    SequenceAnnotation,
    SequenceInterval,
)

from .ontology import Ontology, OntologyTerm  # noqa
