import json
import re
import uuid
from abc import ABCMeta, abstractmethod
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from itertools import chain, groupby
from subprocess import PIPE, Popen
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Tuple, TypeVar, Union
from urllib.parse import unquote
from warnings import warn

import numpy as np
import numpy.typing as npt

from .dag import DAGElement, DirectedAcyclicGraph

SequenceType = TypeVar("SequenceType")  # Used for multiple dispatch


# Characters that must be percent-encoded in gff3 attribute values: percent, control characters (including tab,
# newline and carriage return), and the column 9 separators. No other characters may be encoded.
# https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
GFF3_ENCODE_PATTERN = re.compile(r"[%;=&,\x00-\x1f\x7f]")

# A single GTF attribute: key, followed by a quoted value (that can contain ';') or an unquoted value
GTF_ATTRIBUTE_PATTERN = re.compile(r'\s*([^\s;"]+)\s+("[^"]*"|[^;"]*?)\s*(?:;|$)')

# Translate dna codons to amino acids
TRANSLATION = dict(
    ACC="T",
    ACA="T",
    ACG="T",
    AGG="R",
    AGC="S",
    GTA="V",
    AGA="R",
    ACT="T",
    GTG="V",
    AGT="S",
    CCA="P",
    CCC="P",
    GGT="G",
    CGA="R",
    CGC="R",
    TAT="Y",
    CGG="R",
    CCT="P",
    GGG="G",
    GGA="G",
    GGC="G",
    TAA="*",
    TAC="Y",
    CGT="R",
    TAG="*",
    ATA="I",
    CTT="L",
    ATG="M",
    CTG="L",
    ATT="I",
    CTA="L",
    TTT="F",
    GAA="E",
    TTG="L",
    TTA="L",
    TTC="F",
    GTC="V",
    AAG="K",
    AAA="K",
    AAC="N",
    ATC="I",
    CAT="H",
    AAT="N",
    GTT="V",
    CAC="H",
    CAA="Q",
    CAG="Q",
    CCG="P",
    TCT="S",
    TGC="C",
    TGA="*",
    TGG="W",
    TCG="S",
    TCC="S",
    TCA="S",
    GAG="E",
    GAC="D",
    TGT="C",
    GCA="A",
    GCC="A",
    GCG="A",
    GCT="A",
    CTC="L",
    GAT="D",
)


@dataclass(frozen=True)
class Alphabet(set):
    """Alphabet of arbitrary biological sequences

    Examples:
        >>> DNA = Alphabet('DNA', 'ACGT')
        >>> DNA
        Alphabet(name='DNA', members='ACGT')

        >>> Protein = Alphabet('AminoAcid', '*-?ACDEFGHIKLMNPQRSTVWXY')
        >>> Protein
        Alphabet(name='AminoAcid', members='*-?ACDEFGHIKLMNPQRSTVWXY')


    Args:
        name (str): Alphabet name
        members (Iterable[str]): Letters of the alphabet
    """

    name: str
    members: Iterable[str]

    def __post_init__(self) -> None:
        super().__init__(self.members)

    def __deepcopy__(self, memo) -> "Alphabet":
        return Alphabet(self.name, self.members)

    def score(
        self,
        sequence: str,
        match: float = 1.0,
        mismatch: float = -1.0,
        n_chars: int = 100,
    ) -> float:
        """Scores how well a sequence matches an alphabet by summing \
        (mis)matches of sequence letters that are not in the alphabet \
        and (mis)matches of alphabet letters that are not in the sequence.

        Args:
            sequence (str): Sequence string for which to determine how well \
                it fits the alphabet
            match (float, optional): match score. Defaults to 1.0.
            mismatch (float, optional): mismatch score. Defaults to -1.0.
            n_chars (int, optional): number of sequence characters to use in \
                scoring. Large numbers incur a significant computational cost.

        Returns:
            (float): Score of how well a sequence matches the alphabet
        """
        return sum(match if s in self else mismatch for s in sequence[:n_chars]) + sum(
            match if s in sequence[:n_chars] else mismatch for s in self
        )

    def validate(self, sequence: str) -> bool:
        """Determine whether a sequence strictly fits an alphabet

        Args:
            sequence (str): Sequence string

        Returns:
            bool: true if all characters in sequence are in the alphabet
        """
        return sum(1 if s not in self else 0 for s in sequence) == 0

    def complement(self, sequence: str) -> str:
        """Returns complementary strand of DNA or RNA sequence strings

        Examples:
            >>> DNA = Alphabet('DNA', 'ACGT')
            >>> DNA.complement('AACTACG')
            'TTGATGC'

        Args:
            sequence (str): Sequence string

        Returns:
            str: complementary strand sequence string
        """
        if self.name == "DNA":
            complement = dict(zip("acgtnACGTN-?", "tgcanTGCAN-?", strict=True))
        elif self.name == "RNA":
            complement = dict(zip("acgunACGUN-?", "ugcanUGCAN-?", strict=True))
        else:
            raise TypeError("Cannot complement non-DNA or non-RNA alphabet")
        return "".join(complement[s] for s in sequence)

    def translate(self, sequence: str) -> str:
        """Translate DNA or RNA sequence string to amino acid string

        Examples:
            >>> DNA = Alphabet('DNA', 'ACGT')
            >>> DNA.translate('ATGACGACGTAA')
            'MTT*'

        Args:
            sequence (str): Sequence string (sequence length must be multiple of 3)

        Returns:
            str: Amino acid string
        """
        if self.name not in ("DNA", "RNA"):
            raise TypeError("Cannot translate non-DNA or non-RNA alphabet")
        codons = re.findall("...", sequence.upper())
        return "".join(TRANSLATION.get(codon, "X") for codon in codons)


def alphabet_factory(alphabet):
    """
    Factory function that returns a specific alphabet
    """
    return dict(
        DNA=lambda: Alphabet("DNA", "-?acgtnACGNT"),
        RNA=lambda: Alphabet("RNA", "-?acgunACGNU"),
        AminoAcid=lambda: Alphabet("AminoAcid", "*-?acdefghiklmnpqrstvwxyACDEFGHIKLMNPQRSTVWXY"),
    )[alphabet]


@dataclass(frozen=True)
class Alphabets:
    """
    Immutable container with commonly used biological sequence alphabets
    """

    DNA: Alphabet = field(default_factory=alphabet_factory("DNA"))
    # RNA: Alphabet = field(default_factory=alphabet_factory("RNA"))
    AminoAcid: Alphabet = field(default_factory=alphabet_factory("AminoAcid"))

    def __iter__(self):
        return iter(self.__dict__.values())


alphabets = Alphabets()


def guess_alphabet(sequence: str, n_chars: int = 100) -> Alphabet:
    """Guess the alphabet of a sequence string: the alphabet that contains the most sequence characters. Ties go to
    the smallest alphabet, so a sequence of only nucleotide characters is DNA.

    Examples:
        >>> guess_alphabet('ACGTTGCA').name
        'DNA'
        >>> guess_alphabet('MKVLAAGIVGLL').name
        'AminoAcid'

    Args:
        sequence (str): Sequence string
        n_chars (int, optional): Number of sequence characters to use. Defaults to 100.

    Returns:
        Alphabet: Best matching alphabet
    """
    chars = sequence[:n_chars]
    return max(sorted(alphabets, key=len), key=lambda alphabet: sum(char in alphabet for char in chars))


def quote_gff3(attribute_value: Union[int, str, float]) -> str:
    """Percent-encode characters that have a special meaning in GFF3 attribute values

    Args:
        attribute_value (Union[int, str, float]): Attribute value

    Returns:
        str: Encoded attribute value
    """
    return GFF3_ENCODE_PATTERN.sub(lambda match: f"%{ord(match.group()):02X}", str(attribute_value))


def encode_attribute_value(attribute_value: Iterable[Union[int, str, float]]) -> str:
    """Encode one or more GFF3 attribute values as a comma-separated string

    Args:
        attribute_value (Iterable[Union[int, str, float]]): Single value, or list/tuple of values

    Returns:
        str: Comma-separated, percent-encoded values
    """
    if not isinstance(attribute_value, (list, tuple)):
        attribute_value = [attribute_value]
    return ",".join(quote_gff3(v) for v in attribute_value)


def format_gtf_attribute_string(attributes: Dict[str, Iterable[Union[int, str, float]]]) -> str:
    """Format attributes as a GTF column 9 string (``key "value";``). Attributes with multiple values are written
    as repeated keys.

    Args:
        attributes (Dict[str, Iterable[Union[int, str, float]]]): Attribute names and values

    Returns:
        str: GTF attribute string
    """
    parts = []
    for key, value in attributes.items():
        values = value if isinstance(value, (list, tuple)) else [value]
        parts.extend(f'{key} "{v}";' for v in values)
    return " ".join(parts) if parts else "."


def format_gff_attribute_string(attributes: Dict[str, Iterable[Union[int, str, float]]]) -> str:
    """Format attributes as a GFF3 column 9 string (``key=value;``). Predefined GFF3 attributes are capitalized, and
    ``ID`` and ``Parent`` are written first.

    Args:
        attributes (Dict[str, Iterable[Union[int, str, float]]]): Attribute names and values, must contain ``ID``

    Returns:
        str: GFF3 attribute string
    """
    partially_formatted = {
        (key.capitalize() if key in SequenceInterval._predefined_gff3_attributes else key): encode_attribute_value(
            value
        )
        for key, value in attributes.items()
    }
    partially_formatted["ID"] = partially_formatted.pop("Id")

    def sort_key(key):
        if key == "ID":
            return 0
        elif key == "Parent":
            return 1
        else:
            return 2

    return ";".join(f"{key}={partially_formatted[key]}" for key in sorted(partially_formatted.keys(), key=sort_key))


def unquote_gff3(attribute_value: str) -> str:
    """Decode percent-encoded characters in a GFF3 attribute value

    Args:
        attribute_value (str): Encoded attribute value

    Returns:
        str: Decoded attribute value
    """
    return unquote(attribute_value)


def decode_attribute_value(attribute_value: str) -> List[str]:
    """Split a comma-separated GFF3 attribute value and decode every part

    Args:
        attribute_value (str): Encoded attribute value

    Returns:
        List[str]: Decoded values
    """
    return [unquote_gff3(v) for v in attribute_value.split(",")]


def parse_gtf_attribute_string(gtf_attribute_string: str) -> Dict[str, List[str]]:
    """Parse a GTF column 9 string (``key "value";``) into a dictionary. Values are lists: repeated keys are
    collected in the same list.

    Examples:
        >>> dict(parse_gtf_attribute_string('gene_id "g1"; tag "basic"; tag "a;b"; exon_number 2;'))
        {'gene_id': ['g1'], 'tag': ['basic', 'a;b'], 'exon_number': ['2']}

    Args:
        gtf_attribute_string (str): GTF attribute string

    Returns:
        Dict[str, List[str]]: Attribute names and values

    Raises:
        ValueError: If the attribute string can not be parsed
    """
    attributes = defaultdict(list)
    string = gtf_attribute_string.strip()
    if string == ".":
        return attributes
    position = 0
    while position < len(string):
        match = GTF_ATTRIBUTE_PATTERN.match(string, position)
        if not match:
            raise ValueError(f"Error parsing GTF attribute string at position {position}: {gtf_attribute_string}")
        key, value = match.groups()
        attributes[key].append(value[1:-1] if value.startswith('"') else value)
        position = match.end()
    return attributes


def parse_gff_attribute_string(
    gff_attribute_string: str, case_sensitive_attribute_keys: bool = False
) -> Dict[str, List[str]]:
    """Parse a GFF3 column 9 string (``key=value;``) into a dictionary.

    Keys are lowercased (except ``ID``), and keys that collide with the eight fixed GFF3 columns get an underscore
    prefix. Values are split on commas and percent-decoded, so every value is a list. See "Column 9: Attributes" in
    https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

    Args:
        gff_attribute_string (str): GFF3 attribute string
        case_sensitive_attribute_keys (bool): Currently unused

    Returns:
        Dict[str, List[str]]: Attribute names and values
    """
    attributes = defaultdict(list)
    if gff_attribute_string.strip() == ".":
        return attributes
    for string_part in gff_attribute_string.split(";"):
        if not string_part:
            continue
        try:
            key, value = string_part.split("=", maxsplit=1)
        except Exception as e:
            print(gff_attribute_string, string_part)
            raise Exception(
                f"{e}. Offending string part: {string_part}. " f"Offending attribute string: {gff_attribute_string}"
            ) from e
        # The gff spec lists the predefined attribute fields as starting with
        # a capital letter, but we process in lowercase so we don't miss
        # anything from poorly formatted files. When writing to gff we convert
        # back to a capital
        # EXCEPT FOR THE ID ATTRIBUTE, since lowercase id is reserved in python
        if key != "ID":
            key = key.lower()
        # First eight columns have predefined names which can collide with what
        # is in the 9th column. Solution is to prefix collisions in the 9th
        # column with underscore. E.g. 'score' becomes '_score' because 'score'
        # is the predefined name of the 6th column
        if key in SequenceInterval._fixed_gff3_fields:
            key = f"_{key}"
        for value_part in decode_attribute_value(value):
            attributes[key].append(value_part)
    return attributes


class SequenceAnnotation(DirectedAcyclicGraph):
    """Sequence annotation: a collection of :class:`SequenceInterval` objects (genes, mRNAs, exons, ...) stored by
    ID, and linked into a directed acyclic graph by their ``Parent`` attributes.

    Intervals with duplicate IDs are renamed (``ID_1``, ``ID_2``, ...) with a warning.

    Examples:
        >>> gff = (
        ...     "ctg1\\t.\\tgene\\t1\\t90\\t.\\t+\\t.\\tID=gene1\\n"
        ...     "ctg1\\t.\\tmRNA\\t1\\t90\\t.\\t+\\t.\\tID=mRNA1;Parent=gene1\\n"
        ... )
        >>> annotation = SequenceAnnotation.from_gff(string=gff)
        >>> len(annotation)
        2
        >>> annotation["mRNA1"].parent
        ['gene1']
        >>> [interval.ID for interval in annotation["gene1"].children]
        ['mRNA1']

    Args:
        sequence (Optional[Sequence]): Annotated sequence. If given, its ``annotation`` is set to this annotation.
    """

    def __init__(self, sequence: Optional["Sequence"] = None) -> None:
        super().__init__()
        if sequence:
            sequence.annotation = self
        self.sequence = sequence
        self._gff_headers = list()

    @property
    def intervals(self):
        """List of all intervals"""
        return list(self)

    def _link_parents(self) -> None:
        """
        Add explicit link from parent to child intervals
        GFF/GTF files only contain links of child to parent
        This modifies elements in place
        """
        for interval in self:
            if interval.parent:
                for parent_ID in interval.parent:
                    try:
                        parent = self[parent_ID]
                    except KeyError as exc:
                        raise KeyError(
                            f"Interval {interval.ID} is listing {parent_ID} "
                            "as Parent, but parent could not be found."
                        ) from exc
                    parent._children.append(interval.ID)

    @classmethod
    def from_gtf(
        cls,
        filename: Optional[str] = None,
        string: Optional[str] = None,
        sequence: Optional["Sequence"] = None,
        link_parents: Optional[bool] = True,
    ) -> "SequenceAnnotation":
        """Read a GTF formatted file or string. Exactly one of ``filename`` or ``string`` must be given.

        GTF lines are linked by their ``gene_id`` and ``transcript_id`` attributes, and ``transcript`` lines become
        ``mRNA`` intervals. Genes and transcripts that have no line of their own are created, spanning all their
        child intervals. Other intervals get IDs based on their transcript and type, e.g. ``<transcript_id>.exon_0``.
        Lines without a ``gene_id`` become intervals without parents, with IDs like ``repeat_region_0``.

        Examples:
            >>> gtf = (
            ...     'ctg1\\t.\\texon\\t100\\t200\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "t1";\\n'
            ...     'ctg1\\t.\\texon\\t300\\t400\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "t1";\\n'
            ... )
            >>> annotation = SequenceAnnotation.from_gtf(string=gtf)
            >>> gene = annotation["g1"]
            >>> gene.start, gene.end
            (100, 400)
            >>> [interval.ID for interval in gene.children]
            ['t1', 't1.exon_0', 't1.exon_1']

        Args:
            filename (Optional[str]): GTF filename
            string (Optional[str]): GTF formatted string
            sequence (Optional[Sequence]): Annotated sequence
            link_parents (Optional[bool]): Link parent intervals to their children, so that
                :attr:`~SequenceInterval.children` works

        Returns:
            SequenceAnnotation: Sequence annotation
        """
        assert filename or string
        assert not (filename and string)
        if filename:
            with open(filename) as filehandle:
                string = filehandle.read()
        sequence_annotation = cls(sequence=sequence)

        header = True
        intervals = []
        for line_number, line in enumerate(string.split("\n")):
            line = line.strip()
            if not line:
                continue
            if line[0] == "#":
                if header:
                    sequence_annotation._gff_headers.append(line)
                continue
            header = False
            intervals.append(SequenceInterval.from_gtf_line(gtf_line=line, line_number=line_number))

        def first_value(interval: "SequenceInterval", key: str) -> Optional[str]:
            values = interval.__dict__.get(key)
            return values[0] if values else None

        # GTF lines have no IDs: genes and transcripts are identified by their gene_id and transcript_id attributes,
        # other intervals get an ID based on their transcript (or gene) and type. Genes and transcripts that are not
        # in the file are created from their first child, and span all their children.
        gene_ids = {first_value(interval, "gene_id") for interval in intervals if interval.interval_type == "gene"}
        transcript_ids = {
            first_value(interval, "transcript_id")
            for interval in intervals
            if interval.interval_type in ("transcript", "mRNA")
        }
        new_genes: Dict[str, SequenceInterval] = dict()
        new_transcripts: Dict[str, SequenceInterval] = dict()
        child_counter = Counter()
        ordered_intervals = []
        for interval in intervals:
            gene_id = first_value(interval, "gene_id")
            transcript_id = first_value(interval, "transcript_id")
            if interval.interval_type == "gene" and gene_id:
                interval._ID = gene_id
            elif interval.interval_type in ("transcript", "mRNA") and transcript_id:
                interval._ID = transcript_id
                interval.interval_type = "mRNA"
                interval.parent = [gene_id] if gene_id else []
            elif not gene_id:
                # not part of a gene (GTF requires a gene_id, but GFF3 converted to GTF can lack it)
                child_count = child_counter[(None, interval.interval_type)]
                child_counter[(None, interval.interval_type)] += 1
                interval._ID = f"{interval.interval_type}_{child_count}"
                interval.parent = []
            else:
                if gene_id not in gene_ids and gene_id not in new_genes:
                    new_genes[gene_id] = SequenceInterval._spanning_interval(
                        interval, ID=gene_id, interval_type="gene", parent=None, gene_id=[gene_id]
                    )
                    ordered_intervals.append(new_genes[gene_id])
                if transcript_id and transcript_id not in transcript_ids and transcript_id not in new_transcripts:
                    new_transcripts[transcript_id] = SequenceInterval._spanning_interval(
                        interval,
                        ID=transcript_id,
                        interval_type="mRNA",
                        parent=[gene_id],
                        gene_id=[gene_id],
                        transcript_id=[transcript_id],
                    )
                    ordered_intervals.append(new_transcripts[transcript_id])
                parent_id = transcript_id or gene_id
                child_count = child_counter[(parent_id, interval.interval_type)]
                child_counter[(parent_id, interval.interval_type)] += 1
                interval._ID = f"{parent_id}.{interval.interval_type}_{child_count}"
                interval.parent = [parent_id]
            ordered_intervals.append(interval)

        for interval in intervals:
            for spanning_interval in (
                new_genes.get(first_value(interval, "gene_id")),
                new_transcripts.get(first_value(interval, "transcript_id")),
            ):
                if spanning_interval is not None:
                    spanning_interval.start = min(spanning_interval.start, interval.start)
                    spanning_interval.end = max(spanning_interval.end, interval.end)

        for interval in ordered_intervals:
            interval._container = sequence_annotation
            sequence_annotation[interval.ID] = interval

        if link_parents:
            sequence_annotation._link_parents()

        return sequence_annotation

    def to_gtf(self) -> str:
        """GTF formatted string with all intervals

        Returns:
            str: GTF formatted string
        """
        return "\n".join(interval.to_gtf_line() for interval in self)

    @classmethod
    def from_gff(
        cls,
        filename: Optional[str] = None,
        string: Optional[str] = None,
        sequence: Optional["Sequence"] = None,
        link_parents: bool = True,
    ) -> "SequenceAnnotation":
        """Read a GFF3 formatted file or string. Exactly one of ``filename`` or ``string`` must be given. Comment lines
        are skipped, and reading stops at a ``##FASTA`` line.

        Args:
            filename (Optional[str]): GFF3 filename
            string (Optional[str]): GFF3 formatted string
            sequence (Optional[Sequence]): Annotated sequence
            link_parents (bool): Link parent intervals to their children, so that
                :attr:`~SequenceInterval.children` works

        Returns:
            SequenceAnnotation: Sequence annotation
        """
        assert filename or string
        assert not (filename and string)
        sequence_annotation = cls(sequence=sequence)
        header = True
        if filename:
            with open(filename) as filehandle:
                string = filehandle.read()
        for line_number, line in enumerate(string.split("\n")):
            line = line.strip()
            if not line:
                continue
            if line == "##FASTA":
                break
            if line[0] == "#":
                if header:
                    sequence_annotation._gff_headers.append(line)
                continue
            else:
                header = False

            interval = SequenceInterval.from_gff_line(gff_line=line, line_number=line_number)
            interval._container = sequence_annotation
            sequence_annotation[interval.ID] = interval

        if link_parents:
            sequence_annotation._link_parents()

        return sequence_annotation

    def to_gff(self) -> str:
        """GFF3 formatted string with all intervals (without header lines)

        Returns:
            str: GFF3 formatted string
        """
        return "".join(interval.to_gff_line(trailing_newline=True) for interval in self)

    @classmethod
    def from_json(
        cls,
        filename: Optional[str] = None,
        string: Optional[str] = None,
        sequence: Optional["Sequence"] = None,
    ) -> "SequenceAnnotation":
        """Read a json formatted file or string, as written by :meth:`to_json`. Exactly one of ``filename`` or
        ``string`` must be given.

        The json must be a list of interval dictionaries (see :meth:`SequenceInterval.to_dict`). Each interval
        dictionary can have a ``children`` list with more interval dictionaries.

        Args:
            filename (Optional[str]): json filename
            string (Optional[str]): json formatted string
            sequence (Optional[Sequence]): Annotated sequence

        Returns:
            SequenceAnnotation: Sequence annotation
        """
        assert filename or string
        assert not (filename and string)
        if filename:
            with open(filename) as filehandle:
                string = filehandle.read()

        sequence_annotation = cls(sequence=sequence)

        gene_dicts = json.loads(string)
        assert isinstance(gene_dicts, list)

        for top_dict in gene_dicts:
            child_dicts = top_dict.pop("children", list())
            top_interval = SequenceInterval.from_dict(interval_dict=top_dict)
            top_interval._container = sequence_annotation
            sequence_annotation[top_interval.ID] = top_interval
            for child_dict in child_dicts:
                child_interval = SequenceInterval.from_dict(interval_dict=child_dict)
                child_interval._container = sequence_annotation
                sequence_annotation[child_interval.ID] = child_interval
        sequence_annotation._link_parents()
        return sequence_annotation

    def to_json(self, indent: Optional[int] = None) -> str:
        """json formatted string with a list of interval dictionaries (see :meth:`SequenceInterval.to_dict`)

        Args:
            indent (Optional[int]): Indentation, passed to :func:`json.dumps`

        Returns:
            str: json formatted string
        """
        interval_dicts = [interval.to_dict() for interval in self]
        return json.dumps(interval_dicts, indent=indent)


class SequenceInterval(DAGElement):
    """Single annotated interval on a sequence, such as a gene, mRNA, exon or CDS. Corresponds to one line in a
    GFF3 or GTF file.

    The eight fixed GFF3 columns are stored as attributes (``seqid``, ``source``, ``interval_type``, ``start``,
    ``end``, ``score``, ``strand``, ``phase``). Column 9 attributes are stored as additional attributes with
    lowercase keys (except ``ID``) and list values, and can also be accessed by key (``interval["name"]``).
    Predefined GFF3 attributes that are not set are ``None``.

    Examples:
        >>> interval = SequenceInterval.from_gff_line("ctg1\\t.\\tgene\\t1000\\t9000\\t.\\t+\\t.\\tID=gene1;Name=EDEN")
        >>> interval.interval_type, interval.start, interval.end, interval.strand
        ('gene', 1000, 9000, '+')
        >>> interval.name
        ['EDEN']

    Args:
        ID (Optional[str]): Unique identifier
        seqid (Optional[str]): Name of the sequence the interval is on (e.g. a chromosome)
        source (Optional[str]): Program or database that produced the interval
        interval_type (Optional[str]): Feature type, e.g. ``gene``, ``mRNA``, ``exon`` or ``CDS``
        start (Optional[int]): Start position (as in GFF3: 1-based, inclusive)
        end (Optional[int]): End position (as in GFF3: 1-based, inclusive)
        score (Optional[float]): Score
        strand (Optional[str]): ``+``, ``-`` or ``.``
        phase (Optional[str]): Phase of CDS intervals: ``0``, ``1``, ``2`` or ``.``
        children (Optional[List[str]]): IDs of child intervals
        container (Optional[SequenceAnnotation]): Annotation the interval belongs to
        **kwargs: Additional attributes. ``parent`` is the list of parent interval IDs.
    """

    _predefined_gff3_attributes = (
        "ID",
        "name",
        "alias",
        "parent",
        "target",
        "gap",
        "derives_from",
        "note",
        "dbxref",
        "ontology_term",
        "is_circular",
    )
    _fixed_gff3_fields = (
        "seqid",
        "source",
        "interval_type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
    )
    _gtf_interval_types = dict(mRNA="transcript")

    def __init__(
        self,
        ID: Optional[str] = None,
        seqid: Optional[str] = None,
        source: Optional[str] = None,
        interval_type: Optional[str] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
        score: Optional[float] = None,
        strand: Optional[str] = None,
        phase: Optional[str] = None,
        children: Optional[List[str]] = None,
        container: Optional[SequenceAnnotation] = None,
        **kwargs,
    ):
        parents = kwargs.pop("parent", None)
        super().__init__(ID=ID, children=children, container=container, parents=parents)

        # Standard gff fields
        self.seqid = seqid
        self.source = source
        self.interval_type = interval_type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase

        # Set attributes with predefined meanings in the gff spec to None
        for attr in self._predefined_gff3_attributes:
            # ID and parent are handled separately in DAG superclass
            if attr in {"ID", "parent"}:
                continue
            self[attr] = kwargs.get(attr, None)

        # Any additional attributes
        for key, value in kwargs.items():
            self[key] = value

    def __repr__(self):
        return (
            f"<SequenceInterval type={self.interval_type} "
            f"ID={self.ID} "
            f"loc={self.seqid}..{self.start}..{self.end}..{self.strand} "
            f"at {hex(id(self))}>"
        )

    @property
    def parent(self):
        """IDs of the direct parent intervals (the GFF3 ``Parent`` attribute). Can be set with a single ID or a list of
        IDs. See :attr:`parents` for all ancestors.
        """
        return self._parents

    @parent.setter
    def parent(self, parent_ID: Union[List[str], str]):
        if isinstance(parent_ID, str):
            parent_ID = [parent_ID]
        self._parents = parent_ID

    @property
    def gff_attributes(self) -> Dict[str, str]:
        """Column 9 attributes as a dictionary, including ``ID`` and ``Parent``. Attributes that are ``None`` are left
        out.
        """
        gff_attributes = {
            attr: self[attr]  # dictionary comprehension
            for attr in self.__dict__
            if attr not in self._fixed_gff3_fields  # skip column 1-8 in gff3
            and attr
            not in (
                "_parents",
                "_children",
                "_container",
                "_ID",
                "_original_ID",
            )  # internal use only
            and self[attr] is not None  # no empty attributes
        }

        # Add attributes handled by DAG
        gff_attributes["ID"] = [self.ID]
        if self._parents:
            gff_attributes["Parent"] = self._parents

        return gff_attributes

    @property
    def gtf_attributes(self) -> Dict[str, str]:
        """Attributes for GTF output: :attr:`gff_attributes` without ``ID`` and ``Parent``, plus a ``<type>_id``
        attribute for the interval itself if it is a gene or mRNA, and for every ancestor (e.g. ``transcript_id`` and
        ``gene_id`` for an exon). Intervals with multiple parents of the same type get a single ``<type>_id``.
        """
        def get_gtf_type(gff_interval_type):
            return self._gtf_interval_types.get(gff_interval_type, gff_interval_type)

        if self.parents:
            parent_ids = {f"{get_gtf_type(parent.interval_type)}_id": parent.ID for parent in self.parents}
        else:
            parent_ids = dict()

        attributes = {key: value for key, value in self.gff_attributes.items() if key not in ("ID", "Parent")}
        if self.interval_type == "gene":
            attributes["gene_id"] = self.ID
        elif self.interval_type == "mRNA":
            attributes["transcript_id"] = self.ID
        return {**attributes, **parent_ids}

    @classmethod
    def _spanning_interval(
        cls, child: "SequenceInterval", ID: str, interval_type: str, parent: Optional[List[str]], **attributes
    ) -> "SequenceInterval":
        """Create an interval (e.g. a gene) on the same sequence and strand as ``child``, starting with the span of
        ``child``"""
        return cls(
            ID=ID,
            seqid=child.seqid,
            source=child.source,
            interval_type=interval_type,
            start=child.start,
            end=child.end,
            score=".",
            strand=child.strand,
            phase=".",
            parent=parent,
            **attributes,
        )

    @classmethod
    def from_gtf_line(cls, gtf_line: Optional[str] = None, line_number: Optional[int] = None) -> "SequenceInterval":
        """Create an interval from a single GTF line

        Args:
            gtf_line (Optional[str]): GTF formatted line
            line_number (Optional[int]): Line number, used in error messages

        Returns:
            SequenceInterval: Interval

        Raises:
            ValueError: If the start, end, score, strand, or phase column has an invalid value
        """
        return cls.from_gff_line(gtf_line, line_number, parse_gtf_attribute_string)

    def to_gtf_line(self) -> str:
        """GTF formatted line (without a trailing newline). The ``mRNA`` type is written as ``transcript``.

        Returns:
            str: GTF formatted line
        """
        interval_type = self._gtf_interval_types.get(self.interval_type, self.interval_type)
        return "\t".join(
            [
                self.seqid,
                self.source,
                interval_type,
                str(self.start),
                str(self.end),
                str(self.score),
                self.strand,
                str(self.phase),
                format_gtf_attribute_string(self.gtf_attributes),
            ]
        )

    @classmethod
    def from_gff_line(
        cls,
        gff_line: Optional[str] = None,
        line_number: Optional[int] = None,
        attribute_parser: Callable = parse_gff_attribute_string,
    ) -> "SequenceInterval":
        """Create an interval from a single GFF3 line. Intervals without an ``ID`` attribute get a random UUID.

        Args:
            gff_line (Optional[str]): GFF3 formatted line
            line_number (Optional[int]): Line number, used in error messages
            attribute_parser (Callable): Function that parses column 9 into a dictionary of lists

        Returns:
            SequenceInterval: Interval

        Raises:
            ValueError: If the line does not have 9 columns, or if the start, end, score, strand, or phase column has
                an invalid value
        """
        gff_parts = gff_line.split("\t")
        if len(gff_parts) != 9:
            error = f"GFF and GTF lines must have 9 tab-separated columns, found {len(gff_parts)}"
            if line_number:
                error = f"{error}, line {line_number}"
            raise ValueError(error)
        seqid, source, interval_type, start, end, score, strand, phase = gff_parts[:8]
        try:
            start = int(start)
            end = int(end)
        except ValueError as err:
            error = "GFF start and end fields must be integer"
            if line_number:
                error = f"{error}, gff line {line_number}"
            raise ValueError(error) from err

        if score != ".":
            try:
                score = float(score)
            except ValueError as err:
                error = "GFF score field must be a float"
                if line_number:
                    error = f"{error}, gff line {line_number}"
                raise ValueError(error) from err

        if strand not in ("+", "-", "."):
            error = 'GFF strand must be one of "+", "-" or "."'
            if line_number:
                error = f"{error}, gff line {line_number}"
            raise ValueError(error)

        if phase not in ("0", "1", "2", "."):
            error = 'GFF phase must be one of "0", "1", "2" or "."'
            if line_number:
                error = f"{error}, gff line {line_number}"
            raise ValueError(error)
        elif phase != ".":
            phase = int(phase)

        # Disable phase checking of CDS for now...
        # if interval_type == 'CDS' and phase not in ('0', '1', '2'):
        #     error = 'GFF intervals of type CDS must have phase of\
        #         "0", "1" or "2"'
        #     if line_number:
        #         error = f'{error}, gff line {line_number}'
        #         raise ValueError(error)

        attributes = attribute_parser(gff_parts[8])

        ID = attributes.pop("ID", [str(uuid.uuid4())])[0]

        return cls(
            seqid=seqid,
            source=source,
            interval_type=interval_type,
            start=start,
            end=end,
            score=score,
            strand=strand,
            phase=phase,
            ID=ID,
            **attributes,
        )

    def to_gff_line(self, trailing_newline: bool = False) -> str:
        """GFF3 formatted line

        Args:
            trailing_newline (bool): End the line with a newline

        Returns:
            str: GFF3 formatted line
        """
        # attributes = dict(ID=self.ID, **self.gff_attributes)

        gff_line = "\t".join(
            [
                self.seqid,
                self.source,
                self.interval_type,
                str(self.start),
                str(self.end),
                str(self.score),
                self.strand,
                str(self.phase),
                format_gff_attribute_string(self.gff_attributes),
            ]
        )
        if trailing_newline:
            gff_line = f"{gff_line}\n"
        return gff_line

    @classmethod
    def from_dict(cls, interval_dict: Dict[str, Any]) -> "SequenceInterval":
        """Create an interval from a dictionary, as created by :meth:`to_dict`

        Args:
            interval_dict (Dict[str, Any]): The eight fixed GFF3 fields, ``ID``, and an ``attributes`` dictionary

        Returns:
            SequenceInterval: Interval
        """
        interval_dict = dict(interval_dict)
        attributes = dict(interval_dict.pop("attributes", dict()))
        parent = attributes.pop("Parent", attributes.pop("parent", None))
        if isinstance(parent, str):
            parent = [parent]
        return cls(**interval_dict, **attributes, parent=parent)

    def to_dict(self, include_children: bool = False) -> Dict[str, Any]:
        """Dictionary with the eight fixed GFF3 fields, ``ID``, and an ``attributes`` dictionary with all other
        attributes

        Args:
            include_children (bool): Add a ``children`` list with the dictionaries of all descendant intervals

        Returns:
            Dict[str, Any]: Interval dictionary
        """
        attributes = dict(**self.gff_attributes)
        attributes.pop("ID")
        interval_dict = dict(
            ID=self.ID,
            seqid=self.seqid,
            source=self.source,
            interval_type=self.interval_type,
            start=self.start,
            end=self.end,
            score=self.score,
            strand=self.strand,
            phase=self.phase,
            attributes=attributes,
        )
        if include_children:
            children = [child.to_dict() for child in self.children]
            interval_dict["children"] = children
        return interval_dict

    def to_json(self, include_children: bool = False, indent: Optional[int] = None) -> str:
        """json formatted string of :meth:`to_dict`

        Args:
            include_children (bool): Include all descendant intervals
            indent (Optional[int]): Indentation, passed to :func:`json.dumps`

        Returns:
            str: json formatted string
        """
        return json.dumps(self.to_dict(include_children=include_children), indent=indent)


@dataclass
class Sequence:
    """Single biological sequence with a header

    Examples:
        >>> dna = Sequence('test_dna', 'ACGATCGACTAGCA')
        >>> dna.alphabet.name
        'DNA'
        >>> protein = Sequence('test_aa', 'QAPISAIWPOIWQ*')
        >>> protein.alphabet.name
        'AminoAcid'
        >>> dna.reverse_complement.sequence
        'TGCTAGTCGATCGT'

    Args:
        header (str): Sequence name
        sequence (str): Sequence string
        alphabet (Alphabet): Sequence alphabet. Detected from the sequence when not given (see
            :func:`guess_alphabet`). Sequences derived from this one (slices, complements, etc.) keep the alphabet.
        annotation (Optional[SequenceAnnotation]): Annotation of the sequence. Defaults to an empty annotation.
    """

    header: str = None
    sequence: str = field(repr=False, default=None)
    alphabet: Alphabet = None
    annotation: Optional[SequenceAnnotation] = field(default_factory=SequenceAnnotation, repr=False)

    def __post_init__(self):
        if self.alphabet is not None:
            return
        if self.sequence is None:
            self.alphabet = alphabets.DNA
        else:
            self.alphabet = guess_alphabet(self.sequence)

    def __getitem__(self, key) -> "Sequence":
        """Subset a sequence based on a key (can be int or slice)

        Examples:
            >>> s = Sequence('test_dna', 'ACGTA')
            >>> s[2:]
            Sequence(header='test_dna', \
alphabet=Alphabet(name='DNA', members='-?acgtnACGNT'))
            >>> len(s[2:])
            3
        """
        return Sequence(self.header, self.sequence[key], alphabet=self.alphabet)

    def __len__(self) -> int:
        """Length of the sequence

        Examples:
            >>> s = Sequence('test_dna', 'ACGTA')
            >>> len(s)
            5
        """
        return len(self.sequence)

    @property
    def uppercase(self) -> "Sequence":
        """All sequence characters in uppercase

        Examples:
            >>> s = Sequence('test_dna', 'acgTA')
            >>> s.uppercase.sequence
            'ACGTA'
        """
        return Sequence(self.header, self.sequence.upper(), alphabet=self.alphabet)

    @property
    def lowercase(self) -> "Sequence":
        """All sequence characters in lowercase

        Examples:
            >>> s = Sequence('test_dna', 'acgTA')
            >>> s.lowercase.sequence
            'acgta'
        """
        return Sequence(self.header, self.sequence.lower(), alphabet=self.alphabet)

    @property
    def reverse(self) -> "Sequence":
        """Reverse sequence order

        Examples:
            >>> s = Sequence('test_dna', 'ACGTA')
            >>> s.reverse.sequence
            'ATGCA'
        """
        return Sequence(self.header, self.sequence[::-1], alphabet=self.alphabet)

    @property
    def complement(self) -> "Sequence":
        """Complement DNA sequences based on watson-crick pairing

        Examples:
            >>> s = Sequence('test_dna', 'ACGTA')
            >>> s.complement.sequence
            'TGCAT'
        """
        return Sequence(self.header, self.alphabet.complement(self.sequence), alphabet=self.alphabet)

    @property
    def reverse_complement(self) -> "Sequence":
        """Reverse sequence order and complement nucleotides vased on watson-crick pairing.
        This is the same as accessing the reversed and then complemented sequence (in arbitrary order)

        Examples:
            >>> s = Sequence('test_dna', 'ACGTA')
            >>> s.reverse_complement.sequence
            'TACGT'
            >>> s.reverse_complement.sequence == s.reverse.complement.sequence == s.complement.reverse.sequence
            True
        """
        return Sequence(self.header, self.alphabet.complement(self.sequence[::-1]), alphabet=self.alphabet)

    @property
    def amino_acids(self) -> "Sequence":
        """Translate nucleotide codon triplets into amino acids

        Examples:
            >>> s = Sequence('test_dna', 'ATGATGTAA')
            >>> s.amino_acids.sequence
            'MM*'
        """
        if self.alphabet.name == "AminoAcid":
            return self
        else:
            return Sequence(self.header, self.alphabet.translate(self.sequence), alphabet=alphabets.AminoAcid)

    def to_dict(self) -> Dict[str, str]:
        """Make dictionary with header and sequence elements

        Examples:
            >>> s = Sequence('test', 'ACGTA')
            >>> s.to_dict()
            {'header': 'test', 'sequence': 'ACGTA'}

        Returns:
            Dict[str, str]: sequence dictionary
        """
        return dict(header=self.header, sequence=self.sequence)

    @classmethod
    def from_fasta(cls, string: str) -> "Sequence":
        """Create a sequence from a fasta formatted string (single sequence only)

        Examples:
            >>> s = Sequence.from_fasta('>test\\nACGT')
            >>> s.header, s.sequence
            ('test', 'ACGT')

        Args:
            string (str): fasta formatted string

        Returns:
            Sequence: Sequence
        """
        lines = string.strip().split("\n")
        header = lines[0][1:]
        sequence = "".join(lines[1:])
        return cls(header, sequence)

    def to_fasta(self, linewidth: int = 80) -> str:
        """fasta formatted string

        Examples:
            >>> s = Sequence('test_dna', 'ACGTA')
            >>> s.to_fasta()
            '>test_dna\\nACGTA'

        Args:
            linewidth (int): Maximum number of sequence characters per line

        Returns:
            str: fasta formatted string
        """
        sequence_lines = "\n".join(re.findall(f".{{1,{linewidth}}}", self.sequence))
        return f">{self.header}\n{sequence_lines}"


class FastaParseError(Exception):
    """Raised when a string can not be parsed as fasta"""

    pass


class SequenceReader:
    """Iterator over the sequences in a fasta or json formatted file or string. Exactly one of ``string`` or
    ``filename`` must be given. json input is a list of objects with ``header`` and ``sequence`` keys.

    Examples:
        >>> fasta_string = '>1\\nACGC\\n>2\\nTGTGTA\\n'
        >>> [seq.header for seq in SequenceReader(string=fasta_string, filetype='fasta')]
        ['1', '2']

    Args:
        string (str): fasta or json formatted string
        filename (str): fasta or json filename
        filetype (str): ``"fasta"`` or ``"json"``

    Raises:
        ValueError: If ``filetype`` is not supported
    """

    def __init__(self, string: str = None, filename: str = None, filetype: str = None) -> None:
        assert bool(string) ^ bool(filename), "Must specify exactly one of string or filename"  # exclusive OR
        if filename:
            with open(filename, "r") as filehandle:
                string = filehandle.read().strip()

        self.string = string
        self._iterator = None

        if filetype == "fasta":
            self._iter = self._fasta_iter
        elif filetype == "json":
            self._iter = self._json_iter
        else:
            raise ValueError(f'filetype "{filetype}" is not supported')

    def __iter__(self) -> Iterable[Sequence]:
        """Iterate over all sequences

        Yields:
            Sequence: Next sequence
        """
        yield from self._iter()

    def __next__(self) -> Sequence:
        """Next sequence

        Returns:
            Sequence: Next sequence
        """
        if self._iterator is None:
            self._iterator = self._iter()
        return next(self._iterator)

    def _fasta_iter(self) -> Iterable[Sequence]:
        if self.string[0] != ">":
            raise FastaParseError(
                'First character in fasta format\
                 must be ">"'
            )
        fasta_iter = (x for _, x in groupby(self.string.strip().split("\n"), lambda line: line[:1] == ">"))
        for header in fasta_iter:
            header = next(header)[1:].strip()
            seq = "".join(s.strip() for s in next(fasta_iter))
            yield Sequence(header, seq)

    def _json_iter(self) -> Iterable[Sequence]:
        for entry in json.loads(self.string):
            yield Sequence(entry["header"], entry["sequence"])


class BatchSequenceReader(SequenceReader):
    """Iterator over batches of sequences in a fasta or json formatted file or string. Every batch is a
    :class:`SequenceCollection` of ``batchsize`` sequences, except for the last batch, which can be smaller.

    Examples:
        >>> fasta_string = '>1\\nACGC\\n>2\\nTGTGTA\\n>3\\nAAGT\\n>4\\nCCA\\n>5\\nGGA\\n'
        >>> reader = BatchSequenceReader(string=fasta_string, filetype='fasta', batchsize=2)
        >>> [batch.headers for batch in reader]
        [['1', '2'], ['3', '4'], ['5']]

    Args:
        string (str): fasta or json formatted string
        filename (str): fasta or json filename
        filetype (str): ``"fasta"`` or ``"json"``
        batchsize (int): Number of sequences per batch
    """

    def __init__(
        self,
        string: str = None,
        filename: str = None,
        filetype: str = None,
        batchsize: int = 10,
    ) -> None:
        super().__init__(string, filename, filetype)
        self.batchsize = batchsize

    def __iter__(self) -> Iterable["SequenceCollection"]:
        """Iterate over all batches. The last batch can have fewer than ``batchsize`` sequences.

        Yields:
            SequenceCollection: Next batch
        """
        batch = SequenceCollection()
        for seq in self._iter():
            batch[seq.header] = seq.sequence
            if len(batch) == self.batchsize:
                yield batch
                batch = SequenceCollection()
        if len(batch):
            yield batch

    def __next__(self) -> "SequenceCollection":
        """Next batch

        Returns:
            SequenceCollection: Next batch
        """
        if self._iterator is None:
            self._iterator = iter(self)
        return next(self._iterator)


SequenceIndexKey = Union[int, List[int], slice]


class SequenceIndex:
    """Position based index of a sequence collection, see :attr:`AbstractSequenceCollection.iloc`"""

    def __init__(
        self,
        sequence_collection: Union["SequenceCollection", "MultipleSequenceAlignment"],
    ):
        self.sequence_collection = sequence_collection

    def __getitem__(self, key: SequenceIndexKey):
        if isinstance(key, int):
            key = [key]
        elif isinstance(key, slice):
            key = range(*key.indices(len(self.sequence_collection)))
        elif not isinstance(key, list):
            raise TypeError(f"SequenceIndex key must be of type f{SequenceIndexKey}")
        new_seq_col = self.sequence_collection.__class__()
        for k in key:
            header = self.sequence_collection.headers[k]
            new_seq_col[header] = self.sequence_collection[header].sequence
        return new_seq_col


class AbstractSequenceCollection(metaclass=ABCMeta):
    """(Partially) abstract base class for sequence collections.

    Subclasses implement storage: ``__setitem__``, ``__getitem__``, ``__delitem__``, :meth:`pop`, :attr:`headers`
    and :attr:`n_seqs`. All other methods (reading and writing fasta and json, iteration, indexing with
    :attr:`iloc`, renaming) build on these.

    Sequences are stored by header, and indexing with a header returns a :class:`Sequence`. Setting a header that
    already exists does not overwrite the existing sequence: the new sequence gets a unique header (``header_1``,
    ``header_2``, ...) and a warning is issued.
    """

    @abstractmethod
    def __init__(
        self,
        sequences: Optional[Iterable[Sequence]] = None,
        sequence_annotation: Optional["SequenceAnnotation"] = None,
    ) -> None:
        raise NotImplementedError(
            ("Classes extending from AbstractSequenceCollection should " "implement __init__ method")
        )

    @abstractmethod
    def __setitem__(self, header: str, seq: str) -> None:
        raise NotImplementedError(
            ("Classes extending from AbstractSequenceCollection should " "implement __setitem__ method")
        )

    @abstractmethod
    def __getitem__(self, header: str) -> Sequence:
        raise NotImplementedError(
            ("Classes extending from AbstractSequenceCollection should " "implement __getitem__ method")
        )

    @abstractmethod
    def __delitem__(self, header: str) -> None:
        raise NotImplementedError(
            ("Classes extending from AbstractSequenceCollection should " "implement __delitem__ method")
        )

    def __iter__(self) -> Iterable[Sequence]:
        for header in self.headers:
            yield self[header]

    def __len__(self) -> int:
        return len(self.headers)

    def __add__(self: SequenceType, other: SequenceType) -> SequenceType:
        new_collection = self.__class__()
        for seq in chain(self, other):
            new_collection[seq.header] = seq.sequence
        return new_collection

    @property
    @abstractmethod
    def headers(self) -> List[str]:
        """Sequence headers, in insertion order

        Returns:
            List[str]: Sequence headers
        """
        raise NotImplementedError(
            ("Classes extending from AbstractSequenceCollection should " "implement headers property")
        )

    @property
    def iloc(self) -> SequenceIndex:
        """Position based indexing: index with an int, a list of ints, or a slice to get a new collection of the same
        type with the selected sequences

        Examples:
            >>> seqs = SequenceCollection.from_fasta(string='>a\\nACGT\\n>b\\nGGCC\\n>c\\nTTAA')
            >>> seqs.iloc[1:].headers
            ['b', 'c']
            >>> seqs.iloc[[0, 2]].headers
            ['a', 'c']
        """
        return SequenceIndex(self)

    @property
    def sequences(self) -> List[str]:
        """List of sequences without headers

        Returns:
            List[str]: list of sequences
        """
        return [self[header].sequence for header in self.headers]

    @property
    @abstractmethod
    def n_seqs(self) -> int:
        """Number of sequences in the collection

        Returns:
            int: Number of sequences
        """
        raise NotImplementedError(
            ("Classes extending from AbstractSequenceCollection should " "implement n_seqs property")
        )

    @classmethod
    def from_sequence_iter(cls, sequence_iter: Iterable[Sequence]) -> "SequenceCollection":
        """Create a collection from :class:`Sequence` objects

        Args:
            sequence_iter (Iterable[Sequence]): Sequences

        Returns:
            New collection of the class this method is called on
        """
        sequencecollection = cls()
        for seq in sequence_iter:
            sequencecollection[seq.header] = seq.sequence
        return sequencecollection

    @classmethod
    def from_fasta(
        cls,
        filename: str = None,
        string: str = None,
    ) -> "SequenceCollection":
        """Read a fasta formatted file or string. Exactly one of ``filename`` or ``string`` must be given.

        Examples:
            >>> seqs = SequenceCollection.from_fasta(string='>a\\nACGT\\n>b\\nGGCC')
            >>> seqs.headers
            ['a', 'b']

        Args:
            filename (str): fasta filename
            string (str): fasta formatted string

        Returns:
            New collection of the class this method is called on
        """
        sequencecollection = cls()

        for seq in SequenceReader(string=string, filename=filename, filetype="fasta"):
            sequencecollection[seq.header] = seq.sequence
        return sequencecollection

    def to_fasta(self, linewidth: int = 80) -> str:
        """fasta formatted string of all sequences

        Args:
            linewidth (int): Maximum number of sequence characters per line

        Returns:
            str: Multi-line fasta formatted string
        """
        return "\n".join([seq.to_fasta(linewidth=linewidth) for seq in self])

    @classmethod
    def from_json(cls, filename: Optional[str] = None, string: Optional[str] = None) -> "SequenceCollection":
        """Read a json formatted file or string, as written by :meth:`to_json`. Exactly one of ``filename`` or
        ``string`` must be given.

        Args:
            filename (Optional[str]): json filename
            string (Optional[str]): json formatted string: a list of objects with ``header`` and ``sequence`` keys

        Returns:
            New collection of the class this method is called on
        """
        sequencecollection = cls()

        for seq in SequenceReader(string=string, filename=filename, filetype="json"):
            sequencecollection[seq.header] = seq.sequence

        return sequencecollection

    def to_json(self, indent: Optional[int] = None) -> str:
        """json formatted string: a list of objects with ``header`` and ``sequence`` keys

        Args:
            indent (Optional[int]): Indentation, passed to :func:`json.dumps`

        Returns:
            str: json formatted string
        """
        gene_dicts = [seq.to_dict() for seq in self]
        return json.dumps(gene_dicts, indent=indent)

    @abstractmethod
    def pop(self, header: str) -> Sequence:
        """Remove a sequence from the collection and return it

        Args:
            header (str): Header of the sequence to remove

        Returns:
            Sequence: The removed sequence

        Raises:
            KeyError: If there is no sequence with this header
        """
        raise NotImplementedError("Classes extending from AbstractSequenceCollection should implement pop method")

    def add(self, seq: Sequence) -> None:
        """Add a sequence to the collection (in place)

        Args:
            seq (Sequence): Sequence to add
        """
        self[seq.header] = seq.sequence

    def modify_inplace(self, mod_func: Callable[[str, str], tuple[str, str]]) -> None:
        """Change headers and/or sequences in place by calling ``mod_func`` on every sequence

        Examples:
            >>> seqs = SequenceCollection.from_fasta(string='>a\\nacgt\\n>b\\nggcc')
            >>> seqs.modify_inplace(lambda header, sequence: (header.upper(), sequence.upper()))
            >>> seqs.to_fasta()
            '>A\\nACGT\\n>B\\nGGCC'

        Args:
            mod_func (Callable[[str, str], tuple[str, str]]): Function that takes a header and a sequence string, and
                returns a new header and sequence string
        """
        for header in self.headers:
            s: Sequence = self.pop(header)
            new_header, new_sequence = mod_func(s.header, s.sequence)
            self[new_header] = new_sequence

    def rename_inplace(self, rename_func: Callable[[str], str]) -> None:
        """Rename all headers in place by calling ``rename_func`` on every header

        Args:
            rename_func (Callable[[str], str]): Function that takes a header and returns a new header
        """
        self.modify_inplace(lambda header, sequence: (rename_func(header), sequence))


class SequenceCollection(AbstractSequenceCollection):
    """Collection of (unaligned) DNA or amino acid sequences

    Examples:
        >>> seqs = SequenceCollection([('a', 'ACGT'), ('b', 'GGCCAA')])
        >>> len(seqs)
        2
        >>> seqs['b'].sequence
        'GGCCAA'

    Args:
        sequences (Iterable[Tuple[str, str]]): (header, sequence) tuples
        sequence_annotation (SequenceAnnotation): Annotation of the sequences
    """

    def __init__(
        self: "SequenceCollection",
        sequences: Iterable[Tuple[str, str]] = None,
        sequence_annotation: "SequenceAnnotation" = None,
    ):
        self._collection = dict()
        if sequences:
            for header, sequence in sequences:
                self[header] = sequence
        self.sequence_annotation = sequence_annotation

    def __setitem__(self, header: str, seq: str) -> None:
        if header in self._collection:
            warn(f'Turning duplicate header "{header}" into unique header')
            new_header = header
            modifier = 0
            while new_header in self.headers:
                modifier += 1
                new_header = f"{header}_{modifier}"
            header = new_header
        self._collection[header] = seq

    def __getitem__(self, header: str) -> Sequence:
        sequence = self._collection[header]
        return Sequence(header, sequence)

    def __delitem__(self, header: str) -> None:
        del self._collection[header]

    @property
    def headers(self) -> List[str]:
        return list(self._collection.keys())

    @property
    def n_seqs(self) -> int:
        return len(self._collection.keys())

    def align(
        self, method: Optional[str] = "mafft", method_kwargs: Optional[Mapping[str, str]] = None
    ) -> "MultipleSequenceAlignment":
        """Align the sequences with an external multiple sequence aligner. The aligner is called as
        ``<method> <method_kwargs> -``, and must read fasta from stdin and write aligned fasta to stdout (like
        `MAFFT <https://mafft.cbrc.jp/alignment/software/>`_).

        Args:
            method (Optional[str]): Name of the aligner executable
            method_kwargs (Optional[Mapping[str, str]]): Command line options for the aligner, e.g.
                ``{"--thread": "4"}``

        Returns:
            MultipleSequenceAlignment: Aligned sequences

        Raises:
            RuntimeError: If the aligner exits with an error
        """
        if not method_kwargs:
            method_kwargs = dict()
        fasta = self.to_fasta()
        command = [method, *chain(*method_kwargs.items()), "-"]
        process = Popen(command, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate(input=fasta.encode())
        if process.returncode != 0:
            raise RuntimeError(f"{method} exited with code {process.returncode}: {stderr.decode().strip()}")
        aligned_fasta = stdout.decode().strip()
        return MultipleSequenceAlignment.from_fasta(string=aligned_fasta)

    def pop(self, header: str) -> Sequence:
        sequence = self._collection.pop(header)
        return Sequence(header, sequence)


class MultipleSequenceAlignment(SequenceCollection):
    """Collection of aligned DNA or amino acid sequences, stored as a numpy matrix with one row per sequence.
    Sequences shorter than the alignment are padded with gaps (``-``) at the end.

    Examples:
        >>> msa = MultipleSequenceAlignment.from_fasta(string='>a\\nAC-GT\\n>b\\nACCGT')
        >>> msa.shape
        (2, 5)

    Args:
        sequences (Optional[Iterable[Sequence]]): Sequences
        sequence_annotation (Optional[SequenceAnnotation]): Annotation of the sequences
    """

    def __init__(
        self,
        sequences: Optional[Iterable[Sequence]] = None,
        sequence_annotation: Optional["SequenceAnnotation"] = None,
    ) -> None:
        self._collection = np.empty((0, 0), dtype="uint8")
        self._header_idx = dict()
        if sequences:
            for seq in sequences:
                self[seq.header] = seq.sequence
        # if sequence_annotation:
        #     sequence_annotation.sequence_collection = self
        self.sequence_annotation = sequence_annotation

    def __setitem__(self, header: str, seq: str) -> None:
        seq = seq.encode()
        if header in self._header_idx:
            warn(f'Turning duplicate header "{header}" into unique header')
            new_header = header
            modifier = 0
            while new_header in self._header_idx:
                modifier += 1
                new_header = f"{header}_{modifier}"
            header = new_header
        n_seq, n_char = self._collection.shape
        if n_seq == 0:
            self._collection = np.array([[*seq]], dtype="uint8")
        else:
            len_diff = len(seq) - n_char

            filler1 = np.array([[*b"-"] * len_diff], dtype="uint8")
            arr = np.hstack((self._collection, np.repeat(filler1, n_seq, axis=0)))

            filler2 = np.array([*b"-"] * -len_diff, dtype="uint8")
            new_row = np.array([[*seq, *filler2]], dtype="uint8")

            arr = np.vstack((arr, new_row))
            self._collection = arr
        self._header_idx[header] = n_seq

    def __delitem__(self, header: str) -> None:
        self.pop(header)

    def __getitem__(self, header: str) -> Sequence:
        idx = self._header_idx[header]
        n_chars = self._collection.shape[1]
        sequence = self._collection[idx].view(f"S{n_chars}")[0].decode()
        return Sequence(header, sequence)

    @property
    def headers(self) -> List[str]:
        return list(self._header_idx.keys())

    @property
    def n_seqs(self) -> int:
        return self._collection.shape[0]

    @property
    def n_chars(self) -> int:
        """Number of alignment columns"""
        return self._collection.shape[1]

    @property
    def shape(self) -> Tuple[int, int]:
        """Number of sequences and number of alignment columns"""
        return self._collection.shape

    def to_nexus(self) -> str:
        """Nexus ``data`` block with the alignment. The datatype is ``protein`` if any sequence is an amino acid
        sequence, and ``dna`` otherwise.

        Examples:
            >>> msa = MultipleSequenceAlignment.from_fasta(string='>a\\nAC-GT\\n>b\\nACCGT')
            >>> print(msa.to_nexus())
            begin data;
                dimensions ntax=2 nchar=5;
                format datatype=dna gap=-;
                matrix
                a AC-GT
                b ACCGT
                ;
            end;

        Returns:
            str: Nexus formatted string
        """
        datatype = "protein" if any(seq.alphabet.name == "AminoAcid" for seq in self) else "dna"
        lines = [
            "begin data;",
            f"    dimensions ntax={self.n_seqs} nchar={self.n_chars};",
            f"    format datatype={datatype} gap=-;",
            "    matrix",
            *(f"    {seq.header} {seq.sequence}" for seq in self),
            "    ;",
            "end;",
        ]
        return "\n".join(lines)

    def pop(self, header: str) -> Sequence:
        pop_idx = self._header_idx[header]
        n_chars = self._collection.shape[1]
        sequence = self._collection[pop_idx].view(f"S{n_chars}")[0].decode()
        del self._header_idx[header]
        self._header_idx = {h: (idx if idx < pop_idx else idx - 1) for h, idx in self._header_idx.items()}
        self._collection = np.delete(self._collection, (pop_idx,), axis=0)
        return Sequence(header, sequence)

    def pairwise_distances(self, distance_measure: str = "identity") -> npt.NDArray[np.float64]:
        """Not implemented yet, returns ``None``"""
        pass
