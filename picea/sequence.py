from warnings import warn
from itertools import groupby, chain
from subprocess import Popen, PIPE
from collections import defaultdict, Counter
from typing import Iterable, List, Tuple, Dict, Any, Optional, Callable, Union
from abc import ABCMeta, abstractmethod, abstractproperty
import uuid
import numpy as np
import numpy.typing as npt
import json
from functools import reduce
import re
from dataclasses import dataclass, field
from copy import deepcopy
from .dag import DAGElement, DirectedAcyclicGraph


# (character, code) tuples for encoding special characters in gff3
# percent (%) MUST GO FIRST
ENCODE_SPECIAL_CHARACTERS = (
    ("%", "%25"),
    ("[   |\t]", "%09"),
    ("\n", "%0A"),
    (";", "%3B"),
    ("=", "%3D"),
    ("&", "%26"),
    (",", "%2C"),
)

# (code, character) tuples for decoding special characters in gff3
# percent (%) MUST GO LAST
DECODE_SPECIAL_CHARACTERS = (
    ("%2C", ","),
    ("%26", "&"),
    ("%3D", "="),
    ("%3B", ";"),
    ("%0A", "\n"),
    ("%09", "\t"),
    ("%25", "%"),
)

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

    Args:
        name (str): Alphabet name
        members (Iterable[str]): Letters of the alphabet
    """

    name: str
    members: Iterable[str]

    def __post_init__(self) -> None:
        super().__init__(self.members)

    def __deepcopy__(self, memo) -> "Alphabet":
        return Alphabet(self, self.name)

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
        if self.name != "DNA":
            raise TypeError("Cannot complement non-DNA alphabet")
        dna_complement = dict(zip("acgtnACGTN-?", "tgcanTGCAN-?"))
        return "".join(dna_complement[s] for s in sequence)

    def translate(self, sequence: str) -> str:
        if self.name != "DNA":
            raise TypeError("Cannot translate non-DNA alphabet")
        codons = re.findall("...", sequence.upper())
        return "".join(TRANSLATION.get(codon, "X") for codon in codons)


def alphabet_factory(alphabet):
    return dict(
        DNA=lambda: Alphabet("DNA", "-?ACGNT"),
        AminoAcid=lambda: Alphabet("AminoAcid", "*-?ACDEFGHIKLMNPQRSTVWXY"),
    )[alphabet]


@dataclass(frozen=True)
class Alphabets:
    DNA: Alphabet = field(default_factory=alphabet_factory("DNA"))
    AminoAcid: Alphabet = field(default_factory=alphabet_factory("AminoAcid"))

    def __iter__(self):
        return iter(self.__dict__.values())


alphabets = Alphabets()


def quote_gff3(attribute_value: Union[int, str, float]) -> str:
    """pattern, repl = ENCODE_SPECIAL_CHARACTERS[0]
    quoted_value = re.sub(pattern, repl, attribute_value)
    for pattern, repl in ENCODE_SPECIAL_CHARACTERS[1:]:
        quoted_value = re.sub(pattern, repl, attribute_value)
    return quoted_value"""
    return reduce(
        lambda acc, code: re.sub(code[0], code[1], acc),  # func
        ENCODE_SPECIAL_CHARACTERS,  # iterable
        str(attribute_value),  # initial accumulator
    )


def encode_attribute_value(attribute_value: Iterable[Union[int, str, float]]) -> str:
    """[summary]

    Args:
        attribute_value (Iterable[Union[int, str, float]]): [description]

    Returns:
        str: [description]
    """
    if not isinstance(attribute_value, (list, tuple)):
        attribute_value = [attribute_value]
    return ",".join(quote_gff3(v) for v in attribute_value)


def format_gtf_attribute_string(
    attributes: Dict[str, Iterable[Union[int, str, float]]]
) -> str:
    """[summary]

    Args:
        attributes (Dict[str, Iterable[Union[int, str, float]]]): [description]

    Returns:
        str: [description]
    """
    return "".join(
        f' {key} "{encode_attribute_value(value)}";'
        for key, value in attributes.items()
    ).strip()


def format_gff_attribute_string(
    attributes: Dict[str, Iterable[Union[int, str, float]]]
) -> str:
    """[summary]

    Args:
        attributes (Dict[str, Iterable[Union[int, str, float]]]): [description]

    Returns:
        str: [description]
    """
    partially_formatted = {
        (
            key.capitalize()
            if key in SequenceInterval._predefined_gff3_attributes
            else key
        ): encode_attribute_value(value)
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

    return ";".join(
        f"{key}={partially_formatted[key]}"
        for key in sorted(partially_formatted.keys(), key=sort_key)
    )


def unquote_gff3(attribute_value: str) -> str:
    """[summary]

    Args:
        attribute_value (str): [description]

    Returns:
        str: [description]
    """
    return reduce(
        lambda acc, code: re.sub(code[0], code[1], acc),  # func
        DECODE_SPECIAL_CHARACTERS,  # iterable
        attribute_value,  # initial
    )


def decode_attribute_value(attribute_value: str) -> List[str]:
    """[summary]

    Args:
        attribute_value (str): [description]

    Returns:
        List[str]: [description]
    """
    return [unquote_gff3(v) for v in attribute_value.split(",")]


def parse_gtf_attribute_string(gtf_attribute_string: str) -> Dict[str, List[str]]:
    """[summary]

    Args:
        gtf_attribute_string (str): [description]

    Returns:
        Dict[str, List[str]]: [description]
    """
    attributes = defaultdict(list)
    for string_part in gtf_attribute_string.split(";"):
        string_part = string_part.strip()
        if not string_part:
            continue
        try:
            key, value = string_part.split(" ", maxsplit=1)
        except Exception as e:
            print(gtf_attribute_string, string_part)
            raise Exception(e)
        attributes[key].append(value.strip('"'))
    return attributes


def parse_gff_attribute_string(
    gff_attribute_string: str, case_sensitive_attribute_keys: bool = False
) -> Dict[str, List[str]]:
    """[summary]
    https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
    See "Column 9: Attributes"
    Args:
        gff_attribute_string ([type]): [description]
    """
    attributes = defaultdict(list)
    for string_part in gff_attribute_string.split(";"):
        if not string_part:
            continue
        try:
            key, value = string_part.split("=", maxsplit=1)
        except Exception as e:
            print(gff_attribute_string, string_part)
            raise Exception(e)
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
    def __init__(self, sequence: Optional["Sequence"] = None) -> None:
        """[summary]

        Args:
            sequence (Optional[Sequence], optional): [description]. Defaults\
                 to None.
        """
        super().__init__()
        if sequence:
            sequence.annotation = self
        self.sequence = sequence
        self._gff_headers = list()

    @property
    def intervals(self):
        return list(self)

    @classmethod
    def from_gtf(
        cls,
        filename: Optional[str] = None,
        string: Optional[str] = None,
        sequence: Optional["Sequence"] = None,
    ) -> "SequenceAnnotation":
        """[summary]

        Raises:
            IndexError: [description]
            IndexError: [description]

        Returns:
            [type]: [description]
        """
        assert filename or string
        assert not (filename and string)
        sequence_annotation = cls(sequence=sequence)
        header = True

        # start with just reading all intervals
        if filename:
            with open(filename) as filehandle:
                string = filehandle.read()
        for line_number, line in enumerate(string.split("\n")):
            line = line.strip()
            if not line:
                continue
            if line[0] == "#":
                if header:
                    sequence_annotation._gff_headers.append(line)
                continue
            else:
                header = False
            interval = SequenceInterval.from_gtf_line(
                gtf_line=line, line_number=line_number
            )
            interval._container = sequence_annotation
            sequence_annotation[interval.ID] = interval
        # fix missing gene and transcript intervals
        transcript_child_counter = Counter()
        new_intervals = dict()
        for interval in sequence_annotation:
            gene_id = interval.gff_attributes["gene_id"][0]
            transcript_id = interval.gff_attributes["transcript_id"][0]
            interval_type = interval.interval_type
            id_tuple = (gene_id, transcript_id, interval_type)
            child_count = transcript_child_counter[id_tuple]
            transcript_child_counter.update([id_tuple])
            interval._ID = f"{transcript_id}.{interval_type}_{child_count}"
            if transcript_id not in new_intervals:
                # new transcript interval
                transcript_interval = deepcopy(interval)
                transcript_interval._container = interval._container
                transcript_interval._ID = transcript_id
                transcript_interval.interval_type = "mRNA"
                transcript_interval.parent = [gene_id]
                # new gene interval
                gene_interval = deepcopy(interval)
                gene_interval._container = interval._container
                gene_interval._ID = gene_id
                gene_interval.interval_type = "gene"
                gene_interval.parent = None

                new_intervals[transcript_id] = transcript_interval
                new_intervals[gene_id] = gene_interval

            interval.parent = [transcript_id]
            new_intervals[interval.ID] = interval
        sequence_annotation._intervals = new_intervals

        # set children
        for interval in sequence_annotation:
            if interval.parent:
                for parent_ID in interval.parent:
                    try:
                        parent = sequence_annotation[parent_ID]
                    except IndexError:
                        raise IndexError(
                            "Interval {interval.ID} is listing {parent_ID} "
                            "as Parent, but parent could not be found."
                        )
                    parent._children.append(interval.ID)

        # fix gene and transcript start and stop coordinates
        genes = sequence_annotation.groupby("interval_type")["gene"]
        for gene in genes:
            # fix gene first
            start = 10e9
            end = 0
            for child in gene.children:
                start = min(start, child.start)
                end = max(end, child.end)
            gene.start = start
            gene.end = end

            # fix transcripts
            transcripts = gene.children.groupby("interval_type")["mRNA"]
            for transcript in transcripts:
                start = 10e9
                end = 0
                for child in transcript.children:
                    start = min(start, child.start)
                    end = max(end, child.end)
                transcript.end = end
                transcript.start = start

        return sequence_annotation

    def to_gtf(self) -> str:
        return "\n".join(interval.to_gtf_line() for interval in self)

    @classmethod
    def from_gff(
        cls,
        filename: Optional[str] = None,
        string: Optional[str] = None,
        sequence: Optional["Sequence"] = None,
        link_parents: bool = True,
    ) -> "SequenceAnnotation":
        """[summary]

        Args:
            filename ([type], optional): [description]. Defaults to None.
            string ([type], optional): [description]. Defaults to None.
            sequence ([type], optional): [description].
                Defaults to None.

        Returns:
            [type]: [description]
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

            interval = SequenceInterval.from_gff_line(
                gff_line=line, line_number=line_number
            )
            interval._container = sequence_annotation
            sequence_annotation[interval.ID] = interval

        if link_parents:
            for interval in sequence_annotation:
                if interval.parent:
                    for parent_ID in interval.parent:
                        try:
                            parent = sequence_annotation[parent_ID]
                        except IndexError:
                            raise IndexError(
                                "Interval {interval.ID} is listing {parent_ID}"
                                " as Parent, but parent could not be found."
                            )
                        parent._children.append(interval.ID)

        return sequence_annotation

    def to_gff(self) -> str:
        """[summary]

        Returns:
            str: [description]
        """
        return "".join(interval.to_gff_line(trailing_newline=True) for interval in self)

    @classmethod
    def from_json(
        cls,
        filename: Optional[str] = None,
        string: Optional[str] = None,
        sequence: Optional["Sequence"] = None,
    ) -> "SequenceAnnotation":
        """[summary]"""
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
        for interval in sequence_annotation:
            if interval.parent:
                for parent_ID in interval.parent:
                    try:
                        parent = sequence_annotation[parent_ID]
                    except IndexError:
                        raise IndexError(
                            "Interval {interval.ID} is listing {parent_ID} "
                            "as Parent, but parent could not be found."
                        )
                    parent._children.append(interval.ID)
        return sequence_annotation

    def to_json(self, indent: Optional[int] = None) -> str:
        """[summary]

        Returns:
            str: [description]
        """
        interval_dicts = [interval.to_dict() for interval in self]
        return json.dumps(interval_dicts, indent=indent)

    def add(self, interval: "SequenceInterval"):
        self[interval.ID] = interval

    def pop(self, value: str) -> "SequenceInterval":
        return self._intervals.pop(value)

    def groupby(self, key="seqid") -> Dict[str, "SequenceAnnotation"]:
        grouped = defaultdict(SequenceAnnotation)
        for interval in self:
            grouped[interval[key]][interval.ID] = interval
            interval._container = self
        return grouped

    def filter(self, filter_func: Callable) -> "SequenceAnnotation":
        result = SequenceAnnotation()
        for interval in self:
            if filter_func(interval):
                result[interval.ID] = interval
        return result


class SequenceInterval(DAGElement):
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
        """[summary]

        Args:
            ID (Optional[str], optional): [description]. Defaults to None.
            seqid (Optional[str], optional): [description]. Defaults to None.
            source (Optional[str], optional): [description]. Defaults to None.
            interval_type (Optional[str], optional): [description]. Defaults
                to None.
            start (Optional[int], optional): [description]. Defaults to None.
            end (Optional[int], optional): [description]. Defaults to None.
            score (Optional[float], optional): [description]. Defaults to None.
            strand (Optional[str], optional): [description]. Defaults to None.
            phase (Optional[str], optional): [description]. Defaults to
                None.
            children (Optional[List], optional): [description]. Defaults to
                None.
            container (Optional[SequenceAnnotation], optional): [description].
                Defaults to None.
        """
        # interval ID is a property (see below) with getter and setter
        # self._ID = ID
        # self._original_ID = ID
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

        # Additional fields, used internally
        # self._container = container
        # if children is None:
        #     children = []
        # self._children = children

    def __repr__(self):
        return (
            f"<SequenceInterval type={self.interval_type} "
            f"ID={self.ID} "
            f"loc={self.seqid}..{self.start}..{self.end}..{self.strand} "
            f"at {hex(id(self))}>"
        )

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __deepcopy__(self, memo):
        new_interval = SequenceInterval()
        for key, value in vars(self).items():
            # _container is strictly a reference, do not deepcopy
            if key == "_container":
                new_interval[key] = value
            else:
                new_interval[key] = deepcopy(value)
        return new_interval

    @property
    def parent(self):
        return self._parents

    @property
    def gff_attributes(self) -> Dict[str, str]:
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
        def get_gtf_type(gff_type):
            return self._gtf_interval_types.get(gff_type, gff_type)

        parent_ids = {
            f"{get_gtf_type(parent.interval_type)}_id": parent.ID
            for parent in self._get_parents()
        }
        return {**self.gff_attributes, **parent_ids}

    @classmethod
    def from_gtf_line(
        cls, gtf_line: Optional[str] = None, line_number: Optional[int] = None
    ) -> "SequenceInterval":
        """[summary]

        Returns:
            [type]: [description]

        Yields:
            [type]: [description]
        """
        return cls.from_gff_line(gtf_line, line_number, parse_gtf_attribute_string)

    def to_gtf_line(self) -> str:
        """[summary]

        Returns:
            str: [description]
        """
        interval_type = self._gtf_interval_types.get(
            self.interval_type, self.interval_type
        )
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
        """[summary]

        Args:
            gff_line (Optional[str], optional): [description]. Defaults
                to None.
            line_number (Optional[int], optional): [description]. Defaults
                to None.

        Returns:
            [type]: [description]
        """
        gff_parts = gff_line.split("\t")
        assert len(gff_parts) == 9, gff_parts
        seqid, source, interval_type, start, end, score, strand, phase = gff_parts[:8]
        try:
            start = int(start)
            end = int(end)
        except ValueError:
            error = "GFF start and end fields must be integer"
            if line_number:
                error = f"{error}, gff line {line_number}"
            raise ValueError(error)

        if score != ".":
            try:
                score = float(score)
            except ValueError:
                error = "GFF score field must be a float"
                if line_number:
                    error = f"{error}, gff line {line_number}"
                raise ValueError(error)

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
        """[summary]

        Returns:
            str: [description]
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
        """[summary]
        Args:
            interval_dict

        Returns:
            [type]: [description]
        """
        attributes = interval_dict.pop("attributes", dict())
        return cls(**interval_dict, **attributes)

    def to_dict(self, include_children: bool = False) -> Dict[str, Any]:
        """[summary]

        Returns:
            Dict[str, Any]: [description]
        """
        attributes = dict(**self.attributes)
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
            children = [child.to_dict() for child in self.children[1:]]
            interval_dict["children"] = children
        return interval_dict

    def to_json(
        self, include_children: bool = False, indent: Optional[int] = None
    ) -> str:
        """[summary]

        Args:
            include_children (bool, optional): [description]. Defaults to \
                False.

        Returns:
            str: [description]
        """
        return json.dumps(
            self.to_dict(include_children=include_children), indent=indent
        )


@dataclass
class Sequence:
    """Container for a single biological sequence

    Examples:
        >>> s1 = Sequence('test_dna', 'ACGATCGACTAGCA')
        >>> s1
        Sequence(header='test_dna', \
alphabet=Alphabet(name='DNA', members='-?ACGNT'))
        >>> s2 = Sequence('test_aa', 'QAPISAIWPOIWQ*')
        >>> s2
        Sequence(header='test_aa', \
alphabet=Alphabet(name='AminoAcid', members='*-?ACDEFGHIKLMNPQRSTVWXY'))

    Returns:
        [type]: [description]
    """

    header: str = None
    sequence: str = field(repr=False, default=None)
    alphabet: Alphabet = None
    annotation: Optional[SequenceAnnotation] = field(
        default_factory=SequenceAnnotation, repr=False
    )

    def __post_init__(self):
        if self.alphabet is not None:
            return
        if self.sequence is None:
            self.alphabet = alphabets.DNA
        else:
            self.alphabet = sorted(
                alphabets, key=lambda alphabet: alphabet.score(self.sequence)
            ).pop()

    def __getitem__(self, key):
        return Sequence(self.header, self.sequence[key])

    def __len__(self):
        """Length of the sequence

        Examples:
            >>> s = Sequence('test_dna', 'ACGTA')
            >>> len(s)
            5
        """
        return len(self.sequence)

    @property
    def reverse(self):
        return Sequence(self.header, self.sequence[::-1])

    @property
    def complement(self):
        return Sequence(self.header, self.alphabet.complement(self.sequence))

    @property
    def reverse_complement(self):
        return Sequence(self.header, self.alphabet.complement(self.sequence[::-1]))

    @property
    def amino_acids(self):
        if self.alphabet.name == "AminoAcid":
            return self
        else:
            return Sequence(self.header, self.alphabet.translate(self.sequence))

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
        """Create a sequence object from a fasta formatted file. _single sequence only_

        Examples:
            >>> fasta_string = '>test\\nACGT'
            >>> Sequence.from_fasta(fasta_string)
            Sequence(header='test', \
alphabet=Alphabet(name='DNA', members='-?ACGNT'))

        Arguments:
            string (str)

        Returns:
            Sequence
        """
        lines = string.strip().split("\n")
        header = lines[0][1:]
        sequence = "".join(lines[1:])
        return cls(header, sequence)

    def to_fasta(self, linewidth: int = 80) -> str:
        """Make fasta formatted sequence entry

        Returns:
            str: sequence in fasta format
        """
        sequence_lines = "\n".join(re.findall(f".{{1,{linewidth}}}", self.sequence))
        return f">{self.header}\n{sequence_lines}"


class FastaParseError(Exception):
    pass


class SequenceReader:
    def __init__(
        self, string: str = None, filename: str = None, filetype: str = None
    ) -> None:
        """Iterator over fasta/json formatted sequence strings

        Args:
            string (str): fasta/json formatted string
            filename (str): fasta/json file name
            filetype (str): fasta or json

        Examples:
            >>> fasta_string = '>1\\nACGC\\n>2\\nTGTGTA\\n'
            >>> fasta_reader = SequenceReader(string=fasta_string, \
filetype='fasta')
            >>> next(fasta_reader)
            Sequence(header='1', \
alphabet=Alphabet(name='DNA', members='-?ACGNT'))
        """
        assert bool(string) ^ bool(
            filename
        ), "Must specify exactly one of string or filename"  # exclusive OR
        if filename:
            with open(filename, "r") as filehandle:
                string = filehandle.read().strip()

        self.string = string

        if filetype == "fasta":
            self._iter = self._fasta_iter
        elif filetype == "json":
            self._iter = self._json_iter
        else:
            raise ValueError(f'filetype "{filetype}" is not supported')

    def __iter__(self) -> Iterable[Sequence]:
        """Iterate over header,sequence tuples

        Returns:
            Iterable[Sequence]: [description]

        Yields:
            Iterable[Sequence]: [description]
        """
        try:
            yield from self._iter()
        except StopIteration:
            return

    def __next__(self) -> Sequence:
        """Next header and sequence in the iterator

        Returns:
            Tuple[str, str]: [description]
        """
        return next(self._iter())

    def _fasta_iter(self) -> Iterable[Sequence]:
        if self.string[0] != ">":
            raise FastaParseError(
                'First character in fasta format\
                 must be ">"'
            )
        fasta_iter = (
            x
            for _, x in groupby(
                self.string.strip().split("\n"), lambda line: line[:1] == ">"
            )
        )
        for header in fasta_iter:
            header = next(header)[1:].strip()
            seq = "".join(s.strip() for s in next(fasta_iter))
            yield Sequence(header, seq)

    def _json_iter(self) -> Iterable[Sequence]:
        for entry in json.loads(self.string):
            yield Sequence(entry["header"], entry["sequence"])


class BatchSequenceReader(SequenceReader):
    def __init__(
        self,
        string: str = None,
        filename: str = None,
        filetype: str = None,
        batchsize: int = 10,
    ) -> None:
        """[summary]

        Args:
            string (str, optional): [description]. Defaults to None.
            filename (str, optional): [description]. Defaults to None.
            filetype (str, optional): [description]. Defaults to None.
            batchsize (int, optional): [description]. Defaults to 10.

        Returns:
            [type]: [description]

        Yields:
            [type]: [description]
        """
        super().__init__(string, filename, filetype)
        self.batchsize = batchsize
        self._currentbatch = SequenceCollection()

    def __iter__(self) -> Iterable["SequenceCollection"]:
        for s in self._iter():
            self._currentbatch[s.header] = s
            if len(self._currentbatch) == self.batchsize:
                yield self._currentbatch
                self._currentbatch = SequenceCollection()

    def __next__(self) -> "SequenceCollection":
        currentbatch = self._currentbatch
        self._currentbatch = SequenceCollection()
        if len(currentbatch) == self.batchsize:
            return currentbatch
        for s in self._iter():
            currentbatch[s.header] = s
            if len(currentbatch) == self.batchsize:
                yield currentbatch


SequenceIndexKey = Union[int, List[int], slice]


class SequenceIndex:
    def __init__(
        self,
        sequence_collection: Union["SequenceCollection", "MultipleSequenceAlignment"],
    ):
        self.sequence_collection = sequence_collection

    def __getitem__(self, key: SequenceIndexKey):
        if isinstance(key, int):
            key = [key]
        elif isinstance(key, slice):
            start = key.start if key.start is not None else 0
            stop = key.stop if key.stop is not None else len(self.sequence_collection)
            step = key.step if key.step is not None else 1
            key = range(start, stop, step)
        elif not isinstance(key, list):
            raise TypeError(f"SequenceIndex key must be of type f{SequenceIndexKey}")
        new_seq_col = self.sequence_collection.__class__()
        for k in key:
            header = self.sequence_collection.headers[k]
            new_seq_col[header] = self.sequence_collection[header].sequence
        return new_seq_col


class AbstractSequenceCollection(metaclass=ABCMeta):
    """
    (Partially) Abstract Base Class for sequence collections.
    Classes extending from this baseclass should override
    `__setitem__`, `__getitem__`, `__delitem__`, `headers`, and `n_seqs`.

    If the above methods are implemented, this automatically enables the
    following methods: `from_fasta`, `to_fasta`, `from_json`, `to_json`.

    Args:
        sequences (Optional[Iterable[Tuple[str, str]]], optional):
            Iterable of (header, sequence) tuples. Defaults to None.
        sequence_annotation (Optional[SequenceAnnotation]):
                picea SequenceAnnotation object. Defaults to None.

    Raises:
        NotImplementedError: Abstract Base Class can not be initialized
            and serves as a template only
    """

    @abstractmethod
    def __init__(
        self,
        sequences: Optional[Iterable[Sequence]] = None,
        sequence_annotation: Optional["SequenceAnnotation"] = None,
    ) -> None:
        raise NotImplementedError(
            (
                "Classes extending from AbstractSequenceCollection should "
                "implement __init__ method"
            )
        )

    @abstractmethod
    def __setitem__(self, header: str, seq: str) -> None:
        raise NotImplementedError(
            (
                "Classes extending from AbstractSequenceCollection should "
                "implement __setitem__ method"
            )
        )

    @abstractmethod
    def __getitem__(self, header: str) -> Sequence:
        raise NotImplementedError(
            (
                "Classes extending from AbstractSequenceCollection should "
                "implement __getitem__ method"
            )
        )

    @abstractmethod
    def __delitem__(self, header: str) -> None:
        raise NotImplementedError(
            (
                "Classes extending from AbstractSequenceCollection should "
                "implement __delitem__ method"
            )
        )

    def __iter__(self) -> Iterable[Sequence]:
        for header in self.headers:
            yield self[header]

    def __len__(self) -> int:
        return len(self.headers)

    @abstractproperty
    def headers(self) -> List[str]:
        """List of sequences headers.
        Overridden in subclasses.

        Raises:
            NotImplementedError

        Returns:
            List[str]: List of sequence headers
        """
        raise NotImplementedError(
            (
                "Classes extending from AbstractSequenceCollection should "
                "implement headers property"
            )
        )

    @property
    def iloc(self) -> SequenceIndex:
        """[summary]

        Returns:
            SequenceIndex: [description]
        """
        return SequenceIndex(self)

    @property
    def sequences(self) -> List[str]:
        """List of sequences without headers

        Returns:
            List[str]: list of sequences
        """
        return [self[header].sequence for header in self.headers]

    @abstractproperty
    def n_seqs(self) -> int:
        """Return the number of sequences in the collection.
        Overridden in subclasses

        Raises:
            NotImplementedError

        Returns:
            int: number of sequences
        """
        raise NotImplementedError(
            (
                "Classes extending from AbstractSequenceCollection should "
                "implement n_seqs property"
            )
        )

    @classmethod
    def from_sequence_iter(
        cls, sequence_iter: Iterable[Sequence]
    ) -> "SequenceCollection":
        """[summary]

        Raises:
            NotImplementedError: [description]

        Returns:
            [type]: [description]
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
        """Parse a fasta formatted string into a SequenceCollection object

        Keyword Arguments:
            filename {String} -- filename string (default: {None})
            string {String} -- fasta formatted string (default: {None})

        Returns:
            SequenceCollection -- SequenceCollection instance
        """
        sequencecollection = cls()

        for seq in SequenceReader(string=string, filename=filename, filetype="fasta"):
            sequencecollection[seq.header] = seq.sequence
        return sequencecollection

    def to_fasta(self, linewidth: int = 80) -> str:
        """Get a fasta-formatted string of the sequence collection

        Returns:
            str: Multi-line fasta-formatted string
        """
        return "\n".join([seq.to_fasta(linewidth=linewidth) for seq in self])

    @classmethod
    def from_json(
        cls, filename: Optional[str] = None, string: Optional[str] = None
    ) -> "SequenceCollection":
        """[summary]

        Keyword Arguments:
            string {String} -- JSON formatted string

        Returns:
            SequenceCollection -- SequenceCollection instance
        """
        sequencecollection = cls()

        for seq in SequenceReader(string=string, filename=filename, filetype="json"):
            sequencecollection[seq.header] = seq.sequence

        return sequencecollection

    def to_json(self, indent: Optional[int] = None) -> str:
        """[summary]

        Returns:
            str: [description]
        """
        gene_dicts = [seq.to_dict() for seq in self]
        return json.dumps(gene_dicts, indent=indent)

    @abstractmethod
    def pop(self, header: str) -> Sequence:
        """[summary]

        Args:
            header (str): [description]

        Returns:
            Sequence: [description]
        """
        raise NotImplementedError(
            (
                "Classes extending from AbstractSequenceCollection should "
                "implement pop method"
            )
        )

    def batch_rename(self, rename_func: Callable[[str], str]) -> None:
        """Rename all headers by calling `rename_func` on each header

        Args:
            rename_func (Callable): [description]
        """
        for header in self.headers:
            s: Sequence = self.pop(header)
            s.header = rename_func(s.header)
            self[s.header] = s.sequence


class SequenceCollection(AbstractSequenceCollection):
    """
    A container for multiple (unaligned) DNA or amino acid sequences
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

    def align(self, method: str = "mafft", method_kwargs: Dict[str, str] = dict()):
        """[summary]

        Args:
            method (str, optional): [description]. Defaults to 'mafft'.
            method_kwargs (Dict[str, str], optional): [description]. \
                Defaults to dict().

        Returns:
            [type]: [description]
        """
        fasta = self.to_fasta()
        command = [method, *chain(*method_kwargs.items()), "-"]
        process = Popen(command, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate(input=fasta.encode())
        aligned_fasta = stdout.decode().strip()
        return MultipleSequenceAlignment.from_fasta(string=aligned_fasta)

    def pop(self, header: str) -> Sequence:
        sequence = self._collection.pop(header)
        return Sequence(header, sequence)


class MultipleSequenceAlignment(SequenceCollection):
    def __init__(
        self,
        sequences: Optional[Iterable[Sequence]] = None,
        sequence_annotation: Optional["SequenceAnnotation"] = None,
    ) -> None:
        """
        A container for multiple aligned DNA or amino acid sequences
        """
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

    def __getitem__(self, header: str) -> Sequence:
        idx = self._header_idx[header]
        n_chars = self._collection.shape[1]
        sequence = self._collection[idx].view(f"S{n_chars}")[0].decode()
        return Sequence(header, sequence)

    def __delitem__(self, header: str) -> None:
        """WIP!

        Args:
            header (str): [description]
        """
        idx = self._header_idx[header]
        self._collection = np.delete(self._collection, idx, axis=0)
        # del

    @property
    def headers(self) -> List[str]:
        return list(self._header_idx.keys())

    @property
    def n_seqs(self) -> int:
        return self._collection.shape[0]

    @property
    def n_chars(self) -> int:
        return self._collection.shape[1]

    @property
    def shape(self) -> int:
        return self._collection.shape

    def pop(self, header: str) -> Sequence:
        pop_idx = self._header_idx[header]
        n_chars = self._collection.shape[1]
        sequence = self._collection[pop_idx].view(f"S{n_chars}")[0].decode()
        del self._header_idx[header]
        self._header_idx = {
            h: (idx if idx < pop_idx else idx - 1)
            for h, idx in self._header_idx.items()
        }
        self._collection = np.delete(self._collection, (pop_idx,), axis=0)
        return Sequence(header, sequence)

    def pairwise_distances(
        self, distance_measure: str = "identity"
    ) -> npt.NDArray[np.float64]:
        pass
