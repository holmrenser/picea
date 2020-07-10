from warnings import warn
from itertools import groupby, chain
from subprocess import Popen, PIPE
from collections import defaultdict
from typing import Iterable, List, Tuple, Dict, Any, Optional
from abc import ABCMeta, abstractmethod, abstractproperty
import uuid
import numpy as np
import json
from functools import reduce
import re
from dataclasses import dataclass, field


# (character, code) tuples for encoding special characters in gff3
# percent (%) MUST GO FIRST
ENCODE_SPECIAL_CHARACTERS = (
    ('%', '%25'),
    ('[   |\t]', '%09'),
    ('\n', '%0A'),
    (';', '%3B'),
    ('=', '%3D'),
    ('&', '%26'),
    (',', '%2C')
)

# (code, character) tuples for decoding special characters in gff3
# percent (%) MUST GO LAST
DECODE_SPECIAL_CHARACTERS = (
    ('%2C', ','),
    ('%26', '&'),
    ('%3D', '='),
    ('%3B', ';'),
    ('%0A', '\n'),
    ('%09', '\t'),
    ('%25', '%')
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

    def __post_init__(self):
        super().__init__(self.members)

    def __deepcopy__(self, memo):
        return Alphabet(self, self.name)

    def score(
        self,
        sequence: str,
        match: float = 1.0,
        mismatch: float = -1.0
    ) -> float:
        """Scores how well a sequence matches an alphabet by summing \
        (mis)matches of sequence letters that are not in the alphabet \
        and (mis)matches of alphabet letters that are not in the sequence.

        Args:
            sequence (str): Sequence string for which to determine how well \
                it fits the alphabet
            match (float, optional): match score. Defaults to 1.0.
            mismatch (float, optional): mismatch score. Defaults to -1.0.

        Returns:
            (float): Score of how well a sequence matches the alphabet
        """
        return sum(match if s in self else mismatch for s in sequence) \
            + sum(match if s in sequence else mismatch for s in self)

    def validate(self, sequence: str) -> bool:
        """Determine whether a sequence strictly fits an alphabet

        Args:
            sequence (str): Sequence string

        Returns:
            bool: true if all characters in sequence are in the alphabet
        """
        return sum(1 if s not in self else 0 for s in sequence) == 0


def alphabet_factory(alphabet):
    return dict(
        DNA=lambda: Alphabet('DNA', '-?ACGNT'),
        AminoAcid=lambda: Alphabet('AminoAcid', '*-?ACDEFGHIKLMNPQRSTVWXY')
    )[alphabet]


@dataclass(frozen=True)
class Alphabets:
    DNA: Alphabet = field(default_factory=alphabet_factory('DNA'))
    AminoAcid: Alphabet = field(default_factory=alphabet_factory('AminoAcid'))

    def __iter__(self):
        return iter(self.__dict__.values())


alphabets = Alphabets()


@dataclass
class Sequence:
    """Container for a single biological sequence

    Examples:
        >>> s1 = Sequence('test_dna', 'ACGATCGACTAGCA')
        >>> s1
        Sequence(header='test_dna', sequence='ACGATCGACTAGCA', \
alphabet=Alphabet(name='DNA', members='-?ACGNT'))
        >>> s2 = Sequence('test_aa', 'QAPISAIWPOIWQ*')
        >>> s2
        Sequence(header='test_aa', sequence='QAPISAIWPOIWQ*', \
alphabet=Alphabet(name='AminoAcid', members='*-?ACDEFGHIKLMNPQRSTVWXY'))

    Returns:
        [type]: [description]
    """
    header: str = None
    sequence: str = None
    alphabet: Alphabet = None

    def __post_init__(self):
        if self.alphabet is not None:
            return
        if self.sequence is None:
            self.alphabet = alphabets.DNA
        else:
            self.alphabet = sorted(
                alphabets,
                key=lambda alphabet: alphabet.score(self.sequence)
            ).pop()

    def to_dict(self) -> Dict[str, str]:
        """Make dictionary with header and sequence elements

        Returns:
            Dict[str, str]: sequence dictionary
        """
        return dict(
            header=self.header,
            sequence=self.sequence
        )

    @classmethod
    def from_fasta(cls, string: str) -> 'Sequence':
        """Create a sequence object from a fasta formatted file. _single sequence only_

        Examples:
            >>> fasta_string = '>test\\nACGT'
            >>> Sequence.from_fasta(fasta_string)
            Sequence(header='test', sequence='ACGT', \
alphabet=Alphabet(name='DNA', members='-?ACGNT'))

        Arguments:
            string (str)

        Returns:
            Sequence
        """
        lines = string.strip().split('\n')
        header = lines[0][1:]
        sequence = ''.join(lines[1:])
        return cls(header, sequence)

    def to_fasta(self) -> str:
        """Make fasta formatted sequence entry

        Returns:
            str: sequence in fasta format
        """
        return f'>{self.header}\n{self.sequence}'


class SequenceReader:
    def __init__(
        self,
        string: str = None,
        filename: str = None,
        filetype: str = None
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
            Sequence(header='1', sequence='ACGC', \
alphabet=Alphabet(name='DNA', members='-?ACGNT'))
        """
        assert bool(string) ^ bool(filename), \
            'Must specify exactly one of string or filename'  # exclusive OR
        if filename:
            with open(filename, 'r') as filehandle:
                string = filehandle.read().strip()

        self.string = string

        if filetype == 'fasta':
            self._iter = self._fasta_iter
        elif filetype == 'json':
            self._iter = self._json_iter
        else:
            raise Exception(f'filetype "{filetype}" is not supported')

    def __iter__(self) -> Iterable[Sequence]:
        """Iterate over header,sequence tuples

        Returns:
            Iterable[Tuple[str, str]]: [description]

        Yields:
            Iterable[Tuple[str, str]]: [description]
        """
        yield from self._iter()

    def __next__(self) -> Sequence:
        """Next header and sequence in the iterator

        Returns:
            Tuple[str, str]: [description]
        """
        #header = next(next(self._iter))[1:].strip()
        #seq = ''.join(s.strip() for s in next(self._iter))
        #return Sequence(header, seq)
        return next(self._iter())

    def _fasta_iter(self):
        fasta_iter = (
            x for _, x in groupby(
                self.string.strip().split('\n'),
                lambda line: line[0] == '>'
            )
        )
        for header in fasta_iter:
            header = next(header)[1:].strip()
            seq = ''.join(s.strip() for s in next(fasta_iter))
            yield Sequence(header, seq)

    def _json_iter(self):
        for entry in json.loads(self.string):
            yield Sequence(entry['header'], entry['sequence'])


class SequenceCollection(metaclass=ABCMeta):
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
        sequence_annotation: Optional['SequenceAnnotation'] = None
    ) -> None:
        raise NotImplementedError('Not implemented in base class')

    @abstractmethod
    def __setitem__(self, header: str, seq: str) -> None:
        raise NotImplementedError()

    @abstractmethod
    def __getitem__(self, header: str) -> Sequence:
        raise NotImplementedError()

    @abstractmethod
    def __delitem__(self, header: str) -> None:
        raise NotImplementedError()

    def __iter__(self) -> Iterable[Sequence]:
        for header in self.headers:
            yield self[header]

    def __next__(self) -> Sequence:
        header = next(self.headers)
        return self[header]

    @abstractproperty
    def headers(self) -> List[str]:
        """List of sequences headers.
        Overridden in subclasses.

        Raises:
            NotImplementedError

        Returns:
            List[str]: List of sequence headers
        """
        raise NotImplementedError()

    @property
    def sequences(self) -> List[str]:
        """List of sequences without headers

        Returns:
            List[str]: list of sequences
        """
        return [self[header] for header in self.headers]

    @abstractproperty
    def n_seqs(self) -> int:
        """Return the number of sequences in the collection.
        Overridden in subclasses

        Raises:
            NotImplementedError

        Returns:
            int: number of sequences
        """
        raise NotImplementedError()

    @classmethod
    def from_fasta(
        cls,
        filename: str = None,
        string: str = None,
    ) -> 'SequenceCollection':
        """Parse a fasta formatted string into a SequenceCollection object

        Keyword Arguments:
            filename {String} -- filename string (default: {None})
            string {String} -- fasta formatted string (default: {None})

        Returns:
            SequenceCollection -- SequenceCollection instance
        """
        sequencecollection = cls()

        for seq in SequenceReader(
            string=string, filename=filename, filetype='fasta'
        ):
            sequencecollection[seq.header] = seq.sequence
        return sequencecollection

    def to_fasta(self) -> str:
        """Get a fasta-formatted string of the sequence collection

        Returns:
            str: Multi-line fasta-formatted string
        """
        return '\n'.join([seq.to_fasta() for seq in self])

    @classmethod
    def from_json(
        cls,
        filename: Optional[str] = None,
        string: Optional[str] = None
    ) -> 'SequenceCollection':
        """[summary]

        Keyword Arguments:
            string {String} -- JSON formatted string

        Returns:
            SequenceCollection -- SequenceCollection instance
        """
        sequencecollection = cls()

        for seq in SequenceReader(
            string=string, filename=filename, filetype='json'
        ):
            sequencecollection[seq.header] = seq.sequence

        return sequencecollection

    def to_json(self, indent: Optional[int] = None) -> str:
        """[summary]

        Returns:
            str: [description]
        """
        gene_dicts = [seq.to_dict() for seq in self]
        return json.dumps(gene_dicts, indent=indent)


class SequenceList(SequenceCollection):
    """
    A container for multiple (unaligned) DNA or amino acid sequences
    """
    def __init__(
        self: 'SequenceCollection',
        sequences: Iterable[Tuple[str, str]] = None,
        sequence_annotation: 'SequenceAnnotation' = None
    ):
        self._collection = dict()
        if sequences:
            for header, sequence in sequences:
                self[header] = sequence
        self.sequence_annotation = sequence_annotation

    def __setitem__(self, header: str, seq: str) -> None:
        if header in self.headers:
            warn(f'Turning duplicate header "{header}" into unique header')
            new_header = header
            modifier = 0
            while new_header in self.headers:
                modifier += 1
                new_header = f'{header}_{modifier}'
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
        self,
        method: str = 'mafft',
        method_kwargs: Dict[str, str] = dict()
    ):
        """[summary]

        Args:
            method (str, optional): [description]. Defaults to 'mafft'.
            method_kwargs (Dict[str, str], optional): [description]. \
                Defaults to dict().

        Returns:
            [type]: [description]
        """
        fasta = self.to_fasta()
        command = [method, *chain(*method_kwargs.items()), '-']
        process = Popen(
            command,
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE
        )
        stdout, stderr = process.communicate(input=fasta.encode())
        aligned_fasta = stdout.decode().strip()
        return MultipleSequenceAlignment.from_fasta(string=aligned_fasta)


class MultipleSequenceAlignment(SequenceCollection):
    def __init__(
        self,
        sequences: Optional[Iterable[Sequence]] = None,
        sequence_annotation: Optional['SequenceAnnotation'] = None
    ) -> None:
        """
        A container for multiple aligned DNA or amino acid sequences
        """
        self._collection = np.empty((0, 0), dtype='uint8')
        self._header_idx = dict()
        if sequences:
            for seq in sequences:
                self[seq.header] = seq.sequence
        if sequence_annotation:
            sequence_annotation.sequence_collection = self
        self.sequence_annotation = sequence_annotation

    def __setitem__(self, header: str, seq: str) -> None:
        seq = seq.encode()
        if header in self._header_idx:
            warn(f'Turning duplicate header "{header}" into unique header')
            new_header = header
            modifier = 0
            while new_header in self._header_idx:
                modifier += 1
                new_header = f'{header}_{modifier}'
            header = new_header
        n_seq, n_char = self._collection.shape
        if n_seq == 0:
            self._collection = np.array([[*seq]], dtype='uint8')
        else:
            len_diff = len(seq) - n_char

            filler1 = np.array([[*b'-'] * len_diff], dtype='uint8')
            arr = np.hstack((
                self._collection,
                np.repeat(filler1, n_seq, axis=0)
            ))

            filler2 = np.array([*b'-'] * -len_diff, dtype='uint8')
            new_row = np.array([[*seq, *filler2]], dtype='uint8')

            arr = np.vstack((arr, new_row))
            self._collection = arr
        self._header_idx[header] = n_seq

    def __getitem__(self, header: str) -> Sequence:
        idx = self._header_idx[header]
        n_chars = self._collection.shape[1]
        sequence = self._collection[idx] \
            .view(f'S{n_chars}')[0] \
            .decode()
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


class SequenceAnnotation:
    def __init__(
        self,
        sequence_collection: Optional['SequenceCollection'] = None
    ) -> None:
        """[summary]

        Args:
            sequence_collection (Optional['SequenceCollection']): \
                [description]. Defaults to None.
        """
        if sequence_collection:
            sequence_collection.sequence_annotation = self
        self.sequence_collection = sequence_collection
        self._intervals = dict()
        self._gff_headers = list()

    def __getitem__(self, key):
        return self._intervals[key]

    def __setitem__(self, key, value):
        if key in self._intervals:
            raise Exception('Duplicate ID')
        self._intervals[key] = value

    def __iter__(self):
        yield from self._intervals.values()

    @property
    def intervals(self):
        return list(self._intervals.values())

    @classmethod
    def from_gff(
        cls,
        filename: Optional[str] = None,
        string: Optional[str] = None,
        sequence_collection: Optional['SequenceCollection'] = None
    ) -> 'SequenceAnnotation':
        """[summary]

        Args:
            filename ([type], optional): [description]. Defaults to None.
            string ([type], optional): [description]. Defaults to None.
            sequence_collection ([type], optional): [description].
                Defaults to None.

        Returns:
            [type]: [description]
        """
        assert filename or string
        assert not (filename and string)
        sequence_annotation = cls(sequence_collection=sequence_collection)
        header = True
        if filename:
            with open(filename) as filehandle:
                string = filehandle.read()
        for line_number, line in enumerate(string.split('\n')):
            line = line.strip()
            if not line:
                continue
            if line[0] == '#':
                if header:
                    sequence_annotation._gff_headers.append(line)
                continue
            else:
                header = False
            interval = SequenceInterval.from_gff_line(gff_line=line,
                                                      line_number=line_number)
            interval._container = sequence_annotation
            sequence_annotation._intervals[interval.ID] = interval
        for interval in sequence_annotation:
            if interval.parent:
                for parent_ID in interval.parent:
                    try:
                        parent = sequence_annotation[parent_ID]
                    except IndexError:
                        raise IndexError(
                            'Interval {interval.ID} is listing {parent_ID} '
                            'as Parent, but parent could not be found.'
                        )
                    parent._children.append(interval.ID)

        return sequence_annotation

    def to_gff(self) -> str:
        """[summary]

        Returns:
            str: [description]
        """
        gff_lines = [
            interval.to_gff_line() for interval in self._intervals.values()
        ]
        return '\n'.join(gff_lines)

    @classmethod
    def from_json(
        cls,
        filename: Optional[str] = None,
        string: Optional[str] = None,
        sequence_collection: Optional['SequenceCollection'] = None
    ) -> 'SequenceAnnotation':
        """[summary]
        """
        assert filename or string
        assert not (filename and string)
        if filename:
            with open(filename) as filehandle:
                string = filehandle.read()

        sequence_annotation = cls(sequence_collection=sequence_collection)

        gene_dicts = json.loads(string)
        assert isinstance(gene_dicts, list)

        for top_dict in gene_dicts:
            child_dicts = top_dict.pop('children', list())
            top_interval = SequenceInterval.from_dict(interval_dict=top_dict)
            top_interval._container = sequence_annotation
            sequence_annotation._intervals[top_interval.ID] = top_interval
            for child_dict in child_dicts:
                child_interval = SequenceInterval.from_dict(
                    interval_dict=child_dict)
                child_interval._container = sequence_annotation
                sequence_annotation._intervals[child_interval.ID] = \
                    child_interval
        for interval in sequence_annotation:
            if interval.parent:
                for parent_ID in interval.parent:
                    try:
                        parent = sequence_annotation[parent_ID]
                    except IndexError:
                        raise IndexError(
                            'Interval {interval.ID} is listing {parent_ID} '
                            'as Parent, but parent could not be found.'
                        )
                    parent._children.append(interval.ID)
        return sequence_annotation

    def to_json(self, indent: Optional[int] = None) -> str:
        """[summary]

        Returns:
            str: [description]
        """
        interval_dicts = [
            interval.to_dict() for interval in self._intervals.values()
        ]
        return json.dumps(interval_dicts, indent=indent)


class SequenceInterval:
    _predefined_gff3_attributes = (
        'ID', 'name', 'alias', 'parent', 'target', 'gap', 'derives_from',
        'note', 'dbxref', 'ontology_term', 'is_circular'
    )
    _fixed_gff3_fields = (
        'seqid', 'source', 'interval_type', 'start', 'end', 'score', 'strand',
        'phase'
    )

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
        **kwargs
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
        # interval ID
        self.ID = ID

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
            if attr == 'ID':
                continue
            self[attr] = kwargs.get(attr, None)

        # Any additional attributes
        for key, value in kwargs.items():
            self[key] = value

        # Additional fields, used internally
        self._container = container
        if children is None:
            children = []
        self._children = children

    def __repr__(self):
        return (
            f'<SequenceInterval type={self.interval_type} '
            f'ID={self.ID} '
            f'loc={self.seqid}..{self.start}..{self.end}..{self.strand} '
            f'at {hex(id(self))}>'
        )

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    @property
    def attributes(self):
        return {
            attr: self[attr]
            for attr in self.__dict__
            if attr not in self._fixed_gff3_fields  # skip column 1-8 in gff3
            and attr not in ('_children', '_container')  # internal use only
            and self[attr] is not None  # no empty attributes
        }

    @property
    def children(self) -> List['SequenceInterval']:
        return list(self._get_children())

    @classmethod
    def from_gff_line(
        cls,
        gff_line: Optional[str] = None,
        line_number: Optional[int] = None
    ):
        """[summary]

        Args:
            gff_line (Optional[str], optional): [description]. Defaults
                to None.
            line_number (Optional[int], optional): [description]. Defaults
                to None.

        Raises:
            ValueError: [description]
            ValueError: [description]
            ValueError: [description]
            ValueError: [description]
            ValueError: [description]

        Returns:
            [type]: [description]
        """
        gff_parts = gff_line.split('\t')
        assert len(gff_parts) == 9, gff_parts
        seqid, source, interval_type, start, end,\
            score, strand, phase = gff_parts[:8]
        try:
            start = int(start)
            end = int(end)
        except ValueError:
            error = 'GFF start and end fields must be integer'
            if line_number:
                error = f'{error}, gff line {line_number}'
            raise ValueError(error)

        if score != '.':
            try:
                score = float(score)
            except ValueError:
                error = 'GFF score field must be a float'
                if line_number:
                    error = f'{error}, gff line {line_number}'
                raise ValueError(error)

        if strand not in ('+', '-', '.'):
            error = 'GFF strand must be one of "+", "-" or "."'
            if line_number:
                error = f'{error}, gff line {line_number}'
            raise ValueError(error)

        if phase not in ('0', '1', '2', '.'):
            error = 'GFF phase must be one of "0", "1", "2" or "."'
            if line_number:
                error = f'{error}, gff line {line_number}'
            raise ValueError(error)
        elif phase != '.':
            phase = int(phase)

        # Disable phase checking of CDS for now...
        # if interval_type == 'CDS' and phase not in ('0', '1', '2'):
        #     error = 'GFF intervals of type CDS must have phase of\
        #         "0", "1" or "2"'
        #     if line_number:
        #         error = f'{error}, gff line {line_number}'
        #         raise ValueError(error)

        attributes = parse_gff_attribute_string(gff_parts[8])

        ID = attributes.pop('ID', [str(uuid.uuid4())])[0]

        return cls(seqid=seqid, source=source, interval_type=interval_type,
                   start=start, end=end, score=score, strand=strand,
                   phase=phase, ID=ID, **attributes)

    def to_gff_line(self) -> str:
        """[summary]

        Returns:
            str: [description]
        """
        attributes = defaultdict(list)
        attributes.update(self.attributes)
        for attr in self._predefined_gff3_attributes:
            if self[attr] is not None:
                attributes[attr] = self[attr]
        attributes['ID'] = [attributes['ID']]

        return '\t'.join([
            self.seqid, self.source, self.interval_type, str(self.start),
            str(self.end), str(self.score), self.strand, str(self.phase),
            format_gff_attribute_string(attributes)
        ])

    @classmethod
    def from_dict(cls, interval_dict: Dict[str, Any]) -> 'SequenceInterval':
        """[summary]
        Args:
            interval_dict

        Returns:
            [type]: [description]
        """
        attributes = interval_dict.pop('attributes', dict())
        return cls(**interval_dict, **attributes)

    def to_dict(self, include_children: bool = False) -> Dict[str, Any]:
        """[summary]

        Returns:
            Dict[str, Any]: [description]
        """
        attributes = dict(**self.attributes)
        attributes.pop('ID')
        interval_dict = dict(
            ID=self.ID, seqid=self.seqid, source=self.source,
            interval_type=self.interval_type, start=self.start,
            end=self.end, score=self.score, strand=self.strand,
            phase=self.phase, attributes=attributes
        )
        if include_children:
            children = [child.to_dict() for child in self.children[1:]]
            interval_dict['children'] = children
        return interval_dict

    def to_json(
        self,
        include_children: bool = False,
        indent: Optional[int] = None
    ) -> str:
        """[summary]

        Args:
            include_children (bool, optional): [description]. Defaults to \
                False.

        Returns:
            str: [description]
        """
        return json.dumps(
            self.to_dict(include_children=include_children),
            indent=indent
        )

    def _get_children(self, _visited: Optional[set] = None):
        if _visited is None:
            _visited = set()
        if self not in _visited:
            yield self
            _visited.add(self)
        if not self._children:
            return
        for child_ID in self._children:
            child = self._container[child_ID]
            yield from child._get_children(_visited=_visited)


def quote_gff3(attribute_value: str) -> str:
    '''pattern, repl = ENCODE_SPECIAL_CHARACTERS[0]
    quoted_value = re.sub(pattern, repl, attribute_value)
    for pattern, repl in ENCODE_SPECIAL_CHARACTERS[1:]:
        quoted_value = re.sub(pattern, repl, attribute_value)
    return quoted_value'''
    return reduce(
        lambda acc, code: re.sub(code[0], code[1], acc),  # func
        ENCODE_SPECIAL_CHARACTERS,  # iterable
        attribute_value  # initial
    )


def encode_attribute_value(attribute_value: List[str]) -> str:
    """[summary]

    Args:
        attribute_value (List[str]): [description]

    Returns:
        str: [description]
    """
    return ','.join([quote_gff3(v) for v in attribute_value])


def format_gff_attribute_string(attributes: Dict[str, List[str]]) -> str:
    """[summary]

    Args:
        attributes (Dict[str, List[str]]): [description]

    Returns:
        str: [description]
    """
    return ';'.join([
        f'{key}={encode_attribute_value(value)}'
        for key, value in attributes.items()
    ])


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
        attribute_value  # initial
    )


def decode_attribute_value(attribute_value: str) -> List[str]:
    """[summary]

    Args:
        attribute_value (str): [description]

    Returns:
        List[str]: [description]
    """
    return [unquote_gff3(v) for v in attribute_value.split(',')]


def parse_gff_attribute_string(
    gff_attribute_string: str,
    case_sensitive_attribute_keys: bool = False
) -> Dict[str, List[str]]:
    """[summary]
    https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
    See "Column 9: Attributes"
    Args:
        gff_attribute_string ([type]): [description]
    """
    attributes = defaultdict(list)
    for string_part in gff_attribute_string.split(';'):
        if not string_part:
            continue
        try:
            key, value = string_part.split('=', maxsplit=1)
        except Exception as e:
            print(gff_attribute_string, string_part)
            raise Exception(e)
        # The gff spec lists the predefined attribute fields as starting with
        # a capital letter, but we process in lowercase so we don't miss
        # anything from poorly formatted files. When writing to gff we convert
        # back to a capital
        # EXCEPT FOR THE ID ATTRIBUTE, since lowercase id is reserved in python
        if key != 'ID':
            key = key.lower()
        for value_part in decode_attribute_value(value):
            attributes[key].append(value_part)
    return attributes
