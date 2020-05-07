from warnings import warn
from itertools import groupby, chain
from subprocess import Popen, PIPE
from collections import defaultdict
from typing import Iterable, List, Tuple, Dict, Callable
from abc import ABCMeta, abstractmethod, abstractproperty
import uuid
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm


def msa_plot(seq, ax=None, figsize=None) -> Callable:
    """Multiple sequence alignment plot

    Args:
        seq ([type]): [description]
        ax ([type], optional): [description]. Defaults to None.
        figsize ([type], optional): [description]. Defaults to None.

    Returns:
        [type]: [description]
    """
    codes = np.array([*b'AaCcGgNnTt-'], dtype='uint8')
    colors = np.array([
        *['green'] * 2,
        *['red'] * 2,
        *['blue'] * 2,
        *['gray'] * 2,
        *['darkorange'] * 2,
        'black'
    ])
    order = np.argsort(codes)
    cmap = ListedColormap(colors[order])
    norm = BoundaryNorm([0, *codes[order]], ncolors=codes.size)

    if not figsize:
        figsize = (20, 5)
    if not ax:
        fig, ax = plt.subplots(figsize=figsize)

    ax.imshow(norm(seq._collection), cmap=cmap)
    for i in range(seq.n_chars):
        for j in range(seq.n_seqs):
            nuc = seq._collection.view('S1')[j, i].decode()
            ax.text(i, j, nuc, ha='center', va='center', color='w')
    ax.set_yticks(np.arange(seq.n_seqs))
    ax.set_yticklabels(seq.headers)
    return ax


class FastaIter:
    def __init__(self, string: str) -> None:
        """Iterator over fasta formatted sequence strings

        Args:
            string (str): Fasta formatted string
        """
        self._iter = (
            x for _, x in groupby(
                string.strip().split('\n'),
                lambda line: line[0] == '>'
            )
        )

    def __iter__(self) -> Iterable[Tuple[str, str]]:
        """Iterate over header,sequence tuples

        Returns:
            Iterable[Tuple[str, str]]: [description]

        Yields:
            Iterable[Tuple[str, str]]: [description]
        """
        for header in self._iter:
            header = next(header)[1:].strip()
            seq = ''.join(s.strip() for s in next(self._iter))
            yield header, seq

    def __next__(self) -> Tuple[str, str]:
        """Next header and sequence in the iterator

        Returns:
            Tuple[str, str]: [description]
        """
        header = next(next(self._iter))[1:].strip()
        seq = ''.join(s.strip() for s in next(self._iter))
        return header, seq


class SequenceCollection(object, metaclass=ABCMeta):
    @abstractmethod
    def __init__(
        self: 'SequenceCollection',
        sequences: Iterable[Tuple[str, str]] = None,
        sequence_annotation: 'SequenceAnnotation' = None
    ):
        """[summary]

        Args:
            sequences (Iterable[Tuple[str, str]], optional): [description]. \
                 Defaults to None.
            sequence_annotation ([type], optional): [description]. Defaults \
                to None.

        Raises:
            NotImplementedError: [description]
        """
        raise NotImplementedError()

    @abstractmethod
    def __setitem__(self, header: str, seq: str) -> None:
        raise NotImplementedError()

    @abstractmethod
    def __getitem__(self, header: str) -> str:
        raise NotImplementedError()

    @abstractmethod
    def __delitem__(self, header: str) -> None:
        raise NotImplementedError()

    def __iter__(self) -> Iterable[Tuple[str, str]]:
        for header in self.headers:
            yield header, self[header]

    def __next__(self) -> Tuple[str, str]:
        header = next(self.headers)
        return header, self[header]

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
        assert filename or string
        assert not (filename and string)
        sequencecollection = cls()
        if filename:
            with open(filename) as filehandle:
                string = filehandle.read()
        fasta_iter = FastaIter(string)
        for header, seq in fasta_iter:
            sequencecollection[header] = seq
        return sequencecollection

    def to_fasta(self) -> str:
        """Get a fasta-formatted string of the sequence collection

        Returns:
            str: Multi-line fasta-formatted string
        """
        fasta_lines = []
        for header in self.headers:
            fasta_lines.append(f'>{header}')
            fasta_lines.append(self[header])
        return '\n'.join(fasta_lines)


class SequenceList(SequenceCollection):
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
            while new_header in self:
                modifier += 1
                new_header = f'{header}_{modifier}'
            header = new_header
        self._collection[header] = seq

    def __getitem__(self, header: str) -> str:
        return self._collection[header]

    def __delitem__(self, header: str) -> None:
        del self._collection[header]

    @property
    def headers(self) -> List[str]:
        return list(self._collection.keys())

    @property
    def n_seqs(self) -> int:
        return len(self._collection.keys())


class MultipleSequenceAlignment(SequenceCollection):
    def __init__(
        self,
        sequences: Iterable[Tuple[str, str]] = None,
        sequence_annotation: 'SequenceAnnotation' = None
    ) -> None:
        """[summary]

        Args:
            self ([type]): [description]
            sequences (Iterable, optional): [description]. Defaults to None.
            aligned (bool, optional): [description]. Defaults to False.
            sequence_annotation ([type], optional): [description]. Defaults \
                to None.
        """
        self._collection = np.empty((0, 0), dtype='uint8')
        self._header_idx = dict()
        if sequences:
            for header, sequence in sequences:
                self[header] = sequence
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

    def __getitem__(self, header: str) -> str:
        idx = self._header_idx[header]
        n_chars = self._collection.shape[1]
        seq = self._collection[idx] \
            .view(f'S{n_chars}')[0] \
            .decode()
        if not self.aligned:
            seq = seq.rstrip('-')
        return seq

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

    @classmethod
    def from_fasta(
        cls,
        filename: str = None,
        string: str = None,
        aligned: bool = False
    ) -> 'SequenceCollection':
        """Parse a fasta formatted string into a SequenceCollection object

        Keyword Arguments:
            filename {String} -- filename string (default: {None})
            string {String} -- fasta formatted string (default: {None})

        Returns:
            SequenceCollection -- SequenceCollection instance
        """
        assert filename or string
        assert not (filename and string)
        sequencecollection = cls(aligned=aligned)
        if filename:
            with open(filename) as filehandle:
                string = filehandle.read()
        fasta_iter = FastaIter(string)
        for header, seq in fasta_iter:
            sequencecollection[header] = seq
        return sequencecollection

    def to_fasta(self) -> str:
        """[summary]

        Returns:
            [type] -- [description]
        """
        fasta_lines = []
        for header in self.headers:
            fasta_lines.append(f'>{header}')
            fasta_lines.append(self[header])
        return '\n'.join(fasta_lines)

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
        return SequenceCollection.from_fasta(
            string=aligned_fasta,
            aligned=True
        )


class SequenceAnnotation:
    def __init__(
        self,
        sequence_collection: 'SequenceCollection' = None
    ):
        """[summary]

        Args:
            sequence_collection ([type], optional): [description].
                Defaults to None.
        """
        if sequence_collection:
            sequence_collection.sequence_annotation = self
        self.sequence_collection = sequence_collection
        self._intervals = dict()
        self._index = dict()

    def __getitem__(self, key):
        return self._intervals[key]

    def __setitem__(self, key, value):
        if key in self._intervals:
            raise Exception('Duplicate ID')
        self._intervals[key] = value

    def __iter__(self):
        yield from self._intervals.values()

    @classmethod
    def from_gff(cls, filename=None, string=None,
                 sequence_collection=None):
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
        if filename:
            with open(filename) as filehandle:
                string = filehandle.read()
        for line_number, line in enumerate(string.split('\n')):
            line = line.strip()
            if not line:
                continue
            if line[0] == '#':
                continue
            interval = SequenceInterval.from_gff_line(gff_line=line,
                                                      line_number=line_number)
            sequence_annotation._intervals[interval.ID] = interval

        for interval in sequence_annotation:
            if interval.parent:
                for parent_ID in interval.parent:
                    parent = sequence_annotation[parent_ID]
                    parent.children.append(interval.ID)

        return sequence_annotation


class SequenceInterval:
    _predefined_attributes = ('ID', 'name', 'alias', 'parent', 'target',
                              'gap', 'derives_from', 'note', 'dbxref',
                              'ontology_term', 'is_circular')

    def __init__(self, seqid=None, source=None, interval_type=None,
                 start=None, end=None, score=None, strand=None, phase=None,
                 attributes=None, children=None, **kwargs):
        """[summary]

        Args:
            seqid ([type], optional): [description]. Defaults to None.
            source ([type], optional): [description]. Defaults to None.
            interval_type ([type], optional): [description]. Defaults to None.
            start ([type], optional): [description]. Defaults to None.
            end ([type], optional): [description]. Defaults to None.
            score ([type], optional): [description]. Defaults to None.
            strand ([type], optional): [description]. Defaults to None.
            phase ([type], optional): [description]. Defaults to None.
            attributes ([type], optional): [description]. Defaults to None.
            children ([type], optional): [description]. Defaults to None.

        Raises:
            NotImplementedError: [description]
        """
        # Standard gff fields
        self.seqid = seqid
        self.source = source
        self.interval_type = interval_type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes

        # Attributes with predefined meanings in the gff spec
        for attr in self._predefined_attributes:
            self[attr] = kwargs.get(attr, None)

        for key in kwargs.keys():
            if key not in self._predefined_attributes:
                raise NotImplementedError(f'{key} is not a valid attribute')

        # Additional fields, used internally
        if children is None:
            children = []
        self.children = children

    def __repr__(self):
        return (
            f'<{self.interval_type} SequenceInterval '
            f'{self.ID} '
            f'{self.start}..{self.end}..{self.strand} '
            f'at {hex(id(self))}'
        )

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    @classmethod
    def from_gff_line(cls, gff_line=None, line_number=None):
        """[summary]

        Args:
            gff_line ([type], optional): [description]. Defaults to None.
            line_number ([type], optional): [description]. Defaults to None.

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

        # Disable phase checking of CDS for now...
        # if interval_type == 'CDS' and phase not in ('0', '1', '2'):
        #     error = 'GFF intervals of type CDS must have phase of\
        #         "0", "1" or "2"'
        #     if line_number:
        #         error = f'{error}, gff line {line_number}'
        #         raise ValueError(error)

        interval = cls(seqid=seqid, source=source, interval_type=interval_type,
                       start=start, end=end, score=score, strand=strand,
                       phase=phase)

        attributes = parse_gff_attribute_string(gff_parts[8])

        for attr_key in cls._predefined_attributes:
            attr_value = attributes[attr_key]

            # ID must always be a single value, use random number if not set
            if attr_key == 'ID':
                if not attr_value:
                    attr_value = str(uuid.uuid4())
                else:
                    attr_value = attr_value[0]

            if attr_value:
                interval[attr_key] = attr_value
            del attributes[attr_key]

        interval.attributes = attributes

        return interval


def parse_gff_attribute_string(gff_attribute_string):
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
        for value_part in value.split(','):
            attributes[key].append(value_part)
    return attributes
