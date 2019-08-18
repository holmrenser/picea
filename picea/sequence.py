from warnings import warn
from itertools import groupby, chain
from subprocess import Popen, PIPE
from collections import defaultdict
import tempfile
import uuid


class SequenceCollection(dict):
    def __init__(self, sequences=None, is_aligned=False,
                 sequence_annotation=None):
        """[summary]

        Arguments:
            dict {[type]} -- [description]

        Keyword Arguments:
            sequences {[type]} -- [description] (default: {None})
            is_aligned {bool} -- [description] (default: {False})
        """
        self.is_aligned = is_aligned
        self._collection = dict()
        if sequences:
            for header, sequence in sequences:
                self[header] = sequence
        self.sequence_annotation = sequence_annotation

    def __repr__(self):
        return self._collection.__repr__()

    def __setitem__(self, key, value):
        if key in self:
            warn(f'Turning duplicate key "{key}" into unique key')
            new_key = key
            modifier = 0
            while new_key in self:
                modifier += 1
                new_key = f'{key}_{modifier}'
            key = new_key
        self._collection[key] = value

    @property
    def headers(self):
        return self._collection.keys()

    @property
    def sequences(self):
        return self._collection.values()

    def __getitem__(self, key):
        return self._collection[key]

    @classmethod
    def from_fasta(cls, filename=None, string=None):
        """Parse a fasta formatted string into a SequenceCollection object

        Keyword Arguments:
            filename {String} -- filename string (default: {None})
            string {String} -- fasta formatted string (default: {None})

        Returns:
            SequenceCollection -- SequenceCollection instance
        """
        assert(filename or string)
        assert(not filename and string)
        sequencecollection = cls()
        if filename:
            with open(filename) as filehandle:
                string = filehandle.read()
        fasta_iter = (
            x for _, x in groupby(
                string.split(),
                lambda line: line[0] == '>'
            )
        )
        for header in fasta_iter:
            header = next(header)[1:].strip()
            seq = ''.join(s.strip() for s in next(fasta_iter))
            sequencecollection[header] = seq
        return sequencecollection

    def to_fasta(self):
        """[summary]

        Returns:
            [type] -- [description]
        """
        fasta_lines = []
        for header, seq in self._collection.items():
            fasta_lines.append(f'>{header}')
            fasta_lines.append(seq)
        return '\n'.join(fasta_lines)

    def align(self, method='mafft', method_kwargs=dict()):
        fasta = self.to_fasta()
        command = [method] + list(chain.from_iterable(method_kwargs.items()))
        process = Popen(
            command,
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE
        )
        with tempfile.NamedTemporaryFile() as temp:
            temp.write(fasta)
            stdout, stderr = process.communicate(input=temp.name)
            print(stdout)


class SequenceAnnotation():
    def __init__(self, sequence_collection=None):
        self.sequence_collection = sequence_collection
        self._intervals = dict()
        self._index = dict()

    @classmethod
    def from_gff(cls, filename=None, string=None,
                 sequence_collection=None):
        assert(filename or string)
        assert(not filename and string)
        sequence_annotation = cls(sequence_collection=sequence_collection)
        if filename:
            with open(filename) as filehandle:
                string = filehandle.read()
        for line_number, line in enumerate(string.split('\n')):
            interval = Interval.from_gff_line(gff_line=line,
                                              line_number=line_number)
            sequence_annotation._intervals[interval.ID] = interval

        for interval in sequence_annotation._intervals.values():
            for parentID in interval.parent:
                parent = sequence_annotation._intervals[parentID]
                parent.children.append(interval.ID)

        return sequence_annotation


class Interval():
    _predefined_attributes = ('ID', 'name', 'alias', 'parent', 'target',
                              'gap', 'derives_from', 'note', 'dbxref',
                              'ontology_term', 'is_circular')

    def __init__(self, seqid=None, source=None, interval_type=None,
                 start=None, end=None, score=None, strand=None, phase=None,
                 attributes=None, children=[], **kwargs):
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
        self.children = children

    @classmethod
    def from_gff_line(cls, gff_line=None, line_number=None):
        gff_parts = gff_line.split('\t')
        assert len(gff_parts) == 9
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

        if interval_type == 'CDS' and phase not in ('0', '1', '2'):
            error = 'GFF intervals of type CDS must have phase of\
                "0", "1" or "2"'
            if line_number:
                error = f'{error}, gff line {line_number}'
                raise ValueError(error)

        interval = cls(seqid=seqid, source=source, interval_type=interval_type,
                       start=start, end=end, score=score, strand=strand,
                       phase=phase)

        attributes = parse_gff_attribute_string(gff_parts[8])

        for attr_key in cls._predefined_attributes:
            attr_value = attributes[attr_key]
            if attr_value:
                interval[attr_key] = attr_value
            del attributes[attr_key]

        interval.attributes = attributes

        if not interval.ID:
            interval.ID = uuid.uuid4()

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
        key, value = string_part.split('=')
        # The gff spec lists the predefined attribute fields as starting with
        # a capital letter, but we process in lowercase so we don't miss
        # anything from poorly formatted files. When writing to gff we convert
        # back to a capital
        key = key.lower()
        for value_part in value.split(','):
            attributes[key].append(value_part)
    return attributes
