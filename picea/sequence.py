from warnings import warn
from itertools import groupby

class SequenceCollection(dict):
    def __init__(self, sequences = None, is_aligned = False):
        self.is_aligned = is_aligned
        self._collection = dict()
        if sequences:
            for header,sequence in sequences:
                self[header] = sequence
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
    def __getitem__(self, key):
        return self._collection[key]
    @classmethod
    def from_fasta(cls, filename = None, string = None):
        assert(filename or string)
        assert(not filename and string)
        sequencecollection = cls()
        if filename:
            with open(filename) as filehandle:
                string = filehandle.read()
        fasta_iter = (x for _,x in groupby(string.split(), 
            lambda line: line[0] == '>'))
        for header in fasta_iter:
            header = next(header)[1:].strip()
            seq = ''.join(s.strip() for s in next(fasta_iter))
            sequencecollection[header] = seq
        return sequencecollection
    def to_fasta(self):
        fasta_lines = []
        for header,seq in self._collection.items():
            fasta_lines.append(f'>{header}')
            fasta_lines.append(seq)
        return '\n'.join(fasta_lines)

