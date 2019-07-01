from warnings import warn
from itertools import groupby, chain
from subprocess import Popen, PIPE
import tempfile

class SequenceCollection(dict):
    def __init__(self, sequences = None, is_aligned = False):
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
        """[summary]
        
        Keyword Arguments:
            filename {[type]} -- [description] (default: {None})
            string {[type]} -- [description] (default: {None})
        
        Returns:
            [type] -- [description]
        """
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
        """[summary]
        
        Returns:
            [type] -- [description]
        """
        fasta_lines = []
        for header,seq in self._collection.items():
            fasta_lines.append(f'>{header}')
            fasta_lines.append(seq)
        return '\n'.join(fasta_lines)
    def align(self, method = 'mafft', method_kwargs = dict()):
        fasta = self.to_fasta()
        command = [method] + list(chain.from_iterable(method_kwargs.items()))
        process = Popen(command, stdin = PIPE, stdout = PIPE, stderr = PIPE)
        with tempfile.NamedTemporaryFile() as temp:
            temp.write(fasta)
            stdout, stderr = process.communicate(input = temp.name)
            print(stdout)


