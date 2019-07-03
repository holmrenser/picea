from picea import SequenceCollection


def test_empty_init():
    SequenceCollection()


def test_parsing():
    fasta = '>A\nABC\n>B\nDEF'
    SequenceCollection.from_fasta(string=fasta)


def test_input_output():
    fasta = '>A\nABC\n>B\nDEF'
    seq = SequenceCollection.from_fasta(string=fasta)
    assert seq.to_fasta() == fasta


'''
def test_align():
    fasta = '>A\nABC\n>B\nDEF'
    seq = SequenceCollection.from_fasta(string = fasta)
    seq.align()
'''
