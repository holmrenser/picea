from picea import SequenceList


def test_empty_init():
    SequenceList()


def test_parsing():
    fasta = '>A\nABC\n>B\nDEF'
    SequenceList.from_fasta(string=fasta)


def test_input_output():
    fasta = '>A\nABC\n>B\nDEF'
    seq = SequenceList.from_fasta(string=fasta)
    assert seq.to_fasta() == fasta


'''
def test_align():
    fasta = '>A\nABC\n>B\nDEF'
    seq = SequenceCollection.from_fasta(string = fasta)
    seq.align()
'''
