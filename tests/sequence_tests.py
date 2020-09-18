from unittest import TestCase

from picea import Sequence, SequenceReader, SequenceList,\
    MultipleSequenceAlignment, alphabets


class SequenceTests(TestCase):
    def setUp(self):
        self.fasta = '>A\nABC\n>B\nDEF'
        self.json = (
            '[{"header":"A","sequence":"ABC"},'
            '{"header":"B","sequence":"DEF"}]'
        )
        self.rename_func = lambda x: f'{x}.test'

    def test_empty_init_sequencereader(self):
        self.assertRaises(AssertionError, SequenceReader)

    def test_fasta_sequencereader(self):
        for _ in SequenceReader(string=self.fasta, filetype='fasta'):
            pass

    def test_json_sequencereader(self):
        for _ in SequenceReader(string=self.json, filetype='json'):
            pass

    def test_empty_init_sequence(self):
        Sequence()

    def test_sequence_detect_dna(self):
        s = Sequence('test', 'ACGATCGACTCGAACT')
        self.assertEqual(s.alphabet, alphabets.DNA)

    def test_sequence_detect_aminoacid(self):
        s = Sequence('test', 'KUDHLSKJSPOIJKMSLKM')
        self.assertEqual(s.alphabet, alphabets.AminoAcid)

    def test_empty_init_sequencelist(self):
        SequenceList()

    def test_empty_init_msa(self):
        MultipleSequenceAlignment()

    def test_fasta_parsing_sequencelist(self):
        SequenceList.from_fasta(string=self.fasta)

    def test_fasta_parsing_msa(self):
        MultipleSequenceAlignment.from_fasta(string=self.fasta)

    def test_fasta_input_output_sequencelist(self):
        seq = SequenceList.from_fasta(string=self.fasta)
        self.assertEqual(seq.to_fasta(), self.fasta)

    def test_fasta_input_output_msa(self):
        seq = MultipleSequenceAlignment.from_fasta(string=self.fasta)
        self.assertEqual(seq.to_fasta(), self.fasta)

    def test_json_parsing_sequencelist(self):
        SequenceList.from_json(string=self.json)

    def test_json_parsing_msa(self):
        MultipleSequenceAlignment.from_json(string=self.json)

    def test_trailing_newline_sequencelist(self):
        fasta = f'{self.fasta}\n'
        seq = SequenceList.from_fasta(string=fasta)
        self.assertEqual(seq.to_fasta(), fasta[:-1])

    def test_trailing_newline_msa(self):
        fasta = f'{self.fasta}\n'
        seq = MultipleSequenceAlignment.from_fasta(string=fasta)
        self.assertEqual(seq.to_fasta(), fasta[:-1])

    def test_sequence_iter_sequencelist(self):
        for _ in SequenceList.from_fasta(string=self.fasta):
            pass

    def test_sequence_iter_msa(self):
        for _ in MultipleSequenceAlignment.from_fasta(string=self.fasta):
            pass

    def test_sequencelist_pop(self):
        seq_list = SequenceList.from_fasta(string=self.fasta)
        pop_seq = seq_list.pop('A')
        self.assertEqual(pop_seq.header, 'A')
        self.assertEqual(pop_seq.sequence, 'ABC')
        self.assertNotIn('A', seq_list.headers)
        self.assertNotIn('ABC', seq_list.sequences)

    def test_msa_pop(self):
        msa = MultipleSequenceAlignment.from_fasta(string=self.fasta)
        pop_seq = msa.pop('A')
        self.assertEqual(pop_seq.header, 'A')
        self.assertEqual(pop_seq.sequence, 'ABC')
        self.assertNotIn('A', msa.headers)
        self.assertNotIn('ABC', msa.sequences)

    def test_seqlist_batch_rename(self):
        seq_list = SequenceList.from_fasta(string=self.fasta)
        seq_list.batch_rename(self.rename_func)
        self.assertEqual(seq_list.headers, ['A.test', 'B.test'])

    def test_msa_batch_rename(self):
        msa = MultipleSequenceAlignment.from_fasta(string=self.fasta)
        msa.batch_rename(self.rename_func)
        self.assertEqual(msa.headers, ['A.test', 'B.test'])
