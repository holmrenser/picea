from unittest import TestCase

from picea import (
    Sequence, SequenceReader, BatchSequenceReader,
    SequenceCollection, MultipleSequenceAlignment, alphabets, SequenceAnnotation
)


class AnnotationTests(TestCase):
    maxDiff = None

    def setUp(self):
        # gff string taken from https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
        # modified to have unique CDS IDs
        self.gff3 = (
            # '##gff-version 3.1.26\n'
            # '##sequence-region ctg123 1 1497228\n'
            'ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001;Name=EDEN\n'
            'ctg123\t.\tTF_binding_site\t1000\t1012\t.\t+\t.\tID=tfbs00001;Parent=gene00001\n' # noqa
            'ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00001;Parent=gene00001;Name=EDEN.1\n' # noqa
            'ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00002;Parent=gene00001;Name=EDEN.2\n' # noqa
            'ctg123\t.\tmRNA\t1300\t9000\t.\t+\t.\tID=mRNA00003;Parent=gene00001;Name=EDEN.3\n' # noqa
            'ctg123\t.\texon\t1300\t1500\t.\t+\t.\tID=exon00001;Parent=mRNA00003\n'
            'ctg123\t.\texon\t1050\t1500\t.\t+\t.\tID=exon00002;Parent=mRNA00001,mRNA00002\n' # noqa
            'ctg123\t.\texon\t3000\t3902\t.\t+\t.\tID=exon00003;Parent=mRNA00001,mRNA00003\n' # noqa
            'ctg123\t.\texon\t5000\t5500\t.\t+\t.\tID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003\n' # noqa
            'ctg123\t.\texon\t7000\t9000\t.\t+\t.\tID=exon00005;Parent=mRNA00001,mRNA00002,mRNA00003\n' # noqa
            'ctg123\t.\tCDS\t1201\t1500\t.\t+\t0\tID=cds00001.1;Parent=mRNA00001;Name=edenprotein.1\n' # noqa
            'ctg123\t.\tCDS\t3000\t3902\t.\t+\t0\tID=cds00001.2;Parent=mRNA00001;Name=edenprotein.1\n' # noqa
            'ctg123\t.\tCDS\t5000\t5500\t.\t+\t0\tID=cds00001.3;Parent=mRNA00001;Name=edenprotein.1\n' # noqa
            'ctg123\t.\tCDS\t7000\t7600\t.\t+\t0\tID=cds00001.4;Parent=mRNA00001;Name=edenprotein.1\n' # noqa
            'ctg123\t.\tCDS\t1201\t1500\t.\t+\t0\tID=cds00002.1;Parent=mRNA00002;Name=edenprotein.2\n' # noqa
            'ctg123\t.\tCDS\t5000\t5500\t.\t+\t0\tID=cds00002.2;Parent=mRNA00002;Name=edenprotein.2\n' # noqa
            'ctg123\t.\tCDS\t7000\t7600\t.\t+\t0\tID=cds00002.3;Parent=mRNA00002;Name=edenprotein.2\n' # noqa
            'ctg123\t.\tCDS\t3301\t3902\t.\t+\t0\tID=cds00003.1;Parent=mRNA00003;Name=edenprotein.3\n' # noqa
            'ctg123\t.\tCDS\t5000\t5500\t.\t+\t1\tID=cds00003.2;Parent=mRNA00003;Name=edenprotein.3\n' # noqa
            'ctg123\t.\tCDS\t7000\t7600\t.\t+\t1\tID=cds00003.3;Parent=mRNA00003;Name=edenprotein.3\n' # noqa
            'ctg123\t.\tCDS\t3391\t3902\t.\t+\t0\tID=cds00004.1;Parent=mRNA00003;Name=edenprotein.4\n' # noqa
            'ctg123\t.\tCDS\t5000\t5500\t.\t+\t1\tID=cds00004.2;Parent=mRNA00003;Name=edenprotein.4\n' # noqa
            'ctg123\t.\tCDS\t7000\t7600\t.\t+\t1\tID=cds00004.3;Parent=mRNA00003;Name=edenprotein.4\n' # noqa
        )

    def test_empty_init_sequenceannotation(self):
        SequenceAnnotation()

    def test_read_gff_string(self):
        SequenceAnnotation.from_gff(string=self.gff3)

    def test_get_children(self):
        ann = SequenceAnnotation.from_gff(string=self.gff3)
        ann['gene00001'].children

    def test_input_output_annotation(self):
        ann = SequenceAnnotation.from_gff(string=self.gff3)
        gff3 = ann.to_gff()
        self.assertEqual(gff3, self.gff3)


class SequenceTests(TestCase):
    def setUp(self):
        self.fasta = '>A\nABC\n>B\nDEF'
        self.big_fasta = '>A\nABC\n>B\nDEF\n>C\nGHI\n>D\nJKL'
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

    def test_batchsequencereader(self):
        n_batches = 0
        for batch in BatchSequenceReader(
            string=self.big_fasta,
            filetype='fasta',
            batchsize=2
        ):
            n_batches += 1
            self.assertEqual(len(batch), 2)
        self.assertEqual(n_batches, 2)

    def test_empty_init_sequence(self):
        Sequence()

    def test_sequence_detect_dna(self):
        s = Sequence('test', 'ACGATCGACTCGAACT')
        self.assertEqual(s.alphabet, alphabets.DNA)

    def test_sequence_detect_aminoacid(self):
        s = Sequence('test', 'KUDHLSKJSPOIJKMSLKM')
        self.assertEqual(s.alphabet, alphabets.AminoAcid)

    def test_empty_init_sequencecollection(self):
        SequenceCollection()

    def test_empty_init_msa(self):
        MultipleSequenceAlignment()

    def test_fasta_parsing_sequencecollection(self):
        SequenceCollection.from_fasta(string=self.fasta)

    def test_fasta_parsing_msa(self):
        MultipleSequenceAlignment.from_fasta(string=self.fasta)

    def test_fasta_input_output_sequencecollection(self):
        seq = SequenceCollection.from_fasta(string=self.fasta)
        self.assertEqual(seq.to_fasta(), self.fasta)

    def test_fasta_input_output_msa(self):
        seq = MultipleSequenceAlignment.from_fasta(string=self.fasta)
        self.assertEqual(seq.to_fasta(), self.fasta)

    def test_json_parsing_sequencecollection(self):
        SequenceCollection.from_json(string=self.json)

    def test_json_parsing_msa(self):
        MultipleSequenceAlignment.from_json(string=self.json)

    def test_trailing_newline_sequencecollection(self):
        fasta = f'{self.fasta}\n'
        seq = SequenceCollection.from_fasta(string=fasta)
        self.assertEqual(seq.to_fasta(), fasta[:-1])

    def test_trailing_newline_msa(self):
        fasta = f'{self.fasta}\n'
        seq = MultipleSequenceAlignment.from_fasta(string=fasta)
        self.assertEqual(seq.to_fasta(), fasta[:-1])

    def test_sequence_iter_sequencecollection(self):
        for _ in SequenceCollection.from_fasta(string=self.fasta):
            pass

    def test_sequence_iter_msa(self):
        for _ in MultipleSequenceAlignment.from_fasta(string=self.fasta):
            pass

    def test_sequencecollection_pop(self):
        seq_col = SequenceCollection.from_fasta(string=self.fasta)
        pop_seq = seq_col.pop('A')
        self.assertEqual(pop_seq.header, 'A')
        self.assertEqual(pop_seq.sequence, 'ABC')
        self.assertNotIn('A', seq_col.headers)
        self.assertNotIn('ABC', seq_col.sequences)

    def test_msa_pop(self):
        msa = MultipleSequenceAlignment.from_fasta(string=self.fasta)
        pop_seq = msa.pop('A')
        self.assertEqual(pop_seq.header, 'A')
        self.assertEqual(pop_seq.sequence, 'ABC')
        self.assertNotIn('A', msa.headers)
        self.assertNotIn('ABC', msa.sequences)

    def test_seqcol_batch_rename(self):
        seq_col = SequenceCollection.from_fasta(string=self.fasta)
        seq_col.batch_rename(self.rename_func)
        self.assertEqual(seq_col.headers, ['A.test', 'B.test'])

    def test_msa_batch_rename(self):
        msa = MultipleSequenceAlignment.from_fasta(string=self.fasta)
        msa.batch_rename(self.rename_func)
        self.assertEqual(msa.headers, ['A.test', 'B.test'])

    def test_iloc(self):
        seq_col = SequenceCollection.from_fasta(string=self.fasta)
        sub = seq_col.iloc[0]
        self.assertEqual(sub.headers, ['A'])
        sub_multiple = seq_col.iloc[[0, 1]]
        self.assertEqual(sub_multiple.headers, ['A', 'B'])
        sub_slice1 = seq_col.iloc[0:2]
        self.assertEqual(sub_slice1.headers, ['A', 'B'])
        sub_slice2 = seq_col.iloc[1:]
        self.assertEqual(sub_slice2.headers, ['B'])
        sub_slice3 = seq_col.iloc[:1]
        self.assertEqual(sub_slice3.headers, ['A'])
        with self.assertRaises(TypeError):
            seq_col.iloc['A']

    def test_len(self):
        seq_col = SequenceCollection.from_fasta(string=self.fasta)
        self.assertEqual(len(seq_col), 2)

    def test_from_iter(self):
        seq_col = SequenceCollection.from_fasta(string=self.fasta)
        seq_col2 = SequenceCollection.from_sequence_iter(seq_col)
        self.assertEqual(seq_col.headers, seq_col2.headers)
