from unittest import TestCase

from picea import SequenceAnnotation, SequenceInterval


class IntervalTests(TestCase):
    def setUp(self) -> None:
        self.gff3_line = "ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001;Name=EDEN"

    def test_empty_init(self):
        SequenceInterval()

    def test_from_gff_line(self):
        SequenceInterval.from_gff_line(self.gff3_line)

    def test_input_output(self):
        interval = SequenceInterval.from_gff_line(self.gff3_line)
        gff3_line = interval.to_gff_line()
        self.assertEqual(gff3_line, self.gff3_line)


class AnnotationTests(TestCase):
    maxDiff = None

    def setUp(self):
        # gff string taken from
        # https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
        # modified to have unique CDS IDs
        self.gff3 = (
            # '##gff-version 3.1.26\n'
            # '##sequence-region ctg123 1 1497228\n'
            "ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001;Name=EDEN\n"
            "ctg123\t.\tTF_binding_site\t1000\t1012\t.\t+\t.\tID=tfbs00001;Parent=gene00001\n"  # noqa
            "ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00001;Parent=gene00001;Name=EDEN.1\n"  # noqa
            "ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00002;Parent=gene00001;Name=EDEN.2\n"  # noqa
            "ctg123\t.\tmRNA\t1300\t9000\t.\t+\t.\tID=mRNA00003;Parent=gene00001;Name=EDEN.3\n"  # noqa
            "ctg123\t.\texon\t1300\t1500\t.\t+\t.\tID=exon00001;Parent=mRNA00003\n"
            "ctg123\t.\texon\t1050\t1500\t.\t+\t.\tID=exon00002;Parent=mRNA00001,mRNA00002\n"  # noqa
            "ctg123\t.\texon\t3000\t3902\t.\t+\t.\tID=exon00003;Parent=mRNA00001,mRNA00003\n"  # noqa
            "ctg123\t.\texon\t5000\t5500\t.\t+\t.\tID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003\n"  # noqa
            "ctg123\t.\texon\t7000\t9000\t.\t+\t.\tID=exon00005;Parent=mRNA00001,mRNA00002,mRNA00003\n"  # noqa
            "ctg123\t.\tCDS\t1201\t1500\t.\t+\t0\tID=cds00001.1;Parent=mRNA00001;Name=edenprotein.1\n"  # noqa
            "ctg123\t.\tCDS\t3000\t3902\t.\t+\t0\tID=cds00001.2;Parent=mRNA00001;Name=edenprotein.1\n"  # noqa
            "ctg123\t.\tCDS\t5000\t5500\t.\t+\t0\tID=cds00001.3;Parent=mRNA00001;Name=edenprotein.1\n"  # noqa
            "ctg123\t.\tCDS\t7000\t7600\t.\t+\t0\tID=cds00001.4;Parent=mRNA00001;Name=edenprotein.1\n"  # noqa
            "ctg123\t.\tCDS\t1201\t1500\t.\t+\t0\tID=cds00002.1;Parent=mRNA00002;Name=edenprotein.2\n"  # noqa
            "ctg123\t.\tCDS\t5000\t5500\t.\t+\t0\tID=cds00002.2;Parent=mRNA00002;Name=edenprotein.2\n"  # noqa
            "ctg123\t.\tCDS\t7000\t7600\t.\t+\t0\tID=cds00002.3;Parent=mRNA00002;Name=edenprotein.2\n"  # noqa
            "ctg123\t.\tCDS\t3301\t3902\t.\t+\t0\tID=cds00003.1;Parent=mRNA00003;Name=edenprotein.3\n"  # noqa
            "ctg123\t.\tCDS\t5000\t5500\t.\t+\t1\tID=cds00003.2;Parent=mRNA00003;Name=edenprotein.3\n"  # noqa
            "ctg123\t.\tCDS\t7000\t7600\t.\t+\t1\tID=cds00003.3;Parent=mRNA00003;Name=edenprotein.3\n"  # noqa
            "ctg123\t.\tCDS\t3391\t3902\t.\t+\t0\tID=cds00004.1;Parent=mRNA00003;Name=edenprotein.4\n"  # noqa
            "ctg123\t.\tCDS\t5000\t5500\t.\t+\t1\tID=cds00004.2;Parent=mRNA00003;Name=edenprotein.4\n"  # noqa
            "ctg123\t.\tCDS\t7000\t7600\t.\t+\t1\tID=cds00004.3;Parent=mRNA00003;Name=edenprotein.4\n"  # noqa
        )

    def test_empty_init_sequenceannotation(self):
        SequenceAnnotation()

    def test_read_gff_string(self):
        SequenceAnnotation.from_gff(string=self.gff3)

    def test_get_children(self):
        ann = SequenceAnnotation.from_gff(string=self.gff3)
        self.assertEqual(22, len(ann["gene00001"].children))

    def test_input_output_annotation(self):
        ann = SequenceAnnotation.from_gff(string=self.gff3)
        gff3 = ann.to_gff()
        self.assertEqual(gff3, self.gff3)
