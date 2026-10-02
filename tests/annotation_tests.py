from unittest import TestCase

from picea import SequenceAnnotation, SequenceInterval


class IntervalTests(TestCase):
    def setUp(self) -> None:
        self.gff3_line = "ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001;Name=EDEN"

    def test_empty_init(self):
        SequenceInterval()

    def test_from_gff_line(self):
        SequenceInterval.from_gff_line(self.gff3_line)

    def test_to_gff_line(self):
        interval = SequenceInterval.from_gff_line(self.gff3_line)
        gff3_line = interval.to_gff_line()
        self.assertEqual(gff3_line, self.gff3_line)

    def test_to_gtf_line(self):
        interval = SequenceInterval.from_gff_line(self.gff3_line)
        interval.to_gtf_line()


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

        # gff string from aparport11 A. thaliana annotation
        self.real_gff3 = (
            "1\taraport11\tgene\t3631\t5899\t.\t+\t.\tID=gene:AT1G01010;Name=NAC001;biotype=protein_coding;description=NAC domain containing protein 1 [Source:NCBI gene (formerly Entrezgene) Acc:839580];gene_id=AT1G01010;logic_name=araport11\n"
            "1\taraport11\tmRNA\t3631\t5899\t.\t+\t.\tID=transcript:AT1G01010.1;Parent=gene:AT1G01010;Name=NAC001-201;biotype=protein_coding;tag=Ensembl_canonical;transcript_id=AT1G01010.1\n"
            "1\taraport11\tfive_prime_UTR\t3631\t3759\t.\t+\t.\tID=a0c9aeaa-d982-4272-99c3-8eceec743efd;Parent=transcript:AT1G01010.1\n"
            "1\taraport11\texon\t3631\t3913\t.\t+\t.\tID=930c4181-a5c3-4641-9b37-6208e71f0d8c;Parent=transcript:AT1G01010.1;Name=AT1G01010.1.exon1;constitutive=1;ensembl_end_phase=1;ensembl_phase=-1;exon_id=AT1G01010.1.exon1;rank=1\n"
            "1\taraport11\tCDS\t3760\t3913\t.\t+\t0\tID=CDS:AT1G01010.1;Parent=transcript:AT1G01010.1;protein_id=AT1G01010.1\n"
            "1\taraport11\texon\t3996\t4276\t.\t+\t.\tID=de24200a-5161-48d1-b0fa-be06d5d7fd8e;Parent=transcript:AT1G01010.1;Name=AT1G01010.1.exon2;constitutive=1;ensembl_end_phase=0;ensembl_phase=1;exon_id=AT1G01010.1.exon2;rank=2\n"
            "1\taraport11\tCDS\t3996\t4276\t.\t+\t2\tID=CDS:AT1G01010.1_1;Parent=transcript:AT1G01010.1;protein_id=AT1G01010.1\n"
            "1\taraport11\texon\t4486\t4605\t.\t+\t.\tID=e306d155-dc8e-4ff1-bd13-290b378df902;Parent=transcript:AT1G01010.1;Name=AT1G01010.1.exon3;constitutive=1;ensembl_end_phase=0;ensembl_phase=0;exon_id=AT1G01010.1.exon3;rank=3\n"
            "1\taraport11\tCDS\t4486\t4605\t.\t+\t0\tID=CDS:AT1G01010.1_2;Parent=transcript:AT1G01010.1;protein_id=AT1G01010.1\n"
            "1\taraport11\texon\t4706\t5095\t.\t+\t.\tID=e73a9dc0-896b-470a-960d-1f9ff15e75d3;Parent=transcript:AT1G01010.1;Name=AT1G01010.1.exon4;constitutive=1;ensembl_end_phase=0;ensembl_phase=0;exon_id=AT1G01010.1.exon4;rank=4\n"
            "1\taraport11\tCDS\t4706\t5095\t.\t+\t0\tID=CDS:AT1G01010.1_3;Parent=transcript:AT1G01010.1;protein_id=AT1G01010.1\n"
            "1\taraport11\texon\t5174\t5326\t.\t+\t.\tID=45bf7a3b-ef3d-4217-8104-c11d7702047e;Parent=transcript:AT1G01010.1;Name=AT1G01010.1.exon5;constitutive=1;ensembl_end_phase=0;ensembl_phase=0;exon_id=AT1G01010.1.exon5;rank=5\n"
            "1\taraport11\tCDS\t5174\t5326\t.\t+\t0\tID=CDS:AT1G01010.1_4;Parent=transcript:AT1G01010.1;protein_id=AT1G01010.1\n"
            "1\taraport11\tCDS\t5439\t5630\t.\t+\t0\tID=CDS:AT1G01010.1_5;Parent=transcript:AT1G01010.1;protein_id=AT1G01010.1\n"
            "1\taraport11\texon\t5439\t5899\t.\t+\t.\tID=22a1667f-55ba-49a6-86a4-237d2da57007;Parent=transcript:AT1G01010.1;Name=AT1G01010.1.exon6;constitutive=1;ensembl_end_phase=-1;ensembl_phase=0;exon_id=AT1G01010.1.exon6;rank=6\n"
            "1\taraport11\tthree_prime_UTR\t5631\t5899\t.\t+\t.\tID=8cc9c37e-60c6-43e8-84e2-0a553385687d;Parent=transcript:AT1G01010.1\n"
        )

    def test_empty_init_sequenceannotation(self):
        SequenceAnnotation()

    def test_read_gff_string(self):
        for gff_string in [self.gff3, self.real_gff3]:
            SequenceAnnotation.from_gff(string=gff_string)

    def test_get_children(self):
        ann = SequenceAnnotation.from_gff(string=self.gff3)
        self.assertEqual(22, len(ann["gene00001"].children))

    def test_input_output_annotation(self):
        ann = SequenceAnnotation.from_gff(string=self.gff3)
        gff3 = ann.to_gff()
        self.assertEqual(gff3, self.gff3)
        ann.to_gtf()

    def test_get_number_of_elements(self):
        ann = SequenceAnnotation.from_gff(string=self.gff3)
        self.assertEqual(23, len(ann))

    def test_groupby(self):
        ann = SequenceAnnotation.from_gff(string=self.gff3)
        groups = ann.groupby(lambda el: el.interval_type)
        self.assertEqual(5, len(groups))
        self.assertSetEqual({"gene", "TF_binding_site", "mRNA", "exon", "CDS"}, set(groups.keys()))

    def test_gene_transcript_exon_chain(self):
        ann = SequenceAnnotation.from_gff(string=self.gff3)
        groups = ann.groupby(lambda el: el.interval_type)
        genes = groups["gene"]
        for gene in genes:
            transcripts = gene.children.groupby(lambda el: el.interval_type)["mRNA"]
            for transcript in transcripts:
                cdss = transcript.children.groupby(lambda el: el.interval_type)["CDS"]
                for cds in cdss:
                    cds.to_gff_line()
                    cds.to_gtf_line()


class AnnotationBugfixTests(TestCase):
    def setUp(self):
        self.gff3 = (
            "ctg1\t.\tgene\t100\t400\t.\t+\t.\tID=g1;Name=foo bar\n"
            "ctg1\t.\tmRNA\t100\t400\t.\t+\t.\tID=t1;Parent=g1;Note=a%3Bb,c%2Cd\n"
            "ctg1\t.\texon\t100\t200\t.\t+\t.\tID=e1;Parent=t1\n"
            "ctg1\t.\texon\t300\t400\t.\t+\t.\tID=e2;Parent=t1\n"
        )
        # GTF without gene and transcript lines
        self.gtf = (
            'ctg1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'ctg1\tsrc\tCDS\t150\t200\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
            'ctg1\tsrc\texon\t300\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'ctg1\tsrc\texon\t50\t250\t.\t+\t.\tgene_id "g1"; transcript_id "t2";\n'
        )
        # GTF with gene and transcript lines (Ensembl style)
        self.full_gtf = (
            'ctg1\tsrc\tgene\t100\t400\t.\t+\t.\tgene_id "g1"; gene_name "G1";\n'
            'ctg1\tsrc\ttranscript\t100\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'ctg1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; exon_number "1";\n'
            'ctg1\tsrc\texon\t300\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; exon_number "2";\n'
        )

    def test_gff_attribute_values_with_special_characters(self):
        ann = SequenceAnnotation.from_gff(string=self.gff3)
        self.assertEqual(ann["g1"].name, ["foo bar"])
        self.assertEqual(ann["t1"].note, ["a;b", "c,d"])
        self.assertEqual(ann.to_gff(), self.gff3)

    def test_gff_decodes_any_percent_encoding(self):
        line = "ctg1\t.\tgene\t1\t9\t.\t+\t.\tID=g1;product=x%20y%2Fz"
        interval = SequenceInterval.from_gff_line(line)
        self.assertEqual(interval.product, ["x y/z"])
        self.assertEqual(interval.to_gff_line(), "ctg1\t.\tgene\t1\t9\t.\t+\t.\tID=g1;product=x y/z")

    def test_real_gff_round_trip(self):
        real_gff3 = AnnotationTests("setUp")
        real_gff3.setUp()
        ann = SequenceAnnotation.from_gff(string=real_gff3.real_gff3)
        self.assertEqual(ann.to_gff(), real_gff3.real_gff3)

    def test_gtf_creates_genes_and_transcripts(self):
        ann = SequenceAnnotation.from_gtf(string=self.gtf)
        groups = ann.groupby(lambda interval: interval.interval_type)
        self.assertEqual({key: len(value) for key, value in groups.items()}, {"gene": 1, "mRNA": 2, "exon": 3, "CDS": 1})
        gene = ann["g1"]
        self.assertEqual((gene.start, gene.end), (50, 400))
        self.assertEqual((ann["t1"].start, ann["t1"].end), (100, 400))
        self.assertEqual((ann["t2"].start, ann["t2"].end), (50, 250))
        self.assertEqual(len(gene.children), 6)
        self.assertEqual(len(ann["t1"].children), 3)

    def test_gtf_with_gene_and_transcript_lines(self):
        ann = SequenceAnnotation.from_gtf(string=self.full_gtf)
        self.assertEqual(len(ann), 4)
        self.assertEqual(ann["g1"].gene_name, ["G1"])
        self.assertEqual(ann["t1"].interval_type, "mRNA")
        self.assertEqual([child.ID for child in ann["t1"].children], ["t1.exon_0", "t1.exon_1"])

    def test_gtf_round_trip(self):
        ann = SequenceAnnotation.from_gtf(string=self.gtf)
        ann2 = SequenceAnnotation.from_gtf(string=ann.to_gtf())
        self.assertEqual(sorted(interval.ID for interval in ann), sorted(interval.ID for interval in ann2))
        self.assertEqual(ann2.to_gtf(), ann.to_gtf())

    def test_gtf_from_gff(self):
        ann = SequenceAnnotation.from_gff(string=self.gff3)
        gtf_lines = ann.to_gtf().split("\n")
        self.assertIn('transcript_id "t1";', gtf_lines[2])
        self.assertIn('gene_id "g1";', gtf_lines[2])
        self.assertNotIn("Parent", ann.to_gtf())
        ann2 = SequenceAnnotation.from_gtf(string=ann.to_gtf())
        self.assertEqual(len(ann2["g1"].children), 3)

    def test_json_round_trip(self):
        ann = SequenceAnnotation.from_gff(string=self.gff3)
        ann2 = SequenceAnnotation.from_json(string=ann.to_json())
        self.assertEqual(ann2.to_gff(), self.gff3)
        self.assertEqual(ann2["t1"].parent, ["g1"])
        self.assertEqual(len(ann2["g1"].children), 3)

    def test_json_include_children(self):
        ann = SequenceAnnotation.from_gff(string=self.gff3)
        gene_dict = ann["g1"].to_dict(include_children=True)
        self.assertEqual([child["ID"] for child in gene_dict["children"]], ["t1", "e1", "e2"])
        ann2 = SequenceAnnotation.from_json(string=f'[{ann["g1"].to_json(include_children=True)}]')
        self.assertEqual(len(ann2["g1"].children), 3)

    def test_rename_interval(self):
        ann = SequenceAnnotation.from_gff(string=self.gff3)
        ann["t1"].ID = "transcript1"
        self.assertNotIn("t1", ann)
        self.assertEqual(ann["e1"].parent, ["transcript1"])
        self.assertEqual([child.ID for child in ann["g1"].children], ["transcript1", "e1", "e2"])
