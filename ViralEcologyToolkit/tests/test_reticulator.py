from unittest import TestCase

import ViralEcologyToolkit


class TestReticulator(TestCase):
    def testExtractReadsUnzipped(self):
        s = ViralEcologyToolkit.extract_fastas_from_list("./testfiles/foo.faa",
                                                         "./testfiles/id.list")
        self.assertEquals(len(s), 1)

    def testExtractReadsZipped(self):
        s = ViralEcologyToolkit.extract_fastas_from_list("./testfiles/foo2.faa.gz",
                                                         "./testfiles/id.list")
        self.assertEquals(len(s), 1)

    def testExtractReads(self):
        s = ViralEcologyToolkit.extract_fastas_from_list("./testfiles/foo.faa",
                                                         ['gene'], list_file=False)
        self.assertEquals(len(s), 1)

    def testParseMCLFile(self):
        s = ViralEcologyToolkit.parse_mcl_dump_file("./testfiles/mcl.test")
        self.assertEquals(len(s), 6)
        self.assertEquals(s.iloc[0, 0], 'PC_00001')
        self.assertEquals(s.iloc[0, 1], 'test1')

    def testLoadGeneMap(self):
        s = ViralEcologyToolkit.load_gene_map("./testfiles/gene.map")
        self.assertEquals(s.iloc[0, 0], 'test1')

    def testMergeGeneMap(self):
        s = ViralEcologyToolkit.parse_mcl_dump_file("./testfiles/mcl.test")
        t = ViralEcologyToolkit.load_gene_map("./testfiles/gene.map")
        u = ViralEcologyToolkit.add_pc_labels(s, t)
        self.assertEquals(u.iloc[0, 0], 'test1')
        self.assertEquals(u.iloc[0, 1], 'contig1')
        self.assertEquals(u.iloc[0, 2], 'PC_00001')
        self.assertEquals(u.iloc[5, 2], 'PC_00003')

    def testCreateCompositionMatrix(self):
        s = ViralEcologyToolkit.parse_mcl_dump_file("./testfiles/mcl.test")
        t = ViralEcologyToolkit.load_gene_map("./testfiles/gene.map")
        u = ViralEcologyToolkit.add_pc_labels(s, t)
        v = ViralEcologyToolkit.create_composition_mtx(u)
        self.assertEquals(v.shape, (3, 3))
        self.assertEquals(v.iloc[0, 0], 2)

    def testCreateCompositionMatrixPresenceAbsence(self):
        s = ViralEcologyToolkit.parse_mcl_dump_file("./testfiles/mcl.test")
        t = ViralEcologyToolkit.load_gene_map("./testfiles/gene.map")
        u = ViralEcologyToolkit.add_pc_labels(s, t)
        v = ViralEcologyToolkit.create_composition_mtx(u, presence_absence=True)
        self.assertEquals(v.shape, (3, 3))
        self.assertEquals(v.iloc[0, 0], 1)

    def testCountShared(self):
        s = ViralEcologyToolkit.parse_mcl_dump_file("./testfiles/mcl.test")
        t = ViralEcologyToolkit.load_gene_map("./testfiles/gene.map")
        u = ViralEcologyToolkit.add_pc_labels(s, t)
        v = ViralEcologyToolkit.create_composition_mtx(u, presence_absence=True)
        w = ViralEcologyToolkit.calculate_shared_content('contig1', 'contig2', v)
        self.assertEquals(w, 1)
        w = ViralEcologyToolkit.calculate_shared_content('contig1', 'contig3', v)
        self.assertEquals(w, 1)
        w = ViralEcologyToolkit.calculate_shared_content('contig2', 'contig3', v)
        self.assertEquals(w, 0)

    def testCountShared(self):
        s = ViralEcologyToolkit.parse_mcl_dump_file("./testfiles/mcl.test")
        t = ViralEcologyToolkit.load_gene_map("./testfiles/gene.map")
        u = ViralEcologyToolkit.add_pc_labels(s, t)
        v = ViralEcologyToolkit.create_composition_mtx(u, presence_absence=True)
        w = ViralEcologyToolkit.calculate_shared_matrix(v)
        self.assertEquals(w.shape, (3, 3))
        self.assertEquals(w.loc['contig1', 'contig1'], 1)
        self.assertEquals(w.loc['contig1', 'contig2'], 1)
        self.assertEquals(w.loc['contig1', 'contig3'], 1)
        self.assertEquals(w.loc['contig2', 'contig3'], 0)

    def testCalculateHyperGeometric(self):
        s = ViralEcologyToolkit.parse_mcl_dump_file("./testfiles/mcl.test")
        t = ViralEcologyToolkit.load_gene_map("./testfiles/gene.map")
        u = ViralEcologyToolkit.add_pc_labels(s, t)
        v = ViralEcologyToolkit.create_composition_mtx(u, presence_absence=True)
        w = ViralEcologyToolkit.calculate_shared_matrix(v)
        x = ViralEcologyToolkit.calculate_hypergeometric_survival(u, w)
        print(x)
        self.assertFalse(True)
