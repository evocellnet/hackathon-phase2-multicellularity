"""Toy-data tests for scripts/data/preprocess_proteomes.py helpers (known answers)."""
import importlib.util, os
P = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
spec = importlib.util.spec_from_file_location('pp', f'{P}/scripts/data/preprocess_proteomes.py')
pp = importlib.util.module_from_spec(spec); spec.loader.exec_module(pp)


def test_clean():
    assert pp.clean('mak*ub-z\n*') == 'MAKXXXX'   # internal stop, U, B, Z -> X; gap removed; terminal * trimmed
    assert pp.clean('ACDE.') == 'ACDE'


def test_sanitise():
    assert pp.sanitise('jgi|Ulvmu1|10000|UM059_0049.1') == 'jgi_Ulvmu1_10000_UM059_0049.1'


def test_trinity_keys():
    k = pp.trinity_gene_key
    assert k('X TRINITY_DN0_c0_g1_i1.p1 TRINITY_DN0_c0_g1~~TRINITY_DN0_c0_g1_i1.p1 ORF') == 'TRINITY_DN0_c0_g1'
    assert k('X c0_g1_i1.p1 GENE.c0_g1_i1~~c0_g1_i1.p1 ORF') == 'c0_g1'
    assert k('X m.922546 g.922546  ORF g.922546 m.922546 type:complete len:259 (+) comp100027_c0_seq1:141-917(+)') == 'comp100027_c0'
    assert k('X Unigene10000_Nuclearia_a|m.30543 x') == 'Unigene10000_Nuclearia_a'
    assert k('X Partr_v1_DN10022_c0_g2_i1_m38271 putative') == 'Partr_v1_DN10022_c0_g2'
    assert k('X CAMPEP_0177630248 /NCGR_PEP_ID=MMETSP0447-20121125|1109_1') is None   # no gene info


def test_genome_keys():
    k = pp.eukprot_gene_key
    assert k('X scpid61325| scgid6254| T-complex') == 'scgid6254'
    assert k('X Cvel_145.t1-p1 | transcript=Cvel_145.t1 | gene=Cvel_145 |') == 'Cvel_145'
    assert k('X jgi|Mesvir1|10000|Mv05792-RA.1') == 'Mv05792'
    assert k('X tig00000017_g1.t1') == 'tig00000017_g1'
    assert k('X Em0001g1a') == 'Em0001g1'
    assert k('X CFRG0001T1') == 'CFRG0001'
    assert k('X Mvib_g6437t1') == 'Mvib_g6437'
    assert k('X ANDGO_02304.mRNA.1 hypothetical') == 'ANDGO_02304'
    assert k('X maker-scaffold_101-snap-gene-0.3-mRNA-1 protein') == 'maker-scaffold_101-snap-gene-0.3'
    assert k('X XP_011385989.1 hypothetical protein') is None   # RefSeq accession: no gene info


def test_eukprot_isoform_pairs_collapse():
    # two isoforms of one gene + one other gene -> 2 keys
    hs = ['a scpid1| scgid7|', 'b scpid2| scgid7|', 'c scpid3| scgid8|']
    assert len({pp.eukprot_gene_key(h) for h in hs}) == 2
