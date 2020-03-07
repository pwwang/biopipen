from pathlib import Path
import pytest
from pyppl import PyPPL
from bioprocs.seq import pSeqMutate
from . import assertInfile

def test_seqmutate(rdata):
    # don't do checkref
    pSeqMutate1 = pSeqMutate.copy()
    pSeqMutate1.args.seqbase = 1
    pSeqMutate1.input = rdata.get('seq/seq.fa'), rdata.get('seq/mut.bed')
    PyPPL().start(pSeqMutate1).run()
    assertInfile(pSeqMutate1.workdir.joinpath('1', 'job.stderr'),
                 "Mutation ['chr3', '100', '100', 'SNP3', '0', '+', 'A', 'g'] not found in sequence")
    assertInfile(pSeqMutate1.workdir.joinpath('1', 'job.stderr'),
                 "Mutation ['chr4', '1', '1', 'SNP4', '0', '+', 'A', 'g'] not found in sequence")
    assertInfile(pSeqMutate1.workdir.joinpath('1', 'output', 'seq.mutated.fa'),
                 'TACACTGTGtATC')
    assertInfile(pSeqMutate1.workdir.joinpath('1', 'output', 'seq.mutated.fa'),
                 'CaGAC')
    assertInfile(pSeqMutate1.workdir.joinpath('1', 'output', 'seq.mutated.fa'),
                 'GGCTCTCCGGGtA')
    assertInfile(pSeqMutate1.workdir.joinpath('1', 'output', 'seq.mutated.fa'),
                 'TGGAATGTAAAGAAGTATGTAGAACGGGGTGGTAGT')
    assertInfile(pSeqMutate1.workdir.joinpath('1', 'output', 'seq.mutated.fa'),
                 'GENE::chr2:2-90')


