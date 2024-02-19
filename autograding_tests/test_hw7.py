import numpy as np
from numpy import allclose

import pickle

from rayuela.base.symbol import Sym, ε
from rayuela.base.semiring import Boolean, Real, Tropical, Rational

from rayuela.base.misc import compare_chart

from rayuela.cfg.production import Production
from rayuela.cfg.nonterminal import NT, S
from rayuela.cfg.cfg import CFG
from rayuela.cfg.transformer import Transformer
from rayuela.cfg.misc import *
from rayuela.cfg.treesum import Treesum

pickles_path = "autograding_tests/pickles"
hw_path = pickles_path + "/hw6"

def test_separate_and_binarize():
    with open(f"{hw_path}/cfgs.pkl", 'rb') as f:
        cfgs = pickle.load(f)
    
    T = Transformer()

    for cfg in cfgs:
        scfg = T.separate_terminals(cfg).trim()
        ncfg = T.nullaryremove(scfg).trim()
        ucfg = T.unaryremove(ncfg).trim()
        bcfg = T.binarize(ucfg)

        for (p,w) in scfg.P:
            assert separated(p)

        for (p,w) in bcfg.P:
            assert binarized(p) or preterminal(p) or unary(p) 

        sums = map(lambda x : Treesum(x).sum(), [cfg, scfg, bcfg])
        for s in sums:
            for t in sums:
                assert allclose(float(s), float(t), atol=10e-5)