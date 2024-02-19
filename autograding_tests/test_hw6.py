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


def test_nullary_example():
    
    R= Tropical

    cfg = CFG.from_string(
        """
        A → ε:	17
        A → b:	5
        A → A A:	40
        A → A D:	8
        A → D D:	33
        D → ε:	34
        D → b:	21
        D → A A:	30
        D → A D:	46
        D → D A:	5
        S → ε:	31
        S → b:	37
        S → A A:	28
        S → A D:	43
        S → D A:	14
        """.strip(), R)

    cfg_sum = Treesum(cfg).sum()

    transformer = Transformer()
    rcfg = transformer.nullaryremove(cfg)
    rcfg_sum = Treesum(rcfg).sum()

    for (p,w) in rcfg.P:
        assert not nullary(p, rcfg)

    assert (allclose(float(cfg_sum), float(rcfg_sum), atol=10e-5))


def test_nullary():

    with open(f"{hw_path}/cnfs_nullary.pkl", 'rb') as f:
        cfgs = pickle.load(f)

    for cfg in cfgs:
        cfg_sum = Treesum(cfg).sum()

        transformer = Transformer()
        rcfg = transformer.nullaryremove(cfg)
        rcfg_sum = Treesum(rcfg).sum()

        for (p,w) in rcfg.P:
            assert not nullary(p, rcfg)
        assert (allclose(float(cfg_sum), float(rcfg_sum), atol=10e-5))



def test_unary_example():
    
    R= Tropical

    cfg = CFG.from_string(
        """
        A → c:	42
        A → A A:	  3
        A → B B:	  36
        A → D B:	  13
        A → B D:	  29
        B → x:	42
        B → A B:	  26
        B → A:	  11
        D → A:	41
        D → y:	24
        D → B:	  3
        S → B A:	  53
        S → ε:	33
        """.strip(), R)

    cfg_sum = Treesum(cfg).sum()

    transformer = Transformer()
    ucfg = transformer.unaryremove(cfg)
    ucfg_sum = Treesum(ucfg).sum()

    for (p,w) in ucfg.P:
        assert not unary(p)

    assert (allclose(float(cfg_sum), float(ucfg_sum), atol=10e-5))


def test_unary():

    with open(f"{hw_path}/cnfs_unary.pkl", 'rb') as f:
        cfgs = pickle.load(f)

    for cfg in cfgs:
        cfg_sum = Treesum(cfg).sum()

        transformer = Transformer()
        ucfg = transformer.unaryremove(cfg)
        ucfg_sum = Treesum(ucfg).sum()

        for (p,w) in ucfg.P:
            assert not unary(p)

        assert allclose(float(cfg_sum), float(ucfg_sum), atol=10e-5)



def test_nullary_extra_example():

    R= Tropical

    cfg = CFG.from_string(
        """
        A → ε:  1
        A → c:  42
        A → A A:    3
        A → B B:    36
        A → D B:    13
        A → B D:    29
        A → ε:  24
        B → c:  42
        B → ε:  22
        B → A B:    26
        B → c B:    13
        B → A ε:    11
        B → c A:    39
        D → A:  41
        D → c:  24
        D → ε:  0
        D → D ε:    5
        D → ε B:    3
        D → B B:    2
        """.strip(), R)

    cfg_sum = Treesum(cfg).sum()

    transformer = Transformer()
    rcfg = transformer.nullaryremove(cfg)
    rcfg_sum = Treesum(rcfg).sum()

    for (p,w) in rcfg.P:
        assert not nullary(p, rcfg)

    assert (allclose(float(cfg_sum), float(rcfg_sum), atol=10e-5))


def test_nullary_extra():

    with open(f"{hw_path}/cfgs.pkl", 'rb') as f:
        cfgs = pickle.load(f)

    for cfg in cfgs:
        cfg_sum = Treesum(cfg).sum()

        transformer = Transformer()
        rcfg = transformer.nullaryremove(cfg)
        rcfg_sum = Treesum(rcfg).sum()

        for (p,w) in rcfg.P:
            assert not nullary(p, rcfg)
        assert (allclose(float(cfg_sum), float(rcfg_sum), atol=10e-5))
