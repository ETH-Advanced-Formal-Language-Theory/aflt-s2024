import pickle
from numpy import allclose
from rayuela.base.semiring import Tropical, Real
from rayuela.cfg.parser import Parser, EarleyParser
from rayuela.cfg.cfg import CFG
from rayuela.base.misc import compare_chart, symify

pickles_path = "autograding_tests/pickles"
hw_path = pickles_path + "/hw8"

class ChartedEarley(EarleyParser):
    def __init__(self, cfg):
        super().__init__(cfg)

    def earley_chart(self, input, strategy="earley"):
        β = None
        if strategy == "earley":
            β = self._earley(input)
        elif strategy == "agenda":
            β = self._agenda(input)
        else:
            raise NotImplementedError

        chart = self.R.chart()
        for item, w in β.items():
            if item.end:
                chart[item.head, item.i, item.k] += w

        return chart     

def test_earley_example_tropical():
    
    R = Tropical

    cfg = CFG.from_string(
        """
        A → b:	0.142
        A → c:	0.24
        A → a:	0.448
        A → d:	0.33
        A → A A:	0.14
        A → A B:	0.398
        A → A D:	0.301
        A → C C:	0.125
        B → b:	0.39
        B → c:	0.433
        B → a:	0.435
        B → d:	0.359
        B → A C:	0.121
        B → B B:	0.036
        B → D B:	0.476
        B → D C:	0.352
        C → b:	0.391
        C → c:	0.173
        C → B D:	0.062
        C → C D:	0.099
        D → d:	0.15
        D → A A:	0.361
        D → A D:	0.33
        D → B A:	0.04
        S → A A:	0.207
        S → A C:	0.104
        """.strip(), R)

    input = symify("abc")
    ep = ChartedEarley(cfg)
    parser = Parser(cfg)
    compare_chart(R, ep.earley_chart(input) ,parser._cky(input))
    assert (allclose(float(ep.earley(input)), float(parser.cky(input)), 1e-5))

def test_earley_example_real():
    
    R = Real

    cfg = CFG.from_string(
        """
        A → b:	0.142
        A → c:	0.24
        A → B C:	0.403
        A → D A:	0.136
        A → D C:	0.389
        A → C A:	0.284
        A → C C:	0.125
        B → b:	0.39
        B → a:	0.435
        B → d:	0.359
        B → D C:	0.352
        B → C A:	0.307
        B → C B:	0.044
        C → c:	0.173
        C → a:	0.324
        C → d:	0.456
        C → A B:	0.401
        C → B B:	0.064
        C → C D:	0.099
        D → d:	0.15
        D → C B:	0.274
        S → A A:	0.207
        S → C C:	0.287
        """.strip(), R)

    input = symify("abc")
    ep = ChartedEarley(cfg)
    parser = Parser(cfg)
    compare_chart(R, ep.earley_chart(input),parser._cky(input))
    assert (allclose(float(ep.earley(input)), float(parser.cky(input)), 1e-5))

def test_earley():
    with open(f"{hw_path}/cfgs.pkl", 'rb') as f:
        cfgs = pickle.load(f)
    with open(f"{hw_path}/charts.pkl", 'rb') as f:
        charts = pickle.load(f)
    with open(f"{hw_path}/scores.pkl", 'rb') as f:
        scores = pickle.load(f)

    input = symify("abc")

    for cfg, gold_chart, gold_score in zip(cfgs, charts, scores):
        ep = ChartedEarley(cfg)
        chart = ep.earley_chart(input)
        score = float(ep.earley(input))
        compare_chart(cfg.R, chart, gold_chart)
        assert allclose(score, gold_score, atol=1e-5)

def parse_unambiguous():
    cfg = CFG.from_string("""
        1.0: S → A B
        0.3: B → A B
        0.5: A → a
        0.4: B → b
    """, Real)

    earley = Earley(cfg)
    assert earley.parse("ab") == Real(0.2)
    assert earley.parse("aab") == Real(0.03)
    assert earley.parse("aaab") == Real(0.0045)

    earleyfast = EarleyFast(cfg)
    assert earleyfast.parse("ab") == Real(0.2)
    assert earleyfast.parse("aab") == Real(0.03)
    assert earleyfast.parse("aaab") == Real(0.0045)


def parse_left_recursive():

    cfg = CFG.from_string("""
        1.0: S → A B
        0.3: A → A B
        0.5: A → a
        0.4: B → b
    """, Real)

    earley = Earley(cfg)
    assert earley.parse("ab") == Real(0.2)
    assert earley.parse("abb") == Real(0.024)
    assert earley.parse("abbb") == Real(0.00288)

    earleyfast = EarleyFast(cfg)
    assert earleyfast.parse("ab") == Real(0.2)
    assert earleyfast.parse("abb") == Real(0.024)
    assert earleyfast.parse("abbb") == Real(0.00288)

def parse_unary():
    # grammar contains non-cyclic unary rules
    cfg = CFG.from_string("""
        1.0: S → B
        0.3: B → A B
        0.2: B → A
        0.5: A → a
    """, Real)

    earley = Earley(cfg)
    assert earley.parse("a") == Real(0.1)
    assert earley.parse("aa") == Real(0.015)
    assert earley.parse("aaa") == Real(0.00225)

    earleyfast = EarleyFast(cfg)
    assert earleyfast.parse("a") == Real(0.1)
    assert earleyfast.parse("aa") == Real(0.015)
    assert earleyfast.parse("aaa") == Real(0.00225)

    cfg = CFG.from_string("""
        1.0: S → A
        0.5: S → c A
        0.3: A → B
        0.2: B → C
        0.5: C → c
    """, Real)

    earley = Earley(cfg)
    assert earley.parse("c") == Real(0.03)
    assert earley.parse("cc") == Real(0.015)

    earleyfast = EarleyFast(cfg)
    assert earleyfast.parse("c") == Real(0.03)
    assert earleyfast.parse("cc") == Real(0.015)

def parse_mixed():

    cfg = CFG.from_string("""
        1.0: S → a B c D
        0.4: S → A b
        0.1: B → b b
        0.5: A → a
        0.3: D → d
    """, Real)

    earley = Earley(cfg)
    assert earley.parse("ab") == Real(0.2)
    assert earley.parse("abbcd") == Real(0.03)

    earleyfast = EarleyFast(cfg)
    assert earleyfast.parse("ab") == Real(0.2)
    assert earleyfast.parse("abbcd") == Real(0.03)

def parse_ambiguous():

    cfg = CFG.from_string("""
        1.0: S → A
        0.4: A → A + A
        0.1: A → A - A
        0.5: A → a
    """, Real)

    earley = Earley(cfg)
    assert earley.parse("a") == Real(0.5)
    assert earley.parse("a+a") == Real(0.1)
    assert earley.parse("a+a+a") == Real(0.04)

    earleyfast = EarleyFast(cfg)
    assert earleyfast.parse("a") == Real(0.5)
    assert earleyfast.parse("a+a") == Real(0.1)
    assert earleyfast.parse("a+a+a") == Real(0.04)

    cfg = CFG.from_string("""
        1.0: S → A
        0.4: A → A + A
        0.1: A → A - A
        0.5: A → a
    """, MaxTimes)

    earley = Earley(cfg)
    assert earley.parse("a") == MaxTimes(0.5)
    assert earley.parse("a+a") == MaxTimes(0.1)
    assert earley.parse("a+a+a") == MaxTimes(0.02)

    earleyfast = EarleyFast(cfg)
    assert earleyfast.parse("a") == MaxTimes(0.5)
    assert earleyfast.parse("a+a") == MaxTimes(0.1)
    assert earleyfast.parse("a+a+a") == MaxTimes(0.02)

    cfg = CFG.from_string("""
        0.4: A → A + A
        0.1: A → A - A
        0.5: A → a
    """, Real, start="A")

    earley = Earley(cfg)
    assert earley.parse("a") == Real(0.5)
    assert earley.parse("a+a") == Real(0.1)
    assert earley.parse("a+a+a") == Real(0.04)

    earleyfast = EarleyFast(cfg)
    assert earleyfast.parse("a") == Real(0.5)
    assert earleyfast.parse("a+a") == Real(0.1)
    assert earleyfast.parse("a+a+a") == Real(0.04)
