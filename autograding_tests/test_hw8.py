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
