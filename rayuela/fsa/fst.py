import copy
from email.policy import default
from frozendict import frozendict
from itertools import chain, combinations, product

from collections import deque, Counter
from collections import defaultdict as dd

from rayuela.base.semiring import Boolean, String, ProductSemiring, product_semiring_builder
from rayuela.base.misc import epsilon_filter
from rayuela.base.symbol import Sym, ε, ε_1, ε_2

from rayuela.fsa.fsa import FSA
from rayuela.fsa.state import State, PairState, PowerState, MinimizeState
from rayuela.fsa.pathsum import Pathsum, Strategy
from rayuela.fsa.transformer import Transformer


class FST(FSA):

    def __init__(self, R=Boolean):

        # DEFINITION
        # A weighted finite-state transducer is a 8-tuple <Σ, Δ, Q, F, I, δ, λ, ρ> where
        # • Σ is an alphabet of symbols;
        # • Δ is an alphabet of symbols;
        # • Q is a finite set of states;
        # • I ⊆ Q is a set of initial states;
        # • F ⊆ Q is a set of final states;
        # • δ is a finite relation Q × Σ × Δ × Q × R;
        # • λ is an initial weight function;
        # • ρ is a final weight function.

        # NOTATION CONVENTIONS
        # • single states (elements of Q) are denoted q
        # • multiple states not in sequence are denoted, p, q, r, ...
        # • multiple states in sequence are denoted i, j, k, ...
        # • symbols (elements of Σ and Δ) are denoted lowercase a, b, c, ...
        # • single weights (elements of R) are denoted w
        # • multiple weights (elements of R) are denoted u, v, w, ...

        super().__init__(R=R)

        # alphabet of output symbols
        self.Delta = set()

    def add_arc(self, i, a, b, j, w=None):
        if w is None: w = self.R.one

        if not isinstance(i, State): i = State(i)
        if not isinstance(j, State): j = State(j)
        if not isinstance(a, Sym): a = Sym(a)
        if not isinstance(b, Sym): b = Sym(b)
        if not isinstance(w, self.R): w = self.R(w)

        self.add_states([i, j])
        self.Sigma.add(a)
        self.Delta.add(b)
        self.δ[i][(a, b)][j] += w
        self.δ_inv[j][(a, b)][i] += w

    def set_arc(self, i, a, b, j, w=None):
        if w is None: w = self.R.one

        if not isinstance(i, State): i = State(i)
        if not isinstance(j, State): j = State(j)
        if not isinstance(a, Sym): a = Sym(a)
        if not isinstance(b, Sym): b = Sym(b)
        if not isinstance(w, self.R): w = self.R(w)

        self.add_states([i, j])
        self.Sigma.add(a)
        self.Delta.add(b)
        self.δ[i][(a, b)][j] = w
        self.δ_inv[j][(a, b)][i] = w

    def freeze(self):
        self.Sigma = frozenset(self.Sigma)
        self.Delta = frozenset(self.Delta)
        self.Q = frozenset(self.Q)
        self.δ = frozendict(self.δ)
        self.λ = frozendict(self.λ)
        self.ρ = frozendict(self.ρ)

    def arcs(self, i, no_eps=False):
        for ab, T in self.δ[i].items():
            if no_eps and ab == (ε, ε):
                continue
            for j, w in T.items():
                if w == self.R.zero:
                    continue
                yield ab, j, w

    def accept(self, string1, string2):
        """ determines whether a string is in the language """
        # Requires composition
        raise NotImplementedError

    def top_compose_brute(self, fst):
        # the two machines need to be in the same semiring
        assert self.R == fst.R

        # add initial states
        product_fst = FST(R=self.R)
        for (q1, w1), (q2, w2) in product(self.I, fst.I):
            product_fst.add_I(PairState(q1, q2), w=w1 * w2)

        self_finals = {q: w for q, w in self.F}
        fsa_finals = {q: w for q, w in fst.F}

        # add "body" states
        for q1, q2 in product(self.Q, fst.Q):

            E1 = [((a,b), j, w) for ((a,b), j, w) in self.arcs(q1)]
            E2 = [((a,b), j, w) for ((a,b), j, w) in fst.arcs(q2)]

            M = [((a, j1, w1), (d, j2, w2))
                 for ((a,b), j1, w1), ((c,d), j2, w2) in product(E1, E2)
                 if b == c]

            for (a, j1, w1), (d, j2, w2) in M:
                product_fst.add_arc(
                    PairState(q1, q2), a, d,
                    PairState(j1, j2), w=w1*w2)

        # final state handling
        for (q1,w1), (q2,w2) in product(self.F, fst.F):
            if PairState(q1, q2) in product_fst.Q:
                product_fst.add_F(
                    PairState(q1, q2), w=w1*w2)

        return product_fst

    def bottom_compose_brute(self, fst):
        return fst.top_compose_brute(self)

    def top_compose(self, fst):
        # Assignment 3: Question 3
        raise NotImplementedError

    def bottom_compose(self, fst):
        # Assignment 3: Question 3
        raise NotImplementedError
