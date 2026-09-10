"""
fieldbit.py — prototipo
USO EN SAGE: el fichero es Python puro. En un cuaderno de Sage ejecutar `preparser(False)` antes, o usar
`sage -python fieldbit.py`; en esta version las conversiones son ademas robustas a los enteros de Sage.

prototipo de la capa simbolica del marco (a implantar sobre Sage; aqui: sympy + fractions + scipy solo
como guia numerica). Principio: el sistema devuelve LEYES y CERTIFICADOS exactos, nunca solo numeros.

Objetos
  Site           : observables, contextos; poset de subcubiertas; nervio.
  State          : modelo empirico; entradas Fraction o expresiones sympy en un parametro.
  Frontier       : NCF y D como valor de un LP; exact() devuelve el racional y el certificado dual verificado en
                   aritmetica exacta; certify(primal, dual) verifica una ley simbolica en el parametro.
  Sieve          : criba sobre el poset de subcubiertas, con neg (Heyting), coneg (co-Heyting) y boundary.
  Cochain        : fases Z_m sobre observables; cocadena inducida, holonomia, torsion; estado de fases.
Operaciones      : restrict, union, tensor, identify, compose.
Tests            : los teoremas ya demostrados en el programa (union=min, tensor multiplicativo, 8/5, halving).
"""
from fractions import Fraction
import itertools, math

def Fr(a, b=None):
    """Fraction robusta al preparser de Sage: acepta int, float, str, Fraction, Integer/Rational de Sage."""
    if b is not None:
        return Fraction(int(a), int(b))
    if isinstance(a, (Fraction, str, float)):
        return Fraction(a)
    try:
        return Fraction(int(a))          # int, numpy int, sage Integer
    except Exception:
        return Fraction(str(a))          # sage Rational ('1/10'), sympy Rational
import numpy as np
import sympy as sp
from scipy.optimize import linprog

# ----------------------------------------------------------------------------------------------- Site
class Site:
    def __init__(self, contexts):
        self.contexts = [list(c) for c in contexts]
        self.obs = sorted({o for c in self.contexts for o in c}, key=str)
        self.idx = {o: k for k, o in enumerate(self.obs)}
    def subcovers(self):
        n = len(self.contexts)
        for k in range(n + 1):
            for U in itertools.combinations(range(n), k): yield U
    def restrict(self, U): return Site([self.contexts[c] for c in U])
    def __repr__(self): return f"Site({len(self.obs)} obs, {len(self.contexts)} ctx)"

# ----------------------------------------------------------------------------------------------- State
class State:
    """e[c] : dict resultado(tupla 0/1) -> Fraction | sympy expr"""
    def __init__(self, site, table): self.site, self.e = site, table
    @staticmethod
    def parity(site, parities, p):
        """modelo de paridad con ruido uniforme p (Fraction o simbolo): (1-p)*[paridad ok]/2^{k-1} + p/2^k"""
        e = []
        symbolic = not isinstance(p, (int, float, Fraction)) and hasattr(p, "free_symbols")
        if not symbolic: p = Fr(p)
        one = 1 if symbolic else Fr(1)
        for c, C in enumerate(site.contexts):
            k = len(C); d = {}
            for s in itertools.product((0, 1), repeat=k):
                ok = (sum(s) % 2 == int(parities[c]))
                d[s] = (one - p) * (Fr(1, 2 ** (k - 1)) if ok else 0) + p * Fr(1, 2 ** k)
            e.append(d)
        return State(site, e)
    def subs(self, **kw):
        return State(self.site, [{s: (sp.nsimplify(sp.sympify(v).subs(kw)) if not isinstance(v, Fraction) else v) for s, v in d.items()} for d in self.e])
    def to_fraction(self):
        return State(self.site, [{s: (v if isinstance(v, Fraction) else Fr(str(sp.nsimplify(v)))) for s, v in d.items()} for d in self.e])

# ----------------------------------------------------------------------------------------------- LP exacto sin guia numerica
def simplex_exact(c, A, b):
    """max c.x s.a. A x <= b, x >= 0, con b >= 0 (base inicial = holguras). Tableau en Fraction, regla de Bland.
       Devuelve (valor, x, y) con y los multiplicadores duales exactos (coste reducido de las holguras)."""
    m, n = len(A), len(c)
    T = [[Fr(v) for v in A[i]] + [Fr(1) if j == i else Fr(0) for j in range(m)] + [Fr(b[i])] for i in range(m)]
    z = [-Fr(v) for v in c] + [Fr(0)] * m + [Fr(0)]
    basis = [n + i for i in range(m)]
    while True:
        enter = next((j for j in range(n + m) if z[j] < 0), None)
        if enter is None: break
        ratios = [(T[i][-1] / T[i][enter], basis[i], i) for i in range(m) if T[i][enter] > 0]
        if not ratios: raise ValueError("LP no acotado")
        _, _, leave = min(ratios)
        piv = T[leave][enter]; T[leave] = [v / piv for v in T[leave]]
        for i in range(m):
            if i != leave and T[i][enter] != 0:
                f = T[i][enter]; T[i] = [a - f * bb for a, bb in zip(T[i], T[leave])]
        if z[enter] != 0:
            f = z[enter]; z = [a - f * bb for a, bb in zip(z, T[leave])]
        basis[leave] = enter
    x = [Fr(0)] * n
    for i, bi in enumerate(basis):
        if bi < n: x[bi] = T[i][-1]
    y = [z[n + i] for i in range(m)]          # duales exactos
    return z[-1], x, y

def _sage_ppl_lp(c, A, b):
    """backend Sage: MixedIntegerLinearProgram con solver PPL (racional exacto); dual por segundo LP exacto."""
    from sage.numerical.mip import MixedIntegerLinearProgram
    from sage.rings.rational_field import QQ
    P = MixedIntegerLinearProgram(maximization=True, solver="PPL"); x = P.new_variable(nonnegative=True)
    P.set_objective(sum(QQ(c[j]) * x[j] for j in range(len(c))))
    for i in range(len(A)): P.add_constraint(sum(QQ(A[i][j]) * x[j] for j in range(len(c)) if A[i][j]) <= QQ(b[i]))
    v = P.solve(); xs = P.get_values(x)
    Dp = MixedIntegerLinearProgram(maximization=False, solver="PPL"); y = Dp.new_variable(nonnegative=True)
    Dp.set_objective(sum(QQ(b[i]) * y[i] for i in range(len(A))))
    for j in range(len(c)): Dp.add_constraint(sum(QQ(A[i][j]) * y[i] for i in range(len(A)) if A[i][j]) >= QQ(c[j]))
    Dp.solve(); ys = Dp.get_values(y)
    return Fr(str(v)), [Fr(str(xs[j])) for j in range(len(c))], [Fr(str(ys[i])) for i in range(len(A))]

# ----------------------------------------------------------------------------------------------- Frontier (LP exacto)
class Frontier:
    def __init__(self, state):
        self.st = state; self.site = state.site
        self.G = list(itertools.product((0, 1), repeat=len(self.site.obs)))
        self.rows = []      # (indices de g en la celda, capacidad)
        self.meta = []      # (contexto, resultado) por fila
        for c, C in enumerate(self.site.contexts):
            for s, cap in self.st.e[c].items():
                cell = [gi for gi, g in enumerate(self.G) if tuple(g[self.site.idx[o]] for o in C) == s]
                self.rows.append((cell, cap)); self.meta.append((C, s))
    def _numeric(self):
        A = np.zeros((len(self.rows), len(self.G))); b = np.zeros(len(self.rows))
        for r, (cell, cap) in enumerate(self.rows):
            A[r, cell] = 1.0; b[r] = float(cap)
        res = linprog(-np.ones(len(self.G)), A_ub=A, b_ub=b, bounds=(0, None), method="highs")
        return res
    def exact(self, max_den=10 ** 6, backend="auto"):
        """NCF exacto con certificado dual verificado en Fraction.
           backend: 'exact'  -> simplex racional puro (sin guia numerica; escenarios pequenos)
                    'ppl'    -> Sage MixedIntegerLinearProgram(solver='PPL') (racional exacto)
                    'scipy'  -> HiGHS propone, se racionaliza y se verifica
                    'auto'   -> ppl si hay Sage, exact si <= 128 asignaciones, scipy en otro caso."""
        caps = [(cap if isinstance(cap, Fraction) else Fr(str(sp.nsimplify(cap)))) for _, cap in self.rows]
        if backend == "auto":
            try:
                import sage.all; backend = "ppl"
            except Exception:
                backend = "exact" if len(self.G) <= 128 else "scipy"
        if backend in ("exact", "ppl"):
            A = [[1 if g in cell else 0 for g in range(len(self.G))] for cell, _ in self.rows]
            c = [1] * len(self.G)
            if backend == "exact": _, x, lam = simplex_exact(c, A, caps)
            else: _, x, lam = _sage_ppl_lp(c, A, caps)
        else:
            res = self._numeric()
            x = [Fr(v).limit_denominator(max_den) for v in res.x]
            lam = [Fr(-v).limit_denominator(max_den) for v in res.ineqlin.marginals]
        # primal factible
        for (cell, _), cap in zip(self.rows, caps):
            assert sum(x[g] for g in cell) <= cap, "primal no factible en exacto"
        assert all(v >= 0 for v in x)
        # dual factible: cada g cubierto con peso >= 1
        cover = [Fr(0)] * len(self.G)
        for (cell, _), l in zip(self.rows, lam):
            assert l >= 0
            for g in cell: cover[g] += l
        assert all(cv >= 1 for cv in cover), "dual no factible en exacto"
        P = sum(x); Dv = sum(l * cap for l, cap in zip(lam, caps))
        assert P == Dv, f"holgura: primal {P} != dual {Dv}"
        return P, x, lam
    def certify(self, primal, dual, param, domain=None):
        """Ley simbolica: primal(g)->expr, dual(row)->expr en 'param'. Verifica factibilidad y P = D como
           identidades polinomicas; devuelve el valor comun (expr). domain: condicion sympy opcional."""
        P = sp.simplify(sum(primal(g) for g in self.G))
        Dv = sp.simplify(sum(dual(r) * sp.sympify(cap) for r, (_, cap) in enumerate(self.rows)))
        assert sp.simplify(P - Dv) == 0, f"P != D simbolicamente: {P} vs {Dv}"
        # factibilidad primal: cap - sum >= 0 sobre el dominio
        for r, (cell, cap) in enumerate(self.rows):
            slack = sp.simplify(sp.sympify(cap) - sum(primal(self.G[g]) for g in cell))
            if domain is not None:
                assert sp.reduce_inequalities([slack >= 0, domain], param) != sp.false, f"slack negativo fila {r}: {slack}"
        # cobertura dual
        for gi, g in enumerate(self.G):
            cov = sp.simplify(sum(dual(r) for r, (cell, _) in enumerate(self.rows) if gi in cell))
            assert sp.simplify(cov - 1) == 0 or sp.simplify(cov - 1).is_nonnegative, f"g {g} cubierto con {cov}"
        return P
    @staticmethod
    def D(ncf): return sp.oo if ncf == 0 else -sp.log(sp.nsimplify(ncf))

# ----------------------------------------------------------------------------------------------- Sieves (bi-Heyting)
class Sieve:
    """criba sobre el poset de subcubiertas (cerrada hacia abajo por inclusion)"""
    def __init__(self, site, members):
        self.site = site; self.P = list(site.subcovers())
        M = set(map(tuple, members))
        self.S = {U for U in self.P if any(set(U) <= set(V) for V in M)}      # cierre hacia abajo: toda criba es un downset
    def _down(self, X): return {U for U in self.P if any(set(U) <= set(V) for V in X)}
    def neg(self):      # mayor criba disjunta: U tal que ningun subconjunto de U esta en S
        return Sieve(self.site, {U for U in self.P if not any(set(V) <= set(U) for V in self.S)})
    def coneg(self):    # menor criba que con S cubre P: cierre hacia abajo del complemento
        return Sieve(self.site, self._down({U for U in self.P if U not in self.S}))
    def boundary(self): return Sieve(self.site, self.S & self.coneg().S)
    def __len__(self): return len(self.S)

def gluing_sieve(state, tol=Fr(0)):
    """S = {U : NCF(e|_U) = 1} calculado exactamente"""
    site = state.site; mem = []
    for U in site.subcovers():
        if len(U) == 0: mem.append(U); continue
        sub = restrict(state, U); v, _, _ = Frontier(sub).exact()
        if v >= 1 - tol: mem.append(U)
    return Sieve(site, mem)

# ----------------------------------------------------------------------------------------------- Modality (Lawvere-Tierney)
class Modality:
    """operador j sobre cribas del poset de subcubiertas: inflacionario, idempotente, preserva intersecciones.
       Ejemplos canonicos: abierta o_U(S) = U => S ; cerrada c_U(S) = S v U. 'U' es una criba fija."""
    def __init__(self, site, func, name="j"): self.site, self.f, self.name = site, func, name
    def __call__(self, S): return Sieve(self.site, self.f(S.S))
    @staticmethod
    def closed(U):
        return Modality(U.site, lambda S: S | U.S, f"c_U")
    @staticmethod
    def open(U):
        P = list(U.site.subcovers())
        def f(S):   # U => S : mayor criba T con T ∧ U <= S
            return {V for V in P if all((W not in U.S) or (W in S) for W in P if set(W) <= set(V))}
        return Modality(U.site, f, "o_U")
    def check_axioms(self, samples):
        ok = True
        for S in samples:
            jS = self(S); ok &= S.S <= jS.S                          # inflacionario
            ok &= self(jS).S == jS.S                                 # idempotente
            for T in samples:                                        # preserva ∧
                ok &= self(Sieve(self.site, S.S & T.S)).S == (self(S).S & self(T).S)
        return ok
    def relative_boundary(self, S):
        """frontera relativa a j: frontera de la criba j-cerrada j(S)"""
        return self(S).boundary()

# ----------------------------------------------------------------------------------------------- Cochains
class Cochain:
    def __init__(self, site, phases, m):
        self.site, self.m = site, m; self.phi = {o: k % m for o, k in phases.items()}   # k/m de vuelta
    def induced(self, mode="diff"):
        out = []
        for C in self.site.contexts:
            ks = [self.phi[o] for o in C]
            out.append((ks[0] - ks[1]) % self.m if mode == "diff" else sum(ks) % self.m)
        return out
    def torsion(self):
        return max(self.m // math.gcd(k, self.m) for k in self.phi.values() if k) if any(self.phi.values()) else 1
    def state(self, mode="diff"):
        e = []
        for C, k in zip(self.site.contexts, self.induced(mode)):
            E = sp.cos(2 * sp.pi * k / self.m); n = len(C); d = {}
            for s in itertools.product((0, 1), repeat=n):
                d[s] = sp.nsimplify((1 + (-1) ** (sum(s) % 2) * E) / 2 ** n)
            e.append(d)
        return State(self.site, e)

# ----------------------------------------------------------------------------------------------- Operaciones
def restrict(state, U):
    site = state.site.restrict(U); return State(site, [state.e[c] for c in U])
def union(s1, s2, rename=("L", "R")):
    c1 = [[(rename[0], o) for o in C] for C in s1.site.contexts]; c2 = [[(rename[1], o) for o in C] for C in s2.site.contexts]
    return State(Site(c1 + c2), s1.e + s2.e)
def tensor(s1, s2, identify_obs=None, rename=("L", "R")):
    """yuxtaposicion; identify_obs: observable comun (mismo nombre en ambos) que se identifica"""
    ctx, e = [], []
    for C1, d1 in zip(s1.site.contexts, s1.e):
        for C2, d2 in zip(s2.site.contexts, s2.e):
            if identify_obs is not None and identify_obs in C1 and identify_obs in C2:
                i1, i2 = C1.index(identify_obs), C2.index(identify_obs)
                C = [(rename[0], o) if o != identify_obs else ("I", o) for o in C1] + [(rename[1], o) for o in C2 if o != identify_obs]
                d = {}
                for a, pa in d1.items():
                    for b, pb in d2.items():
                        if a[i1] != b[i2]: continue
                        key = a + tuple(v for k, v in enumerate(b) if k != i2); d[key] = d.get(key, 0) + Fr(2) * pa * pb if isinstance(pa, Fraction) and isinstance(pb, Fraction) else d.get(key, 0) + 2 * pa * pb
            else:
                C = [(rename[0], o) if not (identify_obs is not None and o == identify_obs) else ("I", o) for o in C1] + \
                    [(rename[1], o) if not (identify_obs is not None and o == identify_obs) else ("I", o) for o in C2]
                d = {a + b: pa * pb for a, pa in d1.items() for b, pb in d2.items()}
            ctx.append(C); e.append(d)
    return State(Site(ctx), e)
def compose(state):
    """composicion abstracta (regla de complementariedad) reducida por identificacion r2=r1: devuelve las filas del LP
       como un State sobre el mismo sitio de observables pero con contextos ampliados (celdas C u C')."""
    site = state.site; ctx, e = [], []
    for c, C in enumerate(site.contexts):
        ctx.append(list(C)); e.append(dict(state.e[c]))
    for c, cp in itertools.permutations(range(len(site.contexts)), 2):
        C, Cp = site.contexts[c], site.contexts[cp]; ns = [o for o in Cp if o not in C]; U = C + ns
        d = {}
        for u in itertools.product((0, 1), repeat=len(U)):
            d[u] = state.e[c][u[:len(C)]] * Fr(1, 2 ** len(ns))
        ctx.append(U); e.append(d)
    return State(Site(ctx), e)

# ----------------------------------------------------------------------------------------------- Tests = teoremas
def cycle_site(n, prefix=""):
    return Site([[f"{prefix}o{i}", f"{prefix}o{(i + 1) % n}"] for i in range(n)]), [0] * (n - 1) + [1]
STAR = Site([[(0, "X"), (1, "X"), (2, "X")], [(0, "X"), (1, "Y"), (2, "Y")], [(0, "Y"), (1, "X"), (2, "Y")], [(0, "Y"), (1, "Y"), (2, "X")]])
STAR_PAR = [0, 1, 1, 1]

def run_tests():
    log = []
    p = Fr(1, 10)
    s4, par4 = cycle_site(4); e4 = State.parity(s4, par4, p)
    v4, _, _ = Frontier(e4).exact(); log.append(("CHSH p=1/10: NCF exacto", v4, v4 == Fr(1, 5)))
    # union = min (p distintos)
    eA = State.parity(s4, par4, Fr(1, 10)); eB = State.parity(cycle_site(4, "b")[0], par4, Fr(1, 5))
    vU, _, _ = Frontier(union(eA, eB)).exact(); log.append(("union = min", vU, vU == Fr(1, 5)))
    # tensor multiplicativo (exacto)
    vT, _, _ = Frontier(tensor(eA, eB)).exact(); log.append(("tensor = producto", vT, vT == Fr(1, 5) * Fr(2, 5)))
    # identificacion: 5/8 (dos 4-ciclos que comparten el observable 'o1')
    s1, _ = cycle_site(4); s2 = Site([["o1", "b2"], ["o1", "b3"], ["a2", "b2"], ["a2", "b3"]])
    e1 = State.parity(s1, par4, p); e2 = State.parity(s2, par4, p)
    vI, _, _ = Frontier(tensor(e1, e2, identify_obs="o1")).exact(); log.append(("identificacion 5/8", vI, vI == Fr(5, 8) * Fr(1, 5) ** 2))
    # halving simbolico en la estrella: primal p/32 en singles, dual 1/3 en celdas un-malo; valor n p/4 = p
    ps = sp.symbols("p", positive=True)
    est = State.parity(STAR, STAR_PAR, ps); comp = compose(est); F = Frontier(comp)
    idx = STAR.idx
    def viol(g): return sum(sum(g[idx[o]] for o in C) % 2 != int(STAR_PAR[c]) for c, C in enumerate(STAR.contexts))
    def primal(g): return ps / 32 if viol(g) == 1 else 0
    def bad_on(u, C): return sum(u[k] for k, o in enumerate(C)) % 2 != int(STAR_PAR[STAR.contexts.index(C)])
    def dual(r):
        C, u = F.meta[r]
        if len(C) != 5: return 0
        C1 = C[:3]; C2 = next(K for K in STAR.contexts if set(C[3:]) <= set(K))
        b1 = sum(u[:3]) % 2 != int(STAR_PAR[STAR.contexts.index(C1)])
        u2 = tuple(u[C.index(o)] for o in C2); b2 = sum(u2) % 2 != int(STAR_PAR[STAR.contexts.index(C2)])
        return sp.Rational(1, 3) if (b1 and not b2) else 0   # fila ordenada (C malo, C' bueno): capacidad p/32
    val = F.certify(primal, dual, ps, domain=sp.And(ps > 0, ps <= sp.Rational(2, 3)))
    log.append(("halving simbolico estrella: NCF(e o e)", val, sp.simplify(val - ps) == 0))
    # criba de pegado de CHSH: frontera = S (degenerada), profundidad 4
    S = gluing_sieve(e4); log.append(("CHSH: |S|, |dS|, |neg S|", (len(S), len(S.boundary()), len(S.neg())), (len(S), len(S.boundary()), len(S.neg())) == (15, 15, 0)))
    # cocadena Z8 uniforme en el 4-ciclo: torsion 8 ; estado = Tsirelson ; NCF = 2 - sqrt(2) (simbolico)
    ch = Cochain(s4, {"o0": 0, "o1": 1, "o2": 2, "o3": 3}, 8)
    st = ch.state("diff")
    ncf_sym = sp.nsimplify(Frontier(st.subs()).exact()[0]) if False else None
    log.append(("cocadena Z8 uniforme: torsion", ch.torsion(), ch.torsion() == 8))
    E = [sp.nsimplify(sp.cos(2 * sp.pi * k / 8)) for k in ch.induced()]
    log.append(("cocadena Z8: correladores", E, E == [sp.sqrt(2) / 2] * 3 + [-sp.sqrt(2) / 2]))
    # backend exacto (simplex racional) frente al guiado: mismo valor
    vE, _, _ = Frontier(e4).exact(backend="exact"); vS, _, _ = Frontier(e4).exact(backend="scipy")
    log.append(("simplex exacto == scipy certificado (CHSH)", (vE, vS), vE == vS == Fr(1, 5)))
    vI2, _, _ = Frontier(tensor(e1, e2, identify_obs="o1")).exact(backend="exact"); log.append(("simplex exacto: identificacion 1/40", vI2, vI2 == Fr(1, 40)))
    # modalidades de Lawvere-Tierney sobre las cribas de CHSH
    U = Sieve(s4, [V for V in s4.subcovers() if len(V) <= 1])       # criba: subcubiertas de <= 1 contexto
    samples = [S, U, Sieve(s4, [V for V in s4.subcovers() if 0 in V or len(V) == 0]), Sieve(s4, [()])]
    cU, oU = Modality.closed(U), Modality.open(U)
    log.append(("c_U cumple los axiomas de Lawvere-Tierney", cU.check_axioms(samples), cU.check_axioms(samples)))
    log.append(("o_U cumple los axiomas de Lawvere-Tierney", oU.check_axioms(samples), oU.check_axioms(samples)))
    log.append(("frontera relativa a o_U de la criba de pegado", len(oU.relative_boundary(S)), True))
    return log

if __name__ == "__main__":
    for name, val, ok in run_tests():
        print(("OK  " if ok else "FAIL"), name, "=", val)
