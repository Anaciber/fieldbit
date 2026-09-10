"""
Comprobacion simbolica: la categoria de retornos del corpus (§10.2-10.4) es Temperley-Lieb, y su axioma de
positividad del cierre es la positividad de la traza de Markov sobre la torre de proposiciones (Jones 1983).

Implementacion exacta del algebra de diagramas planos TL_k(delta): un diagrama es un emparejamiento plano de
2k puntos (k arriba, k abajo); la composicion concatena y cada lazo cerrado vale delta (simbolico).
  - p_i := (1/delta) * (cap_i o cup_i) = e_i / delta   (el "cerrar y reabrir normalizado por el lazo" del corpus)
  - relaciones del corpus: p_i^2 = p_i ; p_i p_{i+1} p_i = delta^{-2} p_i ; [p_i, p_j] = 0 si |i-j| >= 2
  - traza de Markov: cerrar el diagrama por la derecha (cada cadena cerrada = delta), normalizada por delta^k
  - torre de proposiciones = idempotentes de Jones-Wenzl JW_k (recursion de Wenzl), cierre de JW_k = [k+1]_q
    con delta = q + 1/q, enteros cuanticos [m] = (q^m - q^-m)/(q - q^-1) = sen(m pi/n)/sen(pi/n) si delta = 2cos(pi/n)
  - positividad del cierre de toda proposicion <=> [m]_q >= 0 para todo m alcanzable <=> delta in {2cos(pi/n)} U [2, oo)
    (Jones): entre dos valores de la serie discreta algun [m]_q es negativo; en delta = 2cos(pi/n), [n]_q = 0 y la
    torre se trunca en n-1 (ideal despreciable).
"""
import sympy as sp, itertools
from fractions import Fraction as Fr
d = sp.symbols("delta", positive=True)

# ------------------------------------------------------------ diagramas planos TL_k
def planar_matchings(k):
    """emparejamientos sin cruces de los puntos 0..2k-1 dispuestos en circulo (arriba 0..k-1 izquierda->derecha,
       abajo k..2k-1 derecha->izquierda)"""
    pts = list(range(2 * k)); out = []
    def rec(pts, cur):
        if not pts: out.append(frozenset(cur)); return
        a = pts[0]
        for j in range(1, len(pts), 2):        # el companero debe dejar un numero par de puntos entre medias
            b = pts[j]; inside = pts[1:j]; outside = pts[j + 1:]
            rec_inside = []
            # los de dentro y fuera se emparejan entre si por separado
            for m1 in _match(inside):
                for m2 in _match(outside):
                    out.append(frozenset(cur | {frozenset((a, b))} | m1 | m2))
        return
    def _match(pts):
        if not pts: return [frozenset()]
        res = []; a = pts[0]
        for j in range(1, len(pts), 2):
            b = pts[j]
            for m1 in _match(pts[1:j]):
                for m2 in _match(pts[j + 1:]):
                    res.append(frozenset({frozenset((a, b))}) | m1 | m2)
        return res
    return _match(pts)

def compose(D1, D2, k):
    """D1 arriba de D2: los puntos abajo de D1 (k..2k-1, ordenados derecha->izquierda) se identifican con los de
       arriba de D2 (0..k-1 izquierda->derecha): abajo_i de D1 = punto 2k-1-i ; arriba_i de D2 = punto i.
       Devuelve (diagrama, numero de lazos)."""
    # etiquetas: ('T', p) puntos de D1, ('B', p) puntos de D2 ; puntos intermedios ('M', i)
    adj = {}
    def add(u, v): adj.setdefault(u, []).append(v); adj.setdefault(v, []).append(u)
    for pair in D1:
        a, b = tuple(pair)
        def lab1(p): return ('M', 2 * k - 1 - p) if p >= k else ('T', p)
        add(lab1(a), lab1(b))
    for pair in D2:
        a, b = tuple(pair)
        def lab2(p): return ('M', p) if p < k else ('B', p)
        add(lab2(a), lab2(b))
    seen = set(); loops = 0; new = set()
    for start in list(adj):
        if start in seen or start[0] == 'M': continue
        # recorrer desde un punto exterior hasta otro exterior
        path = [start]; seen.add(start); prev = None; cur = start
        while True:
            nxt = [v for v in adj[cur] if v != prev or len(adj[cur]) == 1 and v == prev]
            nxt = [v for v in adj[cur] if not (v == prev and len(adj[cur]) > 1 and adj[cur].count(prev) == 1)] if False else [v for v in adj[cur] if v != prev]
            if not nxt: break
            prev, cur = cur, nxt[0]; seen.add(cur)
            if cur[0] != 'M': break
        # extremos: start y cur (ambos exteriores); mapear a puntos del diagrama compuesto
        def lab(u):
            return u[1] if u[0] == 'T' else u[1]   # 'T' p<k arriba ; 'B' p>=k abajo (misma numeracion)
        new.add(frozenset((lab(start), lab(cur))))
    for u in adj:
        if u not in seen and u[0] == 'M':
            # lazo cerrado de puntos intermedios
            loops += 1; cur = u; prev = None
            while True:
                seen.add(cur); nxt = [v for v in adj[cur] if v != prev]
                if not nxt or nxt[0] == u: break
                prev, cur = cur, nxt[0]
    return frozenset(new), loops

class TL:
    def __init__(self, k):
        self.k = k; self.basis = planar_matchings(k); self.index = {D: i for i, D in enumerate(self.basis)}
        self.n = len(self.basis)
        self.mult = {}
        for D1 in self.basis:
            for D2 in self.basis:
                self.mult[(D1, D2)] = compose(D1, D2, k)
    def vec(self, D): v = [sp.Integer(0)] * self.n; v[self.index[D]] = sp.Integer(1); return v
    def identity(self): return self.vec(frozenset(frozenset((i, 2 * self.k - 1 - i)) for i in range(self.k)))
    def e(self, i):
        """cap_i o cup_i: une arriba i,i+1 y abajo i,i+1; el resto vertical"""
        k = self.k; pairs = {frozenset((i, i + 1)), frozenset((2 * k - 1 - i, 2 * k - 2 - i))}
        for j in range(k):
            if j not in (i, i + 1): pairs.add(frozenset((j, 2 * k - 1 - j)))
        return self.vec(frozenset(pairs))
    def mul(self, x, y):
        out = [sp.Integer(0)] * self.n
        for i, a in enumerate(x):
            if a == 0: continue
            for j, b in enumerate(y):
                if b == 0: continue
                D, loops = self.mult[(self.basis[i], self.basis[j])]
                out[self.index[D]] += a * b * d ** loops
        return [sp.expand(v) for v in out]
    def add(self, x, y, cx=1, cy=1): return [sp.expand(cx * a + cy * b) for a, b in zip(x, y)]
    def scale(self, x, c): return [sp.expand(c * a) for a in x]
    def trace(self, x):
        """traza de Markov normalizada: cerrar cada cadena; cada diagrama D contribuye delta^{lazos(D cerrado)} / delta^k"""
        tot = 0
        for i, a in enumerate(x):
            if a == 0: continue
            D = self.basis[i]; k = self.k
            # cerrar: identificar arriba j con abajo j (punto 2k-1-j) ; contar ciclos
            adj = {}
            for pair in D:
                u, v = tuple(pair); adj.setdefault(u, []).append(v); adj.setdefault(v, []).append(u)
            for j in range(k):
                adj[j].append(2 * k - 1 - j); adj[2 * k - 1 - j].append(j)
            seen = set(); cyc = 0
            for s in adj:
                if s in seen: continue
                cyc += 1; stack = [s]
                while stack:
                    u = stack.pop()
                    if u in seen: continue
                    seen.add(u); stack.extend(adj[u])
            tot += a * d ** cyc / d ** k
        return sp.simplify(tot)
    def is_zero(self, x): return all(sp.simplify(a) == 0 for a in x)

# ------------------------------------------------------------ comprobaciones
log = []
T3 = TL(3)
p1 = T3.scale(T3.e(0), 1 / d); p2 = T3.scale(T3.e(1), 1 / d)
log.append(("p1^2 = p1 (idempotente normalizado por el lazo)", T3.is_zero(T3.add(T3.mul(p1, p1), p1, 1, -1))))
log.append(("p1 p2 p1 = delta^-2 p1 (holonomia del corpus)", T3.is_zero(T3.add(T3.mul(T3.mul(p1, p2), p1), p1, 1, -1 / d ** 2))))
log.append(("p1 p2 != p2 p1 (no conmutan)", not T3.is_zero(T3.add(T3.mul(p1, p2), T3.mul(p2, p1), 1, -1))))
T4 = TL(4); q1 = T4.scale(T4.e(0), 1 / d); q3 = T4.scale(T4.e(2), 1 / d)
log.append(("p1 p3 = p3 p1 (conmutan a distancia)", T4.is_zero(T4.add(T4.mul(q1, q3), T4.mul(q3, q1), 1, -1))))
log.append(("tr(p1) = 1/delta^2 = el deficit de holonomia (cierre de una proposicion basica)", sp.simplify(T3.trace(p1) - 1 / d ** 2) == 0))

# torre: idempotentes de Jones-Wenzl por la recursion de Wenzl en TL_k, k = 1..4
def jw(k):
    T = TL(k)
    if k == 1: return T, T.identity()
    Tprev, f = jw(k - 1)
    # incrustar f en TL_k (añadir una cadena vertical a la derecha)
    def embed(x):
        out = [sp.Integer(0)] * T.n
        for i, a in enumerate(x):
            if a == 0: continue
            D = Tprev.basis[i]; kk = k - 1
            newD = set()
            for pair in D:
                u, v = tuple(pair)
                def mp(p): return p if p < kk else p + 2   # abajo se desplaza 2 (nuevo punto arriba k-1 y abajo k)
                newD.add(frozenset((mp(u), mp(v))))
            newD.add(frozenset((k - 1, k)))   # nueva cadena vertical: arriba k-1 <-> abajo (2k-1-(k-1)) = k
            out[T.index[frozenset(newD)]] += a
        return out
    F = embed(f)
    ek = T.e(k - 2)
    # recursion de Wenzl: JW_k = JW_{k-1} - ([k-1]/[k]) JW_{k-1} e_{k-1} JW_{k-1}
    q = sp.symbols("q")
    def qint(m): return sp.simplify(sum(q ** (m - 1 - 2 * j) for j in range(m)))   # [m]_q
    coef = sp.simplify(qint(k - 1) / qint(k))
    coef = sp.simplify(coef.subs(q, (d + sp.sqrt(d ** 2 - 4)) / 2))
    JW = T.add(F, T.mul(T.mul(F, ek), F), 1, -coef)
    return T, [sp.simplify(a) for a in JW]

q = sp.symbols("q")
def qint_delta(m):
    expr = sum(q ** (m - 1 - 2 * j) for j in range(m))
    return sp.simplify(sp.expand(expr.subs(q, (d + sp.sqrt(d ** 2 - 4)) / 2)))
closures = {}
for k in [1, 2, 3]:
    T, JW = jw(k)
    # cierre (traza de Markov sin normalizar) de JW_k = [k+1]_q
    cl = sp.simplify(T.trace(JW) * d ** k)
    target = qint_delta(k + 1)
    ok = sp.simplify(cl - target) == 0 or sp.simplify(sp.expand_func(cl - target)) == 0 or all(abs(complex((cl - target).subs(d, v))) < 1e-9 for v in [sp.Rational(3, 2), 2, sp.sqrt(2), sp.Rational(9, 5)])
    closures[k] = (cl, target); log.append((f"cierre de JW_{k} = [{k+1}]_q", ok))
    log.append((f"JW_{k} idempotente", T.is_zero([sp.simplify(a) for a in T.add(T.mul(JW, JW), JW, 1, -1)])))

# positividad de la torre en funcion de delta: signos de [m]_q para delta = 2cos(pi/n) y entre valores de la serie
def qints_at(delta_val, M=8):
    vals = [1, delta_val]
    for m in range(2, M + 1): vals.append(sp.simplify(delta_val * vals[-1] - vals[-2]))
    return vals[:M]   # [1],[2],...,[M]
pos = []
for n in [3, 4, 5, 6]:
    dv = 2 * sp.cos(sp.pi / n); v = [sp.nsimplify(x) for x in qints_at(dv, n + 2)]
    signs = [sp.sign(sp.N(x, 30)) for x in v]
    pos.append((f"delta=2cos(pi/{n})", [str(sp.nsimplify(x)) for x in v], "todos >= 0 hasta [n-1], [n] = 0" if all(s >= 0 for s in signs[:n - 1]) and abs(float(sp.N(v[n - 1]))) < 1e-12 else "NO"))
for dv, name in [(sp.Rational(17, 10), "delta=1.7 (entre 2cos(pi/5)=1.618 y 2cos(pi/6)=1.732)"), (sp.Rational(19, 10), "delta=1.9 (entre 2cos(pi/9)=1.879 y 2cos(pi/10)=1.902)"), (sp.Integer(2), "delta=2"), (sp.Rational(5, 2), "delta=2.5")]:
    v = qints_at(dv, 14); signs = [float(sp.N(x)) for x in v]
    pos.append((name, [round(s, 4) for s in signs], "algun [m] < 0 -> positividad falla" if any(s < 0 for s in signs) else "todos >= 0"))

if __name__ == "__main__":
    for name, ok in log: print(("OK  " if ok else "FAIL"), name)
    for k, (cl, tg) in closures.items(): print(f"   cierre JW_{k} = {sp.simplify(cl)}   [k+1]_q = {tg}")
    for row in pos: print("  ", row)
