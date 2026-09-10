# fieldbit

Exact symbolic layer for sheaf-theoretic contextuality on finite sites: sites, empirical models (states), the
boundary `D = -ln NCF` computed as a linear program **with an exactly verified certificate**, bi-Heyting sieves
(`¬`, `∼`, `∂`), Lawvere–Tierney modalities, phase cochains, and the operations *restrict, union, tensor,
identify, compose*. The composition laws proved in the accompanying notes are shipped as **tests that must pass
in exact arithmetic**.

Supporting code for: *Contexts as covers: boundary, naturality and non-compositionality in a linguistic topos*
(I. M. Ozcáriz Arraiza, 2026), section "An arithmetic for the boundary".

## What it computes
- `Site`, `State` (entries as `Fraction` or symbolic expressions in a noise parameter).
- `Frontier.exact()` — the non-contextual fraction `NCF` as an exact rational with a verified dual certificate.
  Backends: pure rational simplex (`exact`), Sage `MixedIntegerLinearProgram(solver="PPL")` (`ppl`), or a
  numeric proposal (HiGHS) that is rationalised and verified (`scipy`). `auto` picks PPL inside Sage.
- `Frontier.certify(primal, dual, param, domain)` — proves a *law* in the parameter as a polynomial identity.
- `Sieve` with `neg` (Heyting), `coneg` (co-Heyting), `boundary`; `Modality.open/closed` (Lawvere–Tierney).
- `Cochain` (Z_m phases): induced 1-cochain, torsion, phase state.

## Laws shipped as tests (all exact)
| test | value | law |
|---|---|---|
| CHSH, p = 1/10 | NCF = 1/5 | n·p/2 |
| union of two cycles | 1/5 = min | D(⊔) = max |
| tensor | 2/25 = product | D(⊗) = sum (general) |
| one shared observable | 1/40 = (5/8)(1/25) | identification term ln(8/5) |
| Mermin star composed with itself, symbolic p | NCF(e∘e) = p, certified for 0 < p ≤ 2/3 | halving (D ≥ D + ln 2) |
| CHSH gluing sieve | (|S|, |∂S|, |¬S|) = (15, 15, 0) | degenerate boundary |
| uniform Z₈ cochain on the 4-cycle | torsion 8; Tsirelson correlators | 2n \| m |
| closed / open modalities on CHSH sieves | Lawvere–Tierney axioms | — |

Run: `python tests/run_tests.py` (13 tests). Inside a Sage notebook, run `preparser(False)` first or use
`sage -python`; the code is plain Python and robust to Sage integers.

`tl_check.py` verifies, in exact arithmetic on planar diagrams, that a category of returns with a positivity
axiom is Temperley–Lieb and that positivity on the Jones–Wenzl tower selects δ = 2cos(π/n) ∪ [2,∞) (Jones 1983).

## Install
```
pip install -r requirements.txt
python tests/run_tests.py
```

## Status
Prototype (v0.1). Exact and tested; not optimised. Frontier via the pure rational simplex is fine up to ~128
global assignments; larger scenarios use the verified numeric proposal.

## Licence
MIT. Cite via `CITATION.cff`.
