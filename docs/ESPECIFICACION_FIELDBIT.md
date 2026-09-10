# Especificación de la capa simbólica `fieldbit`

Prototipo funcional: `fieldbit.py` (sympy + fractions; scipy sólo como guía numérica). Ejecutado: 8 tests-teorema, todos exactos. Destino: paquete sobre Sage, donde posets, retículos, complejos simpliciales, cohomología y poliedros exactos ya existen.

## 1. Principio
El sistema devuelve **leyes y certificados**, no números. Toda cantidad se calcula en ℚ o en el anillo simbólico; el solver numérico sólo propone una solución que después se **racionaliza y se verifica en aritmética exacta** (primal factible, dual factible, holgura complementaria). Si no certifica, falla; no devuelve un decimal. Las leyes en un parámetro (el ruido p, la fuerza de un acoplamiento) se demuestran como identidades polinómicas sobre un dominio, no se ajustan.

Lo que esto elimina: el suelo espurio de D por muestreo, los detectores ajustados a mano, los promedios que fabrican niveles. Es lo que nos costó tres experimentos.

## 2. Objetos
| clase | contenido | en Sage |
|---|---|---|
| `Site` | observables, contextos, poset de subcubiertas, nervio | `Poset`, `SimplicialComplex` |
| `State` | modelo empírico; entradas `Fraction` o expresiones simbólicas | `QQ`, `SR` |
| `Frontier` | NCF y D como LP; `exact()` devuelve racional + certificado dual verificado; `certify(primal, dual, param, domain)` demuestra una ley simbólica | `Polyhedron(base_ring=QQ)` (vértices exactos, LP paramétrico sin guía numérica) |
| `Sieve` | criba sobre subcubiertas; `neg` (Heyting), `coneg` (co-Heyting), `boundary` | `LatticePoset` con pseudocomplementos y su dual |
| `Cochain` | fases Z_m sobre observables; cocadena inducida (diferencia o suma), holonomía, torsión, estado de fases | `CyclicPermutationGroup`, cohomología del nervio con coeficientes en Z_m |
| operaciones | `restrict`, `union`, `tensor` (con `identify_obs`), `compose` (regla de complementariedad, reducida por r2 = r1) | constructores de `State` |

Pendiente de definir como objeto: `Modality(j)` (topología de Lawvere–Tierney sobre el poset de subcubiertas; el endofuntor idempotente del corpus), con `Frontier` relativa a j.

## 3. Los teoremas como tests (ejecutados)
| test | valor exacto | teorema |
|---|---|---|
| CHSH, p = 1/10 | NCF = 1/5 | ley n·p/2 |
| unión de dos ciclos con p = 1/10 y 1/5 | 1/5 = min | regla ⊔ |
| yuxtaposición | 2/25 = (1/5)(2/5) | multiplicatividad general |
| identificación de un observable | 1/40 = (5/8)(1/25) | término ln(8/5) |
| composición de la estrella, p simbólico | NCF(e∘e) = p, certificado con primal p/32 en singles y dual 1/3 en celdas (C malo, C′ bueno), como identidad en p con dominio 0 < p ≤ 2/3 | halving |
| criba de pegado de CHSH | (|S|, |∂S|, |¬S|) = (15, 15, 0) | frontera degenerada, ¬S = ∅ |
| cocadena Z₈ uniforme en el 4-ciclo | torsión 8; correladores (√2/2, √2/2, √2/2, −√2/2) | regla 2n \| m |

Una conjetura nueva entra al sistema como test que **debe fallar** hasta que se demuestre; una ley se registra sólo con su certificado.

## 4. Qué haría con Physarum y Del Vecchio
- Physarum: en lugar de Monte Carlo sobre S(t), la reducción de Ott–Antonsen (Chaos 18, 037113, 2008) da la dinámica del parámetro de orden de un Kuramoto con distribución de frecuencias en forma cerrada; la pregunta "¿el bloqueo destruye la memoria de τ?" se responde con el espectro de la variedad de Ott–Antonsen bajo forzamiento periódico, sin detector.
- Del Vecchio: el modelo bien mezclado es un sistema de ecuaciones polinómicas; sus equilibrios y su estabilidad en función de km′ se calculan con bases de Gröbner o resultantes en Sage, y "no hay régimen multiestable" se demuestra, no se observa. El modelo espacial sigue siendo estocástico; la cantidad que sí es simbólica es la probabilidad de absorción como solución de un sistema lineal sobre el espacio de patrones (5^15 estados: demasiado; sobre las clases de simetría, quizá no).

## 5. Límite declarado
Es emulación exacta, no procesamiento analógico: un LP simbólico describe una frontera, no la atraviesa. Su valor es negativo (no puede confirmar lo inexacto) y estructural (convierte medidas en leyes). El procesamiento real es la otra mitad: `ANYONES_Y_RETORNOS.md`.
