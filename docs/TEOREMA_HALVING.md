# Teorema de halving para la composición

Script: `frontera/halving_proof.py`; datos: `halving_proof.json`. LP del corpus.

## Enunciado
Sea un escenario de paridad con n contextos, todos de k observables, ruido uniforme p (e_C(s) = (2−p)/2^k en resultados de paridad correcta, p/2^k en los de paridad errónea), compuesto consigo mismo por la regla de complementariedad (observable compartido repite exacto; no compartido sale uniforme e independiente).

**Teorema.** Si ninguna asignación global viola los n contextos a la vez, entonces

NCF(e ∘ e) ≤ n·p/4.

**Corolario.** Si además NCF(e) = n·p/2 (régimen p ≤ 2/n), entonces D(e ∘ e) ≥ D(e) + ln 2: **componer un recurso consigo mismo cuesta al menos un nat-bit de deuda**, en todo escenario de paridad que cumpla la hipótesis.

**Igualdad (halving).** NCF(e ∘ e) = n·p/4 si y sólo si existe un primal con soporte en asignaciones que violan exactamente 1 o exactamente n−1 contextos y que satura todas las celdas un-malo. Se da en la estrella de Mermin, en CHSH y en Peres–Mermin (primales explícitos abajo); no se da en los ciclos n ≥ 5.

## Prueba de la cota
Tras la identificación g^{r2} = g^{r1} (forzada por la repetición), el LP reducido tiene por variables las asignaciones g de los |X| observables y por restricciones, para cada par ordenado (C, C′), las celdas u sobre C ∪ C′ con capacidad e_C(u|_C)/2^{|C′∖C|}. Para un par no ordenado {C, C′} con j = |C ∩ C′| la celda tiene 2k − j bits y su capacidad efectiva (mínimo de los dos órdenes) es p/2^{2k−j} si u es de paridad errónea en al menos uno de los dos contextos.

*Celdas un-malo.* Las paridades de C y de C′ son dos funcionales lineales distintos sobre F₂^{2k−j}, luego cada patrón (buena/mala, buena/mala) ocupa exactamente un cuarto de las 2^{2k−j} celdas. Las celdas malas en exactamente uno de los dos contextos son 2·2^{2k−j−2} = 2^{2k−j−1}, y su capacidad total por par es 2^{2k−j−1}·p/2^{2k−j} = **p/2, independiente de j y de k**.

*Cubrimiento.* Una asignación g que viola v contextos está en exactamente v(n−v) celdas un-malo (una por cada par con un contexto violado y otro satisfecho). Si 1 ≤ v ≤ n−1, entonces v(n−v) ≥ n−1. Por hipótesis v ≠ n, y v ≥ 1 porque el escenario es contextual (un contexto impar). Luego el peso λ = 1/(n−1) en todas las celdas un-malo da Σ ≥ 1 sobre toda g: es un cubrimiento fraccionario factible, de coste

C(n,2) · (p/2) · 1/(n−1) = n·p/4. ∎

## Necesidad de la hipótesis
En el 3-ciclo la asignación que viola los tres contextos existe (v = 3 = n), no está en ninguna celda un-malo, el dual es infactible y la cota **falla**: NCF(e∘e) = p > 3p/4. Lo mismo ocurre en los ciclos impares n = 5, 7 (v = n posible), donde la cota se cumple numéricamente pero no por este argumento.

## Igualdad: holgura complementaria
Con λ > 0 uniforme sobre todas las celdas un-malo, un primal es óptimo si y sólo si (i) su soporte está en asignaciones con v(n−v) = n−1, es decir v ∈ {1, n−1}, y (ii) satura todas las celdas un-malo. Primales explícitos, verificados restricción a restricción:
- **CHSH** (n = 4, k = 2): masa p/16 en **las 16 asignaciones** (todas tienen v ∈ {1, 3}); total p = 4p/4. Saturación: celdas de par disjunto (1 asignación, capacidad p/16) y de par adyacente (2 asignaciones, p/8).
- **Estrella de Mermin** (n = 4, k = 3): masa p/32 en las 32 asignaciones de violación única; total p (`punto3/DERIVACION_COMPOSICION.md`).
- **Peres–Mermin** (n = 6, k = 3): masa p/64 en las 96 asignaciones de violación única; total 1.5p = 6p/4; factible con exceso máximo 0.

En los ciclos n = 6, 8 (hipótesis cumplida) el óptimo del LP tiene soporte con v = 3 ∉ {1, n−1}: la condición (i) no puede satisfacerse y la cota no es tensa (13p/12 y 11p/9 frente a 1.5p y 2p).

## Lectura
- El halving no es una coincidencia de tres escenarios: es la cota general del teorema alcanzada. Lo que distingue estrella, CHSH y PM de los ciclos largos es una propiedad combinatoria precisa: que las celdas un-malo se puedan saturar con asignaciones de una o de n−1 violaciones.
- La constante ln 2 es estructural: proviene de que la capacidad un-malo por par es p/2 sea cual sea el solape, y del peso 1/(n−1) sobre C(n,2) pares. No depende de la profundidad ni del ruido, lo que explica por qué ambas conjeturas fallaron.
- Para AQC: reutilizar un recurso contextual de paridad sobre sí mismo cuesta *al menos* ln 2 en todo escenario sin violación total, y exactamente ln 2 cuando la estructura permite saturar; más en los ciclos largos. La contabilidad de amortización queda acotada por abajo con teorema.
- Queda abierta la forma cerrada de los ciclos n ≥ 5 (sucesión 17/16, 13/12, 19/16, 11/9 en unidades de p) y la extensión a contextos de tamaño distinto.
