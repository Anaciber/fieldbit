# Ley de la frontera bajo cambio de sitio (§7.2 del programa)

Scripts: `frontera/frontier_law.py`, `compose_general.py`; datos: `frontier_law.json`, `compose_general.json`. Todos sobre `ks_shards_dist.ncf` del corpus.

## 1. La forma general del gauge: un cubrimiento fraccionario
Sea H(e) el hipergrafo cuyos vértices son las asignaciones globales g y cuyas hiperaristas son los eventos locales (C, s), con peso e_C(s); el evento (C, s) contiene las g con g|_C = s. El LP de la fracción no contextual es un empaquetamiento fraccionario de vértices con capacidades en las hiperaristas; su dual es

NCF(e) = min { Σ_{(C,s)} λ_{C,s} · e_C(s) : λ ≥ 0, Σ_{(C,s) ∋ g} λ_{C,s} ≥ 1 ∀ g },

el **número de cubrimiento fraccionario ponderado** de H(e). Verificado primal = dual en CHSH (0.4 = 0.4, `1_primal_vs_dual`). Los certificados duales de las derivaciones anteriores (8/5, composición de la estrella) son cubrimientos fraccionarios de este tipo.

D = −ln NCF es, por tanto, el logaritmo de un número de cubrimiento. Esto fija de golpe su comportamiento bajo cambio de sitio, porque los cambios de sitio son operaciones sobre H(e).

## 2. Leyes generales (demostradas)

**(a) Automorfismo del sitio**: isomorfismo de hipergrafos ponderados; NCF invariante. (Teorema de Restricción (i).)

**(b) Restricción y coarse-graining**: eliminar hiperaristas, o fusionarlas sumando pesos, nunca hace más caro cubrir... con más precisión: toda descomposición e = λ e_NC + (1−λ) e′ se empuja a una descomposición del modelo restringido o agregado, luego NCF no decrece y D no crece. (Teorema de Restricción (ii); ABM 2017.)

**(c) Yuxtaposición (⊗, observables renombrados): NCF(e ⊗ e′) = NCF(e)·NCF(e′), para modelos cualesquiera.** Prueba: el producto de dos empaquetamientos factibles es un empaquetamiento factible del producto (la restricción de una celda producto es el producto de restricciones), luego NCF(⊗) ≥ NCF·NCF′; el producto de dos cubrimientos duales factibles es un cubrimiento factible del producto (λλ′ cubre cada par con Σ ≥ 1·1), de coste producto, luego NCF(⊗) ≤ NCF·NCF′. ∎ Esto generaliza la aditividad del corpus §3.2, que allí estaba comprobada sólo en productos de modelos de paridad, a todos los modelos empíricos: **D es aditiva bajo yuxtaposición sin hipótesis**.

**(d) Unión (⊔, sin contextos conjuntos): NCF(e ⊔ e′) = min(NCF(e), NCF(e′)).** Prueba en `punto2/PUNTO2B_RESULTADOS.md` (marginales de masa común); dual: un cubrimiento de un solo factor cubre todos los pares. Verificada con estructura compartida en nueve escenarios.

**(e) Identificación de un observable en ⊗**: NCF = NCF·NCF′·[1/2 + k₁k₂/(2n₁n₂)] para ciclos de paridad (`punto2/DERIVACION_8_5.md`); en general, es el sub-hipergrafo del producto restringido a la diagonal del observable identificado, con capacidades dobladas en las celdas que lo contienen por ambos lados. Fórmula general para varios observables: abierta.

**(f) Composición (∘, medida secuencial) = identificación + refinamiento.** La repetición exacta fuerza g^{r2} = g^{r1} (identificación total); el cambio de contexto parte cada evento (C, s) en celdas más finas (C ∪ C′, u) de peso e_C(s)/2^{|C′∖C|} (refinamiento con conservación del peso). Refinar no baja el cubrimiento y no sube el empaquetamiento: **D(e ∘ e) ≥ D(e) siempre**. Verificado estricto en once casos (estrella, CHSH, PM, ciclos n = 3…8). Refinar *sin* correlación (un observable ficticio uniforme común a todos los contextos) no cambia NCF (`3_refinement_uniform_extra_observable`): lo que crea deuda no es partir, es partir según un observable que el modelo oculto tiene que fijar.

## 3. Composición: lo que se sabe y lo que falló

### 3.1 Dos conjeturas refutadas, registradas antes de calcular
- **Fracción uniforme**: NCF(e∘e) = máx{λ : e ≥ λ·Uniforme}. Vale en la estrella (p), en CHSH (1 − E) y en los ciclos n = 3, 4. **Falla** en Peres–Mermin (1.5p frente a p) y en los ciclos n ≥ 5.
- **ln(n/2)**: D(e∘e) − D(e) = ln(n/2). Vale en n = 4. **Falla** en n = 6 (1.0186 frente a 1.0986) y en PM (0.693 frente a 1.099).

### 3.2 Lo que se sostiene
- **Halving**: NCF(e∘e) = NCF(e)/2, es decir +ln 2, en tres escenarios de estructura distinta: estrella de Mermin (4 contextos de 3 observables, pares que comparten uno), CHSH con ángulos de Tsirelson (4 contextos de 2, con dos pares disjuntos), Peres–Mermin (6 contextos de 3, con seis pares disjuntos). La composición de CHSH es exacta también en el régimen no contextual: NCF(e∘e) = 1 − E (visibilidad), como la estrella con p_eff. La derivación dual de la estrella (`punto3/DERIVACION_COMPOSICION.md`) no se transfiere a CHSH ni a PM tal cual; el halving en los tres pide una prueba común que no tengo.
- **Ciclos** (p = 0.1): ratio NCF(e∘e)/NCF(e) = 2/3, 1/2, 17/40, 13/36, 19/56, 11/36 para n = 3…8; en p: NCF(e∘e) = p, p, 17p/16, 13p/12, 19p/16, 11p/9. Es una sucesión racional sin forma cerrada identificada. Es el problema combinatorio abierto de la composición, y está bien planteado: cubrimiento fraccionario del hipergrafo refinado de un n-ciclo.

## 4. Relación con la valoración del corpus §16
§16 pide una valoración aditiva bajo composición. Con D = ln τ*:
- aditiva bajo **yuxtaponer** (teorema (c), sin hipótesis);
- tropical bajo **unir** (d);
- no aditiva bajo **componer**: el incremento es el logaritmo de un índice de refinamiento que vale 2 en tres escenarios y sigue una sucesión racional en los ciclos.

Propuesta para el corpus: la valoración de §16 es ln τ* y su aditividad es la de ⊗; la composición no es una operación de la valoración sino un cambio de hipergrafo (identificar y refinar), con su propio coste. Si §16 quiere una valoración aditiva bajo componer, la deuda no es ella.

## 5. Preguntas que pueden fallar, actualizadas
1. ¿Es el halving una ley para todo escenario en que cada observable está en exactamente dos contextos y cada contexto tiene ≥ 3 observables, o hay un contraejemplo pequeño? Test: un escenario de 3 observables por contexto que no sea la estrella ni PM.
2. Forma cerrada de NCF(e∘e)/p para n-ciclos. Los datos: 1, 1, 17/16, 13/12, 19/16, 11/9. Candidato a comprobar antes de creerlo: crece hacia un límite < 2.
3. ¿Se transfiere el certificado dual de la estrella (peso 1/3 en celdas "un-malo") a PM con los pesos apropiados? Si sí, el halving tiene una prueba común para escenarios de tres observables por contexto.
4. Identificación de varios observables: contar celdas malo-malo concordantes en todos los compartidos; predicción antes de calcular para dos observables entre dos 4-ciclos: 1/2 + (términos cruzados), no derivada.

## 6. Balance
Lo demostrado en general: (a), (b), (c), (d) y la monotonía de (f). De ellos, (c) es nuevo en generalidad (el corpus lo tenía verificado, no demostrado) y su prueba es de tres líneas sobre el dual. Lo medido y no derivado: el halving fuera de la estrella y la sucesión de los ciclos. Lo refutado: dos conjeturas mías, registradas con su predicción previa. La "ley de la frontera bajo cambio de sitio" existe como teorema para tres de las cuatro operaciones y como problema abierto bien planteado para la cuarta.
