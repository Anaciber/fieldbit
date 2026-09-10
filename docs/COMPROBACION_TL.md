# Comprobación: la categoría de retornos del corpus es Temperley–Lieb, y su positividad es la de Jones

Script: `tl_check.py` (álgebra de diagramas planos implementada en exacto; sympy). Todas las relaciones se verifican como identidades en δ, no numéricamente.

## 1. Los axiomas del corpus, extraídos de §10.2–10.4
- (A1) Un solo objeto, yuxtaponible consigo mismo; abrir produce dos extremos; composición en sucesión y en yuxtaposición.
- (A2) Zigzag abrir–cerrar = identidad girada por una vuelta de cinta; cerrar contra abrir = lazo cerrado = escalar δ.
- (A3) La estructura libre es la de los enredos planos de una cinta consigo misma con el lazo como único escalar. (El corpus lo declara: "esta estructura existe y está estudiada".)
- (A4) Proposiciones = cerrar-y-reabrir normalizado por el lazo, idempotentes, no conmutativos sobre pares que comparten extremo; holonomía p·q·p = δ⁻²·p; torre de proposiciones sobre bloques crecientes.
- (A5) Positividad del cierre: cerrar un retorno (una proposición) sobre sí mismo nunca da un número negativo.

No hay más axiomas en §10. El §10.5.5 deja explícitamente fuera los cruces de la cinta ("permitir cruces de la cinta añade una dimensión").

## 2. Identificación de la estructura libre (A1–A4) con Temperley–Lieb
La categoría cuyos objetos son n extremos, cuyos morfismos son enredos planos de cintas y cuyo único escalar es el lazo δ es la **categoría de Temperley–Lieb TL(δ)**, y es un teorema que TL(δ) es la categoría pivotal libre generada por un objeto autodual de dimensión δ (Kauffman–Lins 1994; Turaev, *Quantum Invariants*, cap. XII; en lenguaje moderno, Etingof–Gelaki–Nikshych–Ostrik). Las proposiciones del corpus son los generadores normalizados p_i = e_i/δ, con e_i = cap_i ∘ cup_i, y sus relaciones son las relaciones definitorias de TL.

Verificado en exacto en TL₃ y TL₄ (`tl_check.py`):
| relación del corpus | relación de TL | resultado |
|---|---|---|
| p idempotente | e_i² = δ·e_i | OK |
| p·q·p = δ⁻²·p para pares que comparten extremo | e_i e_{i±1} e_i = e_i | OK |
| p, q no conmutan | e_i e_{i+1} ≠ e_{i+1} e_i | OK |
| proposiciones sobre pares disjuntos conmutan | e_i e_j = e_j e_i, \|i − j\| ≥ 2 | OK |
| cierre de p | tr(p_i) = 1/δ² | OK |

Un detalle que el corpus no había escrito: **el cierre de una proposición básica es exactamente su déficit de holonomía**, tr(p) = δ⁻². "El bucle restituye la proposición atenuada por un escalar fijado por la estructura" y "lo que vale la proposición cerrada sobre sí misma" son el mismo número.

**El giro del zigzag (A2).** En una categoría plana no hay curvas cerradas con cruce; "una vuelta de cinta" es un automorfismo escalar del objeto generador, y la única libertad que deja es un signo: es la elección de estructura pivotal, es decir, el indicador de Frobenius–Schur del objeto autodual (±1). Para SU(2) el generador es simpléctico (indicador −1), lo que en el corchete de Kauffman aparece como δ = −A² − A⁻². La "vuelta de cinta" del corpus es ese signo, no una estructura adicional.

**La torre de proposiciones (A4)** es la sucesión de idempotentes de Jones–Wenzl JW_k. Verificado: JW₁, JW₂, JW₃ construidos por la recursión de Wenzl son idempotentes y sus cierres son δ, δ² − 1, δ(δ² − 2), es decir, los enteros cuánticos [2]_q, [3]_q, [4]_q con δ = q + q⁻¹.

## 3. La positividad (A5) es el teorema del índice de Jones
"Cerrar una proposición sobre sí misma nunca da un número negativo" es, sobre la torre, la condición **[m]_q ≥ 0 para todo m**, con [m]_q = sen(mπ/n)/sen(π/n) cuando δ = 2cos(π/n). Verificado:
- δ = 2cos(π/n), n = 3, 4, 5, 6: [1], …, [n−1] ≥ 0 y **[n] = 0**: la proposición de nivel n tiene cierre nulo ("el vacío que retorna y no cuenta"), y a partir de ahí los cierres son negativos. La torre se trunca en n − 1: es el ideal despreciable de TL en raíz de la unidad, y su cociente es la categoría de Temperley–Lieb–Jones TLJ(n).
- δ = 1.7 (entre 2cos(π/5) y 2cos(π/6)): [6]_q < 0. δ = 1.9 (entre 2cos(π/9) y 2cos(π/10)): [10]_q < 0. Positividad violada en algún nivel finito: exactamente lo que el corpus llama "lazo con valor prohibido por positividad".
- δ ≥ 2: todos los [m]_q > 0, torre infinita, sin torsión: el régimen que el corpus llama continuo/clásico.

Es el contenido del teorema de Jones (*Invent. Math.* 72, 1–25, 1983): la traza de Markov sobre la torre de álgebras de Temperley–Lieb es definida positiva si y sólo si δ² ∈ {4cos²(π/n)} ∪ [4, ∞), y la demostración pasa por los signos de [m]_q sobre los idempotentes de Jones–Wenzl (Wenzl 1987). El corpus lo re-deriva desde su axioma y llega al mismo conjunto. El corpus debe citarlo.

## 4. Equivalencia
**Proposición.** La categoría de retornos del corpus con los axiomas A1–A5, cocientada por los morfismos de cierre nulo, es equivalente como categoría pivotal a TLJ(n) cuando δ = 2cos(π/n), y a TL(δ) cuando δ ≥ 2.
*Prueba.* A1–A4 definen la categoría pivotal libre sobre un objeto autodual con dimensión δ, que es TL(δ) por el teorema de libertad (§2); A5 es la positividad de la traza de Markov (§3), que por Jones fuerza δ a la serie discreta o a [2, ∞); en la serie discreta los morfismos despreciables forman el ideal generado por JW_n, y el cociente TL(δ)/(JW_n) es TLJ(n), único (Goodman–de la Harpe–Jones 1989, cap. 2). ∎

Condición de la prueba: que A1–A5 sean *todos* los axiomas. Si el corpus tiene alguno más (§10 no lo muestra), la equivalencia se restringe a un cociente.

## 5. Lo que falta para llegar a los anyones: el trenzado
TLJ(n) es una categoría de fusión esférica; las teorías anyónicas son categorías de fusión **trenzadas** (modulares). El corpus se detiene, por su propia declaración (§10.5.5), antes de permitir cruces. El paso que falta es único y canónico: el corchete de Kauffman define un cruce como A·id + A⁻¹·(cup∘cap) con δ = −A² − A⁻², y ese trenzado hace de TLJ(n) la categoría de representaciones de SU(2)_{n−2}, es decir, los anyones de nivel k = n − 2: Ising para n = 4 (δ = √2), Fibonacci para n = 5 (δ = φ, tomando el sector entero de SU(2)_3). La observación del corpus de que "permitir cruces separa cruce de enmarcado" es exactamente la distinción entre el trenzado y el giro (twist) de una categoría de cintas.

Con ello queda establecido, con la salvedad del §4:
- corpus §10 con positividad = TLJ(n) (plano, esférico) = **la parte de fusión** de los anyones de SU(2)_{n−2};
- corpus §10 + cruces (por hacer) = SU(2)_{n−2} completo = **los anyones**;
- la pregunta abierta del corpus sobre "cuántas maneras de cruzar admite la positividad" tiene respuesta en la teoría de categorías trenzadas: para TLJ, exactamente las del corchete de Kauffman con A raíz primitiva de la unidad de orden 4n, salvo conjugación (Turaev–Wenzl).

## 6. Lo que este resultado cambia en el programa
1. El problema abierto nº 1 del corpus (positividad intrínseca) se divide en dos: *cuándo* hay positividad está cerrado (Jones); *por qué* la reflexividad orientada de §15.1 la produce sigue abierto, y es el contenido propio.
2. La "torsión de orden 2n" no es un resultado del juguete 1+1 sino la raíz q = e^{iπ/n} de TLJ; el juguete la encontró por su cuenta, lo que es una confirmación de consistencia, no un teorema nuevo.
3. La coincidencia δ = φ de §16 es el anyón de Fibonacci; la frontera n = 4 / n = 5 del corpus (Z₄ Clifford frente a orden 8) es la frontera Ising / Fibonacci de la universalidad por trenzado. Sigue siendo correspondencia hasta que se calcule D sobre estados anyónicos.
4. El sector no copiable del marco tiene realización física: es la computación cuántica topológica, con la positividad puesta por la física. Lo que Ω puede añadir ahí es lo que ni Jones ni los anyones tienen: la lógica bi-Heyting de los fragmentos y la deuda estructural como función sobre estados; ninguna de las dos está hecha sobre anyones.

## 7. Registro
| afirmación | estado |
|---|---|
| A1–A4 = Temperley–Lieb | verificado en exacto (relaciones, traza, JW) |
| tr(p) = δ⁻² = déficit de holonomía | verificado; no estaba en el corpus |
| A5 = positividad de [m]_q = Jones 1983 | verificado sobre la torre; el corpus debe citar |
| truncamiento en n−1 = ideal despreciable | verificado ([n]_q = 0) |
| equivalencia con TLJ(n) | demostrada condicionada a que A1–A5 sean todos los axiomas |
| anyones = TLJ + trenzado de Kauffman | estándar; el corpus lo deja fuera explícitamente |
| Z₄/orden 8 ↔ Ising/Fibonacci | correspondencia, no teorema |
