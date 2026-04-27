#set page(flipped: true, margin: 1.5cm)

// ── Macros de variables de estado ────────────────────────────────────────────
#let TGFBRA = $"TGF"beta"RA"$
#let TGFB = $"TGF"beta$
#let TGFBRI = $"TGF"beta"RI"$
#let RSMADI = $"RSMADI"$
#let RSMADA = $"RSMADA"$
#let SMURFSMAD = $"SMURF/SMAD7"$
#let SMURF = $"SMURF"$
#let SMAD = $"SMAD7"$

// ── Función auxiliar: parámetro coloreado en ecuación ────────────────────────
// Uso dentro de math: #pc(color, $simbolo$)
#let pc(c, body) = text(fill: c, body)

// ── Paleta de colores ─────────────────────────────────────────────────────────
#let C = (
  delta1: rgb("#e74c3c"), // rojo
  delta2: rgb("#e67e22"), // naranja
  alpha1: rgb("#2980b9"), // azul
  alpha2: rgb("#8e44ad"), // púrpura
  alpha3: rgb("#16a085"), // verde azulado
  alpha0: rgb("#d4a017"), // amarillo-ocre  (α sin subíndice)
  gamma1: rgb("#27ae60"), // verde
  N1: rgb("#c0392b"), // rojo oscuro
  K1: rgb("#7d3c98"), // morado oscuro
  K2: rgb("#117a65"), // verde oscuro
  Omega1: rgb("#1a5276"), // azul marino
  Omega2: rgb("#784212"), // marrón
  Omega3: rgb("#1f618d"), // azul acero
)

// ── Parámetros coloreados (reutilizables en ecuaciones con #nombre) ───────────
#let d1 = pc(C.delta1, $delta_1$)
#let d2 = pc(C.delta2, $delta_2$)
#let a1 = pc(C.alpha1, $alpha_1$)
#let a2 = pc(C.alpha2, $alpha_2$)
#let a3 = pc(C.alpha3, $alpha_3$)
#let a0 = pc(C.alpha0, $alpha$)
#let g1 = pc(C.gamma1, $gamma_1$)
#let N1 = pc(C.N1, $N_1$)
#let K1 = pc(C.K1, $K_1$)
#let K2 = pc(C.K2, $K_2$)
#let Om1 = pc(C.Omega1, $Omega_1$)
#let Om2 = pc(C.Omega2, $Omega_2$)
#let Om3 = pc(C.Omega3, $Omega_3$)

#let hill = $op("Hill")$

// ── Lista de variables ────────────────────────────────────────────────────────
== Lista de variables de estado

- $TGFBRA$: Receptor de TGF#sym.beta *activo*
- $TGFBRI$: Receptor de TGF#sym.beta *inactivo*
- $RSMADI$: Rsmad inactivo
- $RSMADA$: Rsmad activo
- $SMURFSMAD$: Complejo SMURF/SMAD7

// ── Leyenda de colores ────────────────────────────────────────────────────────
== Leyenda de parámetros

#let legend-entry(sym-display, color, label) = {
  box(fill: color, width: 0.7em, height: 0.7em, radius: 1pt)
  h(0.35em)
  sym-display
  h(0.5em)
  text(fill: luma(80))[#label]
  h(1.2em)
}

#block(
  inset: (x: 0.5em, y: 0.6em),
  stroke: luma(200),
  radius: 4pt,
)[
  #legend-entry($delta_1$, C.delta1, [tasa de degradación mediada por SMURF/SMAD7 sobre TGF#sym.beta RA])
  #legend-entry($delta_2$, C.delta2, [tasa de degradación mediada por SMURF/SMAD7 sobre RSMAD#sub[A]])
  \ #v(0.3em)
  #legend-entry($alpha_1$, C.alpha1, [tasa de producción de SMAD7])
  #legend-entry($alpha_2$, C.alpha2, [tasa de formación de SMURF/SMAD7])
  #legend-entry($alpha_3$, C.alpha3, [tasa de producción de SMURF])
  #legend-entry($alpha$, C.alpha0, [tasa de degradación de SMAD7 por SMURF libre])
  \ #v(0.3em)
  #legend-entry($gamma_1$, C.gamma1, [tasa de inactivación de RSMAD])
  \ #v(0.3em)
  #legend-entry($N_1$, C.N1, [coeficiente de Hill para TGF#sym.beta RA])
  #legend-entry($K_1$, C.K1, [constante de Hill para TGF#sym.beta RA])
  #legend-entry($K_2$, C.K2, [constante de Hill para RSMAD#sub[I]])
  \ #v(0.3em)
  #legend-entry($Omega_1$, C.Omega1, [tasa de degradación de SMAD7])
  #legend-entry($Omega_2$, C.Omega2, [tasa de disociación de SMURF/SMAD7])
  #legend-entry($Omega_3$, C.Omega3, [tasa de degradación de SMURF])
]

// ── Sistema de ecuaciones ─────────────────────────────────────────────────────
== Modelo

$
     TGFBRA & = hill(TGFBRA, #N1, #K1) dot TGFB
              - #d1 dot SMURFSMAD dot TGFBRA \
     RSMADA & = (TGFBRA dot hill(RSMADI, 1, #K2))
              dot (1/(1 + SMAD))
              - #g1 RSMADA
              - #d2 RSMADA dot SMURFSMAD \
     RSMADI & = #g1 dot RSMADA
              - TGFBRA dot hill(RSMADI, 1, #K2) \
       SMAD & = #a1 dot RSMADA
              - #Om1 dot SMAD
              - #a0 dot SMAD dot SMURF
              + #Om2 dot SMURFSMAD \
  SMURFSMAD & = #a2 dot SMAD dot SMURF
              - #Om2 dot SMURFSMAD \
      SMURF & = #a3 dot RSMADA
              - #Om3 dot SMURF
              - #a2 dot SMAD dot SMURF
              + #Om2 dot SMURFSMAD
$

*Recordatorio:* $display(hill(X, n, K) = X^n / (K^n + X))$
