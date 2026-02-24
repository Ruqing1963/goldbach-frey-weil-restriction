# Goldbach–Frey Jacobians as Weil Restrictions

**Paper #15 in the Titan Project series**

Trace vanishing, Asai *L*-functions, and modularity via Q(*i*).

## Key Results

**Theorem 2.1 (Trace Vanishing Law):** For the Goldbach–Frey curve *C*: *y*² = *x*(*x*²−*p*²)(*x*²−*q*²), the Frobenius trace *a*ᵣ = 0 at every good prime *r* ≡ 3 (mod 4).

**Proposition 6.1 (Weil Restriction):** Over Q(*i*), Jac(*C*) is (2,2)-isogenous to *E*₁ × *E*₂, where
- *E*₁: *Y*² = *X*(*X* − *p*²)(*X* − *q*²)
- *E*₂: *Y*² = *X*(*X* + *p*²)(*X* + *q*²)

and *E*₂ ≅ *E*₁^(−1) is the quadratic twist by −1.

**Theorem 7.1 (Modularity):** Jac(*C*)/Q is modular: there exists a weight-2 Siegel paramodular form *F* with *L*(Jac(*C*), *s*) = *L*(*F*, *s*), obtained by Asai transfer from GL₂(Q(*i*)) to GSp₄(Q).

## Logical Chain

```
p + q = 2N  ⟹  f(−x) = −f(x)  ⟹  V_ℓ = Ind(V⁺)  ⟹  Jac(C) ~ Res(E)
```

The Goldbach constraint forces the Jacobians to be Weil restrictions of elliptic curves over Q(*i*), reducing modularity to the known modularity of *E*/Q(*i*) (Allen et al.).

## Repository Structure

```
├── paper/
│   ├── Paramodular_Prospects.pdf    Final paper (7 pages)
│   └── Paramodular_Prospects.tex    LaTeX source
├── figures/
│   ├── fig_trace_vanishing.pdf      Figure 1: trace vanishing law
│   └── fig_cross_validate.pdf       Figure 2: cross-validation
├── scripts/
│   ├── paramodular.py               Paramodularity analysis
│   ├── sato_tate_analysis.py        Trace-zero phenomenon analysis
│   └── fig_paramodular.py           Figure generation
├── README.md
└── LICENSE
```

## Series Context

| # | Paper | Key Result |
|---|-------|------------|
| 12 | True Conductor Validation | Cond_odd = [rad_odd(pqM\|p−q\|)]² |
| 13 | Universal Tame Semistability | f_r = 2 at all odd primes |
| 14 | Conductor Census | 425,082 pairs, bandwidth stability |
| **15** | **Weil Restrictions (this paper)** | **Jac(C) ~ Res(E), Asai lift, modularity** |

## License

MIT
