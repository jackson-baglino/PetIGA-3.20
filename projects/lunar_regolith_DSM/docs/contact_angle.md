# Prescribed contact angle: the wall free-energy term

Full derivation of the boundary term added in `src/assembly.c`. Written down
because the result is surprising — neither ε nor any surface energy survives
into the term — and someone will otherwise "fix" it.

## 1. What the interior form already is

The ice residual (`src/assembly.c`) is

```
R_ice = N·φ_t + 3·M·ε·∇N·∇φ + (3M/ε)·f1·N − source,      M = mob_sub
```

with `f1 = W'(φ)` and `W(φ) = ½φ²(1−φ)²` — the deliberate half-normalisation
documented above `DoubleWellDeriv`, which is what makes the equilibrium profile
come out as `φ = ½(1 + tanh(x/2ε))` with interface-thickness parameter exactly ε.

So the interior form is the gradient flow

```
∂φ/∂t = −(3M/ε) · δ𝔉/δφ,      𝔉[φ] = ∫_Ω [ (ε²/2)|∇φ|² + W(φ) ] dΩ
```

Note 𝔉 is **dimensionless in the energy sense** — no γ appears. γ_ia enters the
solver only through `-d0_sub0 → lambda_sub → tau_sub → mob_sub/alph_sub`.

## 2. The normalisation constant

On the equilibrium profile the two terms of 𝔉 are equal (equipartition, from
the first integral `(ε²/2)φ'² = W`), and

```
∫_{-∞}^{∞} [ (ε²/2)φ'² + W ] dx = ε²∫φ'² dx = ε/6
```

using `φ' = φ(1−φ)/ε` and `∫sech⁴u du = 4/3`. Verified numerically to 3e-11.

So 𝔉 carries an interfacial excess of **ε/6 per unit area**, and the physical
free energy is

```
F_phys = (6 γ_ia / ε) · 𝔉
```

## 3. The wall term

Physically,

```
f_w(φ) = γ_as + (γ_is − γ_as) h(φ) = γ_as − γ_ia cos θ · h(φ)
h(φ) = φ²(3−2φ),   h(0)=0, h(1)=1, h'(0)=h'(1)=0
```

so `f_w(0) = γ_as`, `f_w(1) = γ_is`, and Young's equation
`γ_ia cos θ = γ_as − γ_is` defines θ. Because `h'` vanishes at both pure phases,
the term is inert in the bulk and acts only at the contact line.

To enter 𝔉 it must be divided by the same normalisation:

```
f̃_w = (ε / 6γ_ia) · f_w = (ε γ_as / 6γ_ia) − (ε/6) cos θ · h(φ)
```

The constant drops out under variation and **γ_ia cancels**:

```
δf̃_w/δφ = −(ε/6)·cos θ·h'(φ) = −(ε/6)·cos θ·6φ(1−φ) = −ε·cos θ·φ(1−φ)
```

## 4. The natural boundary condition

Stationarity of the boundary terms of `δ𝔉` gives

```
ε² ∂φ/∂n + δf̃_w/δφ = 0   ⟹   ∂φ/∂n = cos θ · φ(1−φ)/ε
```

Since `|∇φ| = φ(1−φ)/ε` on the equilibrium profile and `∇φ = |∇φ|·m` with `m`
the interface normal into the ice,

```
m · n = cos θ
```

with `n` the outward wall normal. That is the contact angle measured through
the ice. Check the limits: `cos θ = +1` puts `m = n`, i.e. the interface normal
points straight out of the wall — ice sheeted flat onto it, θ = 0.
`cos θ = −1` gives θ = 180°, a point contact.

## 5. The residual term

The Galerkin weak form does not integrate the gradient term by parts, so

```
⟨δ𝔉/δφ, N⟩ = ∫_Ω [ε²∇N·∇φ + W'N] dΩ + ∮_Γ (δf̃_w/δφ) N dΓ
```

Multiplying by `(3M/ε)` reproduces the interior coefficients `3Mε` and `3M/ε`
exactly, and on the boundary **cancels the remaining ε**:

```
R_bnd[a][0] = (3M/ε)·(−ε cos θ φ(1−φ))·N0[a] = −3·M·cos θ·φ(1−φ)·N0[a]

J_bnd[a][0][b][0] = −3·M·cos θ·(1−2φ)·N0[a]·N0[b]
```

No ε. No γ. `cos θ` is the whole of it. The Jacobian has no `shift` term — the
wall energy has no `φ_t` dependence.

**Sign.** `cos θ > 0` (γ_as > γ_is, wetting) makes `R_bnd < 0`; since
`R = N·φ_t + … = 0`, that gives `φ_t > 0` — ice grows at the wall, the contact
line advances, θ decreases. Gated as G4 in
`studies/contact_angle/verification/verify_wall_bc.sh`.

## 6. Cross-check against Metamorph

`demo/Metamorph.c` carries a contact-angle term in this same PetIGA tree:

```c
R[a][0] = -N0[a]*3.0*eps*mob*modgradice*costhet;
```

Since `ε|∇φ| = φ(1−φ)` on the equilibrium profile, that is **identical** to the
expression above. Two independent routes to the same coefficient.

The `h(φ)` form is nevertheless the better one:

- no `1/|∇φ|` regularisation. Metamorph needs `if (modgradice < 1.0e-5)
  modgradice = 1.0e-5;` to avoid dividing by zero in its Jacobian; `h'(φ)` is a
  polynomial and needs nothing.
- its Jacobian is an exact mass-matrix block rather than a normalised-gradient
  derivative.
- it is exact away from the equilibrium profile too, where `ε|∇φ| ≠ φ(1−φ)`.

## 7. Why γ_ia defaults to `-Sigma_i`

`monitoring.c` forms the capillary length as `d0 = Etai·V_m/(R·T)`. Inverting it
at the *physical* `d0_sub0 = 1.0166e-9 m` and −20 °C:

```
γ = 1.0166e-9 · 8.314 · 253.15 / 1.963e-5 = 0.109 J/m²
```

which is exactly `Sigma_i`. So `Sigma_i` is the ice–vapour surface energy this
solver already implies, and it is the self-consistent default for `γ_ia`.

⚠️ `Sigma_a = 0.132` in this project looks like a bug — `enceladus_DSM` sets
both to 0.109 with the comment *"same interface, air side: must equal Sigma_i"*.
It plays no part in the wall term, so it is flagged rather than fixed here.
