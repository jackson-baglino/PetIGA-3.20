# Back-of-the-envelope neck ODE: which α_c and humidity could reproduce Molaro?

A one-equation model of vapour-route neck growth answers three questions cheaply:
- which α_c the Molaro −20 °C data would need;
- whether a local undersaturation is required;
- what exponent a our model should give in each regime.

It is calibrated on our own Fig. 2 run and checked against the `mesh_pair` run. No HPC time is used.

```bash
python studies/molaro_2019/envelope_ode/neck_ode_fit.py    # -> results/envelope_fits.csv, envelope_fit.png
```

![envelope fit](results/envelope_fit.png)

## Model

**Geometry.**
- Neck radius x, effective grain radius R = 84.4 µm, i.e. 2R₀R₁/(R₀+R₁) for the 101/72.5 µm pair.
- Kuczynski fillet ρ = x²/(2(R − x)).
- The concave neck, of mean curvature 1/x − 1/ρ, is fed from the convex grain surface (2/R).

**Two resistances in series:**
- attachment, β_sub ∝ 1/α_c;
- gas diffusion over a length c·ρ, with k_D = ρ_i/(D_v ρ_vs).

The model also allows a local undersaturation s at the neck:

    v_n   = [ d₀ (1/ρ − 1/x + 2/R) − s ] / ( β_sub + c ρ k_D )
    dx/dt = (π/2) v_n            (small-neck volume and area of the fillet)

**Limits:**
- Attachment-limited (small α_c) gives x³ ∝ t, so a = 1/3.
- Diffusion-limited (large α_c) gives x⁵ ∝ t, so a = 1/5.
- s > 0 stalls the neck once d₀/ρ falls to s, which pulls a lower still.

The model has one geometric unknown, c. Thermodynamics come from `preprocess/comp_eps.py`.

## Results (−20 °C)

| step | what | result |
|---|---|---|
| 1. calibrate | fit c and s to our Fig. 2 run (α_c = 0.1 known) | **c = 0.275**, s = 1.6e-4, 0.85 % rms; a 0.147 (ODE) vs 0.154 (run) |
| 2. validate | predict the `mesh_pair` fine run (α_c = 1e-3, sealed h = 1, tangent start); no refit | reaches 32.81 µm at 54,176 s (run: 54,174 s); 51.1 µm at 64.5 h (run: 78.6 h); a = 0.297 (run: 0.283) |
| 3. fit data | Molaro −20 °C neck widths, α_c and s free | α_c ≈ 0.84, s → 0, 5.0 % rms, **a = 0.197 (data 0.20)** |
| 3. profile | best local s at fixed α_c | α_c = 1: 5.0 % · 0.3: 6.0 % · **0.1: 11.6 %** · 0.03: 21 % · 0.01: 28 % · 1e-3: 32 % rms. Best local s is 0 at every α_c; this is the neck, not the ambient (step 4) |
| 4. shrinkage | large-grain slope dR/dt = −6.25e-10 m/s | far-field undersaturation s_∞ = 3.6e-3 (α_c = 0.1); the Fig. 2 walls are at 2.85e-3 |

The 32.81 µm match in step 2 is closer than the model's accuracy justifies. The 18 % miss at the end of that run is the honest error bar.

## What it says

1. **No literature α_c reproduces the Molaro rate.**
   - Even with no undersaturation at the neck, the best fit at α_c = 0.1 is 11.6 % rms, about 70 % of the observed width increase.
   - Matching needs α_c ≈ 0.3–1. Above 0.1 the gain is small, because gas diffusion caps the rate.
   - Choosing 0.1 is therefore the **tuning choice**: the largest literature value, taken to best reproduce the data. The remaining shortfall is consistent with the missing surface-diffusion share.
2. **The data's exponent (0.20) is what diffusion-limited vapour growth gives with no undersaturation at the neck.**
   - The ODE reaches a ≈ 0.20 at α_c ≥ 0.1 with s = 0 at the neck. No neck undersaturation is needed to explain the slope.
3. **The ambient was undersaturated. The open question is how much of that reaches the neck.**
   - Both the ODE and the Fig. 2 run say the chamber was not saturated. The
     measured shrinkage implies a far-field undersaturation s_∞ ≈ 3.6e-3
     (step 4); the run's walls sit at 2.85e-3.
   - They differ only on the **local** undersaturation at the neck:
     - In the run, about 6 % of the wall value reaches the neck (s = 1.6e-4).
       That pulls a from 0.20 down to 0.15.
     - The ODE's fit to the data gains nothing from any local undersaturation
       (s → 0).
   - Likely readings:
     - The neck sits in a crevice shielded by the two grain surfaces, which
       hold the local vapour near their own equilibrium. Our domain may
       shield it less than the experiment did.
     - The ODE's single diffusion length c·ρ cannot represent that shielding.
     - The nine data points identify s only weakly.
   - We report both and explain the difference, without forcing them to agree.
     The ODE fixes the rate and exponent bounds; the simulation carries the
     geometry.
4. **Predictions for the planned saturated runs (equal grains, R = 84.4 µm, r₀ = 14 µm):**
   - α_c = 0.1: a ≈ **0.20**, reaching x/R = 0.35 in ~1.9 h.
   - α_c = 1e-3: a ≈ **0.29–0.30**, inside Demmenie's 0.26–0.33. It needs ~96 h simulated to reach x/R = 0.35, and 7.5 / 25 / 54 h for x/R = 0.20 / 0.25 / 0.30.
   - This is why the Demmenie comparison run uses α_c = 1e-3, which has a theoretical justification: it is the attachment-limited regime in which t^(1/3) is derived.

## Caveats

- **Small-neck geometry.** It is used up to x/R ≈ 0.39. The exact fillet slows growth further at large necks.
- **Grain asymmetry.** The unequal pair is collapsed to one R_eff.
- **c is a single number.** It is calibrated at one α_c and one temperature. Step 2 tests it at α_c = 1e-3, but it is not tested at −5 °C.
- **s is identifiable only weakly** from nine points. "s → 0" means the fit gains nothing from undersaturation, not that s is measured to be zero.
- **The ODE's exponents use the repo's `d_free` protocol** over each curve's own window, in the same way as the tables in `studies/historical_sintering/README.md`.
