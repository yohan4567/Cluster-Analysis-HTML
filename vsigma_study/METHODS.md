# V/σ Analysis of 24 Globular Clusters — Refined Methodology

*Prepared overnight from the team's plan. Every step of the original plan was checked;
two steps contained real methodological errors that are corrected below (§2).*

---

## 1. The original plan, checked

> **Team plan:** make a 3D plot of each GC → analyze → calculate V/σ of the stars →
> compare V/σ across clusters against parameters (metallicity, radius, …).
> Colleague's position: "the large parallax error is OK because we have no better
> method to get the CORRECT 3D plot → calculate the V/σ."

| Step | Verdict | Why |
|---|---|---|
| Make a 3D plot per GC | ✅ Fine **as visualization** | But it is *not an input* to V/σ (see §2.1) |
| 3D plot → calculate V/σ | ❌ **Logically backwards** | V/σ never uses 3D positions. It comes from velocities. |
| "Large error OK, no better method" | ❌ **False premise** | A better method exists and is standard in the literature: V/σ directly from proper motions, with error deconvolution. No 3D reconstruction is needed at any step. |
| V/σ "of every star" | ❌ **Not defined** | A single star has a velocity; V (rotation) and σ (dispersion) are *population* statistics. V/σ is per cluster (or per radial bin), never per star. |
| Compare V/σ vs [Fe/H], radius, … | ✅ Sound | This is a real, publishable-style question (cf. Bianchini et al. 2018). |

**Bottom line:** the pipeline below computes V/σ *directly from Gaia proper motions*
— the same approach as Bianchini et al. (2018) and Sollima et al. (2019), the exact
papers cited in the team's own target list. The 3D plots stay as visualization only.

## 2. Two errors that had to be fixed

### 2.1 Why "correct 3D plot → V/σ" cannot work — and doesn't need to

Per-star distances from parallax at GC distances are noise-dominated (for ω Cen:
median parallax 0.16 mas vs median error 0.30 mas; 30% of members have *negative*
parallax). No filtering can fix this — the information simply is not in the data.

But V/σ is a **ratio of velocities**. Converting proper motion (mas/yr) to km/s
multiplies by 4.74·d, and that factor **cancels between V and σ**. So V/σ is
distance-independent and 3D-position-independent. The "correct 3D plot" the
colleague wanted is unnecessary — which is fortunate, because it is also impossible.

### 2.2 The sinusoid bug in the colleague's rotation fit

The colleague's script fits `v_tan(θ) = a·sin θ + b·cos θ` and calls the amplitude
"rotation". This is wrong, and provably so:

- Consider a cluster with **zero rotation** but a slightly mis-measured systemic
  proper motion (offset vector **δ**). Every star's residual PM contains the constant
  vector δ. Projecting a constant vector onto the tangential direction t̂(θ) gives
  exactly `v_tan(θ) = |δ|·sin(θ − φ)` — a perfect sinusoid.
- So the sinusoidal fit measures the **error in the systemic PM**, not rotation.
  (Also, the script *hardcodes* the axis PA to 100° right after fitting it, so its
  plotted "fitted axis" was never real.)

What actual rotation looks like in proper motions: rotation about the line-of-sight
component of the spin axis produces a mean tangential streaming that is
**independent of θ**: ⟨v_tan⟩(R) ≠ 0, the same value all around the ring.
Rotation about axes lying in the sky plane produces **no mean PM signal at all**
(it lives in the radial velocities). Therefore:

> **Correct PM rotation measure:** the mean v_tan per radial bin — a *constant*
> offset around each annulus, not a sinusoid. (Bianchini et al. 2018 method.)

We keep the sinusoid fit only as a **diagnostic**: after subtracting our systemic PM
its amplitude should be ≈ 0; if it isn't, the mean PM is off.

## 3. Refined pipeline (implemented in `vsigma_pipeline.py`)

For each cluster:

1. **Select members** — membership_prob ≥ 0.90 from the team's filtered catalogs
   (fallback to ≥ 0.5 if < 100 stars survive, e.g. sparse Pal 5).
2. **Project coordinates** — x = ΔRA·cos δ₀ (East+), y = ΔDec (North+), in arcmin;
   convert to parsecs with the *literature* distance from the target list.
3. **Systemic PM** — 3σ-clipped mean of member PMs; residuals dμ per star.
4. **Velocities** — v = 4.74·d·dμ (km/s); decompose into v_rad (outward+) and
   v_tan (counterclockwise+ in (East,North) axes as drawn by the pipeline).
5. **Rotation profile** — ⟨v_tan⟩ in equal-count radial bins ± standard error.
   V_peak = largest |⟨v_tan⟩| among bins with ≥ 200 stars.
6. **Dispersion profile** — per bin, σ²_true = σ²_observed − ⟨(4.74·d·ε_μ)²⟩
   (measurement-error deconvolution — skipping this *inflates σ and crushes V/σ*;
   this is why a naive calculation gave 0.02 for ω Cen instead of ~0.3).
   σ₀ = deconvolved σ inside the half-light radius.
7. **V/σ** — primary statistic V_peak/σ₀; also a global variant. Uncertainties from
   200 bootstrap resamples of the member list.
8. **Diagnostics** — sinusoid amplitude (§2.2, should be ~0; > 1 km/s empirically
   marks contaminated samples); ⟨v_rad⟩(R) profile (contains the *perspective*
   signal, −V_LOS·θ in km/s with θ in radians — expected nonzero for nearby
   clusters with big radial velocities; not an error).
   **Rotation detection** — because V_peak is a maximum over noisy bins it is
   biased high for non-rotators; a χ² test of all bin means against zero gives
   a per-cluster detection p-value, and only p < 0.01 counts as "rotation
   detected". **Quality flag** — each cluster is graded good / error-dominated /
   contaminated / unmeasurable using: deconvolved-variance signal-to-noise,
   whether the membership fallback cut was needed, an outward-*rising* dispersion
   profile (field-star signature), and the dipole diagnostic.
9. **Outputs** — per-cluster `rotation.json` + `rotation.png` (4 panels: rotation
   curve, σ profile, V/σ(R), 2D rotation map) into the app's data folder and this
   folder; cross-cluster `vsigma_summary.csv/.json`.

## 4. Cross-cluster correlations

Spearman rank correlation (robust, no linearity assumed) of V/σ against:
[Fe/H], concentration c, half-light radius (pc), σ₀, N_members, distance.
Small-N caveat: with ~23 clusters, |ρ| ≳ 0.42 is needed for p < 0.05 — reported
alongside every coefficient.

**Validation gate:** the target list's literature "Rotation status" column
(Bianchini/Sollima) is held out as ground truth — the 7 confirmed rotators should
rank near the top of our V/σ table. If they don't, the pipeline is wrong, not the
literature.

## 5. 多旋轉軸 (nested / multiple rotation axes) — what is actually testable

The moon–earth–sun analogy doesn't transfer literally: cluster stars all orbit in the
*collective* gravitational potential — there are no local sub-centers with their own
satellites (two-body attachments are destroyed by relaxation in a few crossing
times). **But** the underlying intuition — "the rotation is not one single rigid
axis" — is a real, actively-studied phenomenon:

- **Radius-dependent spin (testable with our PMs):** compute ⟨v_tan⟩(R). If its
  sign flips or its amplitude is strongly non-monotonic, the line-of-sight spin
  component *changes with radius* — e.g. a counter-rotating core.
- **Axis tilt vs radius (needs radial velocities):** with Gaia RVs for bright members
  we fit the full 3D spin vector Ω per radial shell: Ω_z from the PM rotation curve,
  (Ω_x, Ω_y) from the RV gradient across the face of the cluster — after removing
  *perspective rotation* (apparent RV gradient = v_transverse·θ, up to ±1.6 km/s for
  ω Cen — comparable to real rotation, must be subtracted; we use the literature
  systemic PM and distance for this). The angle between shell axes measures the
  "nested axis" tilt directly. ω Cen is the known best case (literature reports
  kinematic twisting of its rotation axis with radius).
- Pipeline stage: queries Gaia for radial_velocity in each cluster field (bright
  stars only have RVs), crossmatches to our members, runs the shell fits where
  ≥ 50 RV members exist. Full results in `multi_axis_report.md`.

## 6. Honest limitations

- PMs see only the line-of-sight component of the spin vector; a cluster rotating
  exactly in the plane of the sky shows V/σ ≈ 0 here (inclination bias — affects
  every PM-only survey equally; another reason the correlation analysis uses ranks).
- Gaia PM error correlations (pmra_pmdec_corr) were not in the original query;
  deconvolution treats errors as independent per axis.
- Crowding: the innermost arcminute of core-collapsed clusters is incomplete
  (Gaia confusion limit) — σ₀ there reflects the resolvable stars.
- Structural parameters are King-model fits; least reliable for core-collapsed
  clusters (flagged in the target list itself).
- The dispersion estimator measures scatter about the bin mean of each PM
  component, so the rotation field itself inflates σ slightly (≈ (V/σ)²/2
  relative, i.e. ~6% for V/σ = 0.5); the true "ordered vs random" ratio is
  marginally higher than reported. Affects all clusters in the same direction.
- V/σ values for clusters flagged *contaminated* or *error-dominated* should be
  quoted only with that caveat; *unmeasurable* means Gaia PM errors exceed the
  intrinsic dispersion entirely (NGC 2419, Pal 5).
