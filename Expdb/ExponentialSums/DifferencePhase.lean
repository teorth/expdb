module

public import Expdb.ExponentialSums.PhaseFunctions

import Expdb.Mathlib.IteratedDeriv
import Mathlib.Analysis.Calculus.Deriv.MeanValue

/-!
# Normalized difference phases

The A-process replaces a phase by a short difference. An affine contraction keeps
both evaluations inside `[1, 2]`, where the model-phase hypotheses are available.
-/

@[expose] public section

open Filter Topology
open scoped ContDiff

noncomputable section

namespace Expdb

/-- The difference phase, contracted to stay inside the phase interval. -/
def differencePhase (F : ℝ → ℝ) (σ r : ℝ) (u : ℝ) : ℝ :=
  -(σ * r)⁻¹ * (F ((1 - r) * u + 2 * r) - F ((1 - r) * u + r))

private theorem differencePhase_mapsTo {r : ℝ} (hr : 0 ≤ r) (hr1 : r ≤ 1) :
    Set.MapsTo (fun u ↦ (1 - r) * u + r) phaseInterval phaseInterval ∧
      Set.MapsTo (fun u ↦ (1 - r) * u + 2 * r) phaseInterval phaseInterval := by
  constructor <;> intro u hu <;> constructor <;>
    dsimp [phaseInterval] at hu ⊢ <;> nlinarith [hu.1, hu.2]

/-- Contracting the two translates makes a difference phase smooth on the full interval. -/
theorem contDiffOn_differencePhase {F : ℝ → ℝ} {σ r : ℝ}
    (hF : ContDiffOn ℝ ∞ F phaseInterval) (hr : 0 ≤ r) (hr1 : r ≤ 1) :
    ContDiffOn ℝ ∞ (differencePhase F σ r) phaseInterval := by
  obtain ⟨hleft, hright⟩ := differencePhase_mapsTo hr hr1
  unfold differencePhase
  exact contDiffOn_const.mul ((hF.comp (by fun_prop) hright).sub
    (hF.comp (by fun_prop) hleft))

private theorem iteratedDerivWithin_differencePhase {F : ℝ → ℝ} {σ r u : ℝ}
    (hF : ContDiffOn ℝ ∞ F phaseInterval) (hr : 0 ≤ r) (hr1 : r ≤ 1)
    (hu : u ∈ phaseInterval) (p : ℕ) :
    iteratedDerivWithin p (differencePhase F σ r) phaseInterval u =
      -(σ * r)⁻¹ * (1 - r) ^ p *
        (iteratedDerivWithin p F phaseInterval ((1 - r) * u + 2 * r) -
          iteratedDerivWithin p F phaseInterval ((1 - r) * u + r)) := by
  obtain ⟨hleft, hright⟩ := differencePhase_mapsTo hr hr1
  have hsmooth : ContDiffOn ℝ p F phaseInterval :=
    hF.of_le (ENat.natCast_le_of_coe_top_le_withTop le_rfl p)
  have hsub := iteratedDerivWithin_sub (f := fun x ↦ F ((1 - r) * x + 2 * r))
    (g := fun x ↦ F ((1 - r) * x + r)) hu uniqueDiffOn_phaseInterval
    ((hsmooth.comp (by fun_prop) hright) u hu)
    ((hsmooth.comp (by fun_prop) hleft) u hu)
  change iteratedDerivWithin p (fun x ↦ F ((1 - r) * x + 2 * r) -
    F ((1 - r) * x + r)) phaseInterval u = _ at hsub
  unfold differencePhase
  rw [iteratedDerivWithin_const_mul_field, hsub,
    iteratedDerivWithin_comp_affine_of_mapsTo hsmooth uniqueDiffOn_phaseInterval
      uniqueDiffOn_phaseInterval hu hright,
    iteratedDerivWithin_comp_affine_of_mapsTo hsmooth uniqueDiffOn_phaseInterval
      uniqueDiffOn_phaseInterval hu hleft]
  ring

/-- A short derivative difference is a higher derivative at an intermediate point.
Using the mean value theorem avoids differentiating a parameter-dependent integral. -/
private theorem exists_iteratedDerivWithin_difference {F : ℝ → ℝ} {v r : ℝ}
    (hF : ContDiffOn ℝ ∞ F phaseInterval) (hr : 0 < r)
    (hv : v ∈ phaseInterval) (hvr : v + r ∈ phaseInterval) (p : ℕ) :
    ∃ w ∈ Set.Ioo v (v + r),
      iteratedDerivWithin p F phaseInterval (v + r) -
          iteratedDerivWithin p F phaseInterval v =
        r * iteratedDerivWithin (p + 1) F phaseInterval w := by
  have hsub : Set.Icc v (v + r) ⊆ phaseInterval := by
    intro x hx
    exact ⟨hv.1.trans hx.1, hx.2.trans hvr.2⟩
  have hcont := hF.continuousOn_iteratedDerivWithin (m := p)
    (ENat.natCast_le_of_coe_top_le_withTop le_rfl p)
    uniqueDiffOn_phaseInterval
  have hdiff := (hF.of_le
    (ENat.natCast_le_of_coe_top_le_withTop le_rfl (p + 1))).differentiableOn_iteratedDerivWithin
    (m := p) (Nat.cast_lt.mpr p.lt_succ_self)
    uniqueDiffOn_phaseInterval
  obtain ⟨w, hw, heq⟩ := exists_hasDerivAt_eq_slope
    (iteratedDerivWithin p F phaseInterval)
    (iteratedDerivWithin (p + 1) F phaseInterval) (by linarith : v < v + r)
    (hcont.mono hsub) (fun x hx ↦ by
      have hxphase : x ∈ phaseInterval := hsub (Set.Ioo_subset_Icc_self hx)
      have hxinterior : phaseInterval ∈ 𝓝 x := Icc_mem_nhds
        (hv.1.trans_lt hx.1) (hx.2.trans_le hvr.2)
      simpa only [iteratedDerivWithin_succ] using
        (hdiff x hxphase).hasDerivWithinAt.hasDerivAt hxinterior)
  refine ⟨w, hw, ?_⟩
  have heq' : iteratedDerivWithin (p + 1) F phaseInterval w =
      (iteratedDerivWithin p F phaseInterval (v + r) -
        iteratedDerivWithin p F phaseInterval v) / r := by
    simpa using heq
  exact (eq_div_iff hr.ne').mp heq' |>.symm.trans (mul_comm _ _)

private theorem iteratedDerivWithin_modelPhase_succ (σ : ℝ) (p : ℕ)
    {u : ℝ} (hu : u ∈ phaseInterval) :
    iteratedDerivWithin (p + 1) (modelPhase σ) phaseInterval u =
      -σ * iteratedDerivWithin p (modelPhase (σ + 1)) phaseInterval u := by
  rw [iteratedDerivWithin_modelPhase σ (p + 1) hu,
    iteratedDerivWithin_modelPhase (σ + 1) p hu,
    descPochhammer_succ_left]
  simp only [Polynomial.eval_mul, Polynomial.eval_X, Polynomial.eval_comp,
    Polynomial.eval_sub, Polynomial.eval_one]
  rw [show -σ - 1 = -(σ + 1) by ring,
    show -σ - ↑(p + 1) = -(σ + 1) - ↑p by push_cast; ring]
  ring

/-- A difference phase has model exponent `σ + 1`. The hypotheses hold for any
positive shifts tending to zero, so the result applies to shifts chosen by pigeonholing. -/
theorem IsModelPhaseFunctionWith.differencePhase
    {F : VariableFunction (VariableObject.fixed ℝ) ℝ} {σ : ℝ}
    (hF : IsModelPhaseFunctionWith F σ) (hσ : 0 < σ)
    {r : VariableObject ℝ} (hr : ∀ i, 0 < r i) (hr1 : ∀ i, r i < 1)
    (hrzero : Tendsto r atTop (𝓝 0)) :
    IsModelPhaseFunctionWith (fun i ↦ differencePhase (F i) σ (r i)) (σ + 1) := by
  refine ⟨fun i ↦ contDiffOn_differencePhase (hF.isPhaseFunction i)
    (hr i).le (hr1 i).le, ?_⟩
  intro p u
  let v : VariableObject ℝ := fun i ↦ (1 - r i) * (u i : ℝ) + r i
  have hv : ∀ i, v i ∈ phaseInterval := fun i ↦
    (differencePhase_mapsTo (hr i).le (hr1 i).le).1 (u i).property
  have hvr : ∀ i, v i + r i ∈ phaseInterval := by
    intro i
    convert (differencePhase_mapsTo (hr i).le (hr1 i).le).2 (u i).property using 1
    dsimp [v]
    ring
  have hmean := fun i ↦ exists_iteratedDerivWithin_difference (hF.isPhaseFunction i)
    (hr i) (hv i) (hvr i) (p + 1)
  choose w hw hderiv using hmean
  let wphase : ∀ i, phaseInterval := fun i ↦ ⟨w i,
    ⟨(hv i).1.trans (hw i).1.le, (hw i).2.le.trans (hvr i).2⟩⟩
  have hwclose : Tendsto (fun i ↦ w i - (u i : ℝ)) atTop (𝓝 0) := by
    apply squeeze_zero_norm' _ hrzero
    filter_upwards [] with i
    rw [Real.norm_eq_abs, abs_le]
    have hu := (u i).property
    have hi := hw i
    dsimp [v, phaseInterval] at hi hu
    constructor <;> nlinarith [hr i, hu.1, hu.2, hi.1, hi.2]
  let ref : ℝ → ℝ := iteratedDerivWithin p (modelPhase (σ + 1)) phaseInterval
  have href := continuousOn_iteratedDerivWithin_modelPhase (σ + 1) p
  have hrefclose : Tendsto (fun i ↦ ref (w i) - ref (u i)) atTop (𝓝 0) :=
    (isEqUpToInfinitesimal_iff_tendsto_sub _ _).1
      (IsEqUpToInfinitesimal.comp_of_isCompact isCompact_Icc href
        ((isEqUpToInfinitesimal_iff_tendsto_sub _ _).2 hwclose)
        (fun i ↦ (wphase i).property) (fun i ↦ (u i).property))
  have hpow : Tendsto (fun i ↦ (1 - r i) ^ (p + 1)) atTop (𝓝 1) := by
    simpa using ((tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ (1 : ℝ))
      atTop (𝓝 1)).sub hrzero).pow (p + 1)
  obtain ⟨D, hD⟩ := isCompact_Icc.exists_bound_of_continuousOn href
  have hrefBounded : IsBoundedUnder (· ≤ ·) atTop ((‖·‖) ∘ fun i ↦ ref (w i)) :=
    Filter.isBoundedUnder_of_eventually_le <| Filter.Eventually.of_forall
      fun i ↦ hD (w i) (wphase i).property
  have hpowzero : Tendsto (fun i ↦ (1 - r i) ^ (p + 1) - 1) atTop (𝓝 0) := by
    simpa using hpow.sub (tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ (1 : ℝ)) atTop (𝓝 1))
  have hscale : Tendsto (fun i ↦ ((1 - r i) ^ (p + 1) - 1) * ref (w i))
      atTop (𝓝 0) :=
    hpowzero.zero_mul_isBoundedUnder_le hrefBounded
  have herr := hF.error_isChoicewiseInfinitesimal (p + 1) wphase
  rw [VariableObject.IsInfinitesimal] at herr ⊢
  have herrzero : Tendsto (fun i ↦ modelPhaseErrorAt (F i) σ (p + 1) (w i))
      atTop (𝓝 0) := by
    exact tendsto_zero_iff_norm_tendsto_zero.mpr herr
  have herror : Tendsto (fun i ↦
      ((1 - r i) ^ (p + 1) - 1) * ref (w i) + (ref (w i) - ref (u i)) -
        (1 - r i) ^ (p + 1) / σ * modelPhaseErrorAt (F i) σ (p + 1) (w i))
      atTop (𝓝 0) := by
    convert (hscale.add hrefclose).sub ((hpow.div_const σ).mul herrzero) using 1
    simp
  convert herror.norm using 1 <;> try simp only [norm_zero]
  ext i
  change ‖modelPhaseErrorAt (Expdb.differencePhase (F i) σ (r i)) (σ + 1) p (u i)‖ = _
  congr 1
  rw [modelPhaseErrorAt, iteratedDerivWithin_differencePhase (hF.isPhaseFunction i)
    (hr i).le (hr1 i).le (u i).property]
  have hvadd : (1 - r i) * (u i : ℝ) + 2 * r i = v i + r i := by dsimp [v]; ring
  rw [hvadd, hderiv i]
  have hmodel := iteratedDerivWithin_modelPhase_succ σ p (wphase i).property
  dsimp [ref, modelPhaseErrorAt]
  rw [show p + 1 + 1 = p + 2 by omega] at *
  rw [hmodel]
  field_simp [hσ.ne', (hr i).ne']
  ring

end Expdb
