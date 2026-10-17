module

public import Mathlib.Analysis.Calculus.IteratedDeriv.Lemmas
public import Mathlib.Analysis.Fourier.AddCircle
import Mathlib.Analysis.Fourier.FourierTransformDeriv
import Mathlib.Analysis.Calculus.Deriv.Inv
import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Analysis.SpecialFunctions.ExpDeriv

public import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.GCongr
import Mathlib.Tactic.FunProp

/-!
# Two integrations by parts for nonstationary phase

Adapted from `Proof/MRTNonstationaryPhaseTwoIBP.lean` in
[`jconorgrogan/prime-minor-arcs-2-15`](https://github.com/jconorgrogan/prime-minor-arcs-2-15)
at commit `f369f267b4dcfebf010c8e7baa2c9602e2960eba`, licensed under Apache-2.0
(see the root LICENSE).

Adaptations for Expdb: retained the two integrations-by-parts estimates with an explicit
inverse phase derivative `q`, and removed the project-specific import. Namespaces
and proofs were updated for the current Mathlib API.
-/

@[expose] public section

namespace Expdb.NonstationaryPhase

open MeasureTheory Set Filter Topology Complex
open scoped ContDiff FourierTransform
noncomputable section

/-- Exact two-fold integration-by-parts identity.  If `E' = Ed` and
`q * Ed = E`, then the oscillatory amplitude can be divided by the phase
derivative twice.  The endpoint hypotheses are the literal compact-support
boundary conditions used for compactly supported cutoffs. -/
theorem two_integration_by_parts_identity
    {a b : ℝ} {E Ed q qd qdd amplitude amplitude' amplitude'' : ℝ → ℂ}
    (hab : a ≤ b)
    (hE : ∀ w ∈ Set.Icc a b, HasDerivAt E (Ed w) w)
    (hq : ∀ w ∈ Set.Icc a b, HasDerivAt q (qd w) w)
    (hqd : ∀ w ∈ Set.Icc a b, HasDerivAt qd (qdd w) w)
    (ha : ∀ w ∈ Set.Icc a b, HasDerivAt amplitude (amplitude' w) w)
    (ha' : ∀ w ∈ Set.Icc a b, HasDerivAt amplitude' (amplitude'' w) w)
    (hEdCont : ContinuousOn Ed (Set.Icc a b))
    (hqddCont : ContinuousOn qdd (Set.Icc a b))
    (ha''Cont : ContinuousOn amplitude'' (Set.Icc a b))
    (hinverse : ∀ w ∈ Set.Icc a b, q w * Ed w = E w)
    (ha_left : amplitude a = 0) (ha_right : amplitude b = 0)
    (ha'_left : amplitude' a = 0) (ha'_right : amplitude' b = 0) :
    (∫ w : ℝ in a..b, E w * amplitude w) =
      ∫ w : ℝ in a..b,
        (((amplitude'' w * q w + amplitude' w * qd w) +
              (amplitude' w * qd w + amplitude w * qdd w)) * q w +
          (amplitude' w * q w + amplitude w * qd w) * qd w) * E w := by
  let u : ℝ → ℂ := fun w ↦ amplitude w * q w
  let ud : ℝ → ℂ := fun w ↦ amplitude' w * q w + amplitude w * qd w
  let udd : ℝ → ℂ := fun w ↦
    (amplitude'' w * q w + amplitude' w * qd w) +
      (amplitude' w * qd w + amplitude w * qdd w)
  let v : ℝ → ℂ := fun w ↦ ud w * q w
  let vd : ℝ → ℂ := fun w ↦ udd w * q w + ud w * qd w
  have hu : ∀ w ∈ Set.Icc a b, HasDerivAt u (ud w) w := by
    intro w hw
    exact (ha w hw).mul (hq w hw)
  have hud : ∀ w ∈ Set.Icc a b, HasDerivAt ud (udd w) w := by
    intro w hw
    exact (ha' w hw).mul (hq w hw) |>.add ((ha w hw).mul (hqd w hw))
  have hv : ∀ w ∈ Set.Icc a b, HasDerivAt v (vd w) w := by
    intro w hw
    exact (hud w hw).mul (hq w hw)
  have hECont : ContinuousOn E (Set.Icc a b) :=
    fun w hw ↦ (hE w hw).continuousAt.continuousWithinAt
  have hqCont : ContinuousOn q (Set.Icc a b) :=
    fun w hw ↦ (hq w hw).continuousAt.continuousWithinAt
  have hqdCont : ContinuousOn qd (Set.Icc a b) :=
    fun w hw ↦ (hqd w hw).continuousAt.continuousWithinAt
  have haCont : ContinuousOn amplitude (Set.Icc a b) :=
    fun w hw ↦ (ha w hw).continuousAt.continuousWithinAt
  have ha'Cont : ContinuousOn amplitude' (Set.Icc a b) :=
    fun w hw ↦ (ha' w hw).continuousAt.continuousWithinAt
  have huCont : ContinuousOn u (Set.Icc a b) := haCont.mul hqCont
  have hudCont : ContinuousOn ud (Set.Icc a b) :=
    (ha'Cont.mul hqCont).add (haCont.mul hqdCont)
  have huddCont : ContinuousOn udd (Set.Icc a b) :=
    ((ha''Cont.mul hqCont).add (ha'Cont.mul hqdCont)).add
      ((ha'Cont.mul hqdCont).add (haCont.mul hqddCont))
  have hvCont : ContinuousOn v (Set.Icc a b) := hudCont.mul hqCont
  have hvdCont : ContinuousOn vd (Set.Icc a b) :=
    (huddCont.mul hqCont).add (hudCont.mul hqdCont)
  have hiEd : IntervalIntegrable Ed volume a b := by
    apply ContinuousOn.intervalIntegrable
    simpa [Set.uIcc_of_le hab] using hEdCont
  have hiud : IntervalIntegrable ud volume a b := by
    apply ContinuousOn.intervalIntegrable
    simpa [Set.uIcc_of_le hab] using hudCont
  have hivd : IntervalIntegrable vd volume a b := by
    apply ContinuousOn.intervalIntegrable
    simpa [Set.uIcc_of_le hab] using hvdCont
  have hibp1 : (∫ w : ℝ in a..b, u w * Ed w) =
      u b * E b - u a * E a - ∫ w : ℝ in a..b, ud w * E w := by
    exact intervalIntegral.integral_mul_deriv_eq_deriv_mul
      (by simpa [Set.uIcc_of_le hab] using hu)
      (by simpa [Set.uIcc_of_le hab] using hE) hiud hiEd
  have hibp2 : (∫ w : ℝ in a..b, v w * Ed w) =
      v b * E b - v a * E a - ∫ w : ℝ in a..b, vd w * E w := by
    exact intervalIntegral.integral_mul_deriv_eq_deriv_mul
      (by simpa [Set.uIcc_of_le hab] using hv)
      (by simpa [Set.uIcc_of_le hab] using hE) hivd hiEd
  have hfirstIntegrand : ∀ w ∈ Set.uIcc a b,
      E w * amplitude w = u w * Ed w := by
    intro w hw
    have hw' : w ∈ Set.Icc a b := by simpa [Set.uIcc_of_le hab] using hw
    unfold u
    rw [mul_assoc, hinverse w hw']
    ac_rfl
  have hsecondIntegrand : ∀ w ∈ Set.uIcc a b,
      ud w * E w = v w * Ed w := by
    intro w hw
    have hw' : w ∈ Set.Icc a b := by simpa [Set.uIcc_of_le hab] using hw
    unfold v
    rw [mul_assoc, hinverse w hw']
  have huEnds : u a = 0 ∧ u b = 0 := by simp [u, ha_left, ha_right]
  have hudEnds : ud a = 0 ∧ ud b = 0 := by
    simp [ud, ha_left, ha_right, ha'_left, ha'_right]
  have hvEnds : v a = 0 ∧ v b = 0 := by simp [v, hudEnds]
  rw [intervalIntegral.integral_congr hfirstIntegrand, hibp1,
    huEnds.1, huEnds.2]
  simp only [zero_mul, sub_zero, zero_sub]
  rw [intervalIntegral.integral_congr hsecondIntegrand, hibp2,
    hvEnds.1, hvEnds.2]
  simp only [zero_mul, sub_zero, zero_sub, neg_neg]
  rfl

/-- Quantitative consequence of the exact two-IBP identity.  The three
amplitude `L¹` budgets and three inverse-phase derivative sup bounds are kept
separate so the two MRT applications can retain their literal scales. -/
theorem norm_integral_le_two_ibp_budgets
    {a b A0 A1 A2 Q0 Q1 Q2 : ℝ}
    {E Ed q qd qdd amplitude amplitude' amplitude'' : ℝ → ℂ}
    (hab : a ≤ b)
    (hE : ∀ w ∈ Set.Icc a b, HasDerivAt E (Ed w) w)
    (hq : ∀ w ∈ Set.Icc a b, HasDerivAt q (qd w) w)
    (hqd : ∀ w ∈ Set.Icc a b, HasDerivAt qd (qdd w) w)
    (ha : ∀ w ∈ Set.Icc a b, HasDerivAt amplitude (amplitude' w) w)
    (ha' : ∀ w ∈ Set.Icc a b, HasDerivAt amplitude' (amplitude'' w) w)
    (hEdCont : ContinuousOn Ed (Set.Icc a b))
    (hqddCont : ContinuousOn qdd (Set.Icc a b))
    (ha''Cont : ContinuousOn amplitude'' (Set.Icc a b))
    (hinverse : ∀ w ∈ Set.Icc a b, q w * Ed w = E w)
    (ha_left : amplitude a = 0) (ha_right : amplitude b = 0)
    (ha'_left : amplitude' a = 0) (ha'_right : amplitude' b = 0)
    (hA0 : (∫ w : ℝ in a..b, ‖amplitude w‖) ≤ A0)
    (hA1 : (∫ w : ℝ in a..b, ‖amplitude' w‖) ≤ A1)
    (hA2 : (∫ w : ℝ in a..b, ‖amplitude'' w‖) ≤ A2)
    (hQ0 : ∀ w ∈ Set.Icc a b, ‖q w‖ ≤ Q0)
    (hQ1 : ∀ w ∈ Set.Icc a b, ‖qd w‖ ≤ Q1)
    (hQ2 : ∀ w ∈ Set.Icc a b, ‖qdd w‖ ≤ Q2)
    (hEUnit : ∀ w ∈ Set.Icc a b, ‖E w‖ = 1)
    (hQ0nonneg : 0 ≤ Q0) (hQ1nonneg : 0 ≤ Q1) (hQ2nonneg : 0 ≤ Q2) :
    ‖∫ w : ℝ in a..b, E w * amplitude w‖ ≤
      Q0 ^ 2 * A2 + 3 * Q0 * Q1 * A1 +
        (Q0 * Q2 + Q1 ^ 2) * A0 := by
  have hid := two_integration_by_parts_identity hab hE hq hqd ha ha'
    hEdCont hqddCont ha''Cont hinverse ha_left ha_right ha'_left ha'_right
  rw [hid]
  let vd : ℝ → ℂ := fun w ↦
    (((amplitude'' w * q w + amplitude' w * qd w) +
          (amplitude' w * qd w + amplitude w * qdd w)) * q w +
      (amplitude' w * q w + amplitude w * qd w) * qd w) * E w
  have hqCont : ContinuousOn q (Set.Icc a b) :=
    fun w hw ↦ (hq w hw).continuousAt.continuousWithinAt
  have hqdCont : ContinuousOn qd (Set.Icc a b) :=
    fun w hw ↦ (hqd w hw).continuousAt.continuousWithinAt
  have haCont : ContinuousOn amplitude (Set.Icc a b) :=
    fun w hw ↦ (ha w hw).continuousAt.continuousWithinAt
  have ha'Cont : ContinuousOn amplitude' (Set.Icc a b) :=
    fun w hw ↦ (ha' w hw).continuousAt.continuousWithinAt
  have hECont : ContinuousOn E (Set.Icc a b) :=
    fun w hw ↦ (hE w hw).continuousAt.continuousWithinAt
  have hvdCont : ContinuousOn vd (Set.Icc a b) := by
    unfold vd
    fun_prop
  have hivd : IntervalIntegrable vd volume a b := by
    apply ContinuousOn.intervalIntegrable
    simpa [Set.uIcc_of_le hab] using hvdCont
  have hmajorCont : ContinuousOn
      (fun w ↦ Q0 ^ 2 * ‖amplitude'' w‖ +
        3 * Q0 * Q1 * ‖amplitude' w‖ +
          (Q0 * Q2 + Q1 ^ 2) * ‖amplitude w‖) (Set.Icc a b) := by
    fun_prop
  have himajor : IntervalIntegrable
      (fun w ↦ Q0 ^ 2 * ‖amplitude'' w‖ +
        3 * Q0 * Q1 * ‖amplitude' w‖ +
          (Q0 * Q2 + Q1 ^ 2) * ‖amplitude w‖) volume a b := by
    apply ContinuousOn.intervalIntegrable
    simpa [Set.uIcc_of_le hab] using hmajorCont
  have hpoint : ∀ w ∈ Set.Icc a b,
      ‖vd w‖ ≤ Q0 ^ 2 * ‖amplitude'' w‖ +
        3 * Q0 * Q1 * ‖amplitude' w‖ +
          (Q0 * Q2 + Q1 ^ 2) * ‖amplitude w‖ := by
    intro w hw
    have h0 := hQ0 w hw
    have h1 := hQ1 w hw
    have h2 := hQ2 w hw
    have ha2q : ‖amplitude'' w * q w‖ ≤ ‖amplitude'' w‖ * Q0 := by
      rw [norm_mul]
      exact mul_le_mul_of_nonneg_left h0 (norm_nonneg _)
    have ha1qd : ‖amplitude' w * qd w‖ ≤ ‖amplitude' w‖ * Q1 := by
      rw [norm_mul]
      exact mul_le_mul_of_nonneg_left h1 (norm_nonneg _)
    have ha0qdd : ‖amplitude w * qdd w‖ ≤ ‖amplitude w‖ * Q2 := by
      rw [norm_mul]
      exact mul_le_mul_of_nonneg_left h2 (norm_nonneg _)
    have ha1q : ‖amplitude' w * q w‖ ≤ ‖amplitude' w‖ * Q0 := by
      rw [norm_mul]
      exact mul_le_mul_of_nonneg_left h0 (norm_nonneg _)
    have ha0qd : ‖amplitude w * qd w‖ ≤ ‖amplitude w‖ * Q1 := by
      rw [norm_mul]
      exact mul_le_mul_of_nonneg_left h1 (norm_nonneg _)
    have hblock1 :
        ‖(amplitude'' w * q w + amplitude' w * qd w) +
            (amplitude' w * qd w + amplitude w * qdd w)‖ ≤
          (‖amplitude'' w‖ * Q0 + ‖amplitude' w‖ * Q1) +
            (‖amplitude' w‖ * Q1 + ‖amplitude w‖ * Q2) := by
      exact (norm_add_le _ _).trans (add_le_add
        ((norm_add_le _ _).trans (add_le_add ha2q ha1qd))
        ((norm_add_le _ _).trans (add_le_add ha1qd ha0qdd)))
    have hblock2 :
        ‖amplitude' w * q w + amplitude w * qd w‖ ≤
          ‖amplitude' w‖ * Q0 + ‖amplitude w‖ * Q1 := by
      exact (norm_add_le _ _).trans (add_le_add ha1q ha0qd)
    have hraw :
        ‖vd w‖ ≤
          ((‖amplitude'' w‖ * Q0 + ‖amplitude' w‖ * Q1) +
            (‖amplitude' w‖ * Q1 + ‖amplitude w‖ * Q2)) * Q0 +
          (‖amplitude' w‖ * Q0 + ‖amplitude w‖ * Q1) * Q1 := by
      unfold vd
      rw [norm_mul, hEUnit w hw, mul_one]
      calc
        ‖(((amplitude'' w * q w + amplitude' w * qd w) +
              (amplitude' w * qd w + amplitude w * qdd w)) * q w +
            (amplitude' w * q w + amplitude w * qd w) * qd w)‖ ≤
            ‖((amplitude'' w * q w + amplitude' w * qd w) +
              (amplitude' w * qd w + amplitude w * qdd w)) * q w‖ +
              ‖(amplitude' w * q w + amplitude w * qd w) * qd w‖ :=
          norm_add_le _ _
        _ = ‖(amplitude'' w * q w + amplitude' w * qd w +
                (amplitude' w * qd w + amplitude w * qdd w))‖ * ‖q w‖ +
              ‖amplitude' w * q w + amplitude w * qd w‖ * ‖qd w‖ := by
          rw [norm_mul, norm_mul]
        _ ≤ ((‖amplitude'' w‖ * Q0 + ‖amplitude' w‖ * Q1) +
              (‖amplitude' w‖ * Q1 + ‖amplitude w‖ * Q2)) * Q0 +
            (‖amplitude' w‖ * Q0 + ‖amplitude w‖ * Q1) * Q1 := by
          apply add_le_add
          · exact mul_le_mul hblock1 h0 (norm_nonneg _)
              (by positivity)
          · exact mul_le_mul hblock2 h1 (norm_nonneg _)
              (by positivity)
    calc
      ‖vd w‖ ≤ _ := hraw
      _ = Q0 ^ 2 * ‖amplitude'' w‖ +
          3 * Q0 * Q1 * ‖amplitude' w‖ +
            (Q0 * Q2 + Q1 ^ 2) * ‖amplitude w‖ := by ring
  have hnorm : ‖∫ w : ℝ in a..b, vd w‖ ≤
      ∫ w : ℝ in a..b,
        (Q0 ^ 2 * ‖amplitude'' w‖ +
          3 * Q0 * Q1 * ‖amplitude' w‖ +
            (Q0 * Q2 + Q1 ^ 2) * ‖amplitude w‖) := by
    calc
      ‖∫ w : ℝ in a..b, vd w‖ ≤ ∫ w : ℝ in a..b, ‖vd w‖ :=
        intervalIntegral.norm_integral_le_integral_norm hab
      _ ≤ _ := intervalIntegral.integral_mono_on hab hivd.norm himajor hpoint
  have hA0nonneg : 0 ≤ A0 :=
    le_trans (intervalIntegral.integral_nonneg hab (fun _ _ ↦ norm_nonneg _)) hA0
  have hA1nonneg : 0 ≤ A1 :=
    le_trans (intervalIntegral.integral_nonneg hab (fun _ _ ↦ norm_nonneg _)) hA1
  have hA2nonneg : 0 ≤ A2 :=
    le_trans (intervalIntegral.integral_nonneg hab (fun _ _ ↦ norm_nonneg _)) hA2
  have hiA0 : IntervalIntegrable (fun w ↦ ‖amplitude w‖) volume a b := by
    apply ContinuousOn.intervalIntegrable
    simpa [Set.uIcc_of_le hab] using haCont.norm
  have hiA1 : IntervalIntegrable (fun w ↦ ‖amplitude' w‖) volume a b := by
    apply ContinuousOn.intervalIntegrable
    simpa [Set.uIcc_of_le hab] using ha'Cont.norm
  have hiA2 : IntervalIntegrable (fun w ↦ ‖amplitude'' w‖) volume a b := by
    apply ContinuousOn.intervalIntegrable
    simpa [Set.uIcc_of_le hab] using ha''Cont.norm
  calc
    ‖∫ w : ℝ in a..b, vd w‖ ≤ _ := hnorm
    _ = Q0 ^ 2 * (∫ w : ℝ in a..b, ‖amplitude'' w‖) +
        3 * Q0 * Q1 * (∫ w : ℝ in a..b, ‖amplitude' w‖) +
          (Q0 * Q2 + Q1 ^ 2) * (∫ w : ℝ in a..b, ‖amplitude w‖) := by
      rw [show (fun w ↦ Q0 ^ 2 * ‖amplitude'' w‖ +
          3 * Q0 * Q1 * ‖amplitude' w‖ +
            (Q0 * Q2 + Q1 ^ 2) * ‖amplitude w‖) =
          (fun w ↦ Q0 ^ 2 * ‖amplitude'' w‖ +
            (3 * Q0 * Q1 * ‖amplitude' w‖ +
              (Q0 * Q2 + Q1 ^ 2) * ‖amplitude w‖)) by
            funext w; ring,
        intervalIntegral.integral_add (hiA2.const_mul _)
          ((hiA1.const_mul _).add (hiA0.const_mul _)),
        intervalIntegral.integral_add (hiA1.const_mul _) (hiA0.const_mul _),
        intervalIntegral.integral_const_mul,
        intervalIntegral.integral_const_mul,
        intervalIntegral.integral_const_mul]
      ring
    _ ≤ Q0 ^ 2 * A2 + 3 * Q0 * Q1 * A1 +
        (Q0 * Q2 + Q1 ^ 2) * A0 := by
      gcongr


/-! ## A derivative-bound estimate for weighted oscillatory integrals

The following estimate is proved within Expdb. -/

private theorem eventually_zero_of_notMem {w : ℝ → ℝ} {a' b' x : ℝ}
    (hw0 : ∀ y, y ∉ Icc a' b' → w y = 0) (hx' : x < a' ∨ b' < x) :
    w =ᶠ[𝓝 x] fun _ ↦ 0 := by
  rcases hx' with h | h
  · filter_upwards [Iio_mem_nhds h] with y hy
    exact hw0 y fun hm ↦ by linarith [hm.1, (show y < a' from hy)]
  · filter_upwards [Ioi_mem_nhds h] with y hy
    exact hw0 y fun hm ↦ by linarith [hm.2, (show b' < y from hy)]

/-- Two integrations by parts for a weight supported strictly inside an interval on which the
phase is nonstationary. -/
theorem norm_integral_nonstationary_le {H w : ℝ → ℝ} {lam ν d C2 C3 W a a' b' b : ℝ}
    (haa' : a < a') (ha'b' : a' ≤ b') (hb'b : b' < b)
    (hH : ContDiff ℝ ∞ H) (hw : ContDiff ℝ ∞ w)
    (hw0 : ∀ x, x ∉ Icc a' b' → w x = 0)
    (hlam : 0 < lam) (hd : 0 < d)
    (hsep : ∀ x ∈ Icc a b, d ≤ |deriv H x - ν|)
    (h2 : ∀ x ∈ Icc a b, |iteratedDeriv 2 H x| ≤ C2)
    (h3 : ∀ x ∈ Icc a b, |iteratedDeriv 3 H x| ≤ C3)
    (hW : ∀ k ≤ 2, ∀ x, |iteratedDeriv k w x| ≤ W) :
    ‖∫ x : ℝ, (w x : ℂ) * 𝐞 (lam * (H x - ν * x))‖ ≤
      (b - a) * W * (1 + 3 * C2 / d + (C3 / d + 3 * C2 ^ 2 / d ^ 2)) / (lam ^ 2 * d ^ 2) := by
  have hab : a ≤ b := by linarith
  have haI : a ∈ Icc a b := ⟨le_rfl, hab⟩
  have hC2 : 0 ≤ C2 := (abs_nonneg _).trans (h2 a haI)
  have hC3 : 0 ≤ C3 := (abs_nonneg _).trans (h3 a haI)
  have hW0 : 0 ≤ W := (abs_nonneg _).trans (hW 0 (by norm_num) 0)
  set D1 := deriv H with hD1
  set D2 := deriv D1 with hD2
  set D3 := deriv D2 with hD3
  have hD1s : ContDiff ℝ ∞ D1 := (contDiff_infty_iff_deriv.1 hH).2
  have hD2s : ContDiff ℝ ∞ D2 := (contDiff_infty_iff_deriv.1 hD1s).2
  have hD3s : ContDiff ℝ ∞ D3 := (contDiff_infty_iff_deriv.1 hD2s).2
  have hHd : ∀ x, HasDerivAt H (D1 x) x := fun x ↦
    ((hH.differentiable (by simp)) x).hasDerivAt
  have hD1d : ∀ x, HasDerivAt D1 (D2 x) x := fun x ↦
    ((hD1s.differentiable (by simp)) x).hasDerivAt
  have hD2d : ∀ x, HasDerivAt D2 (D3 x) x := fun x ↦
    ((hD2s.differentiable (by simp)) x).hasDerivAt
  have hi2 : iteratedDeriv 2 H = D2 := by
    rw [iteratedDeriv_succ, iteratedDeriv_one]
  have hi3 : iteratedDeriv 3 H = D3 := by
    rw [iteratedDeriv_succ, hi2]
  rw [hi2] at h2
  rw [hi3] at h3
  set w1 := deriv w with hw1
  set w2 := deriv w1 with hw2
  have hw1s : ContDiff ℝ ∞ w1 := (contDiff_infty_iff_deriv.1 hw).2
  have hw2s : ContDiff ℝ ∞ w2 := (contDiff_infty_iff_deriv.1 hw1s).2
  have hwd : ∀ x, HasDerivAt w (w1 x) x := fun x ↦
    ((hw.differentiable (by simp)) x).hasDerivAt
  have hw1d : ∀ x, HasDerivAt w1 (w2 x) x := fun x ↦
    ((hw1s.differentiable (by simp)) x).hasDerivAt
  have hWw : ∀ x, |w x| ≤ W := fun x ↦ by simpa using hW 0 (by norm_num) x
  have hWw1 : ∀ x, |w1 x| ≤ W := fun x ↦ by simpa [iteratedDeriv_one] using hW 1 (by norm_num) x
  have hWw2 : ∀ x, |w2 x| ≤ W := fun x ↦ by
    have := hW 2 (by norm_num) x
    rwa [iteratedDeriv_succ, iteratedDeriv_one] at this
  set L : ℝ := 2 * Real.pi * lam with hL
  have hLpos : 0 < L := by positivity
  have hLlam : lam ≤ L := by
    rw [hL]; nlinarith [Real.two_le_pi]
  -- phase functions
  set E : ℝ → ℂ := fun x ↦ (𝐞 (lam * (H x - ν * x)) : ℂ) with hE
  set g : ℝ → ℂ := fun x ↦ ((L * (D1 x - ν) : ℝ) : ℂ) * I with hg
  set g1 : ℝ → ℂ := fun x ↦ ((L * D2 x : ℝ) : ℂ) * I with hg1
  set g2 : ℝ → ℂ := fun x ↦ ((L * D3 x : ℝ) : ℂ) * I with hg2
  have hnormg : ∀ x, ‖g x‖ = L * |D1 x - ν| := fun x ↦ by
    simp only [g]
    rw [norm_mul, Complex.norm_I, mul_one, Complex.norm_real, Real.norm_eq_abs, abs_mul,
      abs_of_pos hLpos]
  have hnormg1 : ∀ x, ‖g1 x‖ = L * |D2 x| := fun x ↦ by
    simp only [g1]
    rw [norm_mul, Complex.norm_I, mul_one, Complex.norm_real, Real.norm_eq_abs, abs_mul,
      abs_of_pos hLpos]
  have hnormg2 : ∀ x, ‖g2 x‖ = L * |D3 x| := fun x ↦ by
    simp only [g2]
    rw [norm_mul, Complex.norm_I, mul_one, Complex.norm_real, Real.norm_eq_abs, abs_mul,
      abs_of_pos hLpos]
  have hglow : ∀ x ∈ Icc a b, L * d ≤ ‖g x‖ := fun x hx ↦ by
    rw [hnormg]; exact mul_le_mul_of_nonneg_left (hsep x hx) hLpos.le
  have hgne : ∀ x ∈ Icc a b, g x ≠ 0 := fun x hx h ↦ by
    have := hglow x hx
    rw [h, norm_zero] at this
    nlinarith
  have hgd : ∀ x, HasDerivAt g (g1 x) x := fun x ↦ by
    have := (((hD1d x).sub_const ν).const_mul L).ofReal_comp.mul_const I
    simpa [g, g1] using this
  have hg1d : ∀ x, HasDerivAt g1 (g2 x) x := fun x ↦ by
    have := ((hD2d x).const_mul L).ofReal_comp.mul_const I
    simpa [g1, g2] using this
  set Ed : ℝ → ℂ := fun x ↦ g x * E x with hEd
  have hEder : ∀ x, HasDerivAt E (Ed x) x := fun x ↦ by
    have hφ : HasDerivAt (fun y ↦ lam * (H y - ν * y)) (lam * (D1 x - ν)) x := by
      have := ((hHd x).sub ((hasDerivAt_id x).const_mul ν)).const_mul lam
      simpa using this
    have hchar := (Real.hasDerivAt_fourierChar (lam * (H x - ν * x))).scomp x hφ
    have : HasDerivAt (fun y ↦ (𝐞 (lam * (H y - ν * y)) : ℂ))
        (((2 * Real.pi : ℝ) : ℂ) * I * (𝐞 (lam * (H x - ν * x)) : ℂ) *
          (lam * (D1 x - ν) : ℝ)) x := by
      simpa [Function.comp_def, Complex.real_smul, mul_comm] using hchar
    convert this using 1
    simp only [Ed, g, E, L]
    push_cast
    ring
  set q : ℝ → ℂ := fun x ↦ (g x)⁻¹ with hq
  set qd : ℝ → ℂ := fun x ↦ -g1 x / g x ^ 2 with hqd
  set qdd : ℝ → ℂ := fun x ↦ -g2 x / g x ^ 2 + 2 * g1 x ^ 2 / g x ^ 3 with hqdd
  have hqder : ∀ x ∈ Icc a b, HasDerivAt q (qd x) x := fun x hx ↦ by
    have := HasDerivAt.fun_inv (hgd x) (hgne x hx)
    simpa [q, qd, neg_div] using this
  have hqdder : ∀ x ∈ Icc a b, HasDerivAt qd (qdd x) x := fun x hx ↦ by
    have hne := hgne x hx
    have h1 : HasDerivAt (fun y ↦ -g1 y / g y ^ 2) _ x :=
      (HasDerivAt.fun_neg (hg1d x)).fun_div ((hgd x).fun_pow 2) (pow_ne_zero 2 hne)
    refine h1.congr_deriv ?_
    simp only [qdd, Nat.cast_ofNat]
    field_simp
    ring
  -- continuity
  have hEc : Continuous E := continuous_iff_continuousAt.2 fun x ↦ (hEder x).continuousAt
  have hgc : Continuous g := continuous_iff_continuousAt.2 fun x ↦ (hgd x).continuousAt
  have hg1c : Continuous g1 := continuous_iff_continuousAt.2 fun x ↦ (hg1d x).continuousAt
  have hg2c : Continuous g2 := by
    simp only [g2]
    exact (continuous_ofReal.comp (continuous_const.mul hD3s.continuous)).mul continuous_const
  have hEdc : ContinuousOn Ed (Icc a b) := (hgc.mul hEc).continuousOn
  have hqddc : ContinuousOn qdd (Icc a b) := by
    apply ContinuousOn.add
    · exact (hg2c.neg.continuousOn).div (hgc.pow 2).continuousOn
        fun x hx ↦ pow_ne_zero 2 (hgne x hx)
    · exact (continuous_const.mul (hg1c.pow 2)).continuousOn.div (hgc.pow 3).continuousOn
        fun x hx ↦ pow_ne_zero 3 (hgne x hx)
  have hw2c : ContinuousOn (fun x ↦ (w2 x : ℂ)) (Icc a b) :=
    (continuous_ofReal.comp hw2s.continuous).continuousOn
  -- endpoint vanishing
  have haz : w =ᶠ[𝓝 a] fun _ ↦ 0 :=
    eventually_zero_of_notMem hw0 (Or.inl haa')
  have hbz : w =ᶠ[𝓝 b] fun _ ↦ 0 :=
    eventually_zero_of_notMem hw0 (Or.inr hb'b)
  have hw1a : w1 a = 0 := by simp [w1, haz.deriv_eq]
  have hw1b : w1 b = 0 := by simp [w1, hbz.deriv_eq]
  -- budgets
  have hbud : ∀ (u : ℝ → ℝ), (∀ x, |u x| ≤ W) →
      (∫ x in a..b, ‖(u x : ℂ)‖) ≤ (b - a) * W := by
    intro u hu
    have := intervalIntegral.norm_integral_le_of_norm_le_const (a := a) (b := b)
      (f := fun x ↦ ‖(u x : ℂ)‖) (C := W) (fun x _ ↦ by simpa using hu x)
    rw [abs_of_nonneg (sub_nonneg.2 hab)] at this
    refine le_trans (le_abs_self _) ?_
    rw [← Real.norm_eq_abs]
    linarith
  have hQ0 : ∀ x ∈ Icc a b, ‖q x‖ ≤ 1 / (L * d) := fun x hx ↦ by
    simp only [q, norm_inv]
    rw [one_div]
    exact inv_anti₀ (by positivity) (hglow x hx)
  have hQ1 : ∀ x ∈ Icc a b, ‖qd x‖ ≤ C2 / (L * d ^ 2) := fun x hx ↦ by
    simp only [qd, norm_div, norm_neg, norm_pow, hnormg1]
    have hl := hglow x hx
    have hgpos : 0 < ‖g x‖ := lt_of_lt_of_le (by positivity) hl
    rw [div_le_div_iff₀ (pow_pos hgpos 2) (by positivity)]
    have : (L * d) ^ 2 ≤ ‖g x‖ ^ 2 := pow_le_pow_left₀ (by positivity) hl 2
    have h2x := h2 x hx
    calc L * |D2 x| * (L * d ^ 2) = L * (L * d ^ 2) * |D2 x| := by ring
      _ ≤ L * (L * d ^ 2) * C2 := by gcongr
      _ = C2 * (L * d) ^ 2 := by ring
      _ ≤ C2 * ‖g x‖ ^ 2 := by gcongr
  have hQ2 : ∀ x ∈ Icc a b, ‖qdd x‖ ≤ C3 / (L * d ^ 2) + 2 * C2 ^ 2 / (L * d ^ 3) :=
    fun x hx ↦ by
    have hl := hglow x hx
    have hgpos : 0 < ‖g x‖ := lt_of_lt_of_le (by positivity) hl
    refine (norm_add_le _ _).trans (add_le_add ?_ ?_)
    · simp only [norm_div, norm_neg, norm_pow, hnormg2]
      rw [div_le_div_iff₀ (pow_pos hgpos 2) (by positivity)]
      have : (L * d) ^ 2 ≤ ‖g x‖ ^ 2 := pow_le_pow_left₀ (by positivity) hl 2
      have h3x := h3 x hx
      calc L * |D3 x| * (L * d ^ 2) = L * (L * d ^ 2) * |D3 x| := by ring
        _ ≤ L * (L * d ^ 2) * C3 := by gcongr
        _ = C3 * (L * d) ^ 2 := by ring
        _ ≤ C3 * ‖g x‖ ^ 2 := by gcongr
    · simp only [norm_div, norm_mul, norm_pow, hnormg1, Complex.norm_ofNat]
      rw [div_le_div_iff₀ (pow_pos hgpos 3) (by positivity)]
      have : (L * d) ^ 3 ≤ ‖g x‖ ^ 3 := pow_le_pow_left₀ (by positivity) hl 3
      have h2x := h2 x hx
      have : (L * |D2 x|) ^ 2 ≤ (L * C2) ^ 2 :=
        pow_le_pow_left₀ (by positivity) (mul_le_mul_of_nonneg_left h2x hLpos.le) 2
      calc 2 * (L * |D2 x|) ^ 2 * (L * d ^ 3) ≤ 2 * (L * C2) ^ 2 * (L * d ^ 3) := by gcongr
        _ = 2 * C2 ^ 2 * (L * d) ^ 3 := by ring
        _ ≤ 2 * C2 ^ 2 * ‖g x‖ ^ 3 := by gcongr
  have hmain := NonstationaryPhase.norm_integral_le_two_ibp_budgets (a := a) (b := b)
    (amplitude := fun x ↦ (w x : ℂ)) (amplitude' := fun x ↦ (w1 x : ℂ))
    (amplitude'' := fun x ↦ (w2 x : ℂ)) hab (fun x _ ↦ hEder x) hqder hqdder
    (fun x _ ↦ (hwd x).ofReal_comp) (fun x _ ↦ (hw1d x).ofReal_comp) hEdc hqddc hw2c
    (fun x hx ↦ by
      simp only [q, Ed]
      rw [← mul_assoc, inv_mul_cancel₀ (hgne x hx), one_mul])
    (by simp [haz.eq_of_nhds]) (by simp [hbz.eq_of_nhds]) (by simp [hw1a]) (by simp [hw1b])
    (hbud w hWw) (hbud w1 hWw1) (hbud w2 hWw2) hQ0 hQ1 hQ2
    (fun x _ ↦ by simp [E]) (by positivity) (by positivity) (by positivity)
  -- rewrite the integral
  have hint : (∫ x : ℝ, (w x : ℂ) * 𝐞 (lam * (H x - ν * x))) =
      ∫ x in a..b, E x * (w x : ℂ) := by
    rw [intervalIntegral.integral_of_le hab, ← setIntegral_eq_integral_of_forall_compl_eq_zero
      (s := Ioc a b)]
    · exact setIntegral_congr_fun measurableSet_Ioc fun x _ ↦ by simp only [E]; ring
    · intro x hx
      have : x ∉ Icc a' b' := fun h ↦ hx ⟨by linarith [h.1], by linarith [h.2]⟩
      simp [hw0 x this]
  rw [hint]
  refine hmain.trans ?_
  have hd2 : 0 < d ^ 2 := by positivity
  have key : (1 / (L * d)) ^ 2 * ((b - a) * W) +
      3 * (1 / (L * d)) * (C2 / (L * d ^ 2)) * ((b - a) * W) +
      (1 / (L * d) * (C3 / (L * d ^ 2) + 2 * C2 ^ 2 / (L * d ^ 3)) +
        (C2 / (L * d ^ 2)) ^ 2) * ((b - a) * W) =
      (b - a) * W * (1 + 3 * C2 / d + (C3 / d + 3 * C2 ^ 2 / d ^ 2)) / (L ^ 2 * d ^ 2) := by
    field_simp
    ring
  rw [key]
  have hnum : 0 ≤ (b - a) * W * (1 + 3 * C2 / d + (C3 / d + 3 * C2 ^ 2 / d ^ 2)) := by
    have : 0 ≤ b - a := sub_nonneg.2 hab
    positivity
  apply div_le_div_of_nonneg_left hnum (by positivity)
  exact mul_le_mul_of_nonneg_right (pow_le_pow_left₀ hlam.le hLlam 2) hd2.le


end
end Expdb.NonstationaryPhase
