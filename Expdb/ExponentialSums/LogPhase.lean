module

public import Expdb.ExponentialSums.PhaseFunctions
public import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic
public import Mathlib.Order.Interval.Finset.Nat

import Mathlib.Analysis.Fourier.FourierTransformDeriv
import Mathlib.Analysis.Real.Pi.Bounds
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv
import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts

/-!
# The logarithmic phase

Reusable facts about the logarithmic model phase, including its derivatives, separation on
dyadic intervals, and the explicit main term of its oscillatory integral.
-/

@[expose] public section

open Filter Topology
open scoped ContDiff Expdb FourierTransform

noncomputable section

namespace Expdb

/-- The fixed logarithmic phase `u ↦ log u`, regarded as a variable family. -/
def logPhase : VariableFunction (VariableObject.fixed ℝ) ℝ :=
  fun _ ↦ Real.log

/-- The derivatives of `log` on the phase interval agree with the reference derivatives for
model exponent one. -/
theorem iteratedDerivWithin_log_eq_rpow_neg_one
    (p : ℕ) (u : phaseInterval) :
    iteratedDerivWithin (p + 1) Real.log phaseInterval u =
      iteratedDerivWithin p (fun v : ℝ ↦ v ^ (-(1 : ℝ))) phaseInterval u := by
  have hu_pos : 0 < (u : ℝ) := lt_of_lt_of_le zero_lt_one u.property.1
  have hunique : UniqueDiffOn ℝ phaseInterval :=
    uniqueDiffOn_Icc (by norm_num [phaseInterval])
  rw [iteratedDerivWithin_eq_iteratedDeriv hunique
      (Real.contDiffAt_log.2 hu_pos.ne') u.property,
    iteratedDerivWithin_eq_iteratedDeriv hunique
      (Real.contDiffAt_rpow_const_of_ne hu_pos.ne') u.property,
    iteratedDeriv_succ', Real.deriv_log']
  congr 1
  funext v
  exact (Real.rpow_neg_one v).symm

/-- The fixed logarithmic phase `u ↦ log u` is a model phase function, with model exponent
`σ = 1`. -/
theorem isModelPhaseFunction_log : IsModelPhaseFunction logPhase := by
  refine ⟨?_, 1, zero_lt_one, ?_⟩
  · intro i
    change ContDiffOn ℝ ∞ Real.log phaseInterval
    exact Real.contDiffOn_log.mono fun u hu ↦ by
      simp only [Set.mem_compl_iff, Set.mem_singleton_iff]
      exact ne_of_gt (lt_of_lt_of_le zero_lt_one hu.1)
  · intro p u
    rw [VariableObject.IsInfinitesimal]
    simp only [modelPhaseError, logPhase, iteratedDerivWithin_log_eq_rpow_neg_one, sub_self,
      norm_zero]
    exact tendsto_const_nhds

/-- On a dyadic interval, `log (n / N)` is `1 / (2 * N)`-separated. -/
theorem log_div_separated
    {N : ℝ} {a b : ℕ} (hN : 0 < N) (ha : N ≤ a) (hb : (b : ℝ) ≤ 2 * N) :
    IsSeparatedFamily (1 / (2 * N))
      (fun n : ↥(Finset.Icc a b) ↦ Real.log ((n : ℝ) / N)) := by
  have hordered (m n : ↥(Finset.Icc a b)) (hmn : (m : ℕ) < n) :
      1 / (2 * N) ≤ Real.log ((n : ℝ) / N) - Real.log ((m : ℝ) / N) := by
    have hm_mem := Finset.mem_Icc.mp m.property
    have hmpos : 0 < (m : ℝ) := lt_of_lt_of_le hN (ha.trans (by exact_mod_cast hm_mem.1))
    have hnpos : 0 < (n : ℝ) := lt_trans hmpos (by exact_mod_cast hmn)
    have hlog : 1 - (m : ℝ) / n ≤ Real.log ((n : ℝ) / m) := by
      simpa [inv_div] using Real.one_sub_inv_le_log_of_pos (div_pos hnpos hmpos)
    have hdiff : 1 / (2 * N) ≤ 1 - (m : ℝ) / n := by
      have hgapNat : (m : ℕ) + 1 ≤ (n : ℕ) := Nat.add_one_le_iff.mpr hmn
      have hgapCast : ((m : ℕ) : ℝ) + 1 ≤ ((n : ℕ) : ℝ) := by exact_mod_cast hgapNat
      have hgap : (1 : ℝ) ≤ (n : ℝ) - m := by linarith
      have hn_mem := Finset.mem_Icc.mp n.property
      have hnle : (n : ℝ) ≤ 2 * N := le_trans (by exact_mod_cast hn_mem.2) hb
      rw [one_sub_div hnpos.ne']
      calc
        1 / (2 * N) ≤ 1 / (n : ℝ) := one_div_le_one_div_of_le hnpos hnle
        _ ≤ ((n : ℝ) - m) / n := (div_le_div_iff_of_pos_right hnpos).2 hgap
    calc
      1 / (2 * N) ≤ Real.log ((n : ℝ) / m) := hdiff.trans hlog
      _ = Real.log ((n : ℝ) / N) - Real.log ((m : ℝ) / N) := by
        rw [Real.log_div (ne_of_gt hnpos) hN.ne', Real.log_div (ne_of_gt hmpos) hN.ne',
          Real.log_div (ne_of_gt hnpos) (ne_of_gt hmpos)]
        ring
  intro m n hmn
  rw [Real.dist_eq]
  rcases lt_or_gt_of_ne (Subtype.coe_ne_coe.mpr hmn) with hmn' | hnm'
  · have hbound := hordered m n hmn'
    have hsign : Real.log ((m : ℝ) / N) - Real.log ((n : ℝ) / N) ≤ 0 := by
      have : 0 < 1 / (2 * N) := by positivity
      linarith
    rw [abs_of_nonpos hsign, neg_sub]
    exact hbound
  · have hbound := hordered n m hnm'
    have hsign : 0 ≤ Real.log ((m : ℝ) / N) - Real.log ((n : ℝ) / N) := by
      have : 0 < 1 / (2 * N) := by positivity
      linarith
    rw [abs_of_nonneg hsign]
    exact hbound

/-- The explicit main term for the logarithmic phase. With Lean's Fourier-character
normalization, the corresponding sum is a norm-one multiple of `∑ n ^ (2π i T)`. -/
def logPhaseMainTerm (N T : ℝ) : ℂ :=
  (N : ℂ) * (2 * 𝐞 (T * Real.log 2) - 1) /
    (1 + (2 * Real.pi : ℂ) * Complex.I * T)

/-- The logarithmic oscillatory integral is exactly `logPhaseMainTerm`. -/
theorem logPhase_integral_eq_mainTerm
    {N T : ℝ} (hN : 0 < N) :
    (∫ x in N..2 * N, (𝐞 (T * Real.log (x / N)) : ℂ)) =
      logPhaseMainTerm N T := by
  let τ : ℂ := (2 * Real.pi : ℂ) * Complex.I
  let v : ℝ → ℂ := fun x ↦ 𝐞 (T * Real.log (x / N))
  let denom : ℂ := 1 + τ * T
  let G : ℝ → ℂ := fun x ↦ x * v x / denom
  have hdenom : denom ≠ 0 := by
    intro h
    have hre := congrArg Complex.re h
    simp [denom, τ] at hre
  have hv (x : ℝ) (hx : 0 < x) :
      HasDerivAt v (τ * T / x * v x) x := by
    have hdiv : HasDerivAt (fun y : ℝ ↦ y / N) (1 / N) x := by
      simpa using (hasDerivAt_id x).div_const N
    have hlog := (Real.hasDerivAt_log (div_pos hx hN).ne').comp x hdiv
    have hinner : HasDerivAt (fun y ↦ T * Real.log (y / N)) (T / x) x := by
      convert hlog.const_mul T using 1 <;>
        first | rfl | field_simp [hN.ne', hx.ne']
    have hchar' := (Real.hasDerivAt_fourierChar (T * Real.log (x / N))).scomp x hinner
    dsimp only [v]
    change HasDerivAt ((𝐞 · : ℝ → ℂ) ∘ fun y ↦ T * Real.log (y / N))
      (τ * T / x * (𝐞 (T * Real.log (x / N)) : ℂ)) x
    rw [show τ * T / x * (𝐞 (T * Real.log (x / N)) : ℂ) =
        (T / x) • ((2 * Real.pi : ℂ) * Complex.I * 𝐞 (T * Real.log (x / N))) by
      dsimp [τ]
      simp only [Complex.ofReal_div]
      ring]
    exact hchar'
  have hG (x : ℝ) (hx : 0 < x) : HasDerivAt G (v x) x := by
    have hprod := (hasDerivAt_id x).ofReal_comp.mul (hv x hx)
    have hquot := hprod.div_const denom
    have hderivEq :
        ((1 : ℂ) * v x + x * (τ * T / x * v x)) / denom = v x := by
      apply (div_eq_iff hdenom).2
      dsimp [denom, τ]
      field_simp [hx.ne']
    rw [← hderivEq]
    simpa only [G, Pi.mul_apply, id_eq, Complex.ofReal_one] using hquot
  have hvCont : ContinuousOn v (Set.uIcc N (2 * N)) := by
    intro x hx
    have hxIcc : x ∈ Set.Icc N (2 * N) := by
      simpa [Set.uIcc_of_le (by linarith : N ≤ 2 * N)] using hx
    exact (hv x (hN.trans_le hxIcc.1)).continuousAt.continuousWithinAt
  have hFTC := intervalIntegral.integral_eq_sub_of_hasDerivAt
    (a := N) (b := 2 * N) (f := G) (f' := v)
    (fun x hx ↦ by
      have hxIcc : x ∈ Set.Icc N (2 * N) := by
        simpa [Set.uIcc_of_le (by linarith : N ≤ 2 * N)] using hx
      exact hG x (hN.trans_le hxIcc.1)) hvCont.intervalIntegrable
  have hformula : (∫ x in N..2 * N, v x) =
      (N : ℂ) * (2 * 𝐞 (T * Real.log 2) - 1) / denom := by
    rw [hFTC]
    dsimp [G, v]
    rw [show (2 * N) / N = 2 by field_simp [hN.ne'],
      show N / N = 1 by field_simp [hN.ne'], Real.log_one]
    simp only [mul_zero]
    rw [show (𝐞 (0 : ℝ) : ℂ) = 1 by norm_num]
    push_cast
    field_simp [hdenom]
  simpa only [v, logPhaseMainTerm, denom, τ, mul_assoc] using hformula

/-- The logarithmic main term has size comparable to `N / T`. -/
theorem norm_logPhaseMainTerm_bounds
    {N T : ℝ} (hN : 0 < N) (hT : 1 ≤ T) :
    (1 / (1 + 2 * Real.pi)) * (N / T) ≤ ‖logPhaseMainTerm N T‖ ∧
      ‖logPhaseMainTerm N T‖ ≤ (3 / (2 * Real.pi)) * (N / T) := by
  let τ : ℂ := (2 * Real.pi : ℂ) * Complex.I
  let denom : ℂ := 1 + τ * T
  have hdenom : denom ≠ 0 := by
    intro h
    have hre := congrArg Complex.re h
    simp [denom, τ] at hre
  have hnumLower : 1 ≤ ‖(2 : ℂ) * 𝐞 (T * Real.log 2) - 1‖ := by
    calc
      1 = ‖(2 : ℂ) * 𝐞 (T * Real.log 2)‖ - ‖(1 : ℂ)‖ := by norm_num
      _ ≤ ‖(2 : ℂ) * 𝐞 (T * Real.log 2) - 1‖ := norm_sub_norm_le _ _
  have hnumUpper : ‖(2 : ℂ) * 𝐞 (T * Real.log 2) - 1‖ ≤ 3 := by
    calc
      ‖(2 : ℂ) * 𝐞 (T * Real.log 2) - 1‖ ≤
          ‖(2 : ℂ) * 𝐞 (T * Real.log 2)‖ + ‖(1 : ℂ)‖ := norm_sub_le _ _
      _ = 3 := by norm_num
  have hdenUpper : ‖denom‖ ≤ (1 + 2 * Real.pi) * T := by
    calc
      ‖denom‖ ≤ ‖(1 : ℂ)‖ + ‖τ * T‖ := by
        dsimp [denom]
        exact norm_add_le _ _
      _ = 1 + 2 * Real.pi * T := by
        dsimp [τ]
        simp [Complex.norm_real, Real.norm_eq_abs, Real.pi_pos.le,
          abs_of_nonneg (le_trans zero_le_one hT)]
      _ ≤ (1 + 2 * Real.pi) * T := by
        nlinarith [Real.pi_pos]
  have hTpos : 0 < T := zero_lt_one.trans_le hT
  have hdenLower : 2 * Real.pi * T ≤ ‖denom‖ := by
    have him := Complex.abs_im_le_norm denom
    simpa [denom, τ, abs_of_pos Real.pi_pos, abs_of_pos hTpos] using him
  have hdenpos : 0 < ‖denom‖ := norm_pos_iff.mpr hdenom
  have hcoefpos : 0 < 1 + 2 * Real.pi := by positivity
  constructor
  · rw [logPhaseMainTerm, show (1 + (2 * Real.pi : ℂ) * Complex.I * T) = denom by
      simp only [denom, τ, mul_assoc], norm_div, norm_mul,
      Complex.norm_real, Real.norm_eq_abs, abs_of_pos hN]
    apply (le_div_iff₀ hdenpos).2
    calc
      (1 / (1 + 2 * Real.pi)) * (N / T) * ‖denom‖ ≤
          (1 / (1 + 2 * Real.pi)) * (N / T) * ((1 + 2 * Real.pi) * T) := by
        gcongr
      _ = N := by field_simp [hcoefpos.ne', hTpos.ne']
      _ ≤ N * ‖(2 : ℂ) * 𝐞 (T * Real.log 2) - 1‖ := by
        simpa only [mul_one] using mul_le_mul_of_nonneg_left hnumLower hN.le
  · rw [logPhaseMainTerm, show (1 + (2 * Real.pi : ℂ) * Complex.I * T) = denom by
      simp only [denom, τ, mul_assoc], norm_div, norm_mul,
      Complex.norm_real, Real.norm_eq_abs, abs_of_pos hN]
    apply (div_le_iff₀ hdenpos).2
    calc
      N * ‖(2 : ℂ) * 𝐞 (T * Real.log 2) - 1‖ ≤ N * 3 :=
        mul_le_mul_of_nonneg_left hnumUpper hN.le
      _ = (3 / (2 * Real.pi) * (N / T)) * (2 * Real.pi * T) := by
        field_simp [Real.pi_ne_zero, hTpos.ne']
      _ ≤ (3 / (2 * Real.pi) * (N / T)) * ‖denom‖ := by
        exact mul_le_mul_of_nonneg_left hdenLower (by positivity)


end Expdb

end
