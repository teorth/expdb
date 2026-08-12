module

public import Expdb.ExponentialSums.PhaseFunctions
public import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic
public import Mathlib.Order.Interval.Finset.Nat

import Mathlib.Analysis.Fourier.FourierTransformDeriv
import Mathlib.Analysis.Real.Pi.Bounds
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv

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
      iteratedDerivWithin p (modelPhase 1) phaseInterval u := by
  have hu_pos : 0 < (u : ℝ) := lt_of_lt_of_le zero_lt_one u.property.1
  rw [show modelPhase 1 = fun v : ℝ ↦ v ^ (-(1 : ℝ)) from rfl]
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
    simp only [modelPhaseError_apply, modelPhaseErrorAt, logPhase,
      iteratedDerivWithin_log_eq_rpow_neg_one, sub_self, norm_zero]
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

/-- The explicit main term for the logarithmic oscillatory integral, equal to that integral when
`N > 0`. The `2π` factor in the normalization comes from Lean's Fourier character `𝐞`. -/
def logPhaseMainTerm (N T : ℝ) : ℂ :=
  (N : ℂ) * (2 * 𝐞 (T * Real.log 2) - 1) /
    (1 + (2 * Real.pi : ℂ) * Complex.I * T)

/-- The logarithmic oscillatory integral is exactly `logPhaseMainTerm`. -/
theorem logPhase_integral_eq_mainTerm
    {N T : ℝ} (hN : 0 < N) :
    (∫ x in N..2 * N, (𝐞 (T * Real.log (x / N)) : ℂ)) =
      logPhaseMainTerm N T := by
  let r : ℂ := (2 * Real.pi : ℂ) * Complex.I * T
  have hchar (u : ℝ) (hu : 0 < u) :
      (𝐞 (T * Real.log u) : ℂ) = (u : ℂ) ^ r := by
    rw [Real.fourierChar_apply,
      Complex.cpow_def_of_ne_zero (Complex.ofReal_ne_zero.mpr hu.ne')]
    rw [← Complex.ofReal_log hu.le]
    congr 1
    dsimp [r]
    push_cast
    ring
  have hscale :
      (∫ x in N..2 * N, (𝐞 (T * Real.log (x / N)) : ℂ)) =
        (N : ℂ) * ∫ u in (1 : ℝ)..2, (𝐞 (T * Real.log u) : ℂ) := by
    have h := intervalIntegral.smul_integral_comp_mul_left
      (f := fun x : ℝ ↦ (𝐞 (T * Real.log (x / N)) : ℂ))
      (a := (1 : ℝ)) (b := 2) N
    rw [show N * (1 : ℝ) = N by ring, show N * 2 = 2 * N by ring] at h
    rw [← h, Complex.real_smul]
    refine congrArg (fun z : ℂ ↦ (N : ℂ) * z) ?_
    apply intervalIntegral.integral_congr
    intro u _
    exact congrArg (fun y : ℝ ↦ (𝐞 (T * Real.log y) : ℂ))
      (by field_simp [hN.ne'] : N * u / N = u)
  rw [hscale]
  have hint :
      (∫ u in (1 : ℝ)..2, (𝐞 (T * Real.log u) : ℂ)) =
        ((2 : ℂ) ^ (r + 1) - (1 : ℂ) ^ (r + 1)) / (r + 1) := by
    rw [intervalIntegral.integral_congr_Ioo_of_le (by norm_num)
      (fun u hu ↦ hchar u (by linarith [hu.1]))]
    exact integral_cpow (Or.inl (by simp [r]))
  rw [hint, Complex.one_cpow]
  rw [Complex.cpow_add _ _ (by norm_num : (2 : ℂ) ≠ 0), Complex.cpow_one]
  have htwo : (2 : ℂ) ^ r = (𝐞 (T * Real.log 2) : ℂ) :=
    (hchar 2 (by norm_num)).symm
  simp only [htwo]
  dsimp [logPhaseMainTerm, r]
  ring

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
