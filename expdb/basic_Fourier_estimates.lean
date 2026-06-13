/-
  ANTEDB Blueprint -- Chapter 3 : Lemma 3.1 (L² integral estimate)
  =================================================================
_/

import Mathlib

open MeasureTheory Real Complex Filter Topology BigOperators

def SeparatedFreqs (N : ℝ) {R : ℕ} (ξ : Fin R → ℝ) : Prop :=
  ∀ r s : Fin R, r ≠ s → 1 / N ≤ |ξ r - ξ s|

def expSum {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ) (t : ℝ) : ℂ :=
  ∑ r, a r * e (ξ r * t)

@[simp] lemma expSum_def {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ) (t : ℝ) :
    expSum a ξ t = ∑ r, a r * e (ξ r * t) := rfl

def psiHat (ψ : ℝ → ℝ) (u : ℝ) : ℂ := ∫ x : ℝ, (ψ x : ℂ) * e (-(x * u))

structure SieveAuxiliary where
  ψ : ℝ → ℝ
  smooth : ContDiff ℝ ⊤ ψ
  support : ∀ x, ψ x ≠ 0 → |x| ≤ 1 / 4
  nonneg : ∀ x, 0 ≤ ψ x
  l2norm_one : ∫ x, (ψ x) ^ 2 = 1

def SieveAuxiliary.psiHatNorm (B : SieveAuxiliary) (u : ℝ) : ℝ :=
  ‖psiHat B.ψ u‖

def intervalIndicator (a₀ b₀ t : ℝ) : ℝ :=
  Set.indicator (Set.Icc a₀ b₀) (fun _ => (1:ℝ)) t

def distBoundary (a₀ b₀ t : ℝ) : ℝ := min |t - a₀| |t - b₀|

lemma goal1_normalization {R : ℕ} (a : Fin R → ℂ) (M : ℝ)
    (hM : M = ∑ r, ‖a r‖ ^ 2) (hMpos : 0 < M) :
    ∑ r, ‖a r / (Real.sqrt M : ℂ)‖ ^ 2 = 1 := by
  have hsqrt : (Real.sqrt M : ℝ) ≠ 0 := Real.sqrt_ne_zero'.mpr hMpos
  have hnormsqrt : ‖(Real.sqrt M : ℂ)‖ = Real.sqrt M := by
    rw [Complex.norm_ofReal, abs_of_nonneg (Real.sqrt_nonneg M)]
  calc ∑ r, ‖a r / (Real.sqrt M : ℂ)‖ ^ 2
      = ∑ r, ‖a r‖ ^ 2 / M := by
        apply Finset.sum_congr rfl
        intro r _
        rw [norm_div, div_pow, hnormsqrt, Real.sq_sqrt hMpos.le]
    _ = (∑ r, ‖a r‖ ^ 2) / M := by rw [Finset.sum_div]
    _ = 1 := by rw [← hM]; exact div_self hMpos.ne'

lemma goal1_trivial_case {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (h : ∀ r, a r = 0) (I : Set ℝ) :
    ∫ t in I, ‖expSum a ξ t‖ ^ 2 = 0 := by
  have hzero : ∀ t, expSum a ξ t = 0 := by
    intro t
    apply Finset.sum_eq_zero
    intro r _
    rw [h r, zero_mul]
  simp [hzero]

theorem goal2_eq_3_1 (B : SieveAuxiliary) {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (N : ℝ) (hN : 0 < N) (t₀ : ℝ)
    (hnorm : ∑ r, ‖a r‖ ^ 2 = 1)
    (hsep : SeparatedFreqs N ξ) :
    ∫ t : ℝ, ‖expSum a ξ t‖ ^ 2 * B.psiHatNorm ((t - t₀) / N) ^ 2 = N := by
  sorry

theorem goal3_eq_3_2 (B : SieveAuxiliary) {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (N : ℝ) (hN : 0 < N)
    (hnorm : ∑ r, ‖a r‖ ^ 2 = 1)
    (hsep : SeparatedFreqs N ξ)
    (a₀ b₀ : ℝ) (hlen : b₀ - a₀ = N) :
    ∃ C : ℝ, 0 < C ∧ ∫ t in Set.Icc a₀ b₀, ‖expSum a ξ t‖ ^ 2 ≤ C * N := by
  sorry

theorem goal4_identity (B : SieveAuxiliary) {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (N : ℝ) (hN : 0 < N)
    (hnorm : ∑ r, ‖a r‖ ^ 2 = 1)
    (hsep : SeparatedFreqs N ξ)
    (a₀ b₀ T : ℝ) (hT : T = b₀ - a₀) (hab : a₀ ≤ b₀) :
    ∫ t in Set.Icc a₀ b₀, ‖expSum a ξ t‖ ^ 2
      = T - ∫ t : ℝ, ‖expSum a ξ t‖ ^ 2 *
          ((1 / N) * (∫ t₀ in Set.Icc a₀ b₀, B.psiHatNorm ((t - t₀) / N) ^ 2)
            - intervalIndicator a₀ b₀ t) := by
  sorry

theorem goal5_decay (B : SieveAuxiliary)
    (N : ℝ) (hN : 0 < N) (a₀ b₀ : ℝ) (hab : a₀ ≤ b₀) :
    ∃ C : ℝ, 0 < C ∧ ∀ t : ℝ,
      |(1 / N) * (∫ t₀ in Set.Icc a₀ b₀, B.psiHatNorm ((t - t₀) / N) ^ 2)
        - intervalIndicator a₀ b₀ t|
      ≤ C * (1 + distBoundary a₀ b₀ t / N) ^ (-(10:ℝ)) := by
  sorry

theorem goal6_dyadic (B : SieveAuxiliary) {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (N : ℝ) (hN : 0 < N)
    (hnorm : ∑ r, ‖a r‖ ^ 2 = 1)
    (hsep : SeparatedFreqs N ξ)
    (a₀ b₀ : ℝ) (hab : a₀ ≤ b₀) :
    ∃ C : ℝ, 0 < C ∧
      |∫ t : ℝ, ‖expSum a ξ t‖ ^ 2 *
          ((1 / N) * (∫ t₀ in Set.Icc a₀ b₀, B.psiHatNorm ((t - t₀) / N) ^ 2)
            - intervalIndicator a₀ b₀ t)|
      ≤ C * N := by
  sorry

theorem large_sieve_inequality (B : SieveAuxiliary) {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (N : ℝ) (hN : 0 < N) (hsep : SeparatedFreqs N ξ)
    (a₀ b₀ T : ℝ) (hT : T = b₀ - a₀) (hab : a₀ ≤ b₀) :
    ∃ C : ℝ, 0 < C ∧ ∃ θ : ℝ, |θ| ≤ C ∧
      ∫ t in Set.Icc a₀ b₀, ‖expSum a ξ t‖ ^ 2 = (T + θ * N) * ∑ r, ‖a r‖ ^ 2 := by
  set M : ℝ := ∑ r, ‖a r‖ ^ 2 with hM
  by_cases hM0 : M = 0
  · refine ⟨1, one_pos, 0, by simp, ?_⟩
    have hzero : ∀ r, a r = 0 := by
      intro r
      have hle : ‖a r‖ ^ 2 ≤ M := by
        rw [hM]
        exact Finset.single_le_sum (fun i _ => sq_nonneg ‖a i‖) (Finset.mem_univ r)
      have hge : 0 ≤ ‖a r‖ ^ 2 := sq_nonneg _
      have hsq : ‖a r‖ ^ 2 = 0 := le_antisymm (hM0 ▸ hle) hge
      exact norm_eq_zero.mp (pow_eq_zero_iff (n := 2) (by norm_num) |>.mp hsq)
    rw [goal1_trivial_case a ξ hzero (Set.Icc a₀ b₀), hM0]
    simp
  · have hMpos : 0 < M := by
      have hMnn : 0 ≤ M := by rw [hM]; exact Finset.sum_nonneg (fun r _ => sq_nonneg _)
      exact lt_of_le_of_ne hMnn (Ne.symm hM0)
    set A : Fin R → ℂ := fun r => a r / (Real.sqrt M : ℂ) with hA
    have hAnorm : ∑ r, ‖A r‖ ^ 2 = 1 := goal1_normalization a M hM hMpos
    have hexpand : ∀ t : ℝ, ‖expSum a ξ t‖ ^ 2 = M * ‖expSum A ξ t‖ ^ 2 := by
      intro t
      have hrescale : expSum a ξ t = (Real.sqrt M : ℂ) * expSum A ξ t := by
        rw [expSum_def, expSum_def, Finset.mul_sum]
        apply Finset.sum_congr rfl
        intro r _
        rw [hA]
        have hsqrt_ne : (Real.sqrt M : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr (Real.sqrt_ne_zero'.mpr hMpos)
        field_simp
      rw [hrescale, norm_mul, Complex.norm_ofReal, abs_of_nonneg (Real.sqrt_nonneg M), mul_pow, Real.sq_sqrt hMpos.le]
    have h4 := goal4_identity B A ξ N hN hAnorm hsep a₀ b₀ T hT hab
    obtain ⟨C, hCpos, hC⟩ := goal6_dyadic B A ξ N hN hAnorm hsep a₀ b₀ hab
    have hMint : ∫ t in Set.Icc a₀ b₀, ‖expSum a ξ t‖ ^ 2 = M * ∫ t in Set.Icc a₀ b₀, ‖expSum A ξ t‖ ^ 2 := by
      simp_rw [hexpand]
      exact MeasureTheory.integral_const_mul M _
    set err : ℝ := ∫ t : ℝ, ‖expSum A ξ t‖ ^ 2 *
        ((1 / N) * (∫ t₀ in Set.Icc a₀ b₀, B.psiHatNorm ((t - t₀) / N) ^ 2) - intervalIndicator a₀ b₀ t) with herr
    have h4' : ∫ t in Set.Icc a₀ b₀, ‖expSum A ξ t‖ ^ 2 = T - err := h4
    refine ⟨C, hCpos, -err / N, ?_, ?_⟩
    · rw [abs_div, abs_neg, abs_of_pos hN]
      have step : |err| / N ≤ (C * N) / N := (div_le_div_right hN).mpr hC
      rwa [mul_div_cancel_right₀ C (ne_of_gt hN)] at step
    · rw [hMint, h4']
      have hNne : N ≠ 0 := ne_of_gt hN
      field_simp; ring

end
