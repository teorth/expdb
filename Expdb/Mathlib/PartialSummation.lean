module

public import Mathlib.Algebra.BigOperators.Module
import Mathlib.Algebra.Order.BigOperators.Group.Finset
public import Mathlib.Analysis.Complex.Basic

/-!
# Bounds from partial summation

Weighted finite sums are bounded by their unweighted partial sums and the variation
of the weights. The Lipschitz specialization applies to weights sampled at a real scale.
-/

@[expose] public section

namespace Expdb.PartialSummation

noncomputable section

/-- Partial summation with partial sums starting at `P`. -/
theorem norm_sum_Ico_mul_le (a z : ℕ → ℂ) {P Q : ℕ} (hPQ : P ≤ Q) {D : ℝ}
    (hD : ∀ k ∈ Finset.Icc P Q, ‖∑ m ∈ Finset.Ico P k, z m‖ ≤ D) :
    ‖∑ m ∈ Finset.Ico P Q, a m * z m‖ ≤
      (‖a (Q - 1)‖ + ∑ m ∈ Finset.Ico P (Q - 1), ‖a (m + 1) - a m‖) * D := by
  have hD0 : 0 ≤ D := by
    have := hD P (Finset.mem_Icc.2 ⟨le_rfl, hPQ⟩)
    simp only [Finset.Ico_self, Finset.sum_empty, norm_zero] at this
    exact this
  rcases eq_or_lt_of_le hPQ with rfl | hlt
  · simp only [Finset.Ico_self, Finset.sum_empty, norm_zero]
    exact mul_nonneg (add_nonneg (norm_nonneg _) (Finset.sum_nonneg fun _ _ ↦ norm_nonneg _)) hD0
  set g : ℕ → ℂ := fun i ↦ if P ≤ i then z i else 0 with hg
  have hparts := Finset.sum_Ico_by_parts a g hlt
  have hlhs : ∑ i ∈ Finset.Ico P Q, a i • g i = ∑ m ∈ Finset.Ico P Q, a m * z m :=
    Finset.sum_congr rfl fun i hi ↦ by
      rw [smul_eq_mul, hg]; simp only [if_pos (Finset.mem_Ico.1 hi).1]
  rw [hlhs] at hparts
  rw [hparts]
  have hZ : ∀ k, P ≤ k → ∑ i ∈ Finset.range k, g i = ∑ i ∈ Finset.Ico P k, z i :=
    by
      intro k hk
      change (∑ i ∈ Finset.range k, if P ≤ i then z i else 0) = _
      rw [← Finset.sum_filter]
      congr 1
      ext i
      simp only [Finset.mem_filter, Finset.mem_range, Finset.mem_Ico]
      omega
  rw [hZ Q hPQ, hZ P le_rfl, Finset.Ico_self, Finset.sum_empty, smul_zero, sub_zero]
  refine (norm_sub_le _ _).trans ?_
  rw [add_mul, Finset.sum_mul]
  apply add_le_add
  · rw [smul_eq_mul, norm_mul]
    exact mul_le_mul_of_nonneg_left (hD Q (Finset.mem_Icc.2 ⟨hPQ, le_rfl⟩)) (norm_nonneg _)
  · refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun i hi ↦ ?_)
    have hi' := Finset.mem_Ico.1 hi
    rw [smul_eq_mul, norm_mul, hZ (i + 1) (by omega)]
    exact mul_le_mul_of_nonneg_left (hD (i + 1) (Finset.mem_Icc.2 ⟨by omega, by omega⟩))
      (norm_nonneg _)

/-- Partial summation with Lipschitz real weights `A (m / M)`. -/
theorem norm_sum_Ico_lipschitz_mul_le (A : ℝ → ℝ) {Amax LA M : ℝ} (hM : 0 < M)
    (hA : ∀ ν, |A ν| ≤ Amax) (hLip : ∀ ν ν', |A ν - A ν'| ≤ LA * |ν - ν'|)
    (z : ℕ → ℂ) {P Q : ℕ} (hPQ : P ≤ Q) {D : ℝ}
    (hD : ∀ k ∈ Finset.Icc P Q, ‖∑ m ∈ Finset.Ico P k, z m‖ ≤ D) :
    ‖∑ m ∈ Finset.Ico P Q, ((A (m / M) : ℝ) : ℂ) * z m‖ ≤
      (Amax + LA * ((Q : ℝ) - P) / M) * D := by
  have hD0 : 0 ≤ D := by
    have := hD P (Finset.mem_Icc.2 ⟨le_rfl, hPQ⟩)
    simp only [Finset.Ico_self, Finset.sum_empty, norm_zero] at this
    exact this
  have hLA : 0 ≤ LA := by
    have := hLip 1 0
    norm_num at this
    linarith [abs_nonneg (A 1 - A 0)]
  refine (norm_sum_Ico_mul_le _ z hPQ hD).trans (mul_le_mul_of_nonneg_right ?_ hD0)
  apply add_le_add
  · rw [Complex.norm_real, Real.norm_eq_abs]
    exact hA _
  · have hstep : ∀ m ∈ Finset.Ico P (Q - 1),
        ‖((A (((m + 1 : ℕ) : ℝ) / M) : ℝ) : ℂ) - ((A ((m : ℝ) / M) : ℝ) : ℂ)‖ ≤ LA / M := by
      intro m _
      rw [← Complex.ofReal_sub, Complex.norm_real, Real.norm_eq_abs]
      refine (hLip _ _).trans (le_of_eq ?_)
      push_cast
      rw [show ((m : ℝ) + 1) / M - m / M = 1 / M by field_simp; ring,
        abs_of_pos (by positivity)]
      ring
    refine (Finset.sum_le_sum hstep).trans ?_
    rw [Finset.sum_const, Nat.card_Ico, nsmul_eq_mul]
    have hcount : ((Q - 1 - P : ℕ) : ℝ) ≤ (Q : ℝ) - P := by
      have h1 : (Q - 1 - P : ℕ) ≤ Q - P := by omega
      calc ((Q - 1 - P : ℕ) : ℝ) ≤ ((Q - P : ℕ) : ℝ) := by exact_mod_cast h1
        _ = (Q : ℝ) - P := by rw [Nat.cast_sub hPQ]
    calc ((Q - 1 - P : ℕ) : ℝ) * (LA / M) ≤ ((Q : ℝ) - P) * (LA / M) :=
          mul_le_mul_of_nonneg_right hcount (div_nonneg hLA hM.le)
      _ = LA * ((Q : ℝ) - P) / M := by ring

end

end Expdb.PartialSummation
