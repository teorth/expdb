module

public import Expdb.ExponentialSums.PhaseFunctions
public import Mathlib.Algebra.Order.BigOperators.Group.Finset
public import Mathlib.Analysis.Calculus.ContDiff.Defs
public import Mathlib.Order.Interval.Finset.Nat

/-!
# Oscillatory factors and fixed-parameter exponential sums

This module contains the fixed-parameter vocabulary shared by the asymptotic,
non-asymptotic, and analytic exponential-sum estimates.
-/

@[expose] public section

open scoped ContDiff FourierTransform

noncomputable section

namespace Expdb

/-- The oscillatory factor `x ↦ e(T F(x / N))`. -/
def oscillatory (F : ℝ → ℝ) (T N : ℝ) (x : ℝ) : ℂ :=
  𝐞 (T * F (x / N))

@[simp] theorem norm_oscillatory (F : ℝ → ℝ) (T N x : ℝ) :
    ‖oscillatory F T N x‖ = 1 := by
  simp [oscillatory]

/-- The exponential sum `∑ n ∈ [a, b], e(T F(n / N))` at fixed parameters. -/
def exponentialSumAt (F : ℝ → ℝ) (T N : ℝ) (a b : ℕ) : ℂ :=
  ∑ n ∈ Finset.Icc a b, oscillatory F T N n

/-- An exponential sum has at most as many unit-sized terms as its summation range. -/
theorem norm_exponentialSumAt_le_card (F : ℝ → ℝ) (T N : ℝ) (a b : ℕ) :
    ‖exponentialSumAt F T N a b‖ ≤ ((Finset.Icc a b).card : ℝ) := by
  simpa [exponentialSumAt] using norm_sum_le (Finset.Icc a b) fun n ↦ oscillatory F T N n

/-- Bounding the cardinality by the upper endpoint gives a convenient endpoint-only estimate. -/
theorem norm_exponentialSumAt_le_add_one (F : ℝ → ℝ) (T N : ℝ) (a b : ℕ) :
    ‖exponentialSumAt F T N a b‖ ≤ (b : ℝ) + 1 := by
  calc
    ‖exponentialSumAt F T N a b‖ ≤ ((Finset.Icc a b).card : ℝ) :=
      norm_exponentialSumAt_le_card _ _ _ _ _
    _ ≤ (b : ℝ) + 1 := by
      rw [Nat.card_Icc]
      exact_mod_cast Nat.sub_le _ _

@[simp] theorem exponentialSumAt_self (F : ℝ → ℝ) (T N : ℝ) (a : ℕ) :
    exponentialSumAt F T N a a = oscillatory F T N a := by
  simp [exponentialSumAt]

theorem exponentialSumAt_of_lt {a b : ℕ} (hab : b < a) (F : ℝ → ℝ) (T N : ℝ) :
    exponentialSumAt F T N a b = 0 := by
  simp [exponentialSumAt, Finset.Icc_eq_empty (by omega : ¬ a ≤ b)]

/-- The scale conditions making the Euler–Maclaurin remainder of order `s` negligible. -/
structure HasSmallEulerMaclaurinRemainder (s : ℕ) (T N K : ℝ) : Prop where
  /-- The normalized phase derivative is at most one. -/
  ratio_le_one : K * T / N ≤ 1
  /-- The full remainder of order `s` is at most one. -/
  remainder_le_one : N * (K * T / N) ^ (s + 1) ≤ 1

end Expdb

end
