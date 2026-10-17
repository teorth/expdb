module

public import Expdb.ExponentialSums.FixedExponentialSum
public import Expdb.ExponentialSums.DifferencePhase
public import Expdb.Mathlib.VanDerCorput

import Mathlib.Tactic.FieldSimp

/-!
# Differencing at fixed exponential-sum parameters

This is the bridge from the finite van der Corput inequality to the model-phase
sums used in the blueprint. It keeps the sum over shifts explicit for the dyadic
pigeonhole argument in the A-process.
-/

@[expose] public section

open scoped FourierTransform

noncomputable section

namespace Expdb

/-- The correlation at shift `j`, restricted to the overlap of the original interval
and its translate. The half-open interval also handles an empty overlap. -/
def correlationSumAt (F : ℝ → ℝ) (T N : ℝ) (a b j : ℕ) : ℂ :=
  ∑ n ∈ Finset.Ico a (b + 1 - j),
    (𝐞 (T * (F (((n + j : ℕ) : ℝ) / N) - F ((n : ℝ) / N))) : ℂ)

/-- The fixed-parameter differencing inequality, with an absolute constant. -/
theorem norm_exponentialSumAt_sq_le_correlations
    (F : ℝ → ℝ) (T : ℝ) {N : ℝ} {a b H : ℕ}
    (hN : 1 ≤ N) (hNa : N ≤ (a : ℝ)) (hbN : (b : ℝ) ≤ 2 * N)
    (hH : 1 ≤ H) (hHN : (H : ℝ) ≤ N) :
    ‖exponentialSumAt F T N a b‖ ^ 2 ≤
      6 * N / H * (N + ∑ j ∈ Finset.Ico 1 H, ‖correlationSumAt F T N a b j‖) := by
  by_cases hab : a ≤ b
  · have hHa : H ≤ a := by exact_mod_cast hHN.trans hNa
    have hraw := VanDerCorput.norm_sum_Icc_sq_le_correlations
      (φ := fun n ↦ T * F ((n : ℝ) / N)) hH hHa hab
    simp only [← mul_sub] at hraw
    change ‖exponentialSumAt F T N a b‖ ^ 2 ≤ _ at hraw
    let L : ℝ := ((b + 1 - a : ℕ) : ℝ)
    let C : ℝ := ∑ j ∈ Finset.Ico 1 H, ‖correlationSumAt F T N a b j‖
    have hL : L ≤ 2 * N := by
      dsimp [L]
      rw [Nat.cast_sub (by omega : a ≤ b + 1), Nat.cast_add, Nat.cast_one]
      linarith
    have hL0 : 0 ≤ L := Nat.cast_nonneg _
    have hC : 0 ≤ C := Finset.sum_nonneg fun _ _ ↦ norm_nonneg _
    have hHpos : (0 : ℝ) < H := by exact_mod_cast hH
    calc
      _ ≤ (L + H) / H * (L + 2 * C) := hraw
      _ ≤ (3 * N) / H * (2 * (N + C)) := by
        apply mul_le_mul
        · exact div_le_div_of_nonneg_right (by linarith) hHpos.le
        · linarith
        · positivity
        · positivity
      _ = _ := by dsimp [C]; ring
  · rw [exponentialSumAt_of_lt (by omega : b < a), norm_zero, zero_pow (by decide)]
    positivity

/-- Reindexing the overlap identifies a correlation with a normalized model-phase
sum of scale `N - j`, up to conjugation. -/
theorem norm_correlationSumAt_eq_differencePhase
    (F : ℝ → ℝ) {σ T N : ℝ} {a b j : ℕ}
    (hσ : σ ≠ 0) (hN : N ≠ 0) (hj : j ≠ 0) (hNj : N - (j : ℝ) ≠ 0)
    (hja : j ≤ a) (hjb : 2 * j ≤ b) :
    ‖correlationSumAt F T N a b j‖ =
      ‖exponentialSumAt (differencePhase F σ ((j : ℝ) / N))
        (σ * T * (j : ℝ) / N) (N - j) (a - j) (b - 2 * j)‖ := by
  have hr : (j : ℝ) ≠ 0 := by exact_mod_cast hj
  have hmap : Finset.Ico a (b + 1 - j) =
      (Finset.Icc (a - j) (b - 2 * j)).map (addRightEmbedding j) := by
    rw [← Finset.Ico_add_one_right_eq_Icc, Finset.map_add_right_Ico]
    congr 1 <;> omega
  have hphase : ∀ m : ℕ,
      T * (F (((m + j + j : ℕ) : ℝ) / N) - F (((m + j : ℕ) : ℝ) / N)) =
        -(σ * T * (j : ℝ) / N * differencePhase F σ ((j : ℝ) / N)
          ((m : ℝ) / (N - j))) := by
    intro m
    have hleft : (1 - (j : ℝ) / N) * ((m : ℝ) / (N - j)) + (j : ℝ) / N =
        ((m + j : ℕ) : ℝ) / N := by
      push_cast
      field_simp
    have hright : (1 - (j : ℝ) / N) * ((m : ℝ) / (N - j)) + 2 * ((j : ℝ) / N) =
        ((m + j + j : ℕ) : ℝ) / N := by
      push_cast
      field_simp
      ring
    unfold differencePhase
    rw [hleft, hright]
    field_simp
  unfold correlationSumAt exponentialSumAt oscillatory
  rw [hmap, Finset.sum_map]
  have hsum : (∑ m ∈ Finset.Icc (a - j) (b - 2 * j),
      (𝐞 (T * (F (((m + j + j : ℕ) : ℝ) / N) -
        F (((m + j : ℕ) : ℝ) / N))) : ℂ)) =
      starRingEnd ℂ (∑ m ∈ Finset.Icc (a - j) (b - 2 * j),
        (𝐞 (σ * T * (j : ℝ) / N * differencePhase F σ ((j : ℝ) / N)
          ((m : ℝ) / (N - j))) : ℂ)) := by
    rw [map_sum]
    apply Finset.sum_congr rfl
    intro m hm
    rw [hphase]
    simp
  simp only [addRightEmbedding_apply, hsum, RCLike.norm_conj]

end Expdb
