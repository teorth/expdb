module

public import Expdb.ExponentialSums.ExponentPairs

/-!
# Non-asymptotic exponent-pair bounds

This module formalizes the non-asymptotic characterization of exponent pairs from the blueprint
lemma following `exp-pair-def`. The reference exponent `σ > 0` and approximation tolerance
`δ > 0` are quantified explicitly: the constants may depend on `ε` and `σ`, while the estimate
is uniform over all fixed scales, intervals, and phases satisfying the finite-order conditions.
-/

@[expose] public section

open Filter Topology
open scoped Expdb

noncomputable section

namespace Expdb

/-- Fixed-parameter exponent-pair estimates for phases close to a model through finitely many
orders. The constants and derivative cutoff may depend on the exponent loss and reference phase. -/
def IsExponentPairBoundNonAsymptotic (k l : ℝ) : Prop :=
  ∀ ε : ℝ, 0 < ε → ∀ σ : ℝ, 0 < σ →
    ∃ δ : ℝ, 0 < δ ∧ ∃ P : ℕ, 1 ≤ P ∧ ∃ C : ℝ, 1 ≤ C ∧
      ∀ (T N : ℝ) (F : ℝ → ℝ) (a b : ℕ),
        1 ≤ N → N ≤ T → IsApproximateModelPhaseFunction F σ P δ →
        IsDyadicSubinterval N a b →
        ‖exponentialSumAt F T N a b‖ ≤ C * ((T / N) ^ (k + ε) * N ^ (l + ε))

/-- The analytic exponent-pair bound is equivalent to uniform finite-order phase estimates. -/
theorem isExponentPairBound_iff_nonAsymptotic {k l : ℝ} :
    IsExponentPairBound k l ↔ IsExponentPairBoundNonAsymptotic k l := by
  constructor
  · intro h ε hε σ hσ
    by_contra hfailure
    push Not at hfailure
    -- Improving the accuracy and derivative cutoff turns counterexamples into a model family.
    let δ : VariableObject ℝ := fun i ↦ 1 / ((i : ℝ) + 1)
    let P : VariableObject ℕ := fun i ↦ i + 1
    let C : VariableObject ℝ := fun i ↦ (i : ℝ) + 1
    have hδpos : ∀ i, 0 < δ i := fun i ↦ by dsimp [δ]; positivity
    have hPone : ∀ i, 1 ≤ P i := fun i ↦ by dsimp [P]; omega
    have hCone : ∀ i, 1 ≤ C i := fun i ↦ le_add_of_nonneg_left (Nat.cast_nonneg i)
    choose T N F a b hN hNT happrox hab hbad using fun i ↦
      hfailure (δ i) (hδpos i) (P i) (hPone i) (C i) (hCone i)
    have hδinfinitesimal : δ.IsInfinitesimal := by
      rw [VariableObject.IsInfinitesimal]
      have hlim : Tendsto (fun i : ℕ ↦ 1 / ((i : ℝ) + 1)) atTop (𝓝 0) :=
        tendsto_one_div_add_atTop_nhds_zero_nat
      convert hlim using 1
      ext i
      rw [Real.norm_eq_abs, abs_of_pos (hδpos i)]
    have hPtop : ∀ p : ℕ, ∀ᶠ i in atTop, p ≤ P i := by
      intro p
      filter_upwards [eventually_ge_atTop p] with i hi
      dsimp [P]
      omega
    have hF : IsModelPhaseFunction F :=
      isModelPhaseFunction_of_approximations hσ hPtop hδinfinitesimal happrox
    have hO := isExponentPairBound_iff_forall_pos.mp h N T F a b hN hNT hF hab ε hε
    obtain ⟨K, hKpos, hK⟩ := hO.exists_pos
    have hCtop : Tendsto C atTop atTop :=
      tendsto_atTop_mono (fun i ↦ by dsimp [C]; linarith) tendsto_natCast_atTop_atTop
    have hCK : ∀ᶠ i in atTop, K ≤ C i := hCtop.eventually (eventually_ge_atTop K)
    obtain ⟨i, hi, hiC⟩ := (hK.bound.and hCK).exists
    have hNpos : 0 < N i := zero_lt_one.trans_le (hN i)
    have hTpos : 0 < T i := hNpos.trans_le (hNT i)
    have hmajor : 0 < (T i / N i) ^ (k + ε) * N i ^ (l + ε) :=
      mul_pos (Real.rpow_pos_of_pos (div_pos hTpos hNpos) _)
        (Real.rpow_pos_of_pos hNpos _)
    simp only [exponentialSum_apply, Real.norm_eq_abs, abs_of_pos hmajor] at hi
    exact (not_lt_of_ge (hi.trans (mul_le_mul_of_nonneg_right hiC hmajor.le))) (hbad i)
  · intro h
    apply isExponentPairBound_iff_forall_pos.mpr
    intro N T F a b hN hNT hF hab ε hε
    obtain ⟨σ, hσ, hFσ⟩ := hF.exists_with
    obtain ⟨δ, hδ, P, hP, C, hC, hfixed⟩ := h ε hε σ hσ
    refine Asymptotics.IsBigO.of_bound C ?_
    filter_upwards [hFσ.eventually_isApproximate P hδ] with i hi
    have hNpos : 0 < N i := zero_lt_one.trans_le (hN i)
    have hTpos : 0 < T i := hNpos.trans_le (hNT i)
    simpa only [exponentialSum_apply, Real.norm_eq_abs,
      abs_of_pos (mul_pos (Real.rpow_pos_of_pos (div_pos hTpos hNpos) _)
        (Real.rpow_pos_of_pos hNpos _))] using
      hfixed (T i) (N i) (F i) (a i) (b i) (hN i) (hNT i) hi (hab i)

/-- The blueprint's non-asymptotic characterization of exponent pairs (Lemma 5.2),
for fixed coefficients in the exponent-pair triangle. -/
theorem isExponentPair_iff_nonAsymptotic {k l : ℝ}
    (htriangle : 0 ≤ k ∧ k ≤ 1 / 2 ∧ 1 / 2 ≤ l ∧ l ≤ 1 ∧ k + l ≤ 1) :
    IsExponentPair k l ↔ IsExponentPairBoundNonAsymptotic k l := by
  change (0 ≤ k ∧ k ≤ 1 / 2 ∧ 1 / 2 ≤ l ∧ l ≤ 1 ∧ k + l ≤ 1) ∧
    IsExponentPairBound k l ↔ _
  rw [and_iff_right htriangle, isExponentPairBound_iff_nonAsymptotic]

end Expdb
