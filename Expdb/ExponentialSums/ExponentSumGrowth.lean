module

public import Expdb.Basic.PowerAsymptotics
public import Expdb.ExponentialSums.LogPhase
public import Expdb.ExponentialSums.Oscillatory
public import Mathlib.Order.Interval.Finset.Nat

/-!
# Exponential sum growth exponents

This module formalizes the asymptotic definition of the exponential sum growth exponent from the
blueprint's Exponential sum growth exponents chapter (`beta-chapter`). It defines admissible
model-phase bounds at a fixed scale and defines `β(α)` as the least such exponent.

## Overview

For `α : ℝ≥0`, an exponent `β` is admissible if every model-phase exponential sum at scale
`N = T ^ (α + o(1))` is bounded by `T ^ (β + o(1))`. We show that the admissible exponents
form `[β(α), ∞)` and then specialize the bound to the logarithmic model phase.
-/

@[expose] public section

open Filter Topology
open scoped Expdb FourierTransform NNReal

noncomputable section

namespace Expdb

/-! ## Exponential-sum growth definition -/

/-- The variable exponential sum
`∑ n ∈ [a, b], e(T F(n / N))`. -/
def exponentialSum
    (F : VariableFunction (VariableObject.fixed ℝ) ℝ)
    (T N : VariableObject ℝ) (a b : VariableObject ℕ) :
    VariableObject ℂ :=
  fun i ↦ exponentialSumAt (F i) (T i) (N i) (a i) (b i)

@[simp] theorem exponentialSum_apply
    (F : VariableFunction (VariableObject.fixed ℝ) ℝ)
    (T N : VariableObject ℝ) (a b : VariableObject ℕ) (i : ℕ) :
    exponentialSum F T N a b i = exponentialSumAt (F i) (T i) (N i) (a i) (b i) :=
  rfl

/-- The assertion that `β` is an admissible exponential-sum growth exponent at scale `α`. -/
def IsExponentSumBound (α : ℝ≥0) (β : ℝ) : Prop :=
  ∀ (N T : VariableObject ℝ)
    (F : VariableFunction (VariableObject.fixed ℝ) ℝ)
    (a b : VariableObject ℕ),
    (∀ i, 1 ≤ N i) →
    (∀ i, 1 ≤ T i) →
    T.IsUnbounded →
    IsPowerAsymptotic N T (α : ℝ) →
    IsModelPhaseFunction F →
    (∀ i, N i ≤ (a i : ℝ) ∧ (b i : ℝ) ≤ 2 * N i) →
    IsPowerBounded (exponentialSum F T N a b) T β

/-- The set of admissible exponential-sum growth exponents at scale `α`. -/
def exponentSumBounds (α : ℝ≥0) : Set ℝ :=
  {β : ℝ | IsExponentSumBound α β}

@[simp] theorem mem_exponentSumBounds {α : ℝ≥0} {β : ℝ} :
    β ∈ exponentSumBounds α ↔ IsExponentSumBound α β :=
  Iff.rfl

/-- The exponential sum growth exponent `β(α)` from the blueprint definition `beta-def`. -/
def exponentSumGrowthExponent (α : ℝ≥0) : ℝ :=
  sInf (exponentSumBounds α)

/-! ## Admissible exponents

The candidate set is upward closed, nonempty, and bounded below.

### Monotonicity -/

/-- Admissible exponential-sum bounds are monotone in the exponent. -/
theorem IsExponentSumBound.mono {α : ℝ≥0} {β γ : ℝ}
    (hβ : IsExponentSumBound α β) (hβγ : β ≤ γ) :
    IsExponentSumBound α γ := by
  intro N T F a b hN hT hTunbounded hNT hF hab
  exact (hβ N T F a b hN hT hTunbounded hNT hF hab).mono
    (Filter.Eventually.of_forall hT) hβγ

/-! ### Nonemptiness and lower bound -/

private lemma norm_exponentialSum_le_three_mul
    (F : VariableFunction (VariableObject.fixed ℝ) ℝ)
    (T N : VariableObject ℝ) (a b : VariableObject ℕ)
    (hN : ∀ i, 1 ≤ N i)
    (hab : ∀ i, N i ≤ (a i : ℝ) ∧ (b i : ℝ) ≤ 2 * N i)
    (i : ℕ) :
    ‖exponentialSum F T N a b i‖ ≤ 3 * N i := by
  calc
    ‖exponentialSum F T N a b i‖ ≤ ((Finset.Icc (a i) (b i)).card : ℝ) :=
      norm_exponentialSumAt_le_card _ _ _ _ _
    _ ≤ (b i : ℝ) + 1 := by
      rw [Nat.card_Icc]
      exact_mod_cast Nat.sub_le _ _
    _ ≤ 3 * N i := by linarith [hab i |>.2, hN i]

/-- The triangle inequality gives the admissible exponent `α`. -/
theorem isExponentSumBound_self (α : ℝ≥0) :
    IsExponentSumBound α (α : ℝ) := by
  intro N T F a b hN hT _ hNT _ hab
  obtain ⟨exponent, hexponent, hNpower⟩ := hNT
  refine ⟨exponent, hexponent, Asymptotics.IsBigO.of_bound 3 ?_⟩
  filter_upwards [hNpower] with i hi
  calc
    ‖exponentialSum F T N a b i‖ ≤ 3 * N i :=
      norm_exponentialSum_le_three_mul F T N a b hN hab i
    _ = 3 * T i ^ exponent i := by rw [hi]
    _ = 3 * ‖T i ^ exponent i‖ := by
      rw [Real.norm_eq_abs,
        abs_of_nonneg (Real.rpow_nonneg (le_trans zero_le_one (hT i)) _)]

private lemma exponentSumBounds_nonempty (α : ℝ≥0) :
    (exponentSumBounds α).Nonempty :=
  ⟨α, isExponentSumBound_self α⟩

/-- Every admissible exponent at a nonnegative scale is nonnegative. -/
theorem IsExponentSumBound.nonneg
    {α : ℝ≥0} {β : ℝ} (hβ : IsExponentSumBound α β) :
    0 ≤ β := by
  let T : VariableObject ℝ := fun i ↦ (i : ℝ) + 2
  let N : VariableObject ℝ := fun i ↦ T i ^ (α : ℝ)
  let a : VariableObject ℕ := fun i ↦ ⌈N i⌉₊
  have hT : ∀ i, 1 ≤ T i := by
    intro i
    dsimp [T]
    have hi : (0 : ℝ) ≤ i := Nat.cast_nonneg i
    linarith
  have hTunbounded : T.IsUnbounded := by
    apply (VariableObject.isUnbounded_iff_tendsto_atTop
      fun i ↦ zero_le_one.trans (hT i)).2
    exact tendsto_atTop_add_const_right atTop 2 tendsto_natCast_atTop_atTop
  have hN : ∀ i, 1 ≤ N i := fun i ↦ Real.one_le_rpow (hT i) α.property
  have hNT : IsPowerAsymptotic N T (α : ℝ) := by
    refine ⟨VariableObject.fixed (α : ℝ), IsEqUpToInfinitesimal.refl _, ?_⟩
    exact Filter.Eventually.of_forall fun _ ↦ rfl
  have ha : ∀ i, N i ≤ (a i : ℝ) ∧ (a i : ℝ) ≤ 2 * N i := by
    intro i
    have hN0 : 0 ≤ N i := zero_le_one.trans (hN i)
    refine ⟨Nat.le_ceil _, ?_⟩
    exact (Nat.ceil_lt_add_one hN0).le.trans (by linarith [hN i])
  have hbound :=
    hβ N T logPhase a a hN hT hTunbounded hNT
      isModelPhaseFunction_log ha
  apply exponent_le_of_isPowerBounded_of_eventually_norm_ge_rpow
    hT hTunbounded hbound zero_lt_one
  exact Filter.Eventually.of_forall fun i ↦ by simp

private lemma exponentSumBounds_bddBelow {α : ℝ≥0} :
    BddBelow (exponentSumBounds α) :=
  ⟨0, fun _ hβ ↦ hβ.nonneg⟩

/-! ## The least admissible exponent -/

private lemma exists_isExponentSumBound_lt_add
    (α : ℝ≥0) {ε : ℝ} (hε : 0 < ε) :
    ∃ β : ℝ, IsExponentSumBound α β ∧
      β < exponentSumGrowthExponent α + ε := by
  let S : Set ℝ := exponentSumBounds α
  have hS : S.Nonempty := exponentSumBounds_nonempty α
  have hinf_lt : sInf S < sInf S + ε := by linarith
  obtain ⟨β, hβ, hβlt⟩ := exists_lt_of_csInf_lt hS hinf_lt
  change IsExponentSumBound α β at hβ
  refine ⟨β, hβ, ?_⟩
  simpa [S, exponentSumGrowthExponent] using hβlt

/-- By underspill, the infimum of the admissible exponents is itself admissible. -/
theorem isExponentSumBound_exponentSumGrowthExponent (α : ℝ≥0) :
    IsExponentSumBound α (exponentSumGrowthExponent α) := by
  intro N T F a b hN hT hTunbounded hNT hF hab
  apply (isPowerBounded_iff_forall_pos
    (exponentialSum F T N a b) T (exponentSumGrowthExponent α)
    hT hTunbounded).2
  intro ε hε
  -- Choose an admissible exponent within `ε / 2` of the infimum.
  obtain ⟨β, hβ, hβlt⟩ :=
    exists_isExponentSumBound_lt_add α (show 0 < ε / 2 by linarith)
  have hβpower := hβ N T F a b hN hT hTunbounded hNT hF hab
  have hβbound := (isPowerBounded_iff_forall_pos
    (exponentialSum F T N a b) T β hT hTunbounded).1
      hβpower (ε / 2) (by linarith)
  refine hβbound.trans (Asymptotics.IsBigO.of_bound 1 ?_)
  exact Filter.Eventually.of_forall fun i ↦ by
    have hT0 : 0 ≤ T i := le_trans zero_le_one (hT i)
    rw [one_mul, Real.norm_eq_abs, abs_of_nonneg (Real.rpow_nonneg hT0 _),
      Real.norm_eq_abs, abs_of_nonneg (Real.rpow_nonneg hT0 _)]
    apply Real.rpow_le_rpow_of_exponent_le (hT i)
    linarith

/-- The defining universal property of the exponential sum growth exponent. -/
theorem exponentSumGrowthExponent_le_iff {α : ℝ≥0} {β : ℝ} :
    exponentSumGrowthExponent α ≤ β ↔ IsExponentSumBound α β := by
  constructor
  · exact fun h ↦ (isExponentSumBound_exponentSumGrowthExponent α).mono h
  · intro hβ
    exact csInf_le exponentSumBounds_bddBelow hβ

/-- The exponent sum growth exponent is the least admissible exponent. -/
theorem isLeast_exponentSumGrowthExponent (α : ℝ≥0) :
    IsLeast (exponentSumBounds α)
      (exponentSumGrowthExponent α) :=
  ⟨isExponentSumBound_exponentSumGrowthExponent α,
    fun _ hβ ↦ exponentSumGrowthExponent_le_iff.mpr hβ⟩

/-- At scale `α`, the admissible exponents form the interval `[β(α), ∞)`. -/
theorem exponentSumBounds_eq_Ici (α : ℝ≥0) :
    exponentSumBounds α = Set.Ici (exponentSumGrowthExponent α) := by
  ext β
  exact exponentSumGrowthExponent_le_iff.symm

/-- The set of admissible exponential-sum exponents at scale `α` is closed. -/
theorem isClosed_exponentSumBounds (α : ℝ≥0) :
    IsClosed (exponentSumBounds α) := by
  rw [exponentSumBounds_eq_Ici]
  exact isClosed_Ici

/-! ## Logarithmic model phase -/

/-- Logarithmic-phase sums at scale `α` are power-bounded by the least admissible exponent.
This specializes `isExponentSumBound_exponentSumGrowthExponent` to `logPhase`. -/
theorem isPowerBounded_logPhase
    (α : ℝ≥0)
    {N T : VariableObject ℝ} {a b : VariableObject ℕ}
    (hN : ∀ i, 1 ≤ N i)
    (hT : ∀ i, 1 ≤ T i)
    (hTunbounded : T.IsUnbounded)
    (hNT : IsPowerAsymptotic N T (α : ℝ))
    (hab : ∀ i, N i ≤ (a i : ℝ) ∧ (b i : ℝ) ≤ 2 * N i) :
    IsPowerBounded (exponentialSum logPhase T N a b) T
      (exponentSumGrowthExponent α) :=
  isExponentSumBound_exponentSumGrowthExponent α
    N T logPhase a b hN hT hTunbounded hNT isModelPhaseFunction_log hab

/-- A lower bound of size `T ^ γ` for logarithmic-phase sums along an admissible family of
scales forces `γ ≤ β(α)`. -/
theorem le_exponentSumGrowthExponent_of_logPhase_lower_bound
    (α : ℝ≥0)
    {N T : VariableObject ℝ} {a b : VariableObject ℕ} {c γ : ℝ}
    (hN : ∀ i, 1 ≤ N i)
    (hT : ∀ i, 1 ≤ T i)
    (hTunbounded : T.IsUnbounded)
    (hNT : IsPowerAsymptotic N T (α : ℝ))
    (hab : ∀ i, N i ≤ (a i : ℝ) ∧ (b i : ℝ) ≤ 2 * N i)
    (hc : 0 < c)
    (hlower : ∀ᶠ i in atTop, c * T i ^ γ ≤ ‖exponentialSum logPhase T N a b i‖) :
    γ ≤ exponentSumGrowthExponent α :=
  exponent_le_of_isPowerBounded_of_eventually_norm_ge_rpow hT hTunbounded
    (isPowerBounded_logPhase α hN hT hTunbounded hNT hab) hc hlower

end Expdb
