module

public import Expdb.Basic.PowerAsymptotics
public import Expdb.ExponentialSums.LogPhase
public import Expdb.ExponentialSums.FixedExponentialSum
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

/-- Evaluate a variable exponential sum at an ambient index. -/
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
    (∀ i, IsDyadicSubinterval (N i) (a i) (b i)) →
    IsPowerBounded (exponentialSum F T N a b) T β

/-- An exponential-sum bound applies when its setup conditions hold only eventually. -/
theorem IsExponentSumBound.isPowerBounded_of_eventually_setup
    {α : ℝ≥0} {β σ : ℝ} (hβ : IsExponentSumBound α β) (hσ : 0 < σ)
    {N T : VariableObject ℝ} {F : VariableFunction (VariableObject.fixed ℝ) ℝ}
    {a b : VariableObject ℕ} (hTtop : Tendsto T atTop atTop)
    (hNT : IsPowerAsymptotic N T α) (hF : IsModelPhaseFunctionWith F σ)
    (hsetup : ∀ᶠ i in atTop, 1 ≤ N i ∧ 1 ≤ T i ∧
      IsDyadicSubinterval (N i) (a i) (b i)) :
    IsPowerBounded (exponentialSum F T N a b) T β := by
  obtain ⟨i₀, hi₀⟩ := eventually_atTop.1 hsetup
  let φ : ℕ → ℕ := fun i ↦ i + i₀
  have hφ : Tendsto φ atTop atTop := by
    apply tendsto_atTop.2
    intro k
    filter_upwards [eventually_ge_atTop k] with i hi
    dsimp [φ]
    omega
  have hT : ∀ i, 1 ≤ T (φ i) := fun i ↦ (hi₀ _ (by dsimp [φ]; omega)).2.1
  have hTunbounded : VariableObject.IsUnbounded (T ∘ φ) :=
    (VariableObject.isUnbounded_iff_tendsto_atTop
      (fun i ↦ zero_le_one.trans (hT i))).2 (hTtop.comp hφ)
  have htail := hβ (N ∘ φ) (T ∘ φ) (fun i ↦ F (φ i)) (a ∘ φ) (b ∘ φ)
    (fun i ↦ (hi₀ _ (by dsimp [φ]; omega)).1) hT hTunbounded
    (hNT.comp_tendsto hφ) ((hF.comp_tendsto hφ).isModelPhaseFunction hσ)
    (fun i ↦ (hi₀ _ (by dsimp [φ]; omega)).2.2)
  let ψ : ℕ → ℕ := fun i ↦ i - i₀
  have hψ : Tendsto ψ atTop atTop := by
    apply tendsto_atTop.2
    intro k
    filter_upwards [eventually_ge_atTop (k + i₀)] with i hi
    dsimp [ψ]
    omega
  obtain ⟨e, he, hbound⟩ := htail.comp_tendsto hψ
  refine ⟨e, he, hbound.congr' ?_ ?_⟩
  · filter_upwards [eventually_ge_atTop i₀] with i hi
    have hidx : φ (ψ i) = i := Nat.sub_add_cancel hi
    change exponentialSumAt (F (φ (ψ i))) (T (φ (ψ i))) (N (φ (ψ i)))
      (a (φ (ψ i))) (b (φ (ψ i))) = exponentialSumAt (F i) (T i) (N i) (a i) (b i)
    rw [hidx]
  · filter_upwards [eventually_ge_atTop i₀] with i hi
    have hidx : φ (ψ i) = i := Nat.sub_add_cancel hi
    change T (φ (ψ i)) ^ e i = T i ^ e i
    rw [hidx]

/-- The exponential sum bound obtained by testing only the full interval `[N, 2N]`. -/
def IsFullIntervalExponentSumBound (α : ℝ≥0) (β : ℝ) : Prop :=
  ∀ (N T : VariableObject ℝ)
    (F : VariableFunction (VariableObject.fixed ℝ) ℝ),
    (∀ i, 1 ≤ N i) →
    (∀ i, 1 ≤ T i) →
    T.IsUnbounded →
    IsPowerAsymptotic N T (α : ℝ) →
    IsModelPhaseFunction F →
    IsPowerBounded
      (exponentialSum F T N (fun i ↦ ⌈N i⌉₊) (fun i ↦ ⌊2 * N i⌋₊)) T β

/-- Every arbitrary-interval bound applies to the full dyadic interval. -/
theorem IsExponentSumBound.fullInterval {α : ℝ≥0} {β : ℝ}
    (h : IsExponentSumBound α β) : IsFullIntervalExponentSumBound α β := by
  intro N T F hN hT hTunbounded hNT hF
  apply h N T F _ _ hN hT hTunbounded hNT hF
  intro i
  exact ⟨Nat.le_ceil _, Nat.floor_le (by linarith [hN i])⟩

/-- The set of admissible exponential-sum growth exponents at scale `α`. -/
def exponentSumBounds (α : ℝ≥0) : Set ℝ :=
  {β : ℝ | IsExponentSumBound α β}

/-- Membership in the set of admissible exponential sum exponents. -/
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

/-- Bounds at every strictly larger exponent imply a bound at the limiting exponent. -/
theorem isExponentSumBound_of_forall_gt {α : ℝ≥0} {β : ℝ}
    (h : ∀ γ : ℝ, β < γ → IsExponentSumBound α γ) : IsExponentSumBound α β := by
  intro N T F a b hN hT hTunbounded hNT hF hab
  exact IsPowerBounded.of_forall_gt hT hTunbounded
    (fun γ hγ ↦ h γ hγ N T F a b hN hT hTunbounded hNT hF hab)

/-! ### Nonemptiness and lower bound -/

private theorem norm_exponentialSum_le_three_mul
    (F : VariableFunction (VariableObject.fixed ℝ) ℝ)
    (T N : VariableObject ℝ) (a b : VariableObject ℕ)
    (hN : ∀ i, 1 ≤ N i)
    (hab : ∀ i, IsDyadicSubinterval (N i) (a i) (b i))
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

private theorem exponentSumBounds_nonempty (α : ℝ≥0) :
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

private theorem exponentSumBounds_bddBelow {α : ℝ≥0} :
    BddBelow (exponentSumBounds α) :=
  ⟨0, fun _ hβ ↦ hβ.nonneg⟩

/-! ## The least admissible exponent -/

/-- By underspill, the infimum of the admissible exponents is itself admissible. -/
theorem isExponentSumBound_exponentSumGrowthExponent (α : ℝ≥0) :
    IsExponentSumBound α (exponentSumGrowthExponent α) := by
  apply isExponentSumBound_of_forall_gt
  intro γ hγ
  obtain ⟨β, hβ, hβγ⟩ := exists_lt_of_csInf_lt (exponentSumBounds_nonempty α) hγ
  exact IsExponentSumBound.mono hβ hβγ.le

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
    (hab : ∀ i, IsDyadicSubinterval (N i) (a i) (b i)) :
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
    (hab : ∀ i, IsDyadicSubinterval (N i) (a i) (b i))
    (hc : 0 < c)
    (hlower : ∀ᶠ i in atTop, c * T i ^ γ ≤ ‖exponentialSum logPhase T N a b i‖) :
    γ ≤ exponentSumGrowthExponent α :=
  exponent_le_of_isPowerBounded_of_eventually_norm_ge_rpow hT hTunbounded
    (isPowerBounded_logPhase α hN hT hTunbounded hNT hab) hc hlower

end Expdb
