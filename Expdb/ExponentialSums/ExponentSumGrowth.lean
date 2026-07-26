module

public import Expdb.ExponentialSums.PhaseFunctions
public import Mathlib.Analysis.Asymptotics.Defs
public import Mathlib.Order.Interval.Finset.Nat
import Mathlib.Analysis.SpecialFunctions.Pow.Asymptotics
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Algebra.Order.Floor.Semiring
import Mathlib.Topology.Order.Monotone

/-!
# Exponential sum growth exponents

This module formalizes Definition 4.2 of the ANTEDB blueprint. It defines admissible
model-phase bounds at a fixed scale and defines `β(α)` as the least such exponent.

## Definition 4.2 at a glance

For `α : ℝ≥0`, an exponent `β` is admissible if every model-phase exponential sum at scale
`N = T ^ (α + o(1))` is bounded by `T ^ (β + o(1))`. We show that the admissible exponents
form `[β(α), ∞)` and then specialize the bound to the logarithmic model phase.
-/

@[expose] public section

open Filter Topology
open scoped Expdb FourierTransform NNReal

noncomputable section

namespace Expdb

/-! ## Definition 4.2 -/

/-- `N = T ^ (α + o(1))`, expressed using a variable exponent that is equal to `α` up to an
infinitesimal. -/
def IsPowerAsymptotic (N T : VariableObject ℝ) (α : ℝ) : Prop :=
  ∃ exponent : VariableObject ℝ,
    exponent =o VariableObject.fixed α ∧
      ∀ᶠ i in atTop, N i = T i ^ exponent i

/-- `X ≪ T ^ (β + o(1))`, expressed using a variable exponent that is equal to `β` up to
an infinitesimal. -/
def IsPowerBounded {E : Type*} [SeminormedAddCommGroup E]
    (X : VariableObject E) (T : VariableObject ℝ) (β : ℝ) : Prop :=
  ∃ exponent : VariableObject ℝ,
    exponent =o VariableObject.fixed β ∧
      Asymptotics.IsBigO atTop X (fun i ↦ T i ^ exponent i)

/-- The variable exponential sum
`∑ n ∈ [a, b], e(T F(n / N))`. -/
def exponentialSum
    (F : VariableFunction (VariableObject.fixed ℝ) ℝ)
    (T N : VariableObject ℝ) (a b : VariableObject ℕ) :
    VariableObject ℂ :=
  fun i ↦
    ∑ n ∈ Finset.Icc (a i) (b i),
      (𝐞 (T i * F i ((n : ℝ) / N i)) : ℂ)

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

/-- The exponential sum growth exponent `β(α)` from Definition 4.2 of the blueprint. -/
def exponentSumGrowthExponent (α : ℝ≥0) : ℝ :=
  sInf {β : ℝ | IsExponentSumBound α β}

/-! ## Power asymptotics -/

/-- If `N = T ^ (α + o(1))`, then `N` eventually lies between the powers with exponents
`α - ε` and `α + ε`. -/
theorem IsPowerAsymptotic.eventually_between
    {N T : VariableObject ℝ} {α ε : ℝ}
    (hNT : IsPowerAsymptotic N T α)
    (hT : ∀ᶠ i in atTop, 1 ≤ T i)
    (hε : 0 < ε) :
    ∀ᶠ i in atTop, T i ^ (α - ε) ≤ N i ∧ N i ≤ T i ^ (α + ε) := by
  obtain ⟨exponent, hexponent, hN⟩ := hNT
  obtain ⟨hexponent_le, halpha_le⟩ :=
    (isEqUpToInfinitesimal_iff_leUpToInfinitesimal
      exponent (VariableObject.fixed α)).1 hexponent
  have hupper := (isLEUpToInfinitesimal_iff_forall_pos
    exponent (VariableObject.fixed α)).1 hexponent_le ε hε
  have hlower := (isLEUpToInfinitesimal_iff_forall_pos
    (VariableObject.fixed α) exponent).1 halpha_le ε hε
  filter_upwards [hT, hN, hupper, hlower] with i hiT hiN hiupper hilower
  rw [hiN]
  exact
    ⟨Real.rpow_le_rpow_of_exponent_le hiT (by linarith),
      Real.rpow_le_rpow_of_exponent_le hiT (by linarith)⟩

/-! ## Power bounds -/

/-- If `T` is at least one and unbounded, then `X ≪ T ^ (β + o(1))` iff
`X = O(T ^ (β + ε))` for every fixed `ε > 0`. -/
theorem isPowerBounded_iff_forall_pos
    {E : Type*} [SeminormedAddCommGroup E]
    (X : VariableObject E) (T : VariableObject ℝ) (β : ℝ)
    (hT : ∀ i, 1 ≤ T i) (hTunbounded : T.IsUnbounded) :
    IsPowerBounded X T β ↔
      ∀ ε : ℝ, 0 < ε →
        Asymptotics.IsBigO atTop X (fun i ↦ T i ^ (β + ε)) := by
  -- Positivity and unboundedness make powers of `T` eventually monotone in the exponent.
  have hTtop : Tendsto T atTop atTop := by
    rw [VariableObject.IsUnbounded] at hTunbounded
    exact hTunbounded.congr' <| Filter.Eventually.of_forall fun i ↦ by
      rw [Real.norm_eq_abs, abs_of_nonneg (le_trans zero_le_one (hT i))]
  have hTgt : ∀ᶠ i in atTop, 1 < T i :=
    (tendsto_atTop.1 hTtop 2).mono fun _ hi ↦ by linarith
  constructor
  -- An infinitesimal exponent error is eventually smaller than each fixed `ε`.
  · rintro ⟨exponent, hexponent, hX⟩ ε hε
    refine hX.trans (Asymptotics.IsBigO.of_bound 1 ?_)
    have hexponent_upper :=
      (isEqUpToInfinitesimal_iff_forall_pos
        exponent (VariableObject.fixed β)).1 hexponent ε hε
    filter_upwards [hexponent_upper] with i hi
    rw [one_mul, Real.norm_eq_abs,
      abs_of_nonneg (Real.rpow_nonneg (le_trans zero_le_one (hT i)) _),
      Real.norm_eq_abs,
      abs_of_nonneg (Real.rpow_nonneg (le_trans zero_le_one (hT i)) _)]
    apply Real.rpow_le_rpow_of_exponent_le (hT i)
    rw [Real.norm_eq_abs] at hi
    linarith [le_abs_self (exponent i - β)]
  -- Conversely, choose the smallest pointwise exponent that bounds `‖X i‖`.
  · intro hX
    let exponent : VariableObject ℝ := fun i ↦
      if ‖X i‖ = 0 then β else
        max β (Real.log ‖X i‖ / Real.log (T i))
    refine ⟨exponent, ?_, Asymptotics.IsBigO.of_bound 1 ?_⟩
    · apply (isEqUpToInfinitesimal_iff_forall_pos
        exponent (VariableObject.fixed β)).2
      intro δ hδ
      have hδ3 : 0 < δ / 3 := by linarith
      obtain ⟨C, hC⟩ := (hX (δ / 3) hδ3).bound
      let D : ℝ := max C 1
      have hDpos : 0 < D := lt_of_lt_of_le zero_lt_one (le_max_right C 1)
      have hC' : ∀ᶠ i in atTop,
          ‖X i‖ ≤ D * T i ^ (β + δ / 3) := hC.mono fun i hi ↦ by
        calc
          ‖X i‖ ≤ C * ‖T i ^ (β + δ / 3)‖ := hi
          _ ≤ D * ‖T i ^ (β + δ / 3)‖ := by
            gcongr
            exact le_max_left C 1
          _ = D * T i ^ (β + δ / 3) := by
            rw [Real.norm_eq_abs,
              abs_of_nonneg (Real.rpow_nonneg
                (le_trans zero_le_one (hT i)) _)]
      have hpowtop : Tendsto (fun i ↦ T i ^ (δ / 3)) atTop atTop :=
        (tendsto_rpow_atTop hδ3).comp hTtop
      have hD : ∀ᶠ i in atTop, D ≤ T i ^ (δ / 3) :=
        tendsto_atTop.1 hpowtop D
      filter_upwards [hTgt, hC', hD] with i hiT hiX hiD
      by_cases hXi : ‖X i‖ = 0
      · simp [exponent, hXi, hδ]
      · have hXpos : 0 < ‖X i‖ := lt_of_le_of_ne (norm_nonneg _) (Ne.symm hXi)
        have hiTpos : 0 < T i := lt_trans zero_lt_one hiT
        have hXpower : ‖X i‖ ≤ T i ^ (β + 2 * δ / 3) := by
          calc
            ‖X i‖ ≤ D * T i ^ (β + δ / 3) := hiX
            _ ≤ T i ^ (δ / 3) * T i ^ (β + δ / 3) := by
              gcongr
            _ = T i ^ (β + 2 * δ / 3) := by
              rw [← Real.rpow_add hiTpos]
              congr 1
              ring
        have hratio :
            Real.log ‖X i‖ / Real.log (T i) ≤ β + 2 * δ / 3 := by
          apply (div_le_iff₀ (Real.log_pos hiT)).2
          exact (Real.le_rpow_iff_log_le hXpos hiTpos).1 hXpower
        rw [Real.norm_eq_abs, abs_lt]
        constructor
        · simp only [exponent, hXi, ↓reduceIte, VariableObject.fixed]
          linarith [le_max_left β
            (Real.log ‖X i‖ / Real.log (T i))]
        · simp only [exponent, hXi, ↓reduceIte, VariableObject.fixed]
          have := max_lt (show β < β + δ by linarith)
            (lt_of_le_of_lt hratio (show β + 2 * δ / 3 < β + δ by linarith))
          linarith
    · filter_upwards [hTgt] with i hiT
      rw [one_mul, Real.norm_eq_abs,
        abs_of_nonneg (Real.rpow_nonneg (le_trans zero_le_one (hT i)) _)]
      by_cases hXi : ‖X i‖ = 0
      · rw [hXi]
        exact Real.rpow_nonneg (le_trans zero_le_one (hT i)) _
      · have hXpos : 0 < ‖X i‖ := lt_of_le_of_ne (norm_nonneg _) (Ne.symm hXi)
        have hiTpos : 0 < T i := lt_trans zero_lt_one hiT
        apply (Real.le_rpow_iff_log_le hXpos hiTpos).2
        apply (div_le_iff₀ (Real.log_pos hiT)).1
        simp only [exponent, hXi, ↓reduceIte]
        exact le_max_right β (Real.log ‖X i‖ / Real.log (T i))

/-! ## Admissible exponents

The candidate set is upward closed, nonempty, and bounded below.

### Monotonicity -/

private lemma isPowerBounded_mono
    {X : VariableObject ℂ} {T : VariableObject ℝ} {β γ : ℝ}
    (hT : ∀ i, 1 ≤ T i) (hβγ : β ≤ γ)
    (hX : IsPowerBounded X T β) :
    IsPowerBounded X T γ := by
  obtain ⟨exponent, hexponent, hX⟩ := hX
  let shiftedExponent : VariableObject ℝ :=
    fun i ↦ exponent i + (γ - β)
  refine ⟨shiftedExponent, ?_, hX.trans (Asymptotics.IsBigO.of_bound 1 ?_)⟩
  · have hshift := hexponent.add
      (IsEqUpToInfinitesimal.refl (VariableObject.fixed (γ - β)))
    convert hshift using 1
    · ext i
      dsimp [shiftedExponent]
    · ext i
      dsimp [shiftedExponent]
      ring
  exact Filter.Eventually.of_forall fun i ↦ by
    have hT0 : 0 ≤ T i := le_trans zero_le_one (hT i)
    rw [one_mul, Real.norm_eq_abs, abs_of_nonneg (Real.rpow_nonneg hT0 _),
      Real.norm_eq_abs, abs_of_nonneg (Real.rpow_nonneg hT0 _)]
    exact Real.rpow_le_rpow_of_exponent_le (hT i) (by
      dsimp [shiftedExponent]
      linarith)

/-- Admissible exponential-sum bounds are monotone in the exponent. -/
theorem IsExponentSumBound.mono {α : ℝ≥0} {β γ : ℝ}
    (hβ : IsExponentSumBound α β) (hβγ : β ≤ γ) :
    IsExponentSumBound α γ := by
  intro N T F a b hN hT hTunbounded hNT hF hab
  exact isPowerBounded_mono hT hβγ
    (hβ N T F a b hN hT hTunbounded hNT hF hab)

/-! ### Nonemptiness and lower bound -/

private lemma norm_exponentialSum_le_three_mul
    (F : VariableFunction (VariableObject.fixed ℝ) ℝ)
    (T N : VariableObject ℝ) (a b : VariableObject ℕ)
    (hN : ∀ i, 1 ≤ N i)
    (hab : ∀ i, N i ≤ (a i : ℝ) ∧ (b i : ℝ) ≤ 2 * N i)
    (i : ℕ) :
    ‖exponentialSum F T N a b i‖ ≤ 3 * N i := by
  calc
    ‖exponentialSum F T N a b i‖ ≤
        ∑ n ∈ Finset.Icc (a i) (b i),
          ‖(𝐞 (T i * F i ((n : ℝ) / N i)) : ℂ)‖ := norm_sum_le _ _
    _ = ((Finset.Icc (a i) (b i)).card : ℝ) := by simp
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
    {β : ℝ | IsExponentSumBound α β}.Nonempty :=
  ⟨α, isExponentSumBound_self α⟩

/-- Every admissible exponent at a nonnegative scale is nonnegative. -/
theorem IsExponentSumBound.nonneg
    {α : ℝ≥0} {β : ℝ} (hβ : IsExponentSumBound α β) :
    0 ≤ β := by
  by_contra hβneg
  have hβlt : β < 0 := lt_of_not_ge hβneg
  let T : VariableObject ℝ := fun i ↦ (i : ℝ) + 2
  let N : VariableObject ℝ := fun i ↦ T i ^ (α : ℝ)
  let a : VariableObject ℕ := fun i ↦ ⌈N i⌉₊
  have hT : ∀ i, 1 ≤ T i := by
    intro i
    dsimp [T]
    have hi : (0 : ℝ) ≤ i := Nat.cast_nonneg i
    linarith
  have hTtendsto : Tendsto T atTop atTop := by
    exact tendsto_atTop_add_const_right atTop 2 tendsto_natCast_atTop_atTop
  have hTunbounded : T.IsUnbounded := by
    rw [VariableObject.IsUnbounded]
    convert hTtendsto using 1
    ext i
    rw [Real.norm_eq_abs, abs_of_nonneg (le_trans zero_le_one (hT i))]
  have hN : ∀ i, 1 ≤ N i := fun i ↦ Real.one_le_rpow (hT i) α.property
  have hNT : IsPowerAsymptotic N T (α : ℝ) := by
    refine ⟨VariableObject.fixed (α : ℝ), IsEqUpToInfinitesimal.refl _, ?_⟩
    exact Filter.Eventually.of_forall fun _ ↦ rfl
  have ha : ∀ i, N i ≤ (a i : ℝ) ∧ (a i : ℝ) ≤ 2 * N i := by
    intro i
    have hN0 : 0 ≤ N i := le_trans zero_le_one (hN i)
    refine ⟨Nat.le_ceil _, ?_⟩
    exact (Nat.ceil_lt_add_one hN0).le.trans (by linarith [hN i])
  -- A singleton log-phase sum has norm one, so it cannot have negative power growth.
  have hbound :=
    hβ N T logPhase a a hN hT hTunbounded hNT
      isModelPhaseFunction_log ha
  have hO := (isPowerBounded_iff_forall_pos
    (exponentialSum logPhase T N a a) T β hT hTunbounded).1
      hbound (-β / 2) (by linarith)
  have hpower_tendsto :
      Tendsto (fun i ↦ T i ^ (β + -β / 2)) atTop (𝓝 0) := by
    have hneg : 0 < -(β + -β / 2) := by linarith
    apply ((tendsto_rpow_neg_atTop hneg).comp hTtendsto).congr'
    exact Filter.Eventually.of_forall fun i ↦ by
      change T i ^ (- -(β + -β / 2)) = T i ^ (β + -β / 2)
      rw [neg_neg]
  have hsum_tendsto :
      Tendsto (exponentialSum logPhase T N a a) atTop (𝓝 0) :=
    hO.trans_tendsto hpower_tendsto
  have hsmall := (Metric.tendsto_nhds.1 hsum_tendsto) (1 / 2) (by norm_num)
  obtain ⟨i, hi⟩ := hsmall.exists
  have hnorm :
      ‖exponentialSum logPhase T N a a i‖ = 1 := by
    simp [exponentialSum]
  rw [dist_zero_right, hnorm] at hi
  norm_num at hi

private theorem exponentSumBounds_bddBelow {α : ℝ≥0} :
    BddBelow {β : ℝ | IsExponentSumBound α β} :=
  ⟨0, fun _ hβ ↦ hβ.nonneg⟩

/-! ## The least admissible exponent -/

private theorem exists_isExponentSumBound_lt_add
    (α : ℝ≥0) {ε : ℝ} (hε : 0 < ε) :
    ∃ β : ℝ, IsExponentSumBound α β ∧
      β < exponentSumGrowthExponent α + ε := by
  let S : Set ℝ := {β : ℝ | IsExponentSumBound α β}
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
    IsLeast {β : ℝ | IsExponentSumBound α β}
      (exponentSumGrowthExponent α) :=
  ⟨isExponentSumBound_exponentSumGrowthExponent α,
    fun _ hβ ↦ exponentSumGrowthExponent_le_iff.mpr hβ⟩

/-- At scale `α`, the admissible exponents form the interval `[β(α), ∞)`. -/
theorem exponentSumBounds_eq_Ici (α : ℝ≥0) :
    {β : ℝ | IsExponentSumBound α β} =
      Set.Ici (exponentSumGrowthExponent α) := by
  ext β
  exact exponentSumGrowthExponent_le_iff.symm

/-- The set of admissible exponential-sum exponents at scale `α` is closed. -/
theorem isClosed_exponentSumBounds (α : ℝ≥0) :
    IsClosed {β : ℝ | IsExponentSumBound α β} := by
  rw [exponentSumBounds_eq_Ici]
  exact isClosed_Ici

/-! ## Logarithmic model phase -/

/-- The least admissible exponent bounds the logarithmic model-phase sum. -/
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

end Expdb
