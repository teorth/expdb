module

public import Expdb.ExponentialSums.LogPhase
public import Mathlib.Analysis.Asymptotics.Defs
public import Mathlib.Analysis.SpecialFunctions.Log.Base
public import Mathlib.Order.Interval.Finset.Nat

import Mathlib.Analysis.SpecialFunctions.Pow.Asymptotics
import Mathlib.Analysis.Asymptotics.Lemmas

/-!
# Exponential sum growth exponents

This module formalizes the asymptotic definition of the exponential sum growth exponent
from Chapter 4 of the ANTEDB blueprint. It defines admissible
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

/-! ## Asymptotic definition -/

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

/-- Convergence of the logarithmic exponent gives a power asymptotic. -/
theorem isPowerAsymptotic_of_logb_tendsto
    {N T : VariableObject ℝ} {α : ℝ}
    (hT : ∀ᶠ i in atTop, 1 < T i) (hN : ∀ᶠ i in atTop, 0 < N i)
    (hexponent : Tendsto (fun i ↦ Real.logb (T i) (N i)) atTop (nhds α)) :
    IsPowerAsymptotic N T α := by
  let exponent : VariableObject ℝ := fun i ↦ Real.logb (T i) (N i)
  refine ⟨exponent, ?_, ?_⟩
  · rw [IsEqUpToInfinitesimal, VariableObject.IsInfinitesimal]
    have hsub : Tendsto (fun i ↦ Real.logb (T i) (N i) - α) atTop (nhds 0) := by
      simpa using hexponent.sub
        (tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ α) atTop (nhds α))
    convert (continuous_norm.tendsto 0).comp hsub using 1 <;>
      simp [Function.comp_def, exponent, VariableObject.fixed]
  · filter_upwards [hT, hN] with i hiT hiN
    exact (Real.rpow_logb (zero_lt_one.trans hiT) hiT.ne' hiN).symm

/-- A variable power sandwich with an infinitesimal error determines a power asymptotic. -/
theorem isPowerAsymptotic_of_between
    {N T δ : VariableObject ℝ} {α : ℝ}
    (hT : ∀ i, 1 < T i) (hN : ∀ i, 0 < N i)
    (hδnonneg : ∀ i, 0 ≤ δ i) (hδ : δ.IsInfinitesimal)
    (hbetween : ∀ i,
      T i ^ (α - δ i) ≤ N i ∧ N i ≤ T i ^ (α + δ i)) :
    IsPowerAsymptotic N T α := by
  let exponent : VariableObject ℝ := fun i ↦ Real.logb (T i) (N i)
  refine ⟨exponent, ?_, Filter.Eventually.of_forall fun i ↦ ?_⟩
  · apply (isEqUpToInfinitesimal_iff_forall_pos
    exponent (VariableObject.fixed α)).2
    intro ε hε
    have hδsmall :=
      (VariableObject.isInfinitesimal_iff_forall_pos δ).1 hδ ε hε
    filter_upwards [hδsmall] with i hi
    have hlower : α - δ i ≤ exponent i :=
      (Real.le_logb_iff_rpow_le (hT i) (hN i)).2 (hbetween i).1
    have hupper : exponent i ≤ α + δ i :=
      (Real.logb_le_iff_le_rpow (hT i) (hN i)).2 (hbetween i).2
    rw [Real.norm_eq_abs, abs_lt]
    rw [Real.norm_eq_abs, abs_of_nonneg (hδnonneg i)] at hi
    constructor <;> dsimp [exponent, VariableObject.fixed] at * <;> linarith
  · exact (Real.rpow_logb (lt_trans zero_lt_one (hT i))
      (hT i).ne' (hN i)).symm

/-- If `T` lies between fixed positive multiples of `N ^ α⁻¹`, then the logarithm of `N` to
base `T` tends to `α`. -/
private lemma tendsto_logb_of_between_const_rpow
    {N T : VariableObject ℝ} {α A B : ℝ}
    (hα : 0 < α) (hA : 0 < A) (hB : 0 < B)
    (hNtop : Tendsto N atTop atTop)
    (hT : ∀ᶠ i in atTop, A * N i ^ α⁻¹ ≤ T i ∧ T i ≤ B * N i ^ α⁻¹) :
    Tendsto (fun i ↦ Real.logb (T i) (N i)) atTop (nhds α) := by
  have hlogNtop : Tendsto (fun i ↦ Real.log (N i)) atTop atTop :=
    Real.tendsto_log_atTop.comp hNtop
  have hNgt : ∀ᶠ i in atTop, 1 < N i := hNtop.eventually (eventually_gt_atTop 1)
  have hlogNpos : ∀ᶠ i in atTop, 0 < Real.log (N i) := by
    filter_upwards [hNgt] with i hi
    exact Real.log_pos hi
  have hpowtop : Tendsto (fun i ↦ N i ^ α⁻¹) atTop atTop :=
    (tendsto_rpow_atTop (inv_pos.mpr hα)).comp hNtop
  have hpowlarge : ∀ᶠ i in atTop, 1 / A < N i ^ α⁻¹ :=
    hpowtop.eventually (eventually_gt_atTop (1 / A))
  have hTgt : ∀ᶠ i in atTop, 1 < T i := by
    filter_upwards [hpowlarge, hT] with i hi hTi
    have : 1 < A * N i ^ α⁻¹ := by
      simpa [mul_comm] using (div_lt_iff₀ hA).1 hi
    exact this.trans_le hTi.1
  have hratio : Tendsto (fun i ↦ Real.log (T i) / Real.log (N i)) atTop (nhds α⁻¹) := by
    have hlowerLimit : Tendsto
        (fun i ↦ Real.log A / Real.log (N i) + α⁻¹) atTop (nhds α⁻¹) := by
      simpa using (tendsto_const_nhds.div_atTop hlogNtop).add
        (tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ α⁻¹) atTop (nhds α⁻¹))
    have hupperLimit : Tendsto
        (fun i ↦ Real.log B / Real.log (N i) + α⁻¹) atTop (nhds α⁻¹) := by
      simpa using (tendsto_const_nhds.div_atTop hlogNtop).add
        (tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ α⁻¹) atTop (nhds α⁻¹))
    apply tendsto_of_tendsto_of_tendsto_of_le_of_le' hlowerLimit hupperLimit
    · filter_upwards [hNgt, hlogNpos, hT] with i hNi hlogNi hTi
      have hbase : 0 < N i := zero_lt_one.trans hNi
      have hlower := Real.log_le_log (mul_pos hA (Real.rpow_pos_of_pos hbase _)) hTi.1
      rw [Real.log_mul hA.ne' (Real.rpow_pos_of_pos hbase _).ne',
          Real.log_rpow hbase] at hlower
      apply (le_div_iff₀ hlogNi).2
      rw [add_mul, div_mul_cancel₀ _ hlogNi.ne']
      simpa [mul_comm] using hlower
    · filter_upwards [hNgt, hlogNpos, hT] with i hNi hlogNi hTi
      have hbase : 0 < N i := zero_lt_one.trans hNi
      have hupper := Real.log_le_log (lt_of_lt_of_le
        (mul_pos hA (Real.rpow_pos_of_pos hbase _)) hTi.1) hTi.2
      rw [Real.log_mul hB.ne' (Real.rpow_pos_of_pos hbase _).ne',
          Real.log_rpow hbase] at hupper
      apply (div_le_iff₀ hlogNi).2
      rw [add_mul, div_mul_cancel₀ _ hlogNi.ne']
      simpa [mul_comm] using hupper
  have hinv := hratio.inv₀ (inv_ne_zero hα.ne')
  have hinv' : Tendsto (fun i ↦ Real.logb (T i) (N i)) atTop (nhds α⁻¹⁻¹) :=
    hinv.congr' <| by
      filter_upwards [hlogNpos, hTgt] with i hlogNi hTi
      dsimp [Real.logb]
      field_simp [hlogNi.ne', (Real.log_pos hTi).ne']
  simpa using hinv'

/-- Fixed positive multiplicative uncertainty in `T = N ^ α⁻¹` does not affect the
power asymptotic `N = T ^ (α + o(1))`. -/
theorem isPowerAsymptotic_of_between_const_rpow
    {N T : VariableObject ℝ} {α A B : ℝ}
    (hα : 0 < α) (hA : 0 < A) (hB : 0 < B)
    (hNtop : Tendsto N atTop atTop)
    (hT : ∀ᶠ i in atTop, A * N i ^ α⁻¹ ≤ T i ∧ T i ≤ B * N i ^ α⁻¹) :
    IsPowerAsymptotic N T α := by
  have hNpos : ∀ᶠ i in atTop, 0 < N i :=
    hNtop.eventually (eventually_gt_atTop 0)
  have hpowtop : Tendsto (fun i ↦ N i ^ α⁻¹) atTop atTop :=
    (tendsto_rpow_atTop (inv_pos.mpr hα)).comp hNtop
  have hpowlarge : ∀ᶠ i in atTop, 1 / A < N i ^ α⁻¹ :=
    hpowtop.eventually (eventually_gt_atTop (1 / A))
  have hTgt : ∀ᶠ i in atTop, 1 < T i := by
    filter_upwards [hpowlarge, hT] with i hi hTi
    have : 1 < A * N i ^ α⁻¹ := by
      simpa [mul_comm] using (div_lt_iff₀ hA).1 hi
    exact this.trans_le hTi.1
  exact isPowerAsymptotic_of_logb_tendsto hTgt
    hNpos (tendsto_logb_of_between_const_rpow hα hA hB hNtop hT)

/-! ## Power bounds -/

private lemma isLittleO_rpow_comp_tendsto_of_lt
    {T : VariableObject ℝ} {β γ : ℝ}
    (hT : Tendsto T atTop atTop) (hβγ : β < γ) :
    (fun i ↦ T i ^ β) =o[atTop] (fun i ↦ T i ^ γ) := by
  have hzero : ∀ᶠ i in atTop, T i ^ γ = 0 → T i ^ β = 0 := by
    filter_upwards [hT.eventually (eventually_gt_atTop 0)] with i hi
    exact fun hzero ↦ ((Real.rpow_pos_of_pos hi γ).ne' hzero).elim
  rw [Asymptotics.isLittleO_iff_tendsto' hzero]
  have hlimit := (tendsto_rpow_neg_atTop (sub_pos.mpr hβγ)).comp hT
  apply hlimit.congr'
  filter_upwards [hT.eventually (eventually_gt_atTop 0)] with i hi
  change T i ^ (-(γ - β)) = T i ^ β / T i ^ γ
  rw [← Real.rpow_sub hi]
  congr 1
  ring

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

/-- A power lower bound for an unbounded object forces its exponent to be no larger than any
admissible power-bound exponent. -/
theorem exponent_le_of_isPowerBounded_of_eventually_norm_ge_rpow
    {E : Type*} [SeminormedAddCommGroup E]
    {X : VariableObject E} {T : VariableObject ℝ} {β γ c : ℝ}
    (hT : ∀ i, 1 ≤ T i) (hTunbounded : T.IsUnbounded)
    (hbound : IsPowerBounded X T β) (hc : 0 < c)
    (hlower : ∀ᶠ i in atTop, c * T i ^ γ ≤ ‖X i‖) :
    γ ≤ β := by
  by_contra hγβ
  have hβγ : β < γ := lt_of_not_ge hγβ
  let ε := (γ - β) / 2
  have hε : 0 < ε := by dsimp [ε]; linarith
  have hO := (isPowerBounded_iff_forall_pos X T β hT hTunbounded).1 hbound ε hε
  have hTtop : Tendsto T atTop atTop := by
    rw [VariableObject.IsUnbounded] at hTunbounded
    exact hTunbounded.congr' <| Filter.Eventually.of_forall fun i ↦ by
      rw [Real.norm_eq_abs, abs_of_nonneg (le_trans zero_le_one (hT i))]
  have hpowO : (fun i ↦ T i ^ γ) =O[atTop] X := by
    apply Asymptotics.IsBigO.of_bound (1 / c)
    filter_upwards [hlower] with i hi
    rw [Real.norm_eq_abs, abs_of_nonneg
      (Real.rpow_nonneg (zero_le_one.trans (hT i)) _)]
    calc
      T i ^ γ ≤ ‖X i‖ / c := (le_div_iff₀ hc).2 (by simpa [mul_comm] using hi)
      _ = 1 / c * ‖X i‖ := by ring
  have hsmall : (fun i ↦ T i ^ (β + ε)) =o[atTop] (fun i ↦ T i ^ γ) :=
    isLittleO_rpow_comp_tendsto_of_lt hTtop (by dsimp [ε]; linarith)
  have hnonzero : ∃ᶠ i in atTop, T i ^ (β + ε) ≠ 0 :=
    (Filter.Eventually.of_forall fun i ↦
      (Real.rpow_pos_of_pos (zero_lt_one.trans_le (hT i)) _).ne').frequently
  exact hsmall.not_isBigO hnonzero (hpowO.trans hO)

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

private lemma exponentSumBounds_nonempty (α : ℝ≥0) :
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

private lemma exponentSumBounds_bddBelow {α : ℝ≥0} :
    BddBelow {β : ℝ | IsExponentSumBound α β} :=
  ⟨0, fun _ hβ ↦ hβ.nonneg⟩

/-! ## The least admissible exponent -/

private lemma exists_isExponentSumBound_lt_add
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
