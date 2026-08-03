module

public import Expdb.ExponentialSums.ExponentSumGrowth

import Expdb.Basic.AutomaticUniformity
import Mathlib.Analysis.Calculus.ContDiff.Bounds
import Mathlib.Analysis.SpecialFunctions.Log.Base
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv

/-!
# Non-asymptotic exponential sum growth bounds

This module formalizes the non-asymptotic definition of exponential sum exponent from the
blueprint's Exponential sum growth exponents chapter (`beta-chapter`). It gives a fixed-parameter,
epsilon--delta characterization of the exponential sum growth exponent.
-/

@[expose] public section

open Filter Topology
open scoped ContDiff Expdb FourierTransform NNReal

noncomputable section

namespace Expdb

/-! ## Fixed-parameter formulation -/

/-- The exponential sum at fixed parameters. -/
def exponentialSumAt
    (F : ℝ → ℝ) (T N : ℝ) (a b : ℕ) : ℂ :=
  ∑ n ∈ Finset.Icc a b,
    (𝐞 (T * F ((n : ℝ) / N)) : ℂ)

/-- A fixed phase function whose model-phase errors through order `P` are at most `δ`. -/
def IsApproximateModelPhaseFunction
    (F : ℝ → ℝ) (σ : ℝ) (P : ℕ) (δ : ℝ) : Prop :=
  ContDiffOn ℝ ∞ F phaseInterval ∧
    ∀ p : ℕ, p ≤ P →
      ∀ u : phaseInterval,
        ‖iteratedDerivWithin (p + 1) F phaseInterval u -
          iteratedDerivWithin p (fun v : ℝ ↦ v ^ (-σ)) phaseInterval u‖ ≤ δ

/-- The fixed-parameter epsilon--delta bound from the blueprint lemma `beta-asymp`
    (Non-asymptotic definition of `β`). -/
def IsExponentSumBoundNonAsymptotic (α : ℝ≥0) (β : ℝ) : Prop :=
  ∀ ε : ℝ, 0 < ε →
    ∀ σ : ℝ, 0 < σ →
      ∃ δ : ℝ, 0 < δ ∧
        ∃ P : ℕ, 1 ≤ P ∧
          ∃ C : ℝ, 1 ≤ C ∧
            ∀ (T N : ℝ) (F : ℝ → ℝ) (a b : ℕ),
              C ≤ T →
              T ^ ((α : ℝ) - δ) ≤ N →
              N ≤ T ^ ((α : ℝ) + δ) →
              IsApproximateModelPhaseFunction F σ P δ →
              N ≤ (a : ℝ) →
              (b : ℝ) ≤ 2 * N →
              ‖exponentialSumAt F T N a b‖ ≤ C * T ^ (β + ε)

/-! ## Uniform approximation by fixed phases -/

private lemma eventually_isApproximate_of_modelError
    {F : VariableFunction (VariableObject.fixed ℝ) ℝ} {σ : ℝ}
    (hphase : IsPhaseFunction F)
    (herror : ∀ p : ℕ, (modelPhaseError F σ p).IsPointwiseInfinitesimal)
    (P : ℕ) {δ : ℝ} (hδ : 0 < δ) :
    ∀ᶠ i in atTop, IsApproximateModelPhaseFunction (F i) σ P δ := by
  have hp : ∀ p ∈ Finset.range (P + 1), ∀ᶠ i in atTop,
      ∀ u : phaseInterval, ‖modelPhaseError F σ p i u‖ < δ :=
    fun p _ ↦
      (VariableFunction.isPointwiseInfinitesimal_iff_forall_pos_uniform
        (VariableObject.fixed phaseInterval)
        (fun _ ↦ ⟨1, by simp [phaseInterval]⟩)
        (modelPhaseError F σ p)).1 (herror p) δ hδ
  have hall :
      ∀ᶠ i in atTop, ∀ p ∈ Finset.range (P + 1),
        ∀ u : phaseInterval, ‖modelPhaseError F σ p i u‖ < δ :=
    (Finset.range (P + 1)).eventually_all.mpr hp
  filter_upwards [hall] with i hi
  refine ⟨hphase i, ?_⟩
  intro p hp' u
  exact (hi p (Finset.mem_range.mpr (by omega)) u).le

/-! ## From fixed bounds to asymptotic bounds -/

private lemma isExponentSumBound_of_nonAsymptotic
    {α : ℝ≥0} {β : ℝ}
    (hbound : IsExponentSumBoundNonAsymptotic α β) :
    IsExponentSumBound α β := by
  intro N T F a b hN hT hTunbounded hNT hF hab
  apply (isPowerBounded_iff_forall_pos
    (exponentialSum F T N a b) T β hT hTunbounded).2
  intro ε hε
  rcases hF with ⟨hphase, σ, hσ, herror⟩
  obtain ⟨δ, hδ, P, hP, C, hC, hfixed⟩ :=
    hbound ε hε σ hσ
  have hTC : ∀ᶠ i in atTop, C ≤ T i := by
    have hnorm := (VariableObject.isUnbounded_iff_forall_eventually_norm_ge T).1
      hTunbounded C
    filter_upwards [hnorm] with i hi
    simpa [Real.norm_eq_abs, abs_of_nonneg
      (le_trans zero_le_one (hT i))] using hi
  have hNT' :=
    hNT.eventually_between
      (Filter.Eventually.of_forall hT) hδ
  have happrox :=
    eventually_isApproximate_of_modelError hphase herror P hδ
  refine Asymptotics.IsBigO.of_bound C ?_
  filter_upwards [hTC, hNT', happrox] with i hiT hiNT hiF
  simpa only [exponentialSum, exponentialSumAt, Real.norm_eq_abs,
    abs_of_nonneg (Real.rpow_nonneg (le_trans zero_le_one (hT i)) _)] using
      hfixed (T i) (N i) (F i) (a i) (b i)
        hiT hiNT.1 hiNT.2 hiF (hab i).1 (hab i).2

/-! ## Uniform bounds for approximate model phases -/

/-- An approximate model phase with sufficiently small error has a uniform positive first
derivative and uniform bounds on all derivatives through order `P + 1`, with constants depending
only on `σ` and `P`. -/
theorem approximate_model_phase_deriv_bounds
    {σ : ℝ} (hσ : 0 < σ) (P : ℕ) :
    ∃ K : ℝ, 1 ≤ K ∧
      ∀ {δ : ℝ} {F : ℝ → ℝ},
        δ ≤ min ((2 : ℝ) ^ (-σ) / 2) 1 →
        IsApproximateModelPhaseFunction F σ P δ →
        (∀ u ∈ phaseInterval,
          (2 : ℝ) ^ (-σ) / 2 ≤ iteratedDerivWithin 1 F phaseInterval u) ∧
        (∀ k : ℕ, 1 ≤ k → k ≤ P + 1 → ∀ u ∈ phaseInterval,
          ‖iteratedDerivWithin k F phaseInterval u‖ ≤ K) := by
  have href : ContDiffOn ℝ ∞ (fun u : ℝ ↦ u ^ (-σ)) phaseInterval := by
    intro u hu
    exact (Real.contDiffAt_rpow_const_of_ne (by
      exact ne_of_gt (lt_of_lt_of_le zero_lt_one hu.1))).contDiffWithinAt
  have hreference (p : ℕ) :
      ∃ B : ℝ, 0 ≤ B ∧ ∀ u : phaseInterval,
        ‖iteratedDerivWithin p (fun v : ℝ ↦ v ^ (-σ)) phaseInterval u‖ ≤ B := by
    have hcont : ContinuousOn
        (iteratedDerivWithin p (fun v : ℝ ↦ v ^ (-σ)) phaseInterval) phaseInterval :=
      href.continuousOn_iteratedDerivWithin
        (ENat.natCast_le_of_coe_top_le_withTop le_rfl p)
        (uniqueDiffOn_Icc (by norm_num [phaseInterval]))
    obtain ⟨B, hB⟩ := bddAbove_def.mp (isCompact_Icc.bddAbove_image hcont.norm)
    refine ⟨max B 0, le_max_right _ _, fun u ↦ ?_⟩
    exact (hB _ ⟨u, u.property, rfl⟩).trans (le_max_left _ _)
  choose B hBnonneg hB using hreference
  let K : ℝ := 1 + ∑ p ∈ Finset.range (P + 1), B p
  have hK : 1 ≤ K := by
    dsimp [K]
    exact le_add_of_nonneg_right (Finset.sum_nonneg fun p _ ↦ hBnonneg p)
  refine ⟨K, hK, ?_⟩
  intro δ F hδ hF
  have hδc : δ ≤ (2 : ℝ) ^ (-σ) / 2 := hδ.trans (min_le_left _ _)
  have hδone : δ ≤ 1 := hδ.trans (min_le_right _ _)
  constructor
  · intro u hu
    have huPos : 0 < (u : ℝ) := zero_lt_one.trans_le hu.1
    have huLower : (2 : ℝ) ^ (-σ) ≤ (u : ℝ) ^ (-σ) :=
      Real.rpow_le_rpow_of_nonpos huPos hu.2 (neg_nonpos.mpr hσ.le)
    have he := hF.2 0 (Nat.zero_le P) ⟨u, hu⟩
    rw [Real.norm_eq_abs, abs_le] at he
    simp only [zero_add, iteratedDerivWithin_zero] at he
    linarith
  · intro k hk hkP u hu
    obtain ⟨p, rfl⟩ := Nat.exists_eq_add_of_le hk
    have hpP : p ≤ P := by omega
    have hpMem : p ∈ Finset.range (P + 1) := Finset.mem_range.mpr (by omega)
    have he := hF.2 p hpP ⟨u, hu⟩
    have htriangle :
        ‖iteratedDerivWithin (p + 1) F phaseInterval u‖ ≤
          ‖iteratedDerivWithin (p + 1) F phaseInterval u -
            iteratedDerivWithin p (fun v : ℝ ↦ v ^ (-σ)) phaseInterval u‖ +
          ‖iteratedDerivWithin p (fun v : ℝ ↦ v ^ (-σ)) phaseInterval u‖ :=
      norm_le_norm_sub_add _ _
    calc
      ‖iteratedDerivWithin (1 + p) F phaseInterval u‖ ≤ 1 + B p := by
        simpa [Nat.add_comm] using htriangle.trans (add_le_add (he.trans hδone) (hB p ⟨u, hu⟩))
      _ ≤ K := by
        dsimp [K]
        gcongr
        exact Finset.single_le_sum (fun q _ ↦ hBnonneg q) hpMem

/-- The logarithm is an exact approximate model phase at every finite order. -/
theorem isApproximateModelPhaseFunction_log (P : ℕ) :
    IsApproximateModelPhaseFunction Real.log 1 P 0 := by
  refine ⟨isModelPhaseFunction_log.1 0, ?_⟩
  intro p _ u
  rw [iteratedDerivWithin_log_eq_rpow_neg_one, sub_self, norm_zero]

/-! ## Building asymptotic counterexamples -/

private lemma isModelPhaseFunction_of_approximations
    {F : VariableFunction (VariableObject.fixed ℝ) ℝ}
    {σ : ℝ} {P : VariableObject ℕ} {δ : VariableObject ℝ}
    (hσ : 0 < σ)
    (hP : ∀ p : ℕ, ∀ᶠ i in atTop, p ≤ P i)
    (hδnonneg : ∀ i, 0 ≤ δ i) (hδ : δ.IsInfinitesimal)
    (happrox : ∀ i, IsApproximateModelPhaseFunction (F i) σ (P i) (δ i)) :
    IsModelPhaseFunction F := by
  refine ⟨fun i ↦ (happrox i).1, σ, hσ, ?_⟩
  intro p
  apply (VariableFunction.isPointwiseInfinitesimal_iff_forall_pos_uniform
    (VariableObject.fixed phaseInterval)
    (fun _ ↦ ⟨1, by simp [phaseInterval]⟩)
    (modelPhaseError F σ p)).2
  intro ε hε
  have hδsmall :=
    (VariableObject.isInfinitesimal_iff_forall_pos δ).1 hδ ε hε
  filter_upwards [hP p, hδsmall] with i hip hiδ
  intro u
  have herror := (happrox i).2 p hip u
  rw [Real.norm_eq_abs, abs_of_nonneg (hδnonneg i)] at hiδ
  exact lt_of_le_of_lt herror hiδ

private lemma norm_exponentialSumAt_le_add_one
    (F : ℝ → ℝ) (T N : ℝ) (a b : ℕ) :
    ‖exponentialSumAt F T N a b‖ ≤ (b : ℝ) + 1 := by
  calc
    ‖exponentialSumAt F T N a b‖ ≤
        ∑ n ∈ Finset.Icc a b,
          ‖(𝐞 (T * F ((n : ℝ) / N)) : ℂ)‖ := norm_sum_le _ _
    _ = ((Finset.Icc a b).card : ℝ) := by simp
    _ ≤ (b : ℝ) + 1 := by
      rw [Nat.card_Icc]
      exact_mod_cast Nat.sub_le _ _

private lemma nonAsymptotic_of_isExponentSumBound
    {α : ℝ≥0} {β : ℝ} (hbound : IsExponentSumBound α β) :
    IsExponentSumBoundNonAsymptotic α β := by
  intro ε hε σ hσ
  by_contra hfailure
  push Not at hfailure
  let δ : VariableObject ℝ := fun i ↦ 1 / ((i : ℝ) + 1)
  let P : VariableObject ℕ := fun i ↦ i + 1
  let C : VariableObject ℝ := fun i ↦ (i : ℝ) + 3
  have hδpos : ∀ i, 0 < δ i := fun i ↦ by
    dsimp [δ]
    positivity
  have hPone : ∀ i, 1 ≤ P i := fun i ↦ by
    dsimp [P]
    omega
  have hCone : ∀ i, 1 ≤ C i := fun i ↦ by
    dsimp [C]
    have hi : (0 : ℝ) ≤ i := Nat.cast_nonneg i
    linarith
  choose T N F a b hdata using fun i ↦
    hfailure (δ i) (hδpos i) (P i) (hPone i) (C i) (hCone i)
  have hCT : ∀ i, C i ≤ T i := fun i ↦ (hdata i).1
  have hlower : ∀ i, T i ^ ((α : ℝ) - δ i) ≤ N i :=
    fun i ↦ (hdata i).2.1
  have hupper : ∀ i, N i ≤ T i ^ ((α : ℝ) + δ i) :=
    fun i ↦ (hdata i).2.2.1
  have happrox : ∀ i,
      IsApproximateModelPhaseFunction (F i) σ (P i) (δ i) :=
    fun i ↦ (hdata i).2.2.2.1
  have ha : ∀ i, N i ≤ (a i : ℝ) :=
    fun i ↦ (hdata i).2.2.2.2.1
  have hb : ∀ i, (b i : ℝ) ≤ 2 * N i :=
    fun i ↦ (hdata i).2.2.2.2.2.1
  have hviolate : ∀ i,
      C i * T i ^ (β + ε) <
        ‖exponentialSumAt (F i) (T i) (N i) (a i) (b i)‖ :=
    fun i ↦ (hdata i).2.2.2.2.2.2
  have hTgt : ∀ i, 1 < T i := fun i ↦ by
    have hCthree : 3 ≤ C i := by
      dsimp [C]
      exact le_add_of_nonneg_left (Nat.cast_nonneg i)
    linarith [hCT i]
  have hTone : ∀ i, 1 ≤ T i := fun i ↦ (hTgt i).le
  have hTunbounded : VariableObject.IsUnbounded T := by
    apply (VariableObject.isUnbounded_iff_forall_eventually_norm_ge T).2
    intro D
    obtain ⟨k : ℕ, hk : D ≤ (k : ℝ)⟩ := exists_nat_ge D
    filter_upwards [eventually_ge_atTop k] with i hi
    rw [Real.norm_eq_abs, abs_of_nonneg (le_trans zero_le_one (hTone i))]
    calc
      D ≤ (k : ℝ) := hk
      _ ≤ (i : ℝ) := by exact_mod_cast hi
      _ ≤ C i := by
        dsimp [C]
        linarith
      _ ≤ T i := hCT i
  have hδinfinitesimal : δ.IsInfinitesimal := by
    rw [VariableObject.IsInfinitesimal]
    have hlim : Tendsto (fun i : ℕ ↦ 1 / ((i : ℝ) + 1))
        atTop (𝓝 0) :=
      tendsto_one_div_add_atTop_nhds_zero_nat
    convert hlim using 1
    ext i
    rw [Real.norm_eq_abs, abs_of_pos (hδpos i)]
  have hPtop : ∀ p : ℕ, ∀ᶠ i in atTop, p ≤ P i := by
    intro p
    filter_upwards [eventually_ge_atTop p] with i hi
    dsimp [P]
    omega
  have hNpos : ∀ i, 0 < N i := fun i ↦
    lt_of_lt_of_le
      (Real.rpow_pos_of_pos (lt_trans zero_lt_one (hTgt i)) _)
      (hlower i)
  have hβnonneg : 0 ≤ β := hbound.nonneg
  have hNone : ∀ i, 1 ≤ N i := by
    intro i
    by_contra hiN
    have hNlt : N i < 1 := lt_of_not_ge hiN
    have hblt : (b i : ℝ) < 2 := lt_of_le_of_lt (hb i) (by linarith)
    have hbNat : b i < 2 := by exact_mod_cast hblt
    have hsum :
        ‖exponentialSumAt (F i) (T i) (N i) (a i) (b i)‖ ≤ 2 := by
      calc
        ‖exponentialSumAt (F i) (T i) (N i) (a i) (b i)‖ ≤
            (b i : ℝ) + 1 :=
          norm_exponentialSumAt_le_add_one
            (F i) (T i) (N i) (a i) (b i)
        _ ≤ 2 := by
          have : b i ≤ 1 := by omega
          exact_mod_cast Nat.add_le_add_right this 1
    have hpow : 1 ≤ T i ^ (β + ε) :=
      Real.one_le_rpow (hTone i) (by linarith)
    have hCthree : 3 ≤ C i := by
      dsimp [C]
      exact le_add_of_nonneg_left (Nat.cast_nonneg i)
    have hrhs : 3 ≤ C i * T i ^ (β + ε) := by
      calc
        3 ≤ C i := hCthree
        _ = C i * 1 := by ring
        _ ≤ C i * T i ^ (β + ε) :=
          mul_le_mul_of_nonneg_left hpow (le_trans (by norm_num) hCthree)
    linarith [hviolate i]
  have hNT : IsPowerAsymptotic N T (α : ℝ) :=
    isPowerAsymptotic_of_between hTgt hNpos
      (fun i ↦ (hδpos i).le) hδinfinitesimal
      (fun i ↦ ⟨hlower i, hupper i⟩)
  have hF : IsModelPhaseFunction F :=
    isModelPhaseFunction_of_approximations hσ hPtop
      (fun i ↦ (hδpos i).le) hδinfinitesimal happrox
  have hasymptotic :=
    hbound N T F a b hNone hTone hTunbounded hNT hF
      (fun i ↦ ⟨ha i, hb i⟩)
  have hO :=
    (isPowerBounded_iff_forall_pos
      (exponentialSum F T N a b) T β hTone hTunbounded).1
      hasymptotic (ε / 2) (by linarith)
  obtain ⟨K, hK⟩ := hO.bound
  have hCK : ∀ᶠ i in atTop, K ≤ C i := by
    obtain ⟨k : ℕ, hk : K ≤ (k : ℝ)⟩ := exists_nat_ge K
    filter_upwards [eventually_ge_atTop k] with i hi
    calc
      K ≤ (k : ℝ) := hk
      _ ≤ (i : ℝ) := by exact_mod_cast hi
      _ ≤ C i := by
        dsimp [C]
        linarith
  obtain ⟨i, hiK, hiCK⟩ := (hK.and hCK).exists
  have hnorm :
      ‖exponentialSumAt (F i) (T i) (N i) (a i) (b i)‖ ≤
        K * T i ^ (β + ε / 2) := by
    simpa only [exponentialSum, exponentialSumAt, Real.norm_eq_abs,
      abs_of_nonneg (Real.rpow_nonneg (le_trans zero_le_one (hTone i)) _)] using hiK
  have hdominated :
      ‖exponentialSumAt (F i) (T i) (N i) (a i) (b i)‖ ≤
        C i * T i ^ (β + ε) := by
    calc
      ‖exponentialSumAt (F i) (T i) (N i) (a i) (b i)‖ ≤
          K * T i ^ (β + ε / 2) := hnorm
      _ ≤ C i * T i ^ (β + ε / 2) := by
        exact mul_le_mul_of_nonneg_right hiCK
          (Real.rpow_nonneg (le_trans zero_le_one (hTone i)) _)
      _ ≤ C i * T i ^ (β + ε) :=
        mul_le_mul_of_nonneg_left
          (Real.rpow_le_rpow_of_exponent_le (hTone i) (by linarith))
          (le_trans zero_le_one (hCone i))
  exact (not_lt_of_ge hdominated) (hviolate i)

/-! ## Equivalence with the asymptotic definition -/

/-- The non-asymptotic characterization of the exponential sum growth exponent. -/
theorem exponentSumGrowthExponent_le_iff_nonAsymptotic
    {α : ℝ≥0} {β : ℝ} :
    exponentSumGrowthExponent α ≤ β ↔
      IsExponentSumBoundNonAsymptotic α β := by
  rw [exponentSumGrowthExponent_le_iff]
  exact ⟨nonAsymptotic_of_isExponentSumBound,
    isExponentSumBound_of_nonAsymptotic⟩

end Expdb
