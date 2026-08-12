module

public import Expdb.ExponentialSums.TrivialBounds

import Expdb.ExponentialSums.ExponentSumGrowthNonAsymptotic
import Expdb.Basic.AutomaticUniformity
import Mathlib.Algebra.Order.Floor.Div
import Mathlib.Analysis.Calculus.IteratedDeriv.Lemmas
import Mathlib.Analysis.Fourier.ZMod
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Data.Nat.Cast.Order.Field
import Mathlib.Data.Finset.Max
import Mathlib.Topology.MetricSpace.Bounded
import Mathlib.Topology.Semicontinuity.Basic

/-!
# Upper semicontinuity of the exponential-sum growth exponent

This file proves the upper semicontinuity assertion in the blueprint.  The main input is scale
transference: dilation and finite Fourier inversion increase the scale, while decomposition into
residue classes decreases it.  Both operations preserve the fixed model order of the phase.
-/

@[expose] public section

open Filter Topology
open scoped ContDiff Expdb FourierTransform NNReal

noncomputable section

namespace Expdb

/-! ## Finite Fourier orthogonality -/

private lemma sum_stdAddChar_mul_eq_ite (q m : ℕ) [NeZero q] :
    (∑ h : ZMod q, ZMod.stdAddChar (h * (m : ZMod q))) =
      if q ∣ m then (q : ℂ) else 0 := by
  rw [AddChar.sum_mulShift (m : ZMod q) (ZMod.isPrimitive_stdAddChar q)]
  simp only [ZMod.natCast_eq_zero_iff m q, ZMod.card]
  split_ifs <;> norm_num

private lemma stdAddChar_mul_nat_eq_fourierChar (q m : ℕ) [NeZero q] (h : ZMod q) :
    ZMod.stdAddChar (h * (m : ZMod q)) =
      (𝐞 ((h.val : ℝ) * (m : ℝ) / q) : ℂ) := by
  rw [← h.natCast_zmod_val, ← Nat.cast_mul, ZMod.stdAddChar_apply,
    ZMod.toCircle_natCast, Real.fourierChar_apply]
  push_cast
  rw [ZMod.val_natCast_of_lt h.val_lt]
  congr 1
  ring_nf

/-! ## Stability of model phases -/

private lemma iteratedDerivWithin_comp_affine_of_mapsTo
    {f : ℝ → ℝ} {s t : Set ℝ} {x c d : ℝ} {n : ℕ}
    (hf : ContDiffOn ℝ n f t) (hs : UniqueDiffOn ℝ s) (ht : UniqueDiffOn ℝ t)
    (hx : x ∈ s) (hst : Set.MapsTo (fun y ↦ c * y + d) s t) :
    iteratedDerivWithin n (fun y ↦ f (c * y + d)) s x =
      c ^ n * iteratedDerivWithin n f t (c * x + d) := by
  induction n generalizing x with
  | zero => simp
  | succ n ih =>
      have hcx : c * x + d ∈ t := hst hx
      have heq : Set.EqOn
          (iteratedDerivWithin n (fun y ↦ f (c * y + d)) s)
          (fun y ↦ c ^ n * iteratedDerivWithin n f t (c * y + d)) s :=
        fun y hy ↦ ih hf.of_succ hy
      have houter : DifferentiableWithinAt ℝ (iteratedDerivWithin n f t) t (c * x + d) :=
        hf.differentiableOn_iteratedDerivWithin (Nat.cast_lt.mpr n.lt_succ_self) ht _ hcx
      have haffine : HasDerivWithinAt (fun y : ℝ ↦ c * y + d) c s x := by
        simpa using (((hasDerivAt_id x).const_mul c).add_const d |>.hasDerivWithinAt)
      have hcomp : DifferentiableWithinAt ℝ
          (fun y ↦ iteratedDerivWithin n f t (c * y + d)) s x :=
        by
          change DifferentiableWithinAt ℝ
            ((iteratedDerivWithin n f t) ∘ fun y ↦ c * y + d) s x
          exact houter.comp x haffine.differentiableWithinAt hst
      rw [iteratedDerivWithin_succ, derivWithin_congr heq (ih hf.of_succ hx),
        derivWithin_const_mul _ hcomp, iteratedDerivWithin_succ,
        ← Function.comp_def, derivWithin.scomp x houter
          haffine.differentiableWithinAt hst]
      rw [haffine.derivWithin (hs.uniqueDiffWithinAt hx)]
      simp only [pow_succ]
      ring

private theorem modelPhaseData_add_infinitesimal_linear
    {F : VariableFunction (VariableObject.fixed ℝ) ℝ} {σ : ℝ}
    (hphase : IsPhaseFunction F)
    (herror : ∀ p : ℕ, (modelPhaseError F σ p).IsPointwiseInfinitesimal)
    {c : VariableObject ℝ} (hc : c.IsInfinitesimal) :
    IsPhaseFunction (fun i u ↦ F i u + c i * u) ∧
      ∀ p : ℕ, (modelPhaseError (fun i u ↦ F i u + c i * u) σ p).IsPointwiseInfinitesimal := by
  refine ⟨fun i ↦ (hphase i).add (by fun_prop), ?_⟩
  intro p u
  rw [VariableObject.IsInfinitesimal]
  have herr := herror p u
  rw [VariableObject.IsInfinitesimal] at herr
  have hc' : Tendsto (fun i ↦ ‖c i‖) atTop (𝓝 0) := hc
  apply squeeze_zero'
  · exact Filter.Eventually.of_forall fun _ ↦ norm_nonneg _
  · filter_upwards [] with i
    change ‖modelPhaseErrorAt (fun x ↦ F i x + c i * x) σ p (u i)‖ ≤
      ‖modelPhaseErrorAt (F i) σ p (u i)‖ + ‖c i‖
    have hu := (u i).property
    have hFi := hphase i
    change ‖modelPhaseErrorAt (F i + fun x ↦ c i * x) σ p (u i)‖ ≤ _
    rw [modelPhaseErrorAt, modelPhaseErrorAt,
      iteratedDerivWithin_add hu
        (uniqueDiffOn_Icc (by norm_num [phaseInterval]))
        ((hFi (u i) hu).of_le
          (ENat.natCast_le_of_coe_top_le_withTop le_rfl (p + 1)))
        ((show ContDiffOn ℝ ∞ (fun x : ℝ ↦ c i * x) phaseInterval by fun_prop)
          (u i) hu |>.of_le
            (ENat.natCast_le_of_coe_top_le_withTop le_rfl (p + 1)))]
    rw [iteratedDerivWithin_const_mul_field,
      iteratedDerivWithin_fun_id hu (uniqueDiffOn_Icc (by norm_num [phaseInterval]))]
    by_cases hp : p = 0
    · subst p
      calc
        _ = ‖(iteratedDerivWithin 1 (F i) phaseInterval (u i) -
              modelPhase σ (u i)) + c i‖ := by
            congr 1
            simp only [iteratedDerivWithin_zero, if_neg (by omega : (1 : ℕ) ≠ 0),
              if_true]
            ring
        _ ≤ ‖iteratedDerivWithin 1 (F i) phaseInterval (u i) -
              modelPhase σ (u i)‖ + ‖c i‖ := norm_add_le _ _
    · simp only [if_neg (by omega : p + 1 ≠ 0), if_neg (by omega : p + 1 ≠ 1),
        mul_zero, add_zero]
      exact le_add_of_nonneg_right (norm_nonneg _)
  · simpa [modelPhaseError] using herr.add hc'

private theorem modelPhaseData_comp_affine_tendsto_id
    {F : VariableFunction (VariableObject.fixed ℝ) ℝ} {σ : ℝ}
    (hphase : IsPhaseFunction F)
    (herror : ∀ p : ℕ, (modelPhaseError F σ p).IsPointwiseInfinitesimal)
    {c d : VariableObject ℝ}
    (hc : Tendsto c atTop (𝓝 1)) (hd : Tendsto d atTop (𝓝 0))
    (hmap : ∀ i, Set.MapsTo (fun u ↦ c i * u + d i) phaseInterval phaseInterval) :
    IsPhaseFunction (fun i u ↦ F i (c i * u + d i)) ∧
      ∀ p : ℕ,
        (modelPhaseError (fun i u ↦ F i (c i * u + d i)) σ p).IsPointwiseInfinitesimal := by
  refine ⟨fun i ↦ (hphase i).comp (by fun_prop) (hmap i), ?_⟩
  intro p u
  rw [VariableObject.IsInfinitesimal]
  let v : ∀ i, phaseInterval := fun i ↦ ⟨c i * (u i : ℝ) + d i, hmap i (u i).property⟩
  have herr := herror p v
  rw [VariableObject.IsInfinitesimal] at herr
  have hvsub : Tendsto (fun i ↦ (v i : ℝ) - (u i : ℝ)) atTop (𝓝 0) := by
    have hcsub : Tendsto (fun i ↦ c i - 1) atTop (𝓝 0) := by
      simpa using hc.sub (tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ (1 : ℝ)) atTop (𝓝 1))
    have hmul : Tendsto (fun i ↦ (c i - 1) * (u i : ℝ)) atTop (𝓝 0) := by
      rw [tendsto_zero_iff_norm_tendsto_zero]
      apply squeeze_zero' (g := fun i ↦ 2 * ‖c i - 1‖)
      · exact Filter.Eventually.of_forall fun _ ↦ norm_nonneg _
      · filter_upwards [] with i
        rw [norm_mul]
        calc
          _ ≤ ‖c i - 1‖ * 2 := mul_le_mul_of_nonneg_left (by
            rw [Real.norm_eq_abs, abs_of_nonneg]
            · exact (u i).property.2
            · exact zero_le_one.trans (u i).property.1) (norm_nonneg _)
          _ = 2 * ‖c i - 1‖ := by ring
      · simpa using hcsub.norm.const_mul 2
    convert hmul.add hd using 1
    · ext i; dsimp [v]; ring_nf
    · ring_nf
  have hcpow : Tendsto (fun i ↦ c i ^ (p + 1)) atTop (𝓝 1) := by
    simpa using hc.pow (p + 1)
  have href : ContinuousOn
      (iteratedDerivWithin p (modelPhase σ) phaseInterval) phaseInterval := by
    have hsmooth : ContDiffOn ℝ ∞ (modelPhase σ) phaseInterval := by
      intro x hx
      exact (Real.contDiffAt_rpow_const_of_ne
        (ne_of_gt (zero_lt_one.trans_le hx.1))).contDiffWithinAt
    exact hsmooth.continuousOn_iteratedDerivWithin
      (ENat.natCast_le_of_coe_top_le_withTop le_rfl p)
      (uniqueDiffOn_Icc (by norm_num [phaseInterval]))
  have hrefUniform := isCompact_Icc.uniformContinuousOn_of_continuous href
  have hrefclose : Tendsto (fun i ↦
      iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) -
        iteratedDerivWithin p (modelPhase σ) phaseInterval (u i)) atTop (𝓝 0) := by
    rw [Metric.uniformContinuousOn_iff] at hrefUniform
    rw [Metric.tendsto_nhds]
    intro ε hε
    obtain ⟨δ, hδ, hclose⟩ := hrefUniform ε hε
    have hevent : ∀ᶠ i in atTop, ‖(v i : ℝ) - (u i : ℝ)‖ < δ := by
      simpa [Real.dist_eq] using (Metric.tendsto_nhds.1 hvsub) δ hδ
    filter_upwards [hevent] with i hi
    simpa [dist_eq_norm] using
      hclose (v i) (v i).property (u i) (u i).property
        (by simpa [Real.dist_eq] using hi)
  obtain ⟨K, hK⟩ := bddAbove_def.mp (isCompact_Icc.bddAbove_image href.norm)
  let D := max K 0
  have hrefBounded : ∀ i,
      ‖iteratedDerivWithin p (modelPhase σ) phaseInterval (v i)‖ ≤ D :=
    fun i ↦ (hK _ ⟨v i, (v i).property, rfl⟩).trans (le_max_left _ _)
  have hscale : Tendsto (fun i ↦
      (c i ^ (p + 1) - 1) *
        iteratedDerivWithin p (modelPhase σ) phaseInterval (v i)) atTop (𝓝 0) := by
    have hzero : Tendsto (fun i ↦ c i ^ (p + 1) - 1) atTop (𝓝 0) := by
      simpa using hcpow.sub
        (tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ (1 : ℝ)) atTop (𝓝 1))
    rw [tendsto_zero_iff_norm_tendsto_zero]
    apply squeeze_zero' (g := fun i ↦ D * ‖c i ^ (p + 1) - 1‖)
    · exact Filter.Eventually.of_forall fun _ ↦ norm_nonneg _
    · filter_upwards [] with i
      rw [norm_mul]
      calc
        _ ≤ ‖c i ^ (p + 1) - 1‖ * D :=
          mul_le_mul_of_nonneg_left (hrefBounded i) (norm_nonneg _)
        _ = D * ‖c i ^ (p + 1) - 1‖ := by ring
    · simpa [mul_comm] using hzero.norm.const_mul D
  have hmajor : Tendsto (fun i ↦
      ‖c i ^ (p + 1)‖ * ‖modelPhaseErrorAt (F i) σ p (v i)‖ +
        ‖(c i ^ (p + 1) - 1) *
            iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) +
          (iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) -
            iteratedDerivWithin p (modelPhase σ) phaseInterval (u i))‖)
      atTop (𝓝 0) := by
    simpa using (hcpow.norm.mul herr).add (hscale.add hrefclose).norm
  apply squeeze_zero' (g := fun i ↦
    ‖c i ^ (p + 1)‖ * ‖modelPhaseErrorAt (F i) σ p (v i)‖ +
      ‖(c i ^ (p + 1) - 1) *
          iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) +
        (iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) -
          iteratedDerivWithin p (modelPhase σ) phaseInterval (u i))‖)
  · exact Filter.Eventually.of_forall fun _ ↦ norm_nonneg _
  · filter_upwards [] with i
    change ‖modelPhaseErrorAt (fun x ↦ F i (c i * x + d i)) σ p (u i)‖ ≤ _
    rw [modelPhaseErrorAt,
      iteratedDerivWithin_comp_affine_of_mapsTo (s := phaseInterval) (t := phaseInterval)
        ((hphase i).of_le
          (ENat.natCast_le_of_coe_top_le_withTop le_rfl (p + 1)))
        (uniqueDiffOn_Icc (by norm_num [phaseInterval]))
        (uniqueDiffOn_Icc (by norm_num [phaseInterval])) (u i).property (hmap i)]
    change ‖c i ^ (p + 1) *
        iteratedDerivWithin (p + 1) (F i) phaseInterval (v i) - _‖ ≤ _
    rw [show iteratedDerivWithin (p + 1) (F i) phaseInterval (v i) =
        modelPhaseErrorAt (F i) σ p (v i) +
          iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) by
      simp [modelPhaseErrorAt]]
    calc
      _ = ‖c i ^ (p + 1) * modelPhaseErrorAt (F i) σ p (v i) +
          ((c i ^ (p + 1) - 1) *
              iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) +
            (iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) -
              iteratedDerivWithin p (modelPhase σ) phaseInterval (u i)))‖ := by
        congr 1; ring
      _ ≤ ‖c i ^ (p + 1)‖ * ‖modelPhaseErrorAt (F i) σ p (v i)‖ +
          ‖(c i ^ (p + 1) - 1) *
              iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) +
            (iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) -
              iteratedDerivWithin p (modelPhase σ) phaseInterval (u i))‖ := by
        simpa only [norm_mul] using norm_add_le
          (c i ^ (p + 1) * modelPhaseErrorAt (F i) σ p (v i))
          ((c i ^ (p + 1) - 1) *
              iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) +
            (iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) -
              iteratedDerivWithin p (modelPhase σ) phaseInterval (u i)))
      _ = _ := by rfl
  · exact hmajor

/-! ## Reindexing exponential sums -/

private lemma exponentialSumAt_dilate (F : ℝ → ℝ) (T N : ℝ) (a b q : ℕ) [NeZero q]
    (hq : 0 < q) (hT : T ≠ 0) (hN : N ≠ 0) :
    exponentialSumAt F T N a b =
      (q : ℂ)⁻¹ * ∑ h : ZMod q,
        exponentialSumAt (fun u ↦ F u + (h.val : ℝ) * N / T * u)
          T (q * N) (q * a) (q * b) := by
  rw [exponentialSumAt]
  simp_rw [exponentialSumAt, oscillatory]
  rw [Finset.sum_comm]
  classical
  symm
  calc
    (q : ℂ)⁻¹ * ∑ m ∈ Finset.Icc (q * a) (q * b), ∑ h : ZMod q,
        (𝐞 (T * (F ((m : ℝ) / (q * N)) +
          (h.val : ℝ) * N / T * ((m : ℝ) / (q * N)))) : ℂ) =
        ∑ m ∈ Finset.Icc (q * a) (q * b),
          (q : ℂ)⁻¹ * ((𝐞 (T * F ((m : ℝ) / (q * N))) : ℂ) *
            ∑ h : ZMod q, ZMod.stdAddChar (h * (m : ZMod q))) := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro m hm
      apply congrArg (fun z : ℂ ↦ (q : ℂ)⁻¹ * z)
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro h _
      have harg :
          T * (F ((m : ℝ) / (q * N)) +
            (h.val : ℝ) * N / T * ((m : ℝ) / (q * N))) =
            T * F ((m : ℝ) / (q * N)) + (h.val : ℝ) * m / q := by
        field_simp [show (q : ℝ) ≠ 0 by exact_mod_cast hq.ne', hT, hN]
      rw [harg, AddChar.map_add_eq_mul, stdAddChar_mul_nat_eq_fourierChar]
      exact Circle.coe_mul _ _
    _ = ∑ m ∈ Finset.Icc (q * a) (q * b),
          if q ∣ m then (𝐞 (T * F ((m : ℝ) / (q * N))) : ℂ) else 0 := by
      apply Finset.sum_congr rfl
      intro m hm
      rw [sum_stdAddChar_mul_eq_ite]
      split_ifs with hdiv
      · field_simp [show (q : ℂ) ≠ 0 by exact_mod_cast hq.ne']
      · simp
    _ = ∑ n ∈ Finset.Icc a b, (𝐞 (T * F ((n : ℝ) / N)) : ℂ) := by
      symm
      rw [← Finset.sum_filter]
      apply Finset.sum_bij (fun n _ ↦ q * n)
      · intro n hn
        simp only [Finset.mem_Icc] at hn
        simp only [Finset.mem_filter, Finset.mem_Icc]
        exact ⟨⟨Nat.mul_le_mul_left q hn.1, Nat.mul_le_mul_left q hn.2⟩,
          dvd_mul_right q n⟩
      · intro n₁ hn₁ n₂ hn₂ heq
        exact Nat.eq_of_mul_eq_mul_left hq heq
      · intro m hm
        simp only [Finset.mem_filter, Finset.mem_Icc] at hm
        obtain ⟨n, rfl⟩ := hm.2
        refine ⟨n, ?_, rfl⟩
        simp only [Finset.mem_Icc]
        constructor
        · exact Nat.le_of_mul_le_mul_left hm.1.1 hq
        · exact Nat.le_of_mul_le_mul_left hm.1.2 hq
      · intro n hn
        push_cast
        congr 3
        field_simp [show (q : ℝ) ≠ 0 by exact_mod_cast hq.ne']

private lemma norm_exponentialSumAt_dilate_le (F : ℝ → ℝ) (T N : ℝ) (a b q : ℕ)
    [NeZero q]
    (hq : 0 < q) (hT : T ≠ 0) (hN : N ≠ 0) :
    ‖exponentialSumAt F T N a b‖ ≤
      ∑ h : ZMod q,
        ‖exponentialSumAt (fun u ↦ F u + (h.val : ℝ) * N / T * u)
          T (q * N) (q * a) (q * b)‖ / q := by
  rw [exponentialSumAt_dilate F T N a b q hq hT hN, norm_mul, norm_inv,
    Complex.norm_natCast]
  calc
    (q : ℝ)⁻¹ * ‖∑ h : ZMod q,
        exponentialSumAt (fun u ↦ F u + (h.val : ℝ) * N / T * u)
          T (q * N) (q * a) (q * b)‖ ≤
        (q : ℝ)⁻¹ * ∑ h : ZMod q,
          ‖exponentialSumAt (fun u ↦ F u + (h.val : ℝ) * N / T * u)
            T (q * N) (q * a) (q * b)‖ :=
      mul_le_mul_of_nonneg_left (norm_sum_le _ _) (by positivity)
    _ = _ := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro h _
      simp only [div_eq_mul_inv]
      ring

private def residueScale (N : ℝ) (q r : ℕ) : ℝ := (N - r) / q

private def residuePhase (F : ℝ → ℝ) (N : ℝ) (r : ℕ) : ℝ → ℝ :=
  fun u ↦ F ((1 - (r : ℝ) / N) * u + (r : ℝ) / N)

private lemma residuePhase_identity (F : ℝ → ℝ) (N : ℝ) (q r m : ℕ)
    (hN : N ≠ 0) (hrN : (r : ℝ) < N) :
    residuePhase F N r ((m : ℝ) / residueScale N q r) =
      F (((q * m + r : ℕ) : ℝ) / N) := by
  simp only [residuePhase, residueScale]
  push_cast
  field_simp [sub_ne_zero.mpr (ne_of_gt hrN)]

private lemma exponentialSumAt_residue (F : ℝ → ℝ) (T N : ℝ) (a b q r : ℕ)
    (hq : 0 < q) (hr : r < q) (hra : r ≤ a) (hab : a ≤ b)
    (hN : N ≠ 0) (hrN : (r : ℝ) < N) :
    (∑ n ∈ Finset.Icc a b with n % q = r, oscillatory F T N n) =
      exponentialSumAt (residuePhase F N r) T (residueScale N q r)
        ((a - r) ⌈/⌉ q) ((b - r) / q) := by
  rw [exponentialSumAt]
  classical
  symm
  apply Finset.sum_bij (fun m _ ↦ q * m + r)
  · intro m hm
    simp only [Finset.mem_Icc] at hm
    simp only [Finset.mem_filter, Finset.mem_Icc]
    constructor
    · constructor
      · have ha : a - r ≤ q * m :=
          (ceilDiv_le_iff_le_mul hq).mp hm.1
        omega
      · have hb : q * m ≤ b - r := by
          simpa [Nat.mul_comm] using (Nat.le_div_iff_mul_le hq).mp hm.2
        omega
    · simp [Nat.add_mod, Nat.mod_eq_of_lt hr]
  · intro m₁ hm₁ m₂ hm₂ heq
    exact Nat.eq_of_mul_eq_mul_left hq (Nat.add_right_cancel heq)
  · intro n hn
    simp only [Finset.mem_filter, Finset.mem_Icc] at hn
    have hnrepr : n = q * (n / q) + r := by
      have hnmod := hn.2
      have hdecomp := Nat.mod_add_div n q
      omega
    refine ⟨n / q, ?_, ?_⟩
    · simp only [Finset.mem_Icc]
      constructor
      · apply (ceilDiv_le_iff_le_mul hq).mpr
        omega
      · apply (Nat.le_div_iff_mul_le hq).mpr
        rw [Nat.mul_comm]
        omega
    · exact hnrepr.symm
  · intro m hm
    simp only [oscillatory]
    rw [residuePhase_identity F N q r m hN hrN]

private lemma residue_upper_endpoint_le (N : ℝ) (b q r : ℕ)
    (hq : 0 < q) (hr : r < q) (hb : (b : ℝ) ≤ 2 * N) :
    (b - r) / q ≤ ⌊2 * residueScale N q r⌋₊ + 1 := by
  by_cases hrb : r ≤ b
  · have hdiv : (((b - r) / q : ℕ) : ℝ) < (2 * residueScale N q r) + 1 := by
      have hnatdiv : (((b - r) / q : ℕ) : ℝ) ≤ ((b - r : ℕ) : ℝ) / q :=
        Nat.cast_div_le
      refine lt_of_le_of_lt hnatdiv ?_
      simp only [residueScale]
      rw [Nat.cast_sub hrb]
      have hqR : (0 : ℝ) < q := by exact_mod_cast hq
      have hrR : (r : ℝ) < q := by exact_mod_cast hr
      rw [div_lt_iff₀ hqR]
      field_simp [show (q : ℝ) ≠ 0 by positivity]
      nlinarith
    have hfloor := Nat.lt_floor_add_one (2 * residueScale N q r)
    have hcast : (((b - r) / q : ℕ) : ℝ) <
        ((⌊2 * residueScale N q r⌋₊ + 2 : ℕ) : ℝ) := by
      push_cast
      linarith
    have hnat : (b - r) / q < ⌊2 * residueScale N q r⌋₊ + 2 := by
      exact_mod_cast hcast
    omega
  · have hbr : b ≤ r := by omega
    simp [Nat.sub_eq_zero_of_le hbr]

private lemma norm_exponentialSumAt_le_sum_residues
    (F : ℝ → ℝ) (T N : ℝ) (a b q : ℕ)
    (hq : 0 < q) (hqN : (q : ℝ) ≤ N) (hN : N ≠ 0)
    (ha : N ≤ (a : ℝ)) (hb : (b : ℝ) ≤ 2 * N) :
    ‖exponentialSumAt F T N a b‖ ≤
      ∑ r ∈ Finset.range q,
        (‖exponentialSumAt (residuePhase F N r) T (residueScale N q r)
          ((a - r) ⌈/⌉ q)
          (min ((b - r) / q) ⌊2 * residueScale N q r⌋₊)‖ + 1) := by
  by_cases hab : a ≤ b
  swap
  · rw [exponentialSumAt, Finset.Icc_eq_empty hab]
    simp only [Finset.sum_empty, norm_zero]
    apply Finset.sum_nonneg
    intro r _
    exact add_nonneg (norm_nonneg (exponentialSumAt
      (residuePhase F N r) T (residueScale N q r)
      ((a - r) ⌈/⌉ q) (min ((b - r) / q) ⌊2 * residueScale N q r⌋₊))) zero_le_one
  rw [exponentialSumAt]
  classical
  have hpartition :
      ∑ n ∈ Finset.Icc a b, oscillatory F T N n =
        ∑ r ∈ Finset.range q,
          ∑ n ∈ Finset.Icc a b with n % q = r, oscillatory F T N n := by
    rw [Finset.sum_fiberwise_of_maps_to]
    intro n hn
    exact Finset.mem_range.mpr (Nat.mod_lt n hq)
  rw [hpartition]
  calc
    ‖∑ r ∈ Finset.range q, ∑ n ∈ Finset.Icc a b with n % q = r,
        oscillatory F T N n‖ ≤
        ∑ r ∈ Finset.range q,
          ‖∑ n ∈ Finset.Icc a b with n % q = r, oscillatory F T N n‖ :=
      norm_sum_le _ _
    _ ≤ _ := by
      apply Finset.sum_le_sum
      intro r hr
      have hrq := Finset.mem_range.mp hr
      have hrN : (r : ℝ) < N := lt_of_lt_of_le (by exact_mod_cast hrq) hqN
      have hra : r ≤ a := by
        have : (r : ℝ) < a := hrN.trans_le ha
        exact_mod_cast this.le
      rw [exponentialSumAt_residue F T N a b q r hq hrq hra hab hN hrN]
      let A := (a - r) ⌈/⌉ q
      let B := (b - r) / q
      let B' := min B ⌊2 * residueScale N q r⌋₊
      have hBB' : B ≤ B' + 1 := by
        dsimp [B, B']
        have hend := residue_upper_endpoint_le N b q r hq hrq hb
        rw [min_def]
        split_ifs <;> omega
      rw [exponentialSumAt]
      change ‖∑ n ∈ Finset.Icc A B,
          oscillatory (residuePhase F N r) T (residueScale N q r) n‖ ≤
        ‖∑ n ∈ Finset.Icc A B',
          oscillatory (residuePhase F N r) T (residueScale N q r) n‖ + 1
      by_cases hAB0 : A ≤ B
      swap
      · rw [Finset.Icc_eq_empty hAB0]
        simp only [Finset.sum_empty, norm_zero]
        positivity
      by_cases hle : B ≤ B'
      · have : B = B' := le_antisymm hle (min_le_left _ _)
        rw [this]
        exact le_add_of_nonneg_right (by norm_num)
      · have hsucc : B = B' + 1 := by omega
        by_cases hAB : A ≤ B'
        · rw [hsucc, Finset.sum_Icc_succ_top (hAB.trans (Nat.le_succ _))]
          calc
            _ ≤ ‖∑ n ∈ Finset.Icc A B',
                  oscillatory (residuePhase F N r) T (residueScale N q r) n‖ +
                ‖oscillatory (residuePhase F N r) T (residueScale N q r)
                  ((B' + 1 : ℕ) : ℝ)‖ :=
              norm_add_le _ _
            _ = _ := by
              rw [norm_oscillatory]
        · have hAeq : A = B' + 1 := by omega
          rw [hAeq, hsucc]
          simp [norm_oscillatory]

/-! ## Integer power scales -/

private def floorRpow (T : VariableObject ℝ) (κ : ℝ) : VariableObject ℕ :=
  fun i ↦ max 1 ⌊T i ^ κ⌋₊

private lemma floorRpow_pos (T : VariableObject ℝ) (κ : ℝ) (i : ℕ) :
    0 < floorRpow T κ i := by
  simp [floorRpow]

private lemma tendsto_floorRpow_cast_div_rpow_one
    {T : VariableObject ℝ} {κ : ℝ} (hκ : 0 < κ)
    (hT : Tendsto T atTop atTop) :
    Tendsto (fun i ↦ (floorRpow T κ i : ℝ) / T i ^ κ) atTop (𝓝 1) := by
  have hpow : Tendsto (fun i ↦ T i ^ κ) atTop atTop :=
    (tendsto_rpow_atTop hκ).comp hT
  have hlarge : ∀ᶠ i in atTop, 1 ≤ T i ^ κ := hpow.eventually (eventually_ge_atTop 1)
  have hfloor : ∀ᶠ i in atTop,
      floorRpow T κ i = ⌊T i ^ κ⌋₊ := by
    filter_upwards [hlarge] with i hi
    exact max_eq_right ((Nat.one_le_floor_iff _).mpr hi)
  apply (tendsto_nat_floor_div_atTop.comp hpow).congr'
  filter_upwards [hfloor] with i hi
  simp only [Function.comp_apply]
  rw [hi]

private lemma isPowerAsymptotic_floorRpow
    {T : VariableObject ℝ} {κ : ℝ} (hκ : 0 < κ)
    (hT : ∀ i, 1 ≤ T i) (hTunbounded : T.IsUnbounded) :
    IsPowerAsymptotic (fun i ↦ (floorRpow T κ i : ℝ)) T κ := by
  have hTtop : Tendsto T atTop atTop :=
    (VariableObject.isUnbounded_iff_tendsto_atTop
      (fun i ↦ zero_le_one.trans (hT i))).mp hTunbounded
  have hratio := tendsto_floorRpow_cast_div_rpow_one hκ hTtop
  have hTgt : ∀ᶠ i in atTop, 1 < T i := hTtop.eventually (eventually_gt_atTop 1)
  have hqpos : ∀ᶠ i in atTop, 0 < (floorRpow T κ i : ℝ) :=
    Filter.Eventually.of_forall fun i ↦ by exact_mod_cast floorRpow_pos T κ i
  apply isPowerAsymptotic_of_logb_tendsto hTgt hqpos
  have hlogratio : Tendsto (fun i ↦
      Real.log ((floorRpow T κ i : ℝ) / T i ^ κ)) atTop (𝓝 0) := by
    simpa [Function.comp_def] using
      (Real.continuousAt_log one_ne_zero).tendsto.comp hratio
  have hlogTtop := Real.tendsto_log_atTop.comp hTtop
  have hsmall := hlogratio.div_atTop hlogTtop
  have hlimit : Tendsto (fun i ↦ κ +
      Real.log ((floorRpow T κ i : ℝ) / T i ^ κ) / (Real.log ∘ T) i)
      atTop (𝓝 κ) := by
    simpa using (tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ κ) atTop (𝓝 κ)).add hsmall
  apply hlimit.congr'
  filter_upwards [hTgt] with i hi
  rw [Real.logb]
  have hTi : 0 < T i := zero_lt_one.trans hi
  have hqi : 0 < (floorRpow T κ i : ℝ) := by exact_mod_cast floorRpow_pos T κ i
  rw [Real.log_div hqi.ne' (Real.rpow_pos_of_pos hTi κ).ne', Real.log_rpow hTi]
  dsimp only [Function.comp_apply]
  field_simp [(Real.log_pos hi).ne']
  ring

private lemma IsPowerAsymptotic.mul
    {A B T : VariableObject ℝ} {α β : ℝ}
    (hA : IsPowerAsymptotic A T α) (hB : IsPowerAsymptotic B T β)
    (hT : ∀ i, 1 ≤ T i) :
    IsPowerAsymptotic (A * B) T (α + β) := by
  obtain ⟨a, ha, hAeq⟩ := hA
  obtain ⟨b, hb, hBeq⟩ := hB
  refine ⟨a + b, ?_, ?_⟩
  · have hab := ha.add hb
    convert hab using 1
    ext i
    rfl
  · filter_upwards [hAeq, hBeq] with i hiA hiB
    simp only [Pi.mul_apply, hiA, hiB]
    simpa only [Pi.add_apply] using
      (Real.rpow_add (zero_lt_one.trans_le (hT i)) (a i) (b i)).symm

private lemma IsPowerAsymptotic.div
    {A B T : VariableObject ℝ} {α β : ℝ}
    (hA : IsPowerAsymptotic A T α) (hB : IsPowerAsymptotic B T β)
    (hT : ∀ i, 1 ≤ T i) :
    IsPowerAsymptotic (A / B) T (α - β) := by
  obtain ⟨a, ha, hAeq⟩ := hA
  obtain ⟨b, hb, hBeq⟩ := hB
  refine ⟨a - b, ?_, ?_⟩
  · have hab := ha.sub hb
    convert hab using 1
    ext i
    rfl
  · filter_upwards [hAeq, hBeq] with i hiA hiB
    simp only [Pi.div_apply, hiA, hiB]
    simpa only [Pi.sub_apply] using
      (Real.rpow_sub (zero_lt_one.trans_le (hT i)) (a i) (b i)).symm

private lemma IsPowerAsymptotic.tendsto_zero_of_neg
    {X T : VariableObject ℝ} {α : ℝ} (hX : IsPowerAsymptotic X T α)
    (hα : α < 0) (hT : ∀ i, 1 ≤ T i) (hTunbounded : T.IsUnbounded)
    (hXnonneg : ∀ i, 0 ≤ X i) : Tendsto X atTop (𝓝 0) := by
  have hTtop : Tendsto T atTop atTop :=
    (VariableObject.isUnbounded_iff_tendsto_atTop
      (fun i ↦ zero_le_one.trans (hT i))).mp hTunbounded
  have hbetween := hX.eventually_between (Filter.Eventually.of_forall hT)
    (show 0 < -α / 2 by linarith)
  have hpowzero : Tendsto (fun i ↦ T i ^ (α + -α / 2)) atTop (𝓝 0) := by
    have : 0 < -(α + -α / 2) := by linarith
    simpa [Function.comp_def, neg_neg] using
      (tendsto_rpow_neg_atTop this).comp hTtop
  apply squeeze_zero' (Filter.Eventually.of_forall hXnonneg)
  · filter_upwards [hbetween] with i hi
    exact hi.2
  · exact hpowzero

private lemma isPowerAsymptotic_zero_of_tendsto_one
    {X T : VariableObject ℝ} (hX : Tendsto X atTop (𝓝 1))
    (hT : ∀ i, 1 ≤ T i) (hTunbounded : T.IsUnbounded) :
    IsPowerAsymptotic X T 0 := by
  have hTtop : Tendsto T atTop atTop :=
    (VariableObject.isUnbounded_iff_tendsto_atTop
      (fun i ↦ zero_le_one.trans (hT i))).mp hTunbounded
  have hTgt : ∀ᶠ i in atTop, 1 < T i := hTtop.eventually (eventually_gt_atTop 1)
  have hXpos : ∀ᶠ i in atTop, 0 < X i := hX.eventually (eventually_gt_nhds zero_lt_one)
  apply isPowerAsymptotic_of_logb_tendsto hTgt hXpos
  have hlogX : Tendsto (fun i ↦ Real.log (X i)) atTop (𝓝 0) := by
    simpa [Function.comp_def] using
      (Real.continuousAt_log one_ne_zero).tendsto.comp hX
  have hlogTtop := Real.tendsto_log_atTop.comp hTtop
  simpa [Real.logb] using hlogX.div_atTop hlogTtop

private lemma eventually_isApproximateModelPhaseFunction
    {F : VariableFunction (VariableObject.fixed ℝ) ℝ} {σ : ℝ}
    (hphase : IsPhaseFunction F)
    (herror : ∀ p : ℕ, (modelPhaseError F σ p).IsPointwiseInfinitesimal)
    (P : ℕ) {δ : ℝ} (hδ : 0 < δ) :
    ∀ᶠ i in atTop, IsApproximateModelPhaseFunction (F i) σ P δ := by
  have hall : ∀ p ∈ Finset.range (P + 1), ∀ᶠ i in atTop,
      ∀ u : phaseInterval, ‖modelPhaseError F σ p i u‖ < δ :=
    fun p hp ↦
      (VariableFunction.isPointwiseInfinitesimal_iff_forall_pos_uniform
        (VariableObject.fixed phaseInterval)
        (fun _ ↦ ⟨1, by simp [phaseInterval]⟩)
        (modelPhaseError F σ p)).mp (herror p) δ hδ
  have hevent := (Finset.range (P + 1)).eventually_all.mpr hall
  filter_upwards [hevent] with i hi
  refine ⟨hphase i, ?_⟩
  intro p hp u
  exact (hi p (Finset.mem_range.mpr (by omega)) u).le

/-! ## Scale transference -/

/-- On scales at most one, increasing the length scale cannot improve the universal growth
exponent. -/
private theorem exponentSumGrowthExponent_monoOn_Icc :
    MonotoneOn exponentSumGrowthExponent (Set.Icc 0 1) := by
  intro α hα γ hγ hαγ
  by_cases heq : α = γ
  · simp [heq]
  have hlt : (α : ℝ) < (γ : ℝ) := by exact_mod_cast lt_of_le_of_ne hαγ heq
  rw [exponentSumGrowthExponent_le_iff]
  intro N T F a b hN hT hTunbounded hNT hF hab
  apply (isPowerBounded_iff_forall_pos _ T _ hT hTunbounded).mpr
  intro ε hε
  obtain ⟨hphase, σ, hσ, herror⟩ := hF
  obtain ⟨δ, hδ, P, hP, C, hC, hfixed⟩ :=
    (exponentSumGrowthExponent_le_iff_nonAsymptotic.mp
      (le_rfl : exponentSumGrowthExponent γ ≤ exponentSumGrowthExponent γ))
      (ε / 2) (by linarith) σ hσ
  let η : ℝ := min (δ / 4) (((γ : ℝ) - (α : ℝ)) / 2)
  have hη : 0 < η := by dsimp [η]; positivity
  have hηδ : η < δ := lt_of_le_of_lt (min_le_left _ _) (by linarith)
  let κ : ℝ := (γ : ℝ) - (α : ℝ) - η
  have hκ : 0 < κ := by
    dsimp [κ, η]
    have := min_le_right (δ / 4) (((γ : ℝ) - (α : ℝ)) / 2)
    linarith
  let q : VariableObject ℕ := floorRpow T κ
  let Q : VariableObject ℝ := fun i ↦ (q i : ℝ)
  let M : VariableObject ℝ := Q * N
  have hTtop : Tendsto T atTop atTop :=
    (VariableObject.isUnbounded_iff_tendsto_atTop
      (fun i ↦ zero_le_one.trans (hT i))).mp hTunbounded
  have hQasym : IsPowerAsymptotic Q T κ :=
    isPowerAsymptotic_floorRpow hκ hT hTunbounded
  have hMasym0 : IsPowerAsymptotic M T (κ + (α : ℝ)) := hQasym.mul hNT hT
  have hMasym : IsPowerAsymptotic M T ((γ : ℝ) - η) := by
    convert hMasym0 using 1
    dsimp [κ]
    ring
  have hMdivAsym : IsPowerAsymptotic (M / T) T ((γ : ℝ) - η - 1) := by
    apply hMasym.div
    · refine ⟨VariableObject.fixed 1, IsEqUpToInfinitesimal.refl _, ?_⟩
      exact Filter.Eventually.of_forall fun i ↦ by simp [VariableObject.fixed]
    · exact hT
  have hMdivZero : Tendsto (M / T) atTop (𝓝 0) := by
    apply hMdivAsym.tendsto_zero_of_neg
    · have hγone : (γ : ℝ) ≤ 1 := by exact_mod_cast hγ.2
      linarith
    · exact hT
    · exact hTunbounded
    · intro i
      dsimp [M, Q, q]
      exact div_nonneg
        (mul_nonneg (Nat.cast_nonneg _) (zero_le_one.trans (hN i)))
        (zero_le_one.trans (hT i))
  classical
  have hmax (i : ℕ) : ∃ h : ZMod (q i), ∀ k : ZMod (q i),
      ‖exponentialSumAt
          (fun u ↦ F i u + (k.val : ℝ) * N i / T i * u)
          (T i) (M i) (q i * a i) (q i * b i)‖ ≤
      ‖exponentialSumAt
          (fun u ↦ F i u + (h.val : ℝ) * N i / T i * u)
          (T i) (M i) (q i * a i) (q i * b i)‖ := by
    letI : NeZero (q i) := ⟨(floorRpow_pos T κ i).ne'⟩
    obtain ⟨h, hh, hmax⟩ := Finset.exists_max_image Finset.univ
      (fun k : ZMod (q i) ↦
        ‖exponentialSumAt
          (fun u ↦ F i u + (k.val : ℝ) * N i / T i * u)
          (T i) (M i) (q i * a i) (q i * b i)‖) Finset.univ_nonempty
    exact ⟨h, fun k ↦ hmax k (Finset.mem_univ k)⟩
  let h : ∀ i, ZMod (q i) := fun i ↦ (hmax i).choose
  let c : VariableObject ℝ := fun i ↦ (h i).val * N i / T i
  let G : VariableFunction (VariableObject.fixed ℝ) ℝ := fun i u ↦ F i u + c i * u
  have hcBound : ∀ i, ‖c i‖ ≤ M i / T i := by
    intro i
    letI : NeZero (q i) := ⟨(floorRpow_pos T κ i).ne'⟩
    have hval : ((h i).val : ℝ) ≤ q i := by exact_mod_cast (h i).val_lt.le
    dsimp [c, M, Q]
    rw [abs_of_nonneg]
    · exact div_le_div_of_nonneg_right
        (mul_le_mul_of_nonneg_right hval (zero_le_one.trans (hN i)))
        (zero_le_one.trans (hT i))
    · exact div_nonneg
        (mul_nonneg (Nat.cast_nonneg _) (zero_le_one.trans (hN i)))
        (zero_le_one.trans (hT i))
  have hcInf : c.IsInfinitesimal := by
    rw [VariableObject.IsInfinitesimal]
    exact squeeze_zero' (Filter.Eventually.of_forall fun i ↦ norm_nonneg (c i))
      (Filter.Eventually.of_forall hcBound) hMdivZero
  have hGdata : IsPhaseFunction G ∧
      ∀ p : ℕ, (modelPhaseError G σ p).IsPointwiseInfinitesimal := by
    simpa [G] using modelPhaseData_add_infinitesimal_linear hphase herror hcInf
  have hG : IsModelPhaseFunction G := ⟨hGdata.1, σ, hσ, hGdata.2⟩
  have hGapprox : ∀ᶠ i in atTop, IsApproximateModelPhaseFunction (G i) σ P δ :=
    eventually_isApproximateModelPhaseFunction hGdata.1 hGdata.2 P hδ
  have hMbetween := hMasym.eventually_between
    (Filter.Eventually.of_forall hT) (sub_pos.mpr hηδ)
  have hthreshold : ∀ᶠ i in atTop, C ≤ T i :=
    hTtop.eventually (eventually_ge_atTop C)
  refine Asymptotics.IsBigO.of_bound C ?_
  filter_upwards [hGapprox, hMbetween, hthreshold] with i hiG hiM hiC
  have hqi : 0 < q i := floorRpow_pos T κ i
  letI : NeZero (q i) := ⟨hqi.ne'⟩
  have hMi : 1 ≤ M i := by
    dsimp [M, Q]
    have hqone : (1 : ℝ) ≤ q i := by exact_mod_cast (Nat.succ_le_iff.mp hqi)
    nlinarith [hN i]
  have hMbounds : T i ^ ((γ : ℝ) - δ) ≤ M i ∧ M i ≤ T i ^ ((γ : ℝ) + δ) := by
    constructor
    · convert hiM.1 using 1
      congr 1
      ring
    · exact hiM.2.trans (Real.rpow_le_rpow_of_exponent_le (hT i) (by linarith))
  have hinterval : M i ≤ ((q i * a i : ℕ) : ℝ) ∧
      ((q i * b i : ℕ) : ℝ) ≤ 2 * M i := by
    dsimp [M, Q]
    push_cast
    constructor <;> nlinarith [hqi.le, (hab i).1, (hab i).2]
  have hselected := hfixed (T i) (M i) (G i) (q i * a i) (q i * b i)
    ⟨hiC, hMbounds.1, hMbounds.2, hiG, hinterval.1, hinterval.2⟩
  have hdilate := norm_exponentialSumAt_dilate_le (F i) (T i) (N i)
    (a i) (b i) (q i) hqi (ne_of_gt (zero_lt_one.trans_le (hT i)))
      (ne_of_gt (zero_lt_one.trans_le (hN i)))
  calc
    ‖exponentialSum F T N a b i‖ ≤
        ‖exponentialSumAt (G i) (T i) (M i) (q i * a i) (q i * b i)‖ := by
      rw [exponentialSum_apply]
      refine hdilate.trans ?_
      have hall : ∀ k : ZMod (q i),
          ‖exponentialSumAt
              (fun u ↦ F i u + (k.val : ℝ) * N i / T i * u)
              (T i) (M i) (q i * a i) (q i * b i)‖ ≤
            ‖exponentialSumAt (G i) (T i) (M i) (q i * a i) (q i * b i)‖ := by
        intro k
        simpa [G, c, h] using (hmax i).choose_spec k
      calc
        _ ≤ ∑ _k : ZMod (q i),
            ‖exponentialSumAt (G i) (T i) (M i) (q i * a i) (q i * b i)‖ / q i := by
          apply Finset.sum_le_sum
          intro k _
          exact div_le_div_of_nonneg_right (hall k) (Nat.cast_nonneg _)
        _ = _ := by
          rw [Finset.sum_const, show Finset.univ.card = q i from ZMod.card (q i),
            nsmul_eq_mul]
          field_simp [show (q i : ℝ) ≠ 0 by exact_mod_cast hqi.ne']
    _ ≤ C * T i ^ (exponentSumGrowthExponent γ + ε / 2) := hselected
    _ ≤ C * ‖T i ^ (exponentSumGrowthExponent γ + ε)‖ := by
      rw [Real.norm_eq_abs,
        abs_of_nonneg (Real.rpow_nonneg (zero_le_one.trans (hT i)) _)]
      exact mul_le_mul_of_nonneg_left
        (Real.rpow_le_rpow_of_exponent_le (hT i) (by linarith))
        (zero_le_one.trans hC)

/-- Increasing the scale by `γ - α` costs at most that amount in the growth exponent. -/
private theorem exponentSumGrowthExponent_le_add_sub
    {α γ : ℝ≥0} (hαγ : α ≤ γ) :
    exponentSumGrowthExponent γ ≤
      exponentSumGrowthExponent α + (γ : ℝ) - (α : ℝ) := by
  by_cases heq : α = γ
  · subst γ; simp
  by_cases hαzero : α = 0
  · subst α
    simpa using (exponentSumGrowthExponent_le_iff.mpr (isExponentSumBound_self γ))
  have hαpos : 0 < (α : ℝ) := by positivity
  have hgap : 0 < (γ : ℝ) - (α : ℝ) := by
    have hlt : (α : ℝ) < (γ : ℝ) := by
      exact_mod_cast lt_of_le_of_ne hαγ heq
    linarith
  rw [exponentSumGrowthExponent_le_iff]
  intro N T F a b hN hT hTunbounded hNT hF hab
  apply (isPowerBounded_iff_forall_pos _ T _ hT hTunbounded).mpr
  intro ε hε
  obtain ⟨hphase, σ, hσ, herror⟩ := hF
  obtain ⟨δ, hδ, P, hP, C, hC, hfixed⟩ :=
    (exponentSumGrowthExponent_le_iff_nonAsymptotic.mp
      (le_rfl : exponentSumGrowthExponent α ≤ exponentSumGrowthExponent α))
      (ε / 4) (by linarith) σ hσ
  let η : ℝ := min (δ / 4) (min (ε / 8) ((α : ℝ) / 4))
  have hη : 0 < η := by dsimp [η]; positivity
  have hηδ : η < δ := lt_of_le_of_lt (min_le_left _ _) (by linarith)
  have hηε : η ≤ ε / 8 :=
    (min_le_right (δ / 4) _).trans (min_le_left _ _)
  have hηα : η < (α : ℝ) := by
    have hηquarter : η ≤ (α : ℝ) / 4 :=
      (min_le_right (δ / 4) (min (ε / 8) ((α : ℝ) / 4))).trans
        (min_le_right (ε / 8) ((α : ℝ) / 4))
    linarith
  let κ : ℝ := (γ : ℝ) - (α : ℝ) + η
  have hκ : 0 < κ := by dsimp [κ]; linarith
  let q : VariableObject ℕ := floorRpow T κ
  let Q : VariableObject ℝ := fun i ↦ (q i : ℝ)
  have hTtop : Tendsto T atTop atTop :=
    (VariableObject.isUnbounded_iff_tendsto_atTop
      (fun i ↦ zero_le_one.trans (hT i))).mp hTunbounded
  have hQasym : IsPowerAsymptotic Q T κ :=
    isPowerAsymptotic_floorRpow hκ hT hTunbounded
  have hQdivNasym : IsPowerAsymptotic (Q / N) T (η - (α : ℝ)) := by
    have := hQasym.div hNT hT
    convert this using 1
    dsimp [κ]
    ring
  have hQdivNzero : Tendsto (Q / N) atTop (𝓝 0) := by
    apply hQdivNasym.tendsto_zero_of_neg
    · linarith
    · exact hT
    · exact hTunbounded
    · intro i
      exact div_nonneg (Nat.cast_nonneg _) (zero_le_one.trans (hN i))
  have hqNevent : ∀ᶠ i in atTop, (q i : ℝ) ≤ N i := by
    have := (Metric.tendsto_nhds.1 hQdivNzero) 1 zero_lt_one
    filter_upwards [this] with i hi
    rw [dist_zero_right, Real.norm_eq_abs] at hi
    change |(q i : ℝ) / N i| < 1 at hi
    rw [
      abs_of_nonneg (div_nonneg (Nat.cast_nonneg _) (zero_le_one.trans (hN i)))] at hi
    have hNpos : 0 < N i := zero_lt_one.trans_le (hN i)
    exact (div_le_one hNpos).mp hi.le
  classical
  let A : ℕ → ℕ → ℕ := fun i r ↦ (a i - r) ⌈/⌉ q i
  let B : ℕ → ℕ → ℕ := fun i r ↦
    min ((b i - r) / q i) ⌊2 * residueScale (N i) (q i) r⌋₊
  have hmax (i : ℕ) : ∃ r ∈ Finset.range (q i), ∀ s ∈ Finset.range (q i),
      ‖exponentialSumAt (residuePhase (F i) (N i) s) (T i)
          (residueScale (N i) (q i) s) (A i s) (B i s)‖ ≤
      ‖exponentialSumAt (residuePhase (F i) (N i) r) (T i)
          (residueScale (N i) (q i) r) (A i r) (B i r)‖ :=
    Finset.exists_max_image (Finset.range (q i))
      (fun r ↦ ‖exponentialSumAt (residuePhase (F i) (N i) r) (T i)
        (residueScale (N i) (q i) r) (A i r) (B i r)‖)
      ⟨0, Finset.mem_range.mpr (floorRpow_pos T κ i)⟩
  let r : VariableObject ℕ := fun i ↦
    if (q i : ℝ) ≤ N i then (hmax i).choose else 0
  have hrq (i : ℕ) : r i < q i := by
    dsimp [r]
    split_ifs
    · exact Finset.mem_range.mp (hmax i).choose_spec.1
    · exact floorRpow_pos T κ i
  have hrN (i : ℕ) : (r i : ℝ) < N i := by
    dsimp [r]
    split_ifs with hi
    · exact lt_of_lt_of_le (by exact_mod_cast Finset.mem_range.mp (hmax i).choose_spec.1) hi
    · simpa only [Nat.cast_zero] using zero_lt_one.trans_le (hN i)
  let d : VariableObject ℝ := fun i ↦ (r i : ℝ) / N i
  let c : VariableObject ℝ := fun i ↦ 1 - d i
  have hdBound : ∀ i, 0 ≤ d i ∧ d i ≤ Q i / N i := by
    intro i
    constructor
    · dsimp [d]
      exact div_nonneg (Nat.cast_nonneg _) (zero_le_one.trans (hN i))
    · dsimp [d, Q]
      exact div_le_div_of_nonneg_right (by exact_mod_cast (hrq i).le)
        (zero_le_one.trans (hN i))
  have hdZero : Tendsto d atTop (𝓝 0) :=
    squeeze_zero' (Filter.Eventually.of_forall fun i ↦ (hdBound i).1)
      (Filter.Eventually.of_forall fun i ↦ (hdBound i).2) hQdivNzero
  have hcOne : Tendsto c atTop (𝓝 1) := by
    simpa [c] using tendsto_const_nhds.sub hdZero
  let M : VariableObject ℝ := fun i ↦ residueScale (N i) (q i) (r i)
  have hNasymDiv : IsPowerAsymptotic (N / Q) T ((γ : ℝ) - κ) := hNT.div hQasym hT
  have hcasym : IsPowerAsymptotic c T 0 :=
    isPowerAsymptotic_zero_of_tendsto_one hcOne hT hTunbounded
  have hMformula : M = (N / Q) * c := by
    funext i
    dsimp [M, residueScale, Q, c, d]
    field_simp [show (q i : ℝ) ≠ 0 by exact_mod_cast (floorRpow_pos T κ i).ne',
      show N i ≠ 0 by exact (ne_of_gt (zero_lt_one.trans_le (hN i)))]
  have hMasym0 : IsPowerAsymptotic M T (((γ : ℝ) - κ) + 0) := by
    rw [hMformula]
    exact hNasymDiv.mul hcasym hT
  have hMasym : IsPowerAsymptotic M T ((α : ℝ) - η) := by
    convert hMasym0 using 1
    dsimp [κ]
    ring
  let G : VariableFunction (VariableObject.fixed ℝ) ℝ := fun i ↦ residuePhase (F i) (N i) (r i)
  have hmap (i : ℕ) : Set.MapsTo (fun u ↦ c i * u + d i) phaseInterval phaseInterval := by
    intro u hu
    have hdI : d i ∈ Set.Ico (0 : ℝ) 1 := by
      exact ⟨(hdBound i).1,
        (div_lt_one (zero_lt_one.trans_le (hN i))).mpr (hrN i)⟩
    rcases hu with ⟨hu1, hu2⟩
    dsimp [c, phaseInterval] at *
    constructor <;> nlinarith [hdI.1, hdI.2]
  have hGdata : IsPhaseFunction G ∧
      ∀ p : ℕ, (modelPhaseError G σ p).IsPointwiseInfinitesimal := by
    have hdata := modelPhaseData_comp_affine_tendsto_id hphase herror hcOne hdZero hmap
    change IsPhaseFunction
        (fun i u ↦ F i ((1 - (r i : ℝ) / N i) * u + (r i : ℝ) / N i)) ∧
      ∀ p : ℕ, (modelPhaseError
        (fun i u ↦ F i ((1 - (r i : ℝ) / N i) * u + (r i : ℝ) / N i)) σ p
          ).IsPointwiseInfinitesimal
    simpa only [c, d] using hdata
  have hG : IsModelPhaseFunction G := ⟨hGdata.1, σ, hσ, hGdata.2⟩
  have hGapprox : ∀ᶠ i in atTop, IsApproximateModelPhaseFunction (G i) σ P δ :=
    eventually_isApproximateModelPhaseFunction hGdata.1 hGdata.2 P hδ
  have hMbetween := hMasym.eventually_between
    (Filter.Eventually.of_forall hT) (sub_pos.mpr hηδ)
  have hMtop : Tendsto M atTop atTop := by
    have hpowtop : Tendsto (fun i ↦ T i ^ (((α : ℝ) - η) / 2)) atTop atTop :=
      (tendsto_rpow_atTop (by linarith)).comp hTtop
    apply tendsto_atTop_mono' atTop _ hpowtop
    have hb := hMasym.eventually_between (Filter.Eventually.of_forall hT)
      (show 0 < ((α : ℝ) - η) / 2 by linarith)
    filter_upwards [hb] with i hi
    have hexp : ((α : ℝ) - η) / 2 =
        (α : ℝ) - η - ((α : ℝ) - η) / 2 := by ring
    rw [hexp]
    exact hi.1
  have hMone : ∀ᶠ i in atTop, 1 ≤ M i := hMtop.eventually (eventually_ge_atTop 1)
  have hthreshold : ∀ᶠ i in atTop, C ≤ T i := hTtop.eventually (eventually_ge_atTop C)
  have hqUpper : ∀ᶠ i in atTop, (q i : ℝ) ≤ T i ^ κ := by
    have hpowlarge : ∀ᶠ i in atTop, 1 ≤ T i ^ κ :=
      ((tendsto_rpow_atTop hκ).comp hTtop).eventually (eventually_ge_atTop 1)
    filter_upwards [hpowlarge] with i hi
    dsimp [q, floorRpow]
    rw [max_eq_right ((Nat.one_le_floor_iff _).mpr hi)]
    exact Nat.floor_le (Real.rpow_nonneg (zero_le_one.trans (hT i)) _)
  refine Asymptotics.IsBigO.of_bound (2 * C) ?_
  filter_upwards [hqNevent, hGapprox, hMbetween, hMone, hthreshold, hqUpper]
    with i hiqN hiG hiM hiMone hiC hiqUpper
  have hqi := floorRpow_pos T κ i
  have hMbounds : T i ^ ((α : ℝ) - δ) ≤ M i ∧ M i ≤ T i ^ ((α : ℝ) + δ) := by
    constructor
    · have hexp : (α : ℝ) - δ = (α : ℝ) - η - (δ - η) := by ring
      rw [hexp]
      exact hiM.1
    · exact hiM.2.trans (Real.rpow_le_rpow_of_exponent_le (hT i) (by linarith))
  have hA : M i ≤ (A i (r i) : ℝ) := by
    dsimp [M, residueScale, A]
    have hceil : a i - r i ≤ q i * ((a i - r i) ⌈/⌉ q i) :=
      le_smul_ceilDiv hqi
    have hra : r i ≤ a i := by
      exact_mod_cast (hrN i).le.trans (hab i).1
    have hceilR : ((a i - r i : ℕ) : ℝ) ≤
        (q i : ℝ) * (((a i - r i) ⌈/⌉ q i : ℕ) : ℝ) := by
      exact_mod_cast hceil
    rw [Nat.cast_sub hra] at hceilR
    rw [div_le_iff₀ (show (0 : ℝ) < q i by exact_mod_cast hqi)]
    nlinarith [(hab i).1]
  have hB : (B i (r i) : ℝ) ≤ 2 * M i := by
    dsimp [B, M]
    calc
      ((min ((b i - r i) / q i)
          ⌊2 * residueScale (N i) (q i) (r i)⌋₊ : ℕ) : ℝ) ≤
          (⌊2 * residueScale (N i) (q i) (r i)⌋₊ : ℕ) := by
        exact_mod_cast min_le_right ((b i - r i) / q i)
          ⌊2 * residueScale (N i) (q i) (r i)⌋₊
      _ ≤ 2 * residueScale (N i) (q i) (r i) :=
        Nat.floor_le (by linarith [hiMone])
  have hselected := hfixed (T i) (M i) (G i) (A i (r i)) (B i (r i))
    ⟨hiC, hMbounds.1, hMbounds.2, hiG, hA, hB⟩
  have hres := norm_exponentialSumAt_le_sum_residues (F i) (T i) (N i)
    (a i) (b i) (q i) hqi hiqN
      (ne_of_gt (zero_lt_one.trans_le (hN i))) (hab i).1 (hab i).2
  have hmaxEvent : ∀ s ∈ Finset.range (q i),
      ‖exponentialSumAt (residuePhase (F i) (N i) s) (T i)
          (residueScale (N i) (q i) s) (A i s) (B i s)‖ ≤
      ‖exponentialSumAt (G i) (T i) (M i) (A i (r i)) (B i (r i))‖ := by
    intro s hs
    have hri : r i = (hmax i).choose := by simp [r, hiqN]
    change ‖exponentialSumAt (residuePhase (F i) (N i) s) (T i)
        (residueScale (N i) (q i) s) (A i s) (B i s)‖ ≤
      ‖exponentialSumAt (residuePhase (F i) (N i) (r i)) (T i)
        (residueScale (N i) (q i) (r i)) (A i (r i)) (B i (r i))‖
    rw [hri]
    exact (hmax i).choose_spec.2 s hs
  have hsum : ‖exponentialSum F T N a b i‖ ≤
      (q i : ℝ) *
        (‖exponentialSumAt (G i) (T i) (M i) (A i (r i)) (B i (r i))‖ + 1) := by
    rw [exponentialSum_apply]
    refine hres.trans ?_
    calc
      _ ≤ ∑ _s ∈ Finset.range (q i),
          (‖exponentialSumAt (G i) (T i) (M i) (A i (r i)) (B i (r i))‖ + 1) := by
        gcongr with s hs
        simpa only [A, B] using hmaxEvent s hs
      _ = _ := by simp [mul_add, mul_comm]
  have hβnonneg : 0 ≤ exponentSumGrowthExponent α :=
    (isExponentSumBound_exponentSumGrowthExponent α).nonneg
  calc
    ‖exponentialSum F T N a b i‖ ≤
        (q i : ℝ) * (C * T i ^ (exponentSumGrowthExponent α + ε / 4) + 1) := by
      exact hsum.trans (mul_le_mul_of_nonneg_left
        (add_le_add hselected le_rfl) (Nat.cast_nonneg _))
    _ ≤ 2 * C * T i ^
        (exponentSumGrowthExponent α + (γ : ℝ) - (α : ℝ) + ε) := by
      have hone : 1 ≤ C * T i ^ (exponentSumGrowthExponent α + ε / 4) := by
        calc
          1 ≤ C := hC
          _ = C * 1 := by ring
          _ ≤ _ := by gcongr; exact Real.one_le_rpow (hT i) (by linarith)
      calc
        _ ≤ (q i : ℝ) * (2 * (C * T i ^ (exponentSumGrowthExponent α + ε / 4))) := by
          gcongr; linarith
        _ ≤ T i ^ κ * (2 * (C * T i ^ (exponentSumGrowthExponent α + ε / 4))) := by
          gcongr
        _ = 2 * C *
            (T i ^ κ * T i ^ (exponentSumGrowthExponent α + ε / 4)) := by ring
        _ = 2 * C * T i ^ (κ + exponentSumGrowthExponent α + ε / 4) := by
          congr 1
          rw [← Real.rpow_add (zero_lt_one.trans_le (hT i))]
          congr 1
          ring
        _ ≤ _ := by
          have hp : T i ^ (κ + exponentSumGrowthExponent α + ε / 4) ≤
              T i ^ (exponentSumGrowthExponent α + (γ : ℝ) - (α : ℝ) + ε) := by
            apply Real.rpow_le_rpow_of_exponent_le (hT i)
            dsimp [κ]
            linarith
          exact mul_le_mul_of_nonneg_left hp
            (mul_nonneg (by norm_num) (zero_le_one.trans hC))
    _ = 2 * C * ‖T i ^
        (exponentSumGrowthExponent α + (γ : ℝ) - (α : ℝ) + ε)‖ := by
      rw [Real.norm_eq_abs,
        abs_of_nonneg (Real.rpow_nonneg (zero_le_one.trans (hT i)) _)]

/-- The growth exponent is `1`-Lipschitz on the subcritical interval. -/
private theorem lipschitzOnWith_exponentSumGrowthExponent :
    LipschitzOnWith 1 exponentSumGrowthExponent (Set.Icc 0 1) := by
  rw [lipschitzOnWith_iff_dist_le_mul]
  intro α hα γ hγ
  simp only [NNReal.coe_one, one_mul, Real.dist_eq]
  wlog hαγ : α ≤ γ generalizing α γ
  · rw [abs_sub_comm, dist_comm]
    exact this γ hγ α hα (le_of_not_ge hαγ)
  have hmono := exponentSumGrowthExponent_monoOn_Icc hα hγ hαγ
  have hslope := exponentSumGrowthExponent_le_add_sub hαγ
  rw [abs_of_nonpos (sub_nonpos.mpr hmono), NNReal.dist_eq,
    abs_of_nonpos (sub_nonpos.mpr (by exact_mod_cast hαγ))]
  linarith

/-- The exponential-sum growth exponent is upper semicontinuous. -/
theorem upperSemicontinuous_exponentSumGrowthExponent :
    UpperSemicontinuous exponentSumGrowthExponent := by
  rw [upperSemicontinuous_iff]
  intro α y hy
  by_cases hαlt : α < 1
  · let ρ : ℝ := min ((y - exponentSumGrowthExponent α) / 2)
        (((1 : ℝ) - (α : ℝ)) / 2)
    have hρ : 0 < ρ := by
      dsimp [ρ]
      apply lt_min
      · linarith
      · have : (α : ℝ) < 1 := by exact_mod_cast hαlt
        linarith
    rw [Metric.eventually_nhds_iff]
    refine ⟨ρ, hρ, ?_⟩
    intro γ hdist
    have hdist' : |(γ : ℝ) - (α : ℝ)| < ρ := by
      simpa [NNReal.dist_eq, abs_sub_comm] using hdist
    have hγone : γ ≤ 1 := by
      have hρone := min_le_right ((y - exponentSumGrowthExponent α) / 2)
        (((1 : ℝ) - (α : ℝ)) / 2)
      have hdiffUpper := (abs_lt.mp hdist').2
      exact_mod_cast (show (γ : ℝ) ≤ 1 by
        linarith)
    by_cases hγα : γ ≤ α
    · exact lt_of_le_of_lt
        (exponentSumGrowthExponent_monoOn_Icc
          ⟨zero_le, hγone⟩ ⟨zero_le, hαlt.le⟩ hγα) hy
    · have hslope := exponentSumGrowthExponent_le_add_sub (le_of_not_ge hγα)
      have hρy := min_le_left ((y - exponentSumGrowthExponent α) / 2)
        (((1 : ℝ) - (α : ℝ)) / 2)
      have hdiff : (γ : ℝ) - (α : ℝ) < ρ := (abs_lt.mp hdist').2
      linarith
  · have hαone : 1 ≤ α := le_of_not_gt hαlt
    by_cases hαeq : α = 1
    · subst α
      have hβnonneg : 0 ≤ exponentSumGrowthExponent 1 :=
        (isExponentSumBound_exponentSumGrowthExponent 1).nonneg
      let ρ : ℝ := min (y / 2) 1
      have hypos : 0 < y := lt_of_le_of_lt hβnonneg hy
      have hρ : 0 < ρ := by dsimp [ρ]; positivity
      rw [Metric.eventually_nhds_iff]
      refine ⟨ρ, hρ, ?_⟩
      intro γ hdist
      by_cases hγone : γ ≤ 1
      · exact lt_of_le_of_lt
          (exponentSumGrowthExponent_monoOn_Icc
            ⟨zero_le, hγone⟩ ⟨zero_le, le_rfl⟩ hγone) hy
      · have hγgt : 1 < γ := lt_of_not_ge hγone
        rw [exponentSumGrowthExponent_eq_sub_one hγgt]
        have hdist' : |(γ : ℝ) - 1| < ρ := by
          simpa [NNReal.dist_eq, abs_sub_comm] using hdist
        have hρy := min_le_left (y / 2) (1 : ℝ)
        have hdiffUpper := (abs_lt.mp hdist').2
        linarith
    · have hαgt : 1 < α := lt_of_le_of_ne hαone (Ne.symm hαeq)
      rw [exponentSumGrowthExponent_eq_sub_one hαgt] at hy
      let ρ : ℝ := min (((α : ℝ) - 1) / 2)
        ((y - ((α : ℝ) - 1)) / 2)
      have hρ : 0 < ρ := by
        dsimp [ρ]
        have hαgtR : (1 : ℝ) < (α : ℝ) := by exact_mod_cast hαgt
        apply lt_min <;> linarith
      rw [Metric.eventually_nhds_iff]
      refine ⟨ρ, hρ, ?_⟩
      intro γ hdist
      have hdist' : |(γ : ℝ) - (α : ℝ)| < ρ := by
        simpa [NNReal.dist_eq, abs_sub_comm] using hdist
      have hρleft := min_le_left (((α : ℝ) - 1) / 2)
        ((y - ((α : ℝ) - 1)) / 2)
      have hγgt : 1 < γ := by
        have hdiffLower := (abs_lt.mp hdist').1
        exact_mod_cast (show (1 : ℝ) < γ by
          linarith)
      rw [exponentSumGrowthExponent_eq_sub_one hγgt]
      have hρright := min_le_right (((α : ℝ) - 1) / 2)
        ((y - ((α : ℝ) - 1)) / 2)
      have hdiffUpper := (abs_lt.mp hdist').2
      linarith

end Expdb
