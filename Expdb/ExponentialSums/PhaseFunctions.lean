module

public import Expdb.Basic.Asymptotics
public import Mathlib.Analysis.Calculus.IteratedDeriv.Defs
public import Mathlib.Analysis.SpecialFunctions.Pow.Real
public import Mathlib.Order.Interval.Finset.Nat

import Expdb.Basic.AutomaticUniformity
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv

/-!
# Phase functions

This module formalizes the phase function definitions from Chapter 4 of the ANTEDB blueprint.
A phase function is smooth on `[1, 2]`. A model phase function is a variable family
of phase functions whose successive derivatives asymptotically agree with those of
`u ↦ u ^ (-σ)` for some fixed positive exponent `σ`.
-/

@[expose] public section

open Filter Topology
open scoped ContDiff

noncomputable section

namespace Expdb

/-- The interval on which phase functions are considered. -/
def phaseInterval : Set ℝ := Set.Icc 1 2

/-- A phase function is a variable real-valued function that is smooth on `[1, 2]` at every
ambient index. -/
def IsPhaseFunction (F : VariableFunction (VariableObject.fixed ℝ) ℝ) : Prop :=
  ∀ i : ℕ, ContDiffOn ℝ ∞ (F i) phaseInterval

/-- The variable error function in the model phase condition at derivative order `p`. -/
def modelPhaseError
    (F : VariableFunction (VariableObject.fixed ℝ) ℝ) (σ : ℝ) (p : ℕ) :
    VariableFunction (VariableObject.fixed phaseInterval) ℝ :=
  fun i u ↦
    iteratedDerivWithin (p + 1) (F i) phaseInterval u -
      iteratedDerivWithin p (fun v : ℝ => v ^ (-σ)) phaseInterval u

/-- A variable phase function is a model phase function when, for some fixed `σ > 0`,
the error between its `(p + 1)`st derivative and the `p`th derivative of `u ↦ u ^ (-σ)` is
pointwise infinitesimal for every fixed derivative order `p`. -/
def IsModelPhaseFunction (F : VariableFunction (VariableObject.fixed ℝ) ℝ) : Prop :=
  IsPhaseFunction F ∧
    ∃ σ : ℝ, 0 < σ ∧
      ∀ p : ℕ, (modelPhaseError F σ p).IsPointwiseInfinitesimal

/-- The fixed logarithmic phase `u ↦ log u`, regarded as a variable family. -/
def logPhase : VariableFunction (VariableObject.fixed ℝ) ℝ :=
  fun _ ↦ Real.log

/-- The derivatives of `log` on the phase interval agree with the reference derivatives for
model exponent one. -/
theorem iteratedDerivWithin_log_eq_rpow_neg_one
    (p : ℕ) (u : phaseInterval) :
    iteratedDerivWithin (p + 1) Real.log phaseInterval u =
      iteratedDerivWithin p (fun v : ℝ ↦ v ^ (-(1 : ℝ))) phaseInterval u := by
  have hu_pos : 0 < (u : ℝ) := lt_of_lt_of_le zero_lt_one u.property.1
  have hunique : UniqueDiffOn ℝ phaseInterval :=
    uniqueDiffOn_Icc (by norm_num [phaseInterval])
  rw [iteratedDerivWithin_eq_iteratedDeriv hunique
      (Real.contDiffAt_log.2 hu_pos.ne') u.property,
    iteratedDerivWithin_eq_iteratedDeriv hunique
      (Real.contDiffAt_rpow_const_of_ne hu_pos.ne') u.property,
    iteratedDeriv_succ', Real.deriv_log']
  congr 1
  funext v
  exact (Real.rpow_neg_one v).symm

/-- The fixed logarithmic phase `u ↦ log u` is a model phase function, with model exponent
`σ = 1`. -/
theorem isModelPhaseFunction_log : IsModelPhaseFunction logPhase := by
  refine ⟨?_, 1, zero_lt_one, ?_⟩
  · intro i
    change ContDiffOn ℝ ∞ Real.log phaseInterval
    exact Real.contDiffOn_log.mono fun u hu ↦ by
      simp only [Set.mem_compl_iff, Set.mem_singleton_iff]
      exact ne_of_gt (lt_of_lt_of_le zero_lt_one hu.1)
  · intro p u
    rw [VariableObject.IsInfinitesimal]
    simp only [modelPhaseError, logPhase, iteratedDerivWithin_log_eq_rpow_neg_one, sub_self,
      norm_zero]
    exact tendsto_const_nhds

/-- On a dyadic interval, `log (n / N)` is `1 / (2 * N)`-separated. -/
theorem log_div_separated
    {N : ℝ} {a b : ℕ} (hN : 0 < N) (ha : N ≤ a) (hb : (b : ℝ) ≤ 2 * N) :
    IsSeparatedFamily (1 / (2 * N))
      (fun n : ↥(Finset.Icc a b) ↦ Real.log ((n : ℝ) / N)) := by
  have hordered (m n : ↥(Finset.Icc a b)) (hmn : (m : ℕ) < n) :
      1 / (2 * N) ≤ Real.log ((n : ℝ) / N) - Real.log ((m : ℝ) / N) := by
    have hm_mem := Finset.mem_Icc.mp m.property
    have hmpos : 0 < (m : ℝ) := lt_of_lt_of_le hN (ha.trans (by exact_mod_cast hm_mem.1))
    have hnpos : 0 < (n : ℝ) := lt_trans hmpos (by exact_mod_cast hmn)
    have hlog : 1 - (m : ℝ) / n ≤ Real.log ((n : ℝ) / m) := by
      simpa [inv_div] using Real.one_sub_inv_le_log_of_pos (div_pos hnpos hmpos)
    have hdiff : 1 / (2 * N) ≤ 1 - (m : ℝ) / n := by
      have hgapNat : (m : ℕ) + 1 ≤ (n : ℕ) := Nat.add_one_le_iff.mpr hmn
      have hgapCast : ((m : ℕ) : ℝ) + 1 ≤ ((n : ℕ) : ℝ) := by exact_mod_cast hgapNat
      have hgap : (1 : ℝ) ≤ (n : ℝ) - m := by linarith
      have hn_mem := Finset.mem_Icc.mp n.property
      have hnle : (n : ℝ) ≤ 2 * N := le_trans (by exact_mod_cast hn_mem.2) hb
      rw [one_sub_div hnpos.ne']
      calc
        1 / (2 * N) ≤ 1 / (n : ℝ) := one_div_le_one_div_of_le hnpos hnle
        _ ≤ ((n : ℝ) - m) / n := (div_le_div_iff_of_pos_right hnpos).2 hgap
    calc
      1 / (2 * N) ≤ Real.log ((n : ℝ) / m) := hdiff.trans hlog
      _ = Real.log ((n : ℝ) / N) - Real.log ((m : ℝ) / N) := by
        rw [Real.log_div (ne_of_gt hnpos) hN.ne', Real.log_div (ne_of_gt hmpos) hN.ne',
          Real.log_div (ne_of_gt hnpos) (ne_of_gt hmpos)]
        ring
  intro m n hmn
  rw [Real.dist_eq]
  rcases lt_or_gt_of_ne (Subtype.coe_ne_coe.mpr hmn) with hmn' | hnm'
  · have hbound := hordered m n hmn'
    have hsign : Real.log ((m : ℝ) / N) - Real.log ((n : ℝ) / N) ≤ 0 := by
      have : 0 < 1 / (2 * N) := by positivity
      linarith
    rw [abs_of_nonpos hsign, neg_sub]
    exact hbound
  · have hbound := hordered n m hnm'
    have hsign : 0 ≤ Real.log ((m : ℝ) / N) - Real.log ((n : ℝ) / N) := by
      have : 0 < 1 / (2 * N) := by positivity
      linarith
    rw [abs_of_nonneg hsign]
    exact hbound

private lemma modelPhaseError_sum_isPointwiseInfinitesimal
    {F : VariableFunction (VariableObject.fixed ℝ) ℝ} {σ : ℝ}
    (herror : ∀ p : ℕ, (modelPhaseError F σ p).IsPointwiseInfinitesimal)
    (P : ℕ) :
    VariableFunction.IsPointwiseInfinitesimal
      ((fun i (u : phaseInterval) ↦
        ∑ p ∈ Finset.range (P + 1), ‖modelPhaseError F σ p i u‖) :
        VariableFunction (VariableObject.fixed phaseInterval) ℝ) := by
  intro u
  rw [VariableObject.IsInfinitesimal]
  have hsum : Tendsto
      (fun i ↦ ∑ p ∈ Finset.range (P + 1), ‖modelPhaseError F σ p i (u i)‖)
      atTop (nhds 0) := by
    simpa using tendsto_finsetSum (Finset.range (P + 1)) fun p _ ↦ herror p u
  simpa using hsum.norm

/-- For any fixed derivative cutoff, the errors of a model phase function have a common
infinitesimal bound, uniform over the phase interval and all derivative orders up to the cutoff,
after passing to a subsequence. -/
theorem IsModelPhaseFunction.exists_subsequence_uniform_error
    {F : VariableFunction (VariableObject.fixed ℝ) ℝ}
    (hF : IsModelPhaseFunction F) (P : ℕ) :
    ∃ σ : ℝ, 0 < σ ∧
      ∃ φ : ℕ → ℕ, StrictMono φ ∧
        ∃ c : VariableObject ℝ, c.IsInfinitesimal ∧
          ∀ i p, p ≤ P →
            ∀ u : phaseInterval,
              ‖modelPhaseError F σ p (φ i) u‖ ≤ c i := by
  rcases hF with ⟨_, σ, hσ, herror⟩
  let g : VariableFunction (VariableObject.fixed phaseInterval) ℝ :=
    fun i u ↦ ∑ p ∈ Finset.range (P + 1), ‖modelPhaseError F σ p i u‖
  have hg : g.IsPointwiseInfinitesimal :=
    modelPhaseError_sum_isPointwiseInfinitesimal herror P
  obtain ⟨φ, hφ, c, hc, hbound⟩ :=
    automatic_uniformity_of_pointwise_infinitesimal
      (E := VariableObject.fixed phaseInterval)
      (fun _ ↦ ⟨1, by simp [phaseInterval]⟩) g hg
  refine ⟨σ, hσ, φ, hφ, c, hc, ?_⟩
  intro i p hp u
  have hmem : p ∈ Finset.range (P + 1) := Finset.mem_range.mpr (by omega)
  calc
    ‖modelPhaseError F σ p (φ i) u‖ ≤
        ∑ q ∈ Finset.range (P + 1), ‖modelPhaseError F σ q (φ i) u‖ :=
      Finset.single_le_sum
        (fun q _ ↦ norm_nonneg (modelPhaseError F σ q (φ i) u)) hmem
    _ = ‖g (φ i) u‖ := by
      rw [Real.norm_eq_abs, abs_of_nonneg (Finset.sum_nonneg fun _ _ ↦ norm_nonneg _)]
    _ ≤ c i := hbound i u

end Expdb
