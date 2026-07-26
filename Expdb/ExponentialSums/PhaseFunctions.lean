module

public import Expdb.Basic.Asymptotics
public import Mathlib.Analysis.Calculus.IteratedDeriv.Defs
public import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Expdb.Basic.AutomaticUniformity
import Mathlib.Analysis.Calculus.IteratedDeriv.Lemmas
import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv

/-!
# Phase functions

This module formalizes the definitions from Section 4.1 of the ANTEDB blueprint. A phase function
is smooth on `[1, 2]`. A model phase function is a variable family of phase functions whose
successive derivatives asymptotically agree with those of `u ↦ u ^ (-σ)` for some fixed
positive exponent `σ`.
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

private lemma iteratedDerivWithin_log_eq_rpow_neg_one
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
