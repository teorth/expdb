module

public import Expdb.Basic.Asymptotics
public import Mathlib.Analysis.Calculus.IteratedDeriv.Defs
public import Mathlib.Analysis.SpecialFunctions.Pow.Real

/-!
# Phase functions

This module formalizes the definitions from Section 4.1 of the ANTEDB blueprint. A phase function
is smooth on `[1, 2]`. A model phase function is a variable family of phase functions whose
successive derivatives asymptotically agree with those of `u ↦ u ^ (-σ)` for some fixed
positive exponent `σ`.
-/

@[expose] public section

open scoped ContDiff

noncomputable section

namespace Expdb

/-- The interval on which phase functions are considered. -/
def phaseInterval : Set ℝ := Set.Icc 1 2

/-- A phase function is a variable real-valued function that is smooth on `[1, 2]` at every
ambient index. -/
def IsPhaseFunction (F : VariableFunction (fun _ ↦ ℝ) ℝ) : Prop :=
  ∀ i : ℕ, ContDiffOn ℝ ∞ (F i) phaseInterval

/-- The error in the model phase condition at derivative order `p`, family index `i`, and
point `u`. -/
def modelPhaseError
    (F : VariableFunction (fun _ ↦ ℝ) ℝ)
    (σ : ℝ) (p i : ℕ) (u : phaseInterval) : ℝ :=
  iteratedDerivWithin (p + 1) (F i) phaseInterval u -
    iteratedDerivWithin p (fun v : ℝ => v ^ (-σ)) phaseInterval u

/-- A variable phase function is a model phase function when, for some fixed `σ > 0`,
the error between its `(p + 1)`st derivative and the `p`th derivative of `u ↦ u ^ (-σ)` is
pointwise infinitesimal for every fixed derivative order `p`. -/
def IsModelPhaseFunction (F : VariableFunction (fun _ ↦ ℝ) ℝ) : Prop :=
  IsPhaseFunction F ∧
    ∃ σ : ℝ, 0 < σ ∧
      ∀ p : ℕ, VariableFunction.IsPointwiseInfinitesimal (modelPhaseError F σ p)

end Expdb
