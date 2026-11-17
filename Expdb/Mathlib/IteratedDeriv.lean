module

public import Mathlib.Analysis.Calculus.IteratedDeriv.Lemmas

/-!
# Supplementary lemmas for iterated derivatives

This file contains chain rules that are useful in the exponential-sums development but are not
currently available in Mathlib in the required two-set form.
-/

@[expose] public section

namespace Expdb

/-- Iterated derivatives of a function precomposed with an affine map.  Unlike Mathlib's
`iteratedDerivWithin_comp_const_smul`, the domain and codomain sets may be different. -/
theorem iteratedDerivWithin_comp_affine_of_mapsTo
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
          (fun y ↦ iteratedDerivWithin n f t (c * y + d)) s x := by
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

end Expdb
