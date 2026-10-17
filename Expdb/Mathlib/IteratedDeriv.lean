module

public import Mathlib.Analysis.Calculus.IteratedDeriv.Lemmas

/-!
# Supplementary lemmas for iterated derivatives

This file supplements Mathlib with derivative iteration and product estimates, two-set affine
chain rules, and restriction and compact-support bounds for iterated derivatives.
-/

@[expose] public section

namespace Expdb

open scoped ContDiff

/-- Iterating the derivative twice adds the derivative orders. -/
theorem iteratedDeriv_iteratedDeriv (f : ℝ → ℝ) (n m : ℕ) :
    iteratedDeriv n (iteratedDeriv m f) = iteratedDeriv (n + m) f := by
  rw [iteratedDeriv_eq_iterate, iteratedDeriv_eq_iterate, iteratedDeriv_eq_iterate,
    ← Function.iterate_add_apply]

/-- Multiplication by a smooth function bounds derivative differences by a weighted
sum of the derivative bounds of the multiplier. -/
theorem norm_iteratedDeriv_mul_sub_le {f g h : ℝ → ℝ} {x : ℝ} {n : ℕ}
    (hf : ContDiffAt ℝ n f x) (hg : ContDiffAt ℝ n g x) (hh : ContDiffAt ℝ n h x)
    {A : ℕ → ℝ} {ε : ℝ}
    (hA : ∀ m ≤ n, ‖iteratedDeriv m f x‖ ≤ A m)
    (hε : ∀ m ≤ n, ‖iteratedDeriv m g x - iteratedDeriv m h x‖ ≤ ε) :
    ‖iteratedDeriv n (fun y ↦ f y * g y) x -
        iteratedDeriv n (fun y ↦ f y * h y) x‖ ≤
      (∑ m ∈ Finset.range (n + 1), (n.choose m : ℝ) * A m) * ε := by
  rw [iteratedDeriv_fun_mul hf hg, iteratedDeriv_fun_mul hf hh,
    ← Finset.sum_sub_distrib, Finset.sum_mul]
  refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun m hm ↦ ?_)
  have hmn : m ≤ n := Nat.le_of_lt_succ (Finset.mem_range.1 hm)
  have hA0 : 0 ≤ A m := (norm_nonneg _).trans (hA m hmn)
  rw [← mul_sub, norm_mul, norm_mul, Real.norm_of_nonneg (Nat.cast_nonneg _)]
  exact mul_le_mul (mul_le_mul_of_nonneg_left (hA m hmn) (Nat.cast_nonneg _))
    (hε (n - m) (Nat.sub_le _ _)) (norm_nonneg _)
    (mul_nonneg (Nat.cast_nonneg _) hA0)

/-- Iterated derivatives of a function precomposed with an affine map. Unlike Mathlib's
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

/-- Iterated derivatives agree under restriction to a smaller set on which the functions agree. -/
theorem iteratedDerivWithin_eq_of_restriction
    {f g : ℝ → ℝ} {s t : Set ℝ} (hst : s ⊆ t)
    (hs : UniqueDiffOn ℝ s) (ht : UniqueDiffOn ℝ t)
    (hg : ContDiffOn ℝ ∞ g t) (heq : Set.EqOn g f s) {x : ℝ} (hx : x ∈ s) (k : ℕ) :
    iteratedDerivWithin k g t x = iteratedDerivWithin k f s x := by
  have hsub := iteratedFDerivWithin_subset hst hs ht
    (hg.of_le (by exact_mod_cast (le_top : (k : ℕ∞) ≤ ⊤)) : ContDiffOn ℝ k g t) hx
  have hderiv : iteratedDerivWithin k g s x = iteratedDerivWithin k g t x := by
    simp only [iteratedDerivWithin, hsub]
  rw [← hderiv]
  exact iteratedDerivWithin_congr heq hx

/-- Every derivative of a compactly supported smooth function is bounded. -/
theorem exists_bound_iteratedDeriv_of_hasCompactSupport {f : ℝ → ℝ}
    (hf : ContDiff ℝ ∞ f) (hc : HasCompactSupport f) (k : ℕ) :
    ∃ C, ∀ x, ‖iteratedDeriv k f x‖ ≤ C := by
  obtain ⟨C, hC⟩ := Continuous.bounded_above_of_compact_support
    (hf.continuous_iteratedFDeriv (by exact_mod_cast le_top)) (hc.iteratedFDeriv (𝕜 := ℝ) k)
  refine ⟨C, fun x ↦ ?_⟩
  simpa only [norm_iteratedFDeriv_eq_norm_iteratedDeriv] using hC x

end Expdb
