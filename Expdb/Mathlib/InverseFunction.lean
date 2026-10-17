module

public import Mathlib.Analysis.Calculus.ContDiff.Basic
import Mathlib.Analysis.Calculus.ContDiff.Operations
import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.Calculus.InverseFunctionTheorem.Deriv
import Mathlib.Order.Hom.Set
import Mathlib.Topology.Order.IntermediateValue
import Mathlib.Topology.Order.Monotone

/-!
# Smooth inverses under global derivative bounds

A real smooth function whose derivative stays uniformly positive or negative has a smooth
inverse. The same derivative bounds give Lipschitz estimates for the inverse.
-/

@[expose] public section

namespace Expdb

open Filter Topology
open scoped ContDiff
noncomputable section

/-- A positive derivative lower bound gives a two-sided linear lower bound. -/
theorem mul_abs_sub_le_abs_sub_of_le_deriv {f : ℝ → ℝ} (hf : Differentiable ℝ f) {c : ℝ}
    (h : ∀ x, c ≤ deriv f x) (a b : ℝ) : c * |a - b| ≤ |f a - f b| := by
  rcases le_total a b with hab | hab
  · have := mul_sub_le_image_sub_of_le_deriv hf h hab
    rw [abs_sub_comm a b, abs_of_nonneg (sub_nonneg.2 hab), abs_sub_comm]
    exact this.trans (le_abs_self _)
  · have := mul_sub_le_image_sub_of_le_deriv hf h hab
    rw [abs_of_nonneg (sub_nonneg.2 hab)]
    exact this.trans (le_abs_self _)

/-- A smooth function with derivative bounded below by a positive constant is a smooth
diffeomorphism of the real line. -/
theorem exists_smooth_inverse_of_le_deriv {f : ℝ → ℝ} (hf : ContDiff ℝ ∞ f) {c : ℝ}
    (hc : 0 < c) (h : ∀ x, c ≤ deriv f x) :
    ∃ g : ℝ → ℝ, ContDiff ℝ ∞ g ∧ (∀ x, g (f x) = x) ∧ (∀ y, f (g y) = y) := by
  have hd : Differentiable ℝ f := hf.differentiable (by simp)
  have hmono : StrictMono f :=
    strictMono_of_deriv_pos fun x ↦ lt_of_lt_of_le hc (h x)
  have htop : Tendsto f atTop atTop := by
    refine tendsto_atTop_mono' atTop ?_
      (tendsto_atTop_add_const_left atTop (f 0)
        (tendsto_id.const_mul_atTop hc))
    filter_upwards [eventually_ge_atTop 0] with x hx
    have := mul_sub_le_image_sub_of_le_deriv hd h hx
    simp only [sub_zero, id] at this ⊢
    linarith
  have hbot : Tendsto f atBot atBot := by
    refine tendsto_atBot_mono' atBot ?_
      (tendsto_atBot_add_const_left atBot (f 0)
        (tendsto_id.const_mul_atBot hc))
    filter_upwards [eventually_le_atBot 0] with x hx
    have := mul_sub_le_image_sub_of_le_deriv hd h hx
    simp only [zero_sub, id] at this ⊢
    linarith
  have hsurj : Function.Surjective f := hd.continuous.surjective htop hbot
  let e : ℝ ≃o ℝ := StrictMono.orderIsoOfSurjective f hmono hsurj
  let E : ℝ ≃ₜ ℝ := e.toHomeomorph
  have hE : (E : ℝ → ℝ) = f := rfl
  have hsmooth : ContDiff ℝ ∞ (E.symm : ℝ → ℝ) := by
    apply E.contDiff_symm_deriv (f' := deriv f)
    · exact fun x ↦ (lt_of_lt_of_le hc (h x)).ne'
    · exact fun x ↦ (hd x).hasDerivAt
    · exact hf
  refine ⟨E.symm, hsmooth, fun x ↦ ?_, fun y ↦ ?_⟩
  · exact E.symm_apply_apply x
  · exact E.apply_symm_apply y

/-- A smooth function with derivative bounded above by a negative constant is a smooth
diffeomorphism of the real line. -/
theorem exists_smooth_inverse_of_deriv_le_neg {f : ℝ → ℝ} (hf : ContDiff ℝ ∞ f) {c : ℝ}
    (hc : 0 < c) (h : ∀ x, deriv f x ≤ -c) :
    ∃ g : ℝ → ℝ, ContDiff ℝ ∞ g ∧ (∀ x, g (f x) = x) ∧ (∀ y, f (g y) = y) := by
  have hneg : ∀ x, c ≤ deriv (fun y ↦ -f y) x := by
    intro x
    rw [deriv.fun_neg]
    linarith [h x]
  obtain ⟨g, hg, hgf, hfg⟩ := exists_smooth_inverse_of_le_deriv hf.neg hc hneg
  refine ⟨fun y ↦ g (-y), hg.comp contDiff_neg, fun x ↦ ?_, fun y ↦ ?_⟩
  · simpa using hgf x
  · have := hfg (-y)
    linarith

/-- Inverses of functions with derivative at least `c > 0` are `1/c`-Lipschitz. -/
private theorem abs_inverse_sub_le {f g : ℝ → ℝ} (hf : Differentiable ℝ f) {c : ℝ} (hc : 0 < c)
    (h : ∀ x, c ≤ deriv f x) (hfg : ∀ y, f (g y) = y) (y y' : ℝ) :
    |g y - g y'| ≤ |y - y'| / c := by
  have := mul_abs_sub_le_abs_sub_of_le_deriv hf h (g y) (g y')
  rw [hfg, hfg] at this
  rw [le_div_iff₀ hc]
  linarith

/-- Inverses of functions with derivative at most `-c < 0` are `1/c`-Lipschitz. -/
theorem abs_inverse_sub_le_of_neg {f g : ℝ → ℝ} (hf : Differentiable ℝ f) {c : ℝ}
    (hc : 0 < c) (h : ∀ x, deriv f x ≤ -c) (hfg : ∀ y, f (g y) = y) (y y' : ℝ) :
    |g y - g y'| ≤ |y - y'| / c := by
  have hneg : ∀ x, c ≤ deriv (fun y ↦ -f y) x := by
    intro x
    rw [deriv.fun_neg]
    linarith [h x]
  have := abs_inverse_sub_le (f := fun y ↦ -f y) (g := fun y ↦ g (-y)) hf.neg hc hneg
    (fun y ↦ by simp [hfg]) (-y) (-y')
  rw [show -y - -y' = -(y - y') by ring, abs_neg] at this
  simpa using this

end
end Expdb
