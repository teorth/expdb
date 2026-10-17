module

public import Mathlib.Analysis.Fourier.PoissonSummation
import Mathlib.Topology.Algebra.InfiniteSum.Module

/-!
# Poisson summation at a real scale

This specializes Mathlib's Schwartz-space Poisson formula to the normalization
used by exponential sums. For a Schwartz function `f` and a positive real scale
`N`, the dual frequency is `N * m` and the Fourier integral gains a factor `N`.
-/

@[expose] public section

namespace Expdb

open MeasureTheory
open scoped SchwartzMap FourierTransform
noncomputable section

/-- Poisson summation at a positive real scale, with the Fourier integrals
written explicitly in the normalization of the exponential sums. -/
theorem poisson_scaled (f : 𝓢(ℝ, ℂ)) {N : ℝ} (hN : 0 < N) :
    (∑' n : ℤ, f ((n : ℝ) / N)) =
      N • ∑' m : ℤ, ∫ u : ℝ, 𝐞 (-(u * (N * m))) • f u := by
  let e : ℝ ≃L[ℝ] ℝ :=
    (LinearEquiv.smulOfNeZero ℝ ℝ N⁻¹ (inv_ne_zero hN.ne')).toContinuousLinearEquiv
  let g : 𝓢(ℝ, ℂ) := SchwartzMap.compCLMOfContinuousLinearEquiv ℂ e f
  have hg (x : ℝ) : g x = f (x / N) := by
    change f (N⁻¹ * x) = f (x / N)
    rw [div_eq_mul_inv, mul_comm]
  have hfourier (m : ℤ) :
      (𝓕 g) (m : ℝ) = N • ∫ u : ℝ, 𝐞 (-(u * (N * m))) • f u := by
    rw [SchwartzMap.fourier_coe, Real.fourier_real_eq]
    have h := Measure.integral_comp_mul_left
      (fun u : ℝ => 𝐞 (-(u * (N * m))) • f u) N⁻¹
    simp only [inv_inv, abs_of_pos hN] at h
    convert h using 1
    apply integral_congr_ae
    filter_upwards with x
    rw [hg]
    congr 2
    · field_simp
    · rw [div_eq_mul_inv, mul_comm]
  have h := g.tsum_eq_tsum_fourier 0
  simp only [zero_add, hg, AddCircle.coe_zero, fourier_eval_zero, mul_one] at h
  rw [h, tsum_congr hfourier, tsum_const_smul'']

end
end Expdb
