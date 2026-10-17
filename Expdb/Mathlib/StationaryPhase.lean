module

public import Expdb.Mathlib.QuadraticNormalization

import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.FunProp
import Mathlib.Tactic.FieldSimp
public import Mathlib.Analysis.Distribution.SchwartzSpace.Fourier
import Mathlib.Analysis.Calculus.Taylor
import Mathlib.Analysis.SpecialFunctions.Gaussian.FourierTransform
import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.Fourier.PoissonSummation

import Mathlib.Analysis.Calculus.IteratedDeriv.Lemmas
import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals
import Mathlib.MeasureTheory.Integral.Bochner.ContinuousLinearMap
import Mathlib.MeasureTheory.Function.Jacobian

/-!
# Quantitative stationary phase

Adapted from `lean/OAI/Dynamics/ThreeTorus/StationaryPhase.lean` and
`lean/OAI/MeasureTheory/Falconer/Fourier/SchwartzFourierMomentBounds.lean` in
[OpenAI Math](https://github.com/openai/math) at commit
`fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`, licensed under Apache-2.0 (see the root LICENSE).

Adaptations for Expdb: retained the quadratic identity, quantitative remainder,
and Fourier moment estimate. Imports were narrowed, and namespaces and proofs
were updated for the current Mathlib API.
-/

@[expose] public section

namespace Expdb.StationaryPhase

noncomputable section

open MeasureTheory Filter Set Complex
open scoped Topology SchwartzMap FourierTransform Real ContDiff ComplexConjugate

private theorem quadratic_gaussian_pairing (f : 𝓢(ℝ, ℂ)) {b : ℂ} (hb : 0 < b.re) :
    (∫ x : ℝ, f x * Complex.exp (-Real.pi * b * x ^ 2)) =
      (1 / b ^ (1 / 2 : ℂ)) * (∫ ξ : ℝ, (𝓕⁻ f) ξ *
        Complex.exp (-Real.pi / b * ξ ^ 2)) := by
  have hi : Integrable (fun x : ℝ => Complex.exp (-Real.pi * b * x ^ 2)) := by
    have hb' : 0 < ((Real.pi : ℂ) * b).re := by
      simpa only [mul_re, ofReal_re, ofReal_im, zero_mul, sub_zero] using mul_pos Real.pi_pos hb
    simpa only [neg_mul, mul_zero, zero_mul, add_zero] using
      integrable_cexp_quadratic hb' 0 0
  have h := VectorFourier.integral_fourierIntegral_smul_eq_flip
    (L := innerₗ ℝ) (μ := volume) (ν := volume) Real.continuous_fourierChar continuous_inner
      (𝓕⁻ f).integrable hi
  have hinv : (𝓕 ((𝓕⁻ f : 𝓢(ℝ, ℂ)) : ℝ → ℂ)) = (f : ℝ → ℂ) := by
    rw [← SchwartzMap.fourier_coe, FourierTransform.fourier_fourierInv_eq]
  have h' : (∫ x : ℝ, f x * Complex.exp (-Real.pi * b * x ^ 2)) =
      ∫ x : ℝ, (𝓕⁻ f) x * (𝓕 fun y : ℝ => Complex.exp (-Real.pi * b * y ^ 2)) x := by
    rw [flip_innerₗ] at h
    change (∫ x : ℝ, (𝓕 ((𝓕⁻ f : 𝓢(ℝ, ℂ)) : ℝ → ℂ)) x •
      Complex.exp (-Real.pi * b * x ^ 2)) = _ at h
    simpa only [hinv, smul_eq_mul] using! h
  rw [h', fourier_gaussian_pi hb, ← integral_const_mul]
  congr 1
  ext x
  ring

private theorem norm_quadratic_exp_le_one {b : ℂ} (hb : 0 ≤ b.re) (x : ℝ) :
    ‖Complex.exp (-(Real.pi : ℂ) * b * x ^ 2)‖ ≤ 1 := by
  rw [Complex.norm_exp, Real.exp_le_one_iff]
  have he : (-(Real.pi : ℂ) * b * (x : ℂ) ^ 2).re = -Real.pi * b.re * x ^ 2 := by
    simp [mul_re, mul_im, pow_two]
  rw [he]
  have := mul_nonneg (mul_nonneg Real.pi_pos.le hb) (sq_nonneg x)
  linarith

/-- The Gaussian identity extends to the imaginary boundary after pairing with
an integrable Schwartz amplitude. This uses dominated convergence, never a
non-absolutely convergent Bochner integral for a bare Fresnel kernel. -/
private theorem quadratic_boundary_pairing (f : 𝓢(ℝ, ℂ)) {b : ℂ}
    (hb : 0 ≤ b.re) (hbslit : b ∈ Complex.slitPlane) :
    (∫ x : ℝ, f x * Complex.exp (-Real.pi * b * x ^ 2)) =
      (1 / b ^ (1 / 2 : ℂ)) * (∫ ξ : ℝ, (𝓕⁻ f) ξ *
        Complex.exp (-Real.pi / b * ξ ^ 2)) := by
  have hb0 : b ≠ 0 := Complex.slitPlane_ne_zero hbslit
  have hparam : Tendsto (fun ε : ℝ => b + (ε : ℂ)) (𝓝[>] 0) (𝓝 b) := by
    have hcast : Tendsto (fun ε : ℝ => (ε : ℂ)) (𝓝[>] 0) (𝓝 (0 : ℂ)) :=
      (Complex.continuous_ofReal.tendsto (0 : ℝ)).mono_left nhdsWithin_le_nhds
    simpa only [add_zero] using hcast.const_add b
  have heps : ∀ᶠ ε : ℝ in 𝓝[>] 0, 0 < ε := self_mem_nhdsWithin
  have hre : ∀ᶠ ε : ℝ in 𝓝[>] 0, 0 < (b + (ε : ℂ)).re := by
    filter_upwards [heps] with ε hε
    simp only [add_re, ofReal_re]
    linarith
  have hleft : Tendsto (fun ε : ℝ => ∫ x : ℝ,
      f x * Complex.exp (-Real.pi * (b + (ε : ℂ)) * x ^ 2)) (𝓝[>] 0)
        (𝓝 (∫ x : ℝ, f x * Complex.exp (-Real.pi * b * x ^ 2))) := by
    refine tendsto_integral_filter_of_dominated_convergence (fun x => ‖f x‖)
      ?_ ?_ f.integrable.norm ?_
    · exact Eventually.of_forall (fun ε => by apply Continuous.aestronglyMeasurable; fun_prop)
    · filter_upwards [hre] with ε hε
      filter_upwards [] with x
      rw [norm_mul]
      exact (mul_le_mul_of_nonneg_left (norm_quadratic_exp_le_one hε.le x)
        (norm_nonneg (f x))).trans_eq (mul_one _)
    · filter_upwards [] with x
      exact tendsto_const_nhds.mul (((hparam.const_mul _).mul_const _).cexp)
  have hintegral : Tendsto (fun ε : ℝ => ∫ x : ℝ,
      (𝓕⁻ f) x * Complex.exp (-Real.pi / (b + (ε : ℂ)) * x ^ 2)) (𝓝[>] 0)
        (𝓝 (∫ x : ℝ, (𝓕⁻ f) x * Complex.exp (-Real.pi / b * x ^ 2))) := by
    refine tendsto_integral_filter_of_dominated_convergence (fun x => ‖(𝓕⁻ f) x‖)
      ?_ ?_ (𝓕⁻ f).integrable.norm ?_
    · exact Eventually.of_forall (fun ε => by apply Continuous.aestronglyMeasurable; fun_prop)
    · filter_upwards [hre] with ε hε
      filter_upwards [] with x
      rw [norm_mul, div_eq_mul_inv]
      have hb' : 0 ≤ (b + (ε : ℂ))⁻¹.re := by
        rw [Complex.inv_re]
        exact div_nonneg hε.le (Complex.normSq_nonneg _)
      exact (mul_le_mul_of_nonneg_left (norm_quadratic_exp_le_one hb' x)
        (norm_nonneg ((𝓕⁻ f) x))).trans_eq (mul_one _)
    · filter_upwards [] with x
      exact tendsto_const_nhds.mul (((tendsto_const_nhds.div hparam hb0).mul_const _).cexp)
  have hpow : Tendsto (fun ε : ℝ => 1 / (b + (ε : ℂ)) ^ (1 / 2 : ℂ))
      (𝓝[>] 0) (𝓝 (1 / b ^ (1 / 2 : ℂ))) := by
    apply tendsto_const_nhds.div
    · exact (continuousAt_cpow_const hbslit).tendsto.comp hparam
    · exact cpow_ne_zero_iff.mpr (Or.inl hb0)
  have hr := hpow.mul hintegral
  have hright := hr.congr' (hre.mono (fun ε hε => (quadratic_gaussian_pairing f hε).symm))
  exact tendsto_nhds_unique hleft hright

private theorem cpow_pos_mul {t : ℝ} (ht : 0 < t) {z : ℂ} (hz : z ≠ 0) (w : ℂ) :
    ((t : ℂ) * z) ^ w = (t : ℂ) ^ w * z ^ w := by
  rw [Complex.cpow_def_of_ne_zero (mul_ne_zero (ofReal_ne_zero.mpr ht.ne') hz),
    Complex.cpow_def_of_ne_zero (ofReal_ne_zero.mpr ht.ne'),
    Complex.cpow_def_of_ne_zero hz, Complex.log_ofReal_mul ht hz,
    ← Complex.ofReal_log ht.le, add_mul, Complex.exp_add]

private theorem quadratic_factor {t : ℝ} (ht : 0 < t) :
    (Real.sqrt t : ℂ) / ((t : ℂ) * (-I)) ^ (1 / 2 : ℂ) =
      Complex.exp ((Real.pi : ℂ) / 4 * I) := by
  have hs : (Real.sqrt t : ℂ) = (t : ℂ) ^ (1 / 2 : ℂ) := by
    rw [Real.sqrt_eq_rpow, Complex.ofReal_cpow ht.le]
    norm_num
  rw [cpow_pos_mul ht (neg_ne_zero.mpr I_ne_zero), hs, div_mul_eq_div_div,
    div_self (cpow_ne_zero_iff.mpr (Or.inl (ofReal_ne_zero.mpr ht.ne'))),
    Complex.cpow_def_of_ne_zero (neg_ne_zero.mpr I_ne_zero), Complex.log_neg_I,
    one_div, ← Complex.exp_neg]
  congr 1
  ring

private theorem integral_fourierInv (f : 𝓢(ℝ, ℂ)) : (∫ x : ℝ, (𝓕⁻ f) x) = f 0 := by
  have he : (𝓕 ((𝓕⁻ f : 𝓢(ℝ, ℂ)) : ℝ → ℂ)) 0 = f 0 := by
    rw [← SchwartzMap.fourier_coe, FourierTransform.fourier_fourierInv_eq]
  simpa [Real.fourier_eq] using he

lemma quadratic_stationary_identity (f : 𝓢(ℝ, ℂ)) {t : ℝ} (ht : 0 < t) :
    (Real.sqrt t : ℂ) * (∫ x : ℝ, f x * Complex.exp (I * (Real.pi * t * x^2 : ℝ))) =
      Complex.exp ((Real.pi : ℂ) / 4 * I) *
        (∫ ξ : ℝ, (𝓕⁻ f) ξ * Complex.exp (I * (-Real.pi / t * ξ^2 : ℝ))) := by
  have hb : 0 ≤ ((t : ℂ) * (-I)).re := by simp
  have hslit : (t : ℂ) * (-I) ∈ Complex.slitPlane := by
    rw [Complex.mem_slitPlane_iff]
    right
    simpa using ht.ne'
  have he := quadratic_boundary_pairing f hb hslit
  have hleft : ∀ x : ℝ, Complex.exp (-(Real.pi : ℂ) * ((t : ℂ) * (-I)) * x^2) =
      Complex.exp (I * (Real.pi * t * x^2 : ℝ)) := by
    intro x
    congr 1
    push_cast
    ring
  have hright : ∀ x : ℝ, Complex.exp (-(Real.pi : ℂ) / ((t : ℂ) * (-I)) * x^2) =
      Complex.exp (I * (-Real.pi / t * x^2 : ℝ)) := by
    intro x
    congr 1
    push_cast
    field_simp
    simp [I_sq]
    ring
  simp only [hleft, hright] at he
  rw [he, ← mul_assoc, ← div_eq_mul_one_div, quadratic_factor ht]

lemma quadratic_stationary_error (f : 𝓢(ℝ, ℂ)) {t : ℝ} (ht : 0 < t) :
    ‖(Real.sqrt t : ℂ) * (∫ x : ℝ, f x * Complex.exp (I * (Real.pi * t * x^2 : ℝ))) -
      Complex.exp ((Real.pi : ℂ) / 4 * I) * f 0‖ ≤
        (Real.pi / t) * (∫ ξ : ℝ, ξ^2 * ‖(𝓕⁻ f) ξ‖) := by
  rw [quadratic_stationary_identity f ht, ← mul_sub, norm_mul]
  have hnorm : ‖Complex.exp ((Real.pi : ℂ) / 4 * I)‖ = 1 := by
    rw [Complex.norm_exp]
    simp
  rw [hnorm, one_mul, ← integral_fourierInv f]
  have hi : Integrable (fun ξ : ℝ => (𝓕⁻ f) ξ *
      Complex.exp (I * (-Real.pi / t * ξ^2 : ℝ))) := by
    apply ((𝓕⁻ f).integrable (μ := volume)).mul_bdd
    · apply Continuous.aestronglyMeasurable
      fun_prop
    · filter_upwards [] with ξ
      have he : ‖Complex.exp (I * (-Real.pi / t * ξ^2 : ℝ))‖ = 1 := by
        rw [mul_comm I, Complex.norm_exp_ofReal_mul_I]
      exact he.le
  rw [← integral_sub hi (𝓕⁻ f).integrable]
  have hd : Integrable (fun ξ : ℝ => ξ^2 * ‖(𝓕⁻ f) ξ‖) := by
    have hh := (𝓕⁻ f).integrable_pow_mul volume 2
    simpa only [Real.norm_eq_abs, sq_abs] using hh
  rw [← integral_const_mul]
  apply norm_integral_le_of_norm_le (hd.const_mul (Real.pi / t))
  filter_upwards [] with ξ
  have heq : (𝓕⁻ f) ξ * Complex.exp (I * (-Real.pi / t * ξ^2 : ℝ)) - (𝓕⁻ f) ξ =
      (𝓕⁻ f) ξ * (Complex.exp (I * (-Real.pi / t * ξ^2 : ℝ)) - 1) := by ring
  rw [heq, norm_mul]
  have he := Real.norm_exp_I_mul_ofReal_sub_one_le (x := -Real.pi / t * ξ^2)
  have he' : ‖Complex.exp (I * (-Real.pi / t * ξ^2 : ℝ)) - 1‖ ≤
      Real.pi / t * ξ^2 := by
    convert he using 1
    rw [Real.norm_eq_abs, abs_mul, abs_div, abs_neg, abs_of_pos Real.pi_pos,
      abs_of_pos ht, abs_of_nonneg (sq_nonneg ξ)]
  calc
    ‖(𝓕⁻ f) ξ‖ * ‖Complex.exp (I * (-Real.pi / t * ξ^2 : ℝ)) - 1‖ ≤
        ‖(𝓕⁻ f) ξ‖ * (Real.pi / t * ξ^2) :=
      mul_le_mul_of_nonneg_left he' (norm_nonneg _)
    _ = _ := by ring


private theorem schwartz_iteratedDeriv_integrable (φ : 𝓢(ℝ, ℂ)) (n : ℕ) :
    Integrable (iteratedDeriv n (φ : ℝ → ℂ)) := by
  apply (φ.integrable_pow_mul_iteratedFDeriv volume 0 n).mono'
    ((φ.smooth ⊤).continuous_iteratedDeriv n
      (by exact_mod_cast (le_top : (n : ℕ∞) ≤ ⊤))).aestronglyMeasurable
  filter_upwards with x
  simp [norm_iteratedFDeriv_eq_norm_iteratedDeriv]

/-- The L¹ norm of a derivative, normalized by the Fourier derivative factor. -/
noncomputable def fourierDerivativeL1 (φ : 𝓢(ℝ, ℂ)) (n : ℕ) : ℝ :=
  (∫ x : ℝ, ‖iteratedDeriv n (φ : ℝ → ℂ) x‖)/(2*Real.pi)^n

private theorem schwartz_fourier_weighted_pointwise (φ : 𝓢(ℝ, ℂ)) (n : ℕ) (ξ : ℝ) :
    ‖ξ‖^n * ‖𝓕 φ ξ‖ ≤ fourierDerivativeL1 φ n := by
  have he := congrFun (Real.fourier_iteratedDeriv (N := (⊤ : ℕ∞)) (φ.smooth ⊤)
    (fun k _ => schwartz_iteratedDeriv_integrable φ k) (n := n) (by simp)) ξ
  have hh := VectorFourier.norm_fourierIntegral_le_integral_norm Real.fourierChar volume
    (innerₗ ℝ) (iteratedDeriv n (φ : ℝ → ℂ)) ξ
  have hh' : ‖𝓕 (iteratedDeriv n (φ : ℝ → ℂ)) ξ‖ ≤ ∫ x, ‖iteratedDeriv n (φ : ℝ → ℂ) x‖ := by
    simpa using! hh
  rw [he] at hh'
  simp only [norm_smul,norm_pow,norm_mul,Complex.norm_ofNat,Complex.norm_real,
    Real.norm_eq_abs,abs_of_pos Real.pi_pos,Complex.norm_I,mul_one] at hh'
  have hh'' : (2*Real.pi)^n * (‖ξ‖^n*‖𝓕 φ ξ‖) ≤ ∫ x, ‖iteratedDeriv n (φ : ℝ → ℂ) x‖ := by
    convert hh' using 1
    simp only [mul_pow,Real.norm_eq_abs,SchwartzMap.fourier_coe]
    ring
  exact (le_div_iff₀ (by positivity : 0 < (2*Real.pi)^n)).mpr (by simpa [mul_comm] using hh'')

/-- An `n`th Fourier moment is bounded by the `n`th and `(n + 2)`nd derivative norms. -/
theorem schwartz_fourier_moment_bound (φ : 𝓢(ℝ, ℂ)) (n : ℕ) :
    (∫ ξ : ℝ, ‖ξ‖^n*‖𝓕 φ ξ‖) ≤
      (fourierDerivativeL1 φ n + fourierDerivativeL1 φ (n+2)) *
        ∫ ξ : ℝ, (1+ξ^2)⁻¹ := by
  have hi : Integrable (fun ξ : ℝ => ‖ξ‖^n*‖𝓕 φ ξ‖) := (𝓕 φ).integrable_pow_mul volume n
  rw [← integral_const_mul]
  apply integral_mono hi (integrable_inv_one_add_sq.const_mul _) _
  intro ξ
  have h₀ := schwartz_fourier_weighted_pointwise φ n ξ
  have h₂ := schwartz_fourier_weighted_pointwise φ (n+2) ξ
  have hp : 0 < 1+ξ^2 := by positivity
  dsimp only
  rw [← div_eq_mul_inv]
  apply (le_div_iff₀ hp).mpr
  have he : ‖ξ‖^(n+2)*‖𝓕 φ ξ‖ = (‖ξ‖^n*‖𝓕 φ ξ‖)*ξ^2 := by
    rw [pow_add,Real.norm_eq_abs,sq_abs]
    ring
  rw [he] at h₂
  nlinarith


/-- The inverse Fourier transform has the same even weighted moment as the
Fourier transform. This is the moment needed by the quadratic remainder. -/
lemma integral_sq_norm_fourierInv (f : 𝓢(ℝ, ℂ)) :
    (∫ ξ : ℝ, ξ ^ 2 * ‖(𝓕⁻ f) ξ‖) = ∫ ξ : ℝ, ξ ^ 2 * ‖(𝓕 f) ξ‖ := by
  have h := MeasureTheory.integral_neg_eq_self
    (fun ξ : ℝ => ξ ^ 2 * ‖(𝓕 f) ξ‖) volume
  simpa only [Real.fourierInv_eq_fourier_neg, SchwartzMap.fourierInv_coe,
    SchwartzMap.fourier_coe, neg_sq] using h

/-- A quantitative remainder controlled by ordinary derivative norms of the
amplitude. Uniform bounds for these norms give a uniform `t⁻³ᐟ²` remainder. -/
lemma quadratic_stationary_error_of_derivative_bounds (f : 𝓢(ℝ, ℂ))
    {t B₂ B₄ : ℝ} (ht : 0 < t)
    (h₂ : fourierDerivativeL1 f 2 ≤ B₂) (h₄ : fourierDerivativeL1 f 4 ≤ B₄) :
    ‖(∫ x : ℝ, f x * Complex.exp (I * (Real.pi * t * x ^ 2 : ℝ))) -
        Complex.exp ((Real.pi : ℂ) / 4 * I) * f 0 / (Real.sqrt t : ℂ)‖ ≤
      (Real.pi * ((B₂ + B₄) * ∫ ξ : ℝ, (1 + ξ ^ 2)⁻¹)) /
        (t * Real.sqrt t) := by
  have hs : 0 < Real.sqrt t := Real.sqrt_pos.2 ht
  have hm := schwartz_fourier_moment_bound f 2
  simp only [Real.norm_eq_abs, sq_abs] at hm
  have hmoment : (∫ ξ : ℝ, ξ ^ 2 * ‖(𝓕⁻ f) ξ‖) ≤
      (B₂ + B₄) * ∫ ξ : ℝ, (1 + ξ ^ 2)⁻¹ := by
    rw [integral_sq_norm_fourierInv]
    exact hm.trans (mul_le_mul_of_nonneg_right (add_le_add h₂ h₄)
      (integral_nonneg (fun _ => by positivity)))
  have he := (quadratic_stationary_error f ht).trans
    (mul_le_mul_of_nonneg_left hmoment (by positivity))
  have hfactor :
      (Real.sqrt t : ℂ) *
          ((∫ x : ℝ, f x * Complex.exp (I * (Real.pi * t * x ^ 2 : ℝ))) -
            Complex.exp ((Real.pi : ℂ) / 4 * I) * f 0 / (Real.sqrt t : ℂ)) =
        (Real.sqrt t : ℂ) *
            (∫ x : ℝ, f x * Complex.exp (I * (Real.pi * t * x ^ 2 : ℝ))) -
          Complex.exp ((Real.pi : ℂ) / 4 * I) * f 0 := by
    have hs0 : (Real.sqrt t : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr hs.ne'
    field_simp
  rw [← hfactor, norm_mul, Complex.norm_real, Real.norm_of_nonneg hs.le] at he
  calc
    _ ≤ ((Real.pi / t) *
        ((B₂ + B₄) * ∫ ξ : ℝ, (1 + ξ ^ 2)⁻¹)) / Real.sqrt t :=
      (le_div_iff₀ hs).2 (by simpa only [mul_comm] using he)
    _ = _ := by ring


/-- The negative-curvature version of the derivative-controlled remainder.
Model phases have negative second derivative, so this is the sign used in the
B-process. It follows by conjugation, with the same amplitude bounds. -/
lemma quadratic_stationary_error_negative_of_derivative_bounds (f : 𝓢(ℝ, ℂ))
    {t B₂ B₄ : ℝ} (ht : 0 < t)
    (h₂ : fourierDerivativeL1 f 2 ≤ B₂) (h₄ : fourierDerivativeL1 f 4 ≤ B₄) :
    ‖(∫ x : ℝ, f x * Complex.exp (-I * (Real.pi * t * x ^ 2 : ℝ))) -
        Complex.exp (-(Real.pi : ℂ) / 4 * I) * f 0 / (Real.sqrt t : ℂ)‖ ≤
      (Real.pi * ((B₂ + B₄) * ∫ ξ : ℝ, (1 + ξ ^ 2)⁻¹)) /
        (t * Real.sqrt t) := by
  let g : 𝓢(ℝ, ℂ) := SchwartzMap.postcompCLM (𝕜 := ℝ)
    Complex.conjCLE.toContinuousLinearMap f
  have hg (x : ℝ) : g x = conj (f x) := rfl
  have hderiv (k : ℕ) (x : ℝ) :
      iteratedDeriv k (g : ℝ → ℂ) x = conj (iteratedDeriv k (f : ℝ → ℂ) x) := by
    change iteratedFDeriv ℝ k (Complex.conjCLE ∘ (f : ℝ → ℂ)) x (fun _ => 1) = _
    rw [Complex.conjCLE.iteratedFDeriv_comp_left]
    rfl
  have hL1 (k : ℕ) : fourierDerivativeL1 g k = fourierDerivativeL1 f k := by
    unfold fourierDerivativeL1
    simp only [hderiv, Complex.norm_conj]
  have he := quadratic_stationary_error_of_derivative_bounds g ht
    (by simpa only [hL1] using h₂) (by simpa only [hL1] using h₄)
  have hphase (a : ℝ) : conj (Complex.exp (-I * (a : ℂ))) =
      Complex.exp (I * (a : ℂ)) := by
    rw [← Complex.exp_conj]
    simp only [map_mul, map_neg, Complex.conj_I, Complex.conj_ofReal, neg_neg]
  have hmain : conj (Complex.exp (-(Real.pi : ℂ) / 4 * I)) =
      Complex.exp ((Real.pi : ℂ) / 4 * I) := by
    rw [← Complex.exp_conj]
    congr 1
    simp only [map_mul, map_div₀, map_neg, Complex.conj_ofReal, Complex.conj_I, map_ofNat]
    ring
  have hconj :
      conj ((∫ x : ℝ, f x * Complex.exp (-I * (Real.pi * t * x ^ 2 : ℝ))) -
        Complex.exp (-(Real.pi : ℂ) / 4 * I) * f 0 / (Real.sqrt t : ℂ)) =
      (∫ x : ℝ, g x * Complex.exp (I * (Real.pi * t * x ^ 2 : ℝ))) -
        Complex.exp ((Real.pi : ℂ) / 4 * I) * g 0 / (Real.sqrt t : ℂ) := by
    rw [map_sub, ← integral_conj, map_div₀, map_mul, hmain, Complex.conj_ofReal]
    simp only [map_mul, hphase, hg]
  rw [← Complex.norm_conj, hconj]
  exact he

/-! ## Stationary phase after quadratic normalization

The following estimates are proved within Expdb. -/

private theorem iteratedDeriv_ofReal_comp {g : ℝ → ℝ} (hg : ContDiff ℝ ∞ g) (k : ℕ) (x : ℝ) :
    iteratedDeriv k (fun y ↦ (g y : ℂ)) x = ((iteratedDeriv k g x : ℝ) : ℂ) := by
  change iteratedFDeriv ℝ k (Complex.ofRealCLM ∘ g) x (fun _ ↦ 1) = _
  rw [Complex.ofRealCLM.iteratedFDeriv_comp_left hg.contDiffAt (by exact_mod_cast le_top)]
  rfl

/-- The integral norm budget of a function supported in an interval. -/
theorem integral_norm_le_of_support_Icc {g : ℝ → ℝ} {a b B : ℝ} (hab : a ≤ b)
    (hzero : ∀ y, y ∉ Icc a b → g y = 0) (hB : ∀ y ∈ Icc a b, ‖g y‖ ≤ B) :
    ∫ y, ‖g y‖ ≤ B * (b - a) := by
  rw [← setIntegral_eq_integral_of_forall_compl_eq_zero (s := Icc a b)]
  · have h := norm_setIntegral_le_of_norm_le_const (μ := volume) (s := Icc a b)
      (f := fun y ↦ ‖g y‖) (C := B) (by simp) (fun y hy ↦ by simpa using hB y hy)
    rw [Real.norm_of_nonneg (integral_nonneg fun _ ↦ norm_nonneg _)] at h
    simpa only [Measure.real, Real.volume_Icc, ENNReal.toReal_ofReal (sub_nonneg.2 hab)] using h
  · intro y hy
    simp [hzero y hy]

/-- An integral of a function supported in `[-R, R]` and bounded by `B` is at most
`2 R B`. -/
private theorem integral_norm_le_of_support {g : ℝ → ℝ} {R B : ℝ} (hR : 0 ≤ R)
    (hzero : ∀ y, R < |y| → g y = 0) (hB : ∀ y, ‖g y‖ ≤ B) :
    ∫ y, ‖g y‖ ≤ 2 * R * B := by
  have h := integral_norm_le_of_support_Icc (g := g) (a := -R) (b := R) (by linarith)
    (fun y hy ↦ hzero y (by
      rw [mem_Icc, not_and_or, not_le, not_le] at hy
      rcases hy with h | h
      · rw [abs_of_neg (by linarith)]; linarith
      · rw [abs_of_pos (by linarith)]; linarith)) (fun y _ ↦ hB y)
  nlinarith

/-- Stationary phase for a fixed globally concave phase. A quadratic change of variables
reduces the integral to the quadratic stationary-phase estimate. -/
theorem norm_integral_stationary_sub_le {H w X : ℝ → ℝ} {c C lam R B v : ℝ}
    (hH : ContDiff ℝ ∞ H) (hc : 0 < c)
    (hcurv : ∀ x, -C ≤ iteratedDeriv 2 H x ∧ iteratedDeriv 2 H x ≤ -c)
    (hw : ContDiff ℝ ∞ w) (hX : ContDiff ℝ ∞ X)
    (hYX : ∀ y, phaseQuadraticCoordinate H v (X y) = y)
    (hlam : 0 < lam) (hR : 0 ≤ R)
    (hsupp : ∀ y, R < |y| → w (X y) = 0)
    (hB2 : ∀ y, ‖iteratedDeriv 2 (fun y ↦ w (X y) * deriv X y) y‖ ≤ B)
    (hB4 : ∀ y, ‖iteratedDeriv 4 (fun y ↦ w (X y) * deriv X y) y‖ ≤ B) :
    ‖(∫ x : ℝ, (w x : ℂ) * 𝐞 (lam * (H x - deriv H v * x))) -
        𝐞 (lam * (H v - deriv H v * v)) *
          (Complex.exp (-(Real.pi : ℂ) / 4 * I) *
            ((w v / Real.sqrt (-iteratedDeriv 2 H v) : ℝ) : ℂ) / (Real.sqrt lam : ℂ))‖ ≤
      4 * Real.pi ^ 2 * R * B / (lam * Real.sqrt lam) := by
  set Y := phaseQuadraticCoordinate H v with hYdef
  have hC : 0 < C := by linarith [(hcurv 0).1, (hcurv 0).2]
  have hY : ContDiff ℝ ∞ Y := contDiff_phaseQuadraticCoordinate hH hc hcurv v
  have hYd : Differentiable ℝ Y := hY.differentiable (by simp)
  have hXd : Differentiable ℝ X := hX.differentiable (by simp)
  have hYpos : ∀ x, 0 < deriv Y x := fun x ↦
    lt_of_lt_of_le (div_pos hc (Real.sqrt_pos.2 hC)) (le_deriv_phaseQuadraticCoordinate hH hc
      hcurv v x)
  -- derivative of the inverse
  have hchain : ∀ y, deriv Y (X y) * deriv X y = 1 := by
    intro y
    have h1 : deriv (Y ∘ X) y = 1 := by
      have : Y ∘ X = id := funext hYX
      rw [this, deriv_id]
    rwa [deriv_comp y (hYd _) (hXd _)] at h1
  have hXpos : ∀ y, 0 < deriv X y := by
    intro y
    have := hchain y
    have := hYpos (X y)
    by_contra h
    push Not at h
    nlinarith
  -- the amplitude
  set amp : ℝ → ℝ := fun y ↦ w (X y) * deriv X y with hamp
  have hampS : ContDiff ℝ ∞ amp := (hw.comp hX).mul (hX.deriv' (n := ∞))
  have hampz : ∀ y, R < |y| → amp y = 0 := fun y hy ↦ by simp [amp, hsupp y hy]
  have hampzd : ∀ k y, R < |y| → iteratedDeriv k amp y = 0 := by
    intro k y hy
    have hev : amp =ᶠ[𝓝 y] fun _ ↦ 0 := by
      have hopen : IsOpen {z : ℝ | R < |z|} := isOpen_lt continuous_const continuous_abs
      filter_upwards [hopen.mem_nhds hy] with z hz
      exact hampz z hz
    rw [hev.iteratedDeriv_eq]
    simp
  set ampC : ℝ → ℂ := fun y ↦ (amp y : ℂ) with hampC
  have hampCS : ContDiff ℝ ∞ ampC := Complex.ofRealCLM.contDiff.comp hampS
  have hampCc : HasCompactSupport ampC := by
    apply HasCompactSupport.intro (isCompact_closedBall (0 : ℝ) R)
    intro y hy
    rw [Metric.mem_closedBall, Real.dist_eq, sub_zero, not_le] at hy
    simp [ampC, hampz y hy]
  let f : 𝓢(ℝ, ℂ) := hampCc.toSchwartzMap hampCS
  have hf : ∀ y, f y = ampC y := fun y ↦ by simp [f]
  have hfcoe : (f : ℝ → ℂ) = ampC := funext hf
  -- derivative bounds
  have hL1 : ∀ k, (∀ y, ‖iteratedDeriv k amp y‖ ≤ B) →
      StationaryPhase.fourierDerivativeL1 f k ≤ 2 * R * B := by
    intro k hk
    unfold StationaryPhase.fourierDerivativeL1
    rw [hfcoe]
    have hint : ∫ y, ‖iteratedDeriv k ampC y‖ ≤ 2 * R * B := by
      have : (fun y ↦ ‖iteratedDeriv k ampC y‖) = fun y ↦ ‖iteratedDeriv k amp y‖ := by
        funext y
        rw [hampC, iteratedDeriv_ofReal_comp hampS, Complex.norm_real]
      rw [this]
      exact integral_norm_le_of_support hR (fun y hy ↦ hampzd k y hy) hk
    have hpow : 1 ≤ (2 * Real.pi) ^ k :=
      one_le_pow₀ (by nlinarith [Real.two_le_pi])
    have hnn : 0 ≤ ∫ y, ‖iteratedDeriv k ampC y‖ := integral_nonneg fun _ ↦ norm_nonneg _
    calc (∫ y, ‖iteratedDeriv k ampC y‖) / (2 * Real.pi) ^ k
        ≤ ∫ y, ‖iteratedDeriv k ampC y‖ := div_le_self hnn hpow
      _ ≤ 2 * R * B := hint
  have herr := StationaryPhase.quadratic_stationary_error_negative_of_derivative_bounds f hlam
    (hL1 2 hB2) (hL1 4 hB4)
  rw [integral_univ_inv_one_add_sq] at herr
  -- value at zero
  have hX0 : X 0 = v := by
    have hYv : Y v = 0 := phaseQuadraticCoordinate_self H v
    have hinj : Function.Injective Y := (strictMono_of_deriv_pos hYpos).injective
    apply hinj
    rw [hYX, hYv]
  have hcurvv : iteratedDeriv 2 H v < 0 := by linarith [(hcurv v).2]
  have hYv' : deriv Y v = Real.sqrt (-iteratedDeriv 2 H v) :=
    (hasDerivAt_phaseQuadraticCoordinate_self H v
      (hH.of_le (ENat.natCast_le_of_coe_top_le_withTop le_rfl 3))
      hcurvv).deriv
  have hf0 : f 0 = ((w v / Real.sqrt (-iteratedDeriv 2 H v) : ℝ) : ℂ) := by
    rw [hf]
    simp only [ampC, amp, hX0]
    have h := hchain 0
    rw [hX0, hYv'] at h
    have hs : 0 < Real.sqrt (-iteratedDeriv 2 H v) := Real.sqrt_pos.2 (by linarith)
    congr 1
    field_simp
    rw [mul_assoc, mul_comm (deriv X 0), h, mul_one]
  -- change of variables
  set ν := deriv H v
  set A := H v - ν * v
  have hcov : (∫ x : ℝ, (w x : ℂ) * 𝐞 (lam * (H x - ν * x))) =
      𝐞 (lam * A) * ∫ y : ℝ, f y * Complex.exp (-I * (Real.pi * lam * y ^ 2 : ℝ)) := by
    have hsurj : X '' univ = univ := by
      rw [image_univ, range_eq_univ]
      intro x
      exact ⟨Y x, by
        have hinj : Function.Injective Y := (strictMono_of_deriv_pos hYpos).injective
        apply hinj
        rw [hYX]⟩
    have hinjX : InjOn X univ := fun a _ b _ hab ↦ by
      have := congrArg Y hab
      rwa [hYX, hYX] at this
    have hcv := integral_image_eq_integral_abs_deriv_smul (s := univ) (f := X) (f' := deriv X)
      MeasurableSet.univ (fun y _ ↦ (hXd y).hasDerivAt.hasDerivWithinAt) hinjX
      (fun x ↦ (w x : ℂ) * 𝐞 (lam * (H x - ν * x)))
    rw [hsurj, setIntegral_univ, setIntegral_univ] at hcv
    rw [hcv, ← integral_const_mul]
    congr 1
    funext y
    have hq := phaseQuadraticFactor_mem hH hcurv v (X y)
    have hphase := phase_sub_tangent_eq_neg_coordinate_sq H v (X y)
      (fun _ _ ↦ (hH.of_le (ENat.natCast_le_of_coe_top_le_withTop le_rfl 2)).contDiffAt) (by
        linarith [hq.2])
    rw [← hYdef, hYX] at hphase
    rw [abs_of_pos (hXpos y), hf]
    simp only [ampC, amp, Real.fourierChar_apply, real_smul]
    have : (2 * Real.pi * (lam * (H (X y) - ν * X y)) : ℝ) =
        2 * Real.pi * (lam * A) + (-(Real.pi * lam * y ^ 2)) := by
      simp only [A]
      linear_combination (2 * Real.pi * lam) * hphase
    rw [this]
    push_cast
    rw [show ((2 : ℂ) * Real.pi * (lam * A) + -(Real.pi * lam * y ^ 2)) * I =
      (2 * Real.pi * (lam * A)) * I + -I * (Real.pi * lam * y ^ 2) by ring, Complex.exp_add]
    ring
  rw [hcov, ← mul_sub, norm_mul, Circle.norm_coe, one_mul]
  rw [hf0] at herr
  calc _ ≤ Real.pi * ((2 * R * B + 2 * R * B) * Real.pi) / (lam * Real.sqrt lam) := herr
    _ = _ := by ring

/-! ## Scalar amplitude estimates (proved for Expdb) -/

/-- Reciprocal square roots are Lipschitz on an interval bounded away from zero. -/
theorem abs_inv_sqrt_sub_le {c s s' : ℝ} (hc : 0 < c) (hs : c ≤ s) (hs' : c ≤ s') :
    |1 / Real.sqrt s - 1 / Real.sqrt s'| ≤ |s - s'| / (2 * c * Real.sqrt c) := by
  have hsc : 0 < Real.sqrt c := Real.sqrt_pos.2 hc
  have hrs : Real.sqrt c ≤ Real.sqrt s := Real.sqrt_le_sqrt hs
  have hrs' : Real.sqrt c ≤ Real.sqrt s' := Real.sqrt_le_sqrt hs'
  have hps : 0 < Real.sqrt s := lt_of_lt_of_le hsc hrs
  have hps' : 0 < Real.sqrt s' := lt_of_lt_of_le hsc hrs'
  have h1 : 1 / Real.sqrt s - 1 / Real.sqrt s' =
      (s' - s) / (Real.sqrt s * Real.sqrt s' * (Real.sqrt s + Real.sqrt s')) := by
    have e1 := Real.sq_sqrt (hc.le.trans hs)
    have e2 := Real.sq_sqrt (hc.le.trans hs')
    field_simp
    nlinarith [e1, e2]
  have hdpos : 0 < Real.sqrt s * Real.sqrt s' * (Real.sqrt s + Real.sqrt s') :=
    mul_pos (mul_pos hps hps') (add_pos hps hps')
  have h2c : 0 < 2 * c * Real.sqrt c := mul_pos (mul_pos two_pos hc) hsc
  rw [h1, abs_div, abs_of_pos hdpos, abs_sub_comm s' s]
  have hden : 2 * c * Real.sqrt c ≤
      Real.sqrt s * Real.sqrt s' * (Real.sqrt s + Real.sqrt s') := by
    have e : Real.sqrt c * Real.sqrt c = c := Real.mul_self_sqrt hc.le
    have h2 : Real.sqrt c * Real.sqrt c ≤ Real.sqrt s * Real.sqrt s' :=
      mul_le_mul hrs hrs' hsc.le hps.le
    have h3 : 2 * Real.sqrt c ≤ Real.sqrt s + Real.sqrt s' := by linarith
    calc 2 * c * Real.sqrt c = (Real.sqrt c * Real.sqrt c) * (2 * Real.sqrt c) := by
          rw [e]; ring
      _ ≤ (Real.sqrt s * Real.sqrt s') * (Real.sqrt s + Real.sqrt s') :=
          mul_le_mul h2 h3 (by positivity) (by positivity)
  calc |s - s'| / (Real.sqrt s * Real.sqrt s' * (Real.sqrt s + Real.sqrt s'))
      ≤ |s - s'| / (2 * c * Real.sqrt c) :=
        div_le_div_of_nonneg_left (abs_nonneg _) h2c hden

/-! ## Stationary points and the Legendre transform (proved for Expdb) -/

/-- The inverse of the derivative of a phase: `phaseInverse H ν` is the stationary point of
`x ↦ H x - ν x`. -/
def phaseInverse (H : ℝ → ℝ) : ℝ → ℝ := Function.invFun (deriv H)

theorem phaseInverse_spec {H : ℝ → ℝ} (hH : ContDiff ℝ ∞ H) {c : ℝ} (hc : 0 < c)
    (hcurv : ∀ x, iteratedDeriv 2 H x ≤ -c) :
    ContDiff ℝ ∞ (phaseInverse H) ∧ (∀ x, phaseInverse H (deriv H x) = x) ∧
      ∀ y, deriv H (phaseInverse H y) = y := by
  have hH1 : ContDiff ℝ ∞ (deriv H) := hH.deriv' (n := ∞)
  have hd : ∀ x, deriv (deriv H) x ≤ -c := fun x ↦ by
    simpa [iteratedDeriv_succ, iteratedDeriv_one] using hcurv x
  obtain ⟨g, hg, hgf, hfg⟩ := exists_smooth_inverse_of_deriv_le_neg hH1 hc hd
  have heq : phaseInverse H = g :=
    Function.invFun_eq_of_injective_of_rightInverse
      (show Function.LeftInverse g _ from hgf).injective hfg
  rw [heq]
  exact ⟨hg, hgf, hfg⟩

/-- The stationary value `ν ↦ H(x(ν)) - ν x(ν)`, using the negative Legendre-transform
convention. Here `x(ν)` is the stationary point of `x ↦ H x - ν x`. -/
def legendre (H : ℝ → ℝ) (ν : ℝ) : ℝ := H (phaseInverse H ν) - ν * phaseInverse H ν

theorem contDiff_legendre {H : ℝ → ℝ} (hH : ContDiff ℝ ∞ H) {c : ℝ} (hc : 0 < c)
    (hcurv : ∀ x, iteratedDeriv 2 H x ≤ -c) : ContDiff ℝ ∞ (legendre H) := by
  have hψ := (phaseInverse_spec hH hc hcurv).1
  exact (hH.comp hψ).sub (contDiff_id.mul hψ)

theorem hasDerivAt_legendre {H : ℝ → ℝ} (hH : ContDiff ℝ ∞ H) {c : ℝ} (hc : 0 < c)
    (hcurv : ∀ x, iteratedDeriv 2 H x ≤ -c) (ν : ℝ) :
    HasDerivAt (legendre H) (-phaseInverse H ν) ν := by
  obtain ⟨hψ, -, hfg⟩ := phaseInverse_spec hH hc hcurv
  have hψd : HasDerivAt (phaseInverse H) (deriv (phaseInverse H) ν) ν :=
    ((hψ.differentiable (by simp)) ν).hasDerivAt
  have hHd : HasDerivAt H (deriv H (phaseInverse H ν)) (phaseInverse H ν) :=
    ((hH.differentiable (by simp)) _).hasDerivAt
  have h := (hHd.comp ν hψd).sub ((hasDerivAt_id ν).mul hψd)
  refine h.congr_deriv ?_
  rw [hfg]
  simp only [id]
  ring

end

end Expdb.StationaryPhase
