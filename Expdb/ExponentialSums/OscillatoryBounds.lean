module

public import Expdb.ExponentialSums.FixedExponentialSum
public import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic
public import Mathlib.Order.Interval.Finset.Nat

import Expdb.Mathlib.EulerMaclaurin
import Expdb.Mathlib.IteratedDeriv
import Mathlib.Analysis.Calculus.ContDiff.Bounds
import Mathlib.Analysis.Fourier.FourierTransformDeriv
import Mathlib.Analysis.Real.Pi.Bounds
import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts

/-!
# Oscillatory integral and exponential-sum bounds

Calculus for the Fourier character and uniform estimates for `oscillatory F T N`, including a
first-derivative integral bound and Euler–Maclaurin comparisons between oscillatory sums and
integrals.
-/

@[expose] public section

open Filter Topology
open scoped ContDiff Expdb FourierTransform

noncomputable section

namespace Expdb

/-! ## Fourier-character calculus -/

/-- The real Fourier character `x ↦ 𝐞 x`, viewed as complex-valued, is smooth. -/
theorem contDiff_fourierChar : ContDiff ℝ ∞ (𝐞 · : ℝ → ℂ) := by
  rw [show (𝐞 · : ℝ → ℂ) = fun x : ℝ ↦
      Complex.exp (((2 * Real.pi * x : ℝ) : ℂ) * Complex.I) by
    funext x
    exact Real.fourierChar_apply x]
  have hreal : ContDiff ℝ ∞ (fun x : ℝ ↦ 2 * Real.pi * x) := by fun_prop
  exact ((Complex.ofRealCLM.contDiff.comp hreal).mul contDiff_const).cexp

/-- The `n`th derivative of the real Fourier character is
`(2πi)ⁿ 𝐞 x`. -/
theorem iteratedDeriv_fourierChar (n : ℕ) (x : ℝ) :
    iteratedDeriv n (𝐞 · : ℝ → ℂ) x =
      ((2 * Real.pi : ℂ) * Complex.I) ^ n * 𝐞 x := by
  induction n with
  | zero => simp
  | succ n ih =>
      rw [iteratedDeriv_succ']
      have hderiv : deriv (𝐞 · : ℝ → ℂ) =
          fun y ↦ ((2 * Real.pi : ℂ) * Complex.I) * 𝐞 y := by
        funext y
        exact Real.deriv_fourierChar y
      rw [hderiv, iteratedDeriv_const_mul _
        (contDiff_fourierChar.contDiffAt.of_le
          (ENat.natCast_le_of_coe_top_le_withTop le_rfl n)), ih]
      ring

/-- The norm of the `n`th derivative of the real Fourier character is exactly `(2π)ⁿ`. -/
@[simp] theorem norm_iteratedDeriv_fourierChar (n : ℕ) (x : ℝ) :
    ‖iteratedDeriv n (𝐞 · : ℝ → ℂ) x‖ = (2 * Real.pi) ^ n := by
  rw [iteratedDeriv_fourierChar, norm_mul, norm_pow]
  simp [Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg Real.pi_pos.le]

private lemma norm_iteratedDeriv_fourierChar_le
    {n k : ℕ} (hk : k ≤ n) (x : ℝ) :
    ‖iteratedDeriv k (𝐞 · : ℝ → ℂ) x‖ ≤ (2 * Real.pi + 1) ^ n := by
  rw [norm_iteratedDeriv_fourierChar]
  exact (pow_le_pow_left₀ (by positivity) (by linarith) k).trans
    (pow_le_pow_right₀ (by linarith [Real.pi_pos]) hk)

/-! ## Derivative bounds for oscillatory phases -/

/-- Uniform derivative control for a rescaled oscillatory phase. -/
theorem norm_iteratedDerivWithin_oscillatory_le
    {F : ℝ → ℝ} {s : Set ℝ} {x T N K : ℝ} {n : ℕ}
    (hF : ContDiffOn ℝ n F phaseInterval) (hs : UniqueDiffOn ℝ s)
    (hx : x ∈ s) (hmap : Set.MapsTo (N⁻¹ * ·) s phaseInterval)
    (hT : 1 ≤ T) (hN : 1 ≤ N) (hK : 1 ≤ K)
    (hderiv : HasPhaseDerivBound F n K) :
    ‖iteratedDerivWithin n (oscillatory F T N) s x‖ ≤
      (n.factorial : ℝ) * (2 * Real.pi + 1) ^ n * (K * T / N) ^ n := by
  let q : ℝ → ℝ := fun y ↦ T * F (N⁻¹ * y)
  have hscale : ContDiffOn ℝ n (N⁻¹ * ·) s := by fun_prop
  have hFscale : ContDiffOn ℝ n (fun y ↦ F (N⁻¹ * y)) s :=
    by
      change ContDiffOn ℝ n (F ∘ fun y ↦ N⁻¹ * y) s
      exact hF.comp hscale hmap
  have hq : ContDiffOn ℝ n q s := by
    dsimp [q]
    exact contDiffOn_const.mul hFscale
  have houter : ContDiffOn ℝ n (𝐞 · : ℝ → ℂ) Set.univ :=
    (contDiff_fourierChar.of_le
      (ENat.natCast_le_of_coe_top_le_withTop le_rfl n)).contDiffOn
  have hC : ∀ k, k ≤ n →
      ‖iteratedFDerivWithin ℝ k (𝐞 · : ℝ → ℂ) Set.univ (q x)‖ ≤
        (2 * Real.pi + 1) ^ n := by
    intro k hk
    rw [norm_iteratedFDerivWithin_eq_norm_iteratedDerivWithin,
      iteratedDerivWithin_eq_iteratedDeriv uniqueDiffOn_univ
        (contDiff_fourierChar.contDiffAt.of_le
          (ENat.natCast_le_of_coe_top_le_withTop le_rfl k)) (Set.mem_univ _)]
    exact norm_iteratedDeriv_fourierChar_le hk _
  have hD : ∀ k, 1 ≤ k → k ≤ n →
      ‖iteratedFDerivWithin ℝ k q s x‖ ≤ (K * T / N) ^ k := by
    intro k hk hkn
    rw [norm_iteratedFDerivWithin_eq_norm_iteratedDerivWithin]
    have hrescale : iteratedDerivWithin k (fun y ↦ F (N⁻¹ * y)) s x =
        N⁻¹ ^ k * iteratedDerivWithin k F phaseInterval (N⁻¹ * x) := by
      simpa only [add_zero] using iteratedDerivWithin_comp_affine_of_mapsTo
        (f := F) (c := N⁻¹) (d := 0) (n := k)
        (hF.of_le (Nat.cast_le.mpr hkn)) hs
        (uniqueDiffOn_Icc (by norm_num [phaseInterval])) hx
        (by simpa only [add_zero] using hmap)
    have hqderiv : iteratedDerivWithin k q s x =
        T * N⁻¹ ^ k * iteratedDerivWithin k F phaseInterval (N⁻¹ * x) := by
      rw [show q = fun y ↦ T * F (N⁻¹ * y) by rfl,
        iteratedDerivWithin_const_mul_field, hrescale]
      ring
    rw [hqderiv, norm_mul, norm_mul, Real.norm_eq_abs,
      abs_of_nonneg (le_trans zero_le_one hT), Real.norm_eq_abs,
      abs_of_nonneg (pow_nonneg (inv_nonneg.mpr (le_trans zero_le_one hN)) _)]
    have hFx : ‖iteratedDerivWithin k F phaseInterval (N⁻¹ * x)‖ ≤ K :=
      hderiv k hk hkn _ (hmap hx)
    calc
      T * N⁻¹ ^ k * ‖iteratedDerivWithin k F phaseInterval (N⁻¹ * x)‖ ≤
          T * N⁻¹ ^ k * K := by gcongr
      _ = (K * T) * N⁻¹ ^ k := by ring
      _ ≤ (K * T) ^ k * N⁻¹ ^ k := by
        gcongr
        simpa using pow_le_pow_right₀ (one_le_mul_of_one_le_of_one_le hK hT) hk
      _ = (K * T / N) ^ k := by simp only [div_eq_mul_inv, mul_pow]
  have hcomp := norm_iteratedFDerivWithin_comp_le
    houter hq (by rfl) uniqueDiffOn_univ hs (Set.mapsTo_univ q s) hx hC hD
  rw [norm_iteratedFDerivWithin_eq_norm_iteratedDerivWithin] at hcomp
  change ‖iteratedDerivWithin n (fun y ↦ (𝐞 (q y) : ℂ)) s x‖ ≤
    (n.factorial : ℝ) * (2 * Real.pi + 1) ^ n * (K * T / N) ^ n at hcomp
  rw [show oscillatory F T N = fun y ↦ (𝐞 (T * F (y / N)) : ℂ) from rfl]
  simpa only [q, div_eq_mul_inv, mul_comm N⁻¹] using hcomp

/-- A smooth phase remains smooth after rescaling and composition with the Fourier character. -/
theorem contDiffOn_oscillatory
    {F : ℝ → ℝ} {s : Set ℝ} {T N : ℝ} {n : ℕ}
    (hF : ContDiffOn ℝ n F phaseInterval)
    (hmap : Set.MapsTo (N⁻¹ * ·) s phaseInterval) :
    ContDiffOn ℝ n (oscillatory F T N) s := by
  have hscale : ContDiffOn ℝ n (N⁻¹ * ·) s := by fun_prop
  have hFscale : ContDiffOn ℝ n (fun x ↦ F (N⁻¹ * x)) s := by
    change ContDiffOn ℝ n (F ∘ fun x ↦ N⁻¹ * x) s
    exact hF.comp hscale hmap
  have hinner : ContDiffOn ℝ n (fun x ↦ T * F (N⁻¹ * x)) s :=
    contDiffOn_const.mul hFscale
  have heq : oscillatory F T N =
      (𝐞 · : ℝ → ℂ) ∘ fun x ↦ T * F (N⁻¹ * x) := by
    funext x
    simp only [oscillatory, Function.comp_apply, div_eq_mul_inv]
    rw [mul_comm x N⁻¹]
  rw [heq]
  exact (contDiff_fourierChar.of_le
    (ENat.natCast_le_of_coe_top_le_withTop le_rfl n)).comp_contDiffOn hinner

/-! ## A first-derivative estimate for the phase integral -/

/-- A phase with first derivative at least `c` and second derivative bounded by `K` has an
    oscillatory integral bounded by `(2 / c + K / c ^ 2) / T`. -/
theorem norm_integral_oscillatory_le
    {F : ℝ → ℝ} {A B T c K : ℝ}
    (hAB : A ≤ B) (hsub : Set.Icc A B ⊆ phaseInterval)
    (hF : ContDiffOn ℝ 2 F phaseInterval) (hT : 0 < T) (hc : 0 < c)
    (hfirst : ∀ x ∈ Set.Icc A B,
      c ≤ iteratedDerivWithin 1 F phaseInterval x)
    (hsecond : ∀ x ∈ Set.Icc A B,
      ‖iteratedDerivWithin 2 F phaseInterval x‖ ≤ K) :
    ‖∫ x in A..B, (𝐞 (T * F x) : ℂ)‖ ≤ (2 / c + K / c ^ 2) / T := by
  let d : ℝ → ℝ := iteratedDerivWithin 1 F phaseInterval
  let dd : ℝ → ℝ := iteratedDerivWithin 2 F phaseInterval
  let τ : ℂ := (2 * Real.pi : ℂ) * Complex.I
  let v : ℝ → ℂ := fun x ↦ 𝐞 (T * F x)
  let vp : ℝ → ℂ := fun x ↦ τ * T * d x * v x
  let z : ℝ → ℂ := fun x ↦ τ * T * d x
  let u : ℝ → ℂ := z⁻¹
  let up : ℝ → ℂ := fun x ↦ -(τ * T * dd x) / (z x) ^ 2
  have hunique : UniqueDiffOn ℝ phaseInterval := uniqueDiffOn_Icc (by
    norm_num [phaseInterval])
  have hτnorm : ‖τ‖ = 2 * Real.pi := by
    dsimp [τ]
    simp [Complex.norm_real, Real.norm_eq_abs, Real.pi_pos.le]
  have hτone : 1 ≤ ‖τ‖ := by rw [hτnorm]; linarith [Real.pi_gt_three]
  have hderivF (x : ℝ) (hx : x ∈ Set.Icc A B) :
      HasDerivWithinAt F (d x) (Set.Icc A B) x := by
    dsimp [d]
    rw [iteratedDerivWithin_one]
    exact ((hF.differentiableOn (by norm_num) x (hsub hx)).hasDerivWithinAt).mono hsub
  have hderivD (x : ℝ) (hx : x ∈ Set.Icc A B) :
      HasDerivWithinAt d (dd x) (Set.Icc A B) x := by
    simpa only [d, dd, show (2 : ℕ) = 1 + 1 by omega,
      iteratedDerivWithin_succ] using
        ((hF.differentiableOn_iteratedDerivWithin (m := 1) (by norm_num) hunique
          x (hsub hx)).hasDerivWithinAt).mono hsub
  have hv (x : ℝ) (hx : x ∈ Set.Icc A B) :
      HasDerivWithinAt v (vp x) (Set.Icc A B) x := by
    have hinner := (hderivF x hx).const_mul T
    have h := (Real.hasDerivAt_fourierChar (T * F x)).hasDerivWithinAt.scomp
      x hinner (Set.mapsTo_univ _ _)
    dsimp only [v]
    change HasDerivWithinAt ((𝐞 · : ℝ → ℂ) ∘ fun y ↦ T * F y)
      (vp x) (Set.Icc A B) x
    rw [show vp x = (T * d x : ℝ) •
        ((2 * Real.pi : ℂ) * Complex.I * 𝐞 (T * F x)) by
      dsimp [vp, τ]
      simp only [Complex.ofReal_mul]
      ring]
    exact h
  have hz (x : ℝ) (hx : x ∈ Set.Icc A B) :
      HasDerivWithinAt z (τ * T * dd x) (Set.Icc A B) x := by
    have h := (hderivD x hx).ofReal_comp.const_mul (τ * T)
    dsimp [z]
    simpa only [mul_assoc] using h
  have hz_ne (x : ℝ) (hx : x ∈ Set.Icc A B) : z x ≠ 0 := by
    have hdpos : 0 < d x := hc.trans_le (hfirst x hx)
    dsimp [z]
    exact mul_ne_zero (mul_ne_zero (norm_ne_zero_iff.mp
      (ne_of_gt (lt_of_lt_of_le zero_lt_one hτone)))
        (Complex.ofReal_ne_zero.mpr hT.ne'))
      (Complex.ofReal_ne_zero.mpr hdpos.ne')
  have hu (x : ℝ) (hx : x ∈ Set.Icc A B) :
      HasDerivWithinAt u (up x) (Set.Icc A B) x := by
    have h := (hz x hx).inv (hz_ne x hx)
    simpa only [u, up, Pi.inv_apply, div_eq_mul_inv] using h
  have hvCont : ContinuousOn v (Set.Icc A B) :=
    fun x hx ↦ (hv x hx).continuousWithinAt
  have huCont : ContinuousOn u (Set.Icc A B) :=
    fun x hx ↦ (hu x hx).continuousWithinAt
  have hvpCont : ContinuousOn vp (Set.Icc A B) := by
    have hdCont : ContinuousOn d (Set.Icc A B) := by
      dsimp [d]
      exact (hF.continuousOn_iteratedDerivWithin (by norm_num) hunique).mono hsub
    dsimp [vp]
    fun_prop
  have hupCont : ContinuousOn up (Set.Icc A B) := by
    have hdCont : ContinuousOn d (Set.Icc A B) := by
      dsimp [d]
      exact (hF.continuousOn_iteratedDerivWithin (by norm_num) hunique).mono hsub
    have hddCont : ContinuousOn dd (Set.Icc A B) := by
      dsimp [dd]
      exact (hF.continuousOn_iteratedDerivWithin (by norm_num) hunique).mono hsub
    have hdC : ContinuousOn (fun x ↦ (d x : ℂ)) (Set.Icc A B) :=
      Complex.continuous_ofReal.comp_continuousOn hdCont
    have hddC : ContinuousOn (fun x ↦ (dd x : ℂ)) (Set.Icc A B) :=
      Complex.continuous_ofReal.comp_continuousOn hddCont
    dsimp [up, z]
    exact ((continuousOn_const.mul continuousOn_const).mul hddC).neg.div
      (((continuousOn_const.mul continuousOn_const).mul hdC).pow 2)
      (fun x hx ↦ pow_ne_zero _ (hz_ne x hx))
  have hu' : ∀ x ∈ Set.uIcc A B, HasDerivWithinAt u (up x) (Set.uIcc A B) x := by
    simpa [Set.uIcc_of_le hAB] using hu
  have hv' : ∀ x ∈ Set.uIcc A B, HasDerivWithinAt v (vp x) (Set.uIcc A B) x := by
    simpa [Set.uIcc_of_le hAB] using hv
  have hupInt : IntervalIntegrable up MeasureTheory.volume A B := by
    have hc : ContinuousOn up (Set.uIcc A B) := by
      simpa [Set.uIcc_of_le hAB] using hupCont
    exact hc.intervalIntegrable
  have hvpInt : IntervalIntegrable vp MeasureTheory.volume A B := by
    have hc : ContinuousOn vp (Set.uIcc A B) := by
      simpa [Set.uIcc_of_le hAB] using hvpCont
    exact hc.intervalIntegrable
  have hparts := intervalIntegral.integral_mul_deriv_eq_deriv_mul_of_hasDerivWithinAt
    hu' hv' hupInt hvpInt
  have huv (x : ℝ) (hx : x ∈ Set.Icc A B) : u x * vp x = v x := by
    change (z x)⁻¹ * (z x * v x) = v x
    rw [← mul_assoc, inv_mul_cancel₀ (hz_ne x hx), one_mul]
  have hIntegral : (∫ x in A..B, v x) =
      u B * v B - u A * v A - ∫ x in A..B, up x * v x := by
    rw [← hparts]
    apply intervalIntegral.integral_congr
    intro x hx
    exact (huv x (by simpa [Set.uIcc_of_le hAB] using hx)).symm
  have huBound (x : ℝ) (hx : x ∈ Set.Icc A B) : ‖u x‖ ≤ 1 / (T * c) := by
    rw [show ‖u x‖ = (‖τ‖ * T * |d x|)⁻¹ by
      dsimp [u, z]
      rw [norm_inv, norm_mul, norm_mul]
      simp only [Complex.norm_real, Real.norm_eq_abs, abs_of_pos hT]]
    have hd : c ≤ |d x| := (hfirst x hx).trans (le_abs_self _)
    rw [one_div]
    apply (inv_le_inv₀
      (mul_pos (mul_pos (lt_of_lt_of_le zero_lt_one hτone) hT)
        (lt_of_lt_of_le hc hd)) (mul_pos hT hc)).2
    calc
      T * c ≤ ‖τ‖ * T * c := by
        calc
          T * c = 1 * (T * c) := by ring
          _ ≤ ‖τ‖ * (T * c) := mul_le_mul_of_nonneg_right hτone (by positivity)
          _ = ‖τ‖ * T * c := by ring
      _ ≤ ‖τ‖ * T * |d x| := by gcongr
  have hupBound (x : ℝ) (hx : x ∈ Set.Icc A B) :
      ‖up x‖ ≤ K / (T * c ^ 2) := by
    have hdpos : 0 < d x := hc.trans_le (hfirst x hx)
    have hzNorm : ‖z x‖ = ‖τ‖ * T * |d x| := by
      dsimp [z]
      rw [norm_mul, norm_mul]
      simp only [Complex.norm_real, Real.norm_eq_abs, abs_of_pos hT]
    have hnumNorm : ‖-(τ * T * dd x)‖ = ‖τ‖ * T * |dd x| := by
      rw [norm_neg, norm_mul, norm_mul]
      simp only [Complex.norm_real, Real.norm_eq_abs, abs_of_pos hT]
    rw [show ‖up x‖ = ‖-(τ * T * dd x)‖ / ‖z x‖ ^ 2 by
      simp only [up, norm_div, norm_pow], hnumNorm, hzNorm]
    rw [show ‖τ‖ * T * |dd x| / (‖τ‖ * T * |d x|) ^ 2 =
        |dd x| / (‖τ‖ * T * |d x| ^ 2) by
      field_simp [hT.ne', (lt_of_lt_of_le zero_lt_one hτone).ne', hdpos.ne']]
    have hden : T * c ^ 2 ≤ ‖τ‖ * T * |d x| ^ 2 := by
      have hd : c ≤ |d x| := (hfirst x hx).trans (le_abs_self _)
      calc
        T * c ^ 2 ≤ T * |d x| ^ 2 := by gcongr
        _ ≤ ‖τ‖ * T * |d x| ^ 2 := by
          have : 0 ≤ T * |d x| ^ 2 := by positivity
          calc
            T * |d x| ^ 2 = 1 * (T * |d x| ^ 2) := by ring
            _ ≤ ‖τ‖ * (T * |d x| ^ 2) := mul_le_mul_of_nonneg_right hτone this
            _ = ‖τ‖ * T * |d x| ^ 2 := by ring
    have hK0 : 0 ≤ K := (norm_nonneg (dd x)).trans (hsecond x hx)
    calc
      |dd x| / (‖τ‖ * T * |d x| ^ 2) ≤
          K / (‖τ‖ * T * |d x| ^ 2) :=
        div_le_div_of_nonneg_right (by simpa [Real.norm_eq_abs] using hsecond x hx) (by positivity)
      _ ≤ K / (T * c ^ 2) :=
        div_le_div_of_nonneg_left hK0 (by positivity) hden
  rw [hIntegral]
  calc
    ‖u B * v B - u A * v A - ∫ x in A..B, up x * v x‖ ≤
        ‖u B‖ + ‖u A‖ + ‖∫ x in A..B, up x * v x‖ := by
      calc
        _ ≤ ‖u B * v B - u A * v A‖ + ‖∫ x in A..B, up x * v x‖ :=
          norm_sub_le _ _
        _ ≤ (‖u B * v B‖ + ‖u A * v A‖) + ‖∫ x in A..B, up x * v x‖ := by
          gcongr
          exact norm_sub_le _ _
        _ = _ := by simp [v]
    _ ≤ 2 / (T * c) + (B - A) * (K / (T * c ^ 2)) := by
      have hA : A ∈ Set.Icc A B := ⟨le_rfl, hAB⟩
      have hB : B ∈ Set.Icc A B := ⟨hAB, le_rfl⟩
      have hint := intervalIntegral.norm_integral_le_of_norm_le_const
        (a := A) (b := B) (C := K / (T * c ^ 2)) (fun x hx ↦ by
          have hx' : x ∈ Set.Icc A B := by
            rw [Set.uIoc_of_le hAB] at hx
            exact ⟨hx.1.le, hx.2⟩
          calc
            ‖up x * v x‖ = ‖up x‖ := by simp [v]
            _ ≤ K / (T * c ^ 2) := hupBound x hx')
      rw [abs_of_nonneg (sub_nonneg.mpr hAB)] at hint
      calc
        ‖u B‖ + ‖u A‖ + ‖∫ x in A..B, up x * v x‖ ≤
            1 / (T * c) + 1 / (T * c) +
              K / (T * c ^ 2) * (B - A) := by
          exact add_le_add (add_le_add (huBound B hB) (huBound A hA)) hint
        _ = 2 / (T * c) + (B - A) * (K / (T * c ^ 2)) := by ring
    _ ≤ (2 / c + K / c ^ 2) / T := by
      have hlen : B - A ≤ 1 := by
        have hA := (hsub ⟨le_rfl, hAB⟩).1
        have hB := (hsub ⟨hAB, le_rfl⟩).2
        linarith
      have hK0 : 0 ≤ K := (norm_nonneg (dd A)).trans (hsecond A ⟨le_rfl, hAB⟩)
      field_simp
      nlinarith

/-! ## Euler–Maclaurin bounds for oscillatory sums -/

private def oscillatoryErrorConstant (s : ℕ) : ℝ :=
  1 +
    (∑ m ∈ Finset.range s, |EulerMaclaurin.saw (m + 2) 0| *
      (2 * (((m + 1).factorial : ℝ) * (2 * Real.pi + 1) ^ (m + 1)))) +
    EulerMaclaurin.sawBound (s + 1) *
      (((s + 1).factorial : ℝ) * (2 * Real.pi + 1) ^ (s + 1))

private lemma one_le_oscillatoryErrorConstant (s : ℕ) :
    1 ≤ oscillatoryErrorConstant s := by
  have hsum : 0 ≤ ∑ m ∈ Finset.range s, |EulerMaclaurin.saw (m + 2) 0| *
      (2 * (((m + 1).factorial : ℝ) * (2 * Real.pi + 1) ^ (m + 1))) :=
    Finset.sum_nonneg fun _ _ ↦ by positivity
  have hrem : 0 ≤ EulerMaclaurin.sawBound (s + 1) *
      (((s + 1).factorial : ℝ) * (2 * Real.pi + 1) ^ (s + 1)) :=
    mul_nonneg EulerMaclaurin.sawBound_nonneg (by positivity)
  dsimp only [oscillatoryErrorConstant]
  linarith

private lemma norm_oscillatory_sum_sub_integral_le
    (s : ℕ) {F : ℝ → ℝ} {T N K : ℝ} {a b : ℕ}
    (hab : a < b) (hN : 1 ≤ N) (hT : 1 ≤ T) (hK : 1 ≤ K)
    (hF : ContDiffOn ℝ ∞ F phaseInterval)
    (hderiv : HasPhaseDerivBound F (s + 1) K)
    (ha : N ≤ a) (hb : (b : ℝ) ≤ 2 * N)
    (hscale : HasSmallEulerMaclaurinRemainder s T N K) :
    ‖exponentialSumAt F T N a b -
      ∫ x in (a : ℝ)..b, oscillatory F T N x‖ ≤
        oscillatoryErrorConstant s := by
  obtain ⟨hratio, hremainder⟩ := hscale
  rw [exponentialSumAt]
  have hNpos : 0 < N := zero_lt_one.trans_le hN
  have hrnonneg : 0 ≤ K * T / N := by positivity
  let interval : Set ℝ := Set.Icc (a : ℝ) (b : ℝ)
  have hmap : Set.MapsTo (N⁻¹ * ·) interval phaseInterval := by
    intro x hx
    simp only [interval, phaseInterval, Set.mem_Icc] at hx ⊢
    constructor
    · rw [inv_mul_eq_div, le_div_iff₀ hNpos]
      simpa using ha.trans hx.1
    · rw [inv_mul_eq_div, div_le_iff₀ hNpos]
      exact hx.2.trans hb
  have hunique : UniqueDiffOn ℝ interval := uniqueDiffOn_Icc (by exact_mod_cast hab)
  let A : ℕ → ℝ := fun k ↦ (k.factorial : ℝ) * (2 * Real.pi + 1) ^ k
  have hosc (k : ℕ) (hk : 1 ≤ k) (hkS : k ≤ s + 1)
      (x : ℝ) (hx : x ∈ interval) :
      ‖iteratedDerivWithin k (oscillatory F T N) interval x‖ ≤
        A k * (K * T / N) ^ k := by
    simpa only [A] using norm_iteratedDerivWithin_oscillatory_le
      (hF.of_le (ENat.natCast_le_of_coe_top_le_withTop le_rfl k))
      hunique hx hmap hT hN hK
      (fun q hq hqk u hu ↦ hderiv q hq (hqk.trans hkS) u hu)
  have hraw := EulerMaclaurin.norm_sum_Icc_nat_sub_integral_le
    (f := oscillatory F T N) (t := interval)
    hab.le
    (contDiffOn_oscillatory
      (hF.of_le (ENat.natCast_le_of_coe_top_le_withTop le_rfl (s + 1))) hmap)
    hunique Set.Subset.rfl
    (C₀ := 1) (fun _ _ ↦ by simp)
    (fun k ↦ A k * (K * T / N) ^ k) hosc
  let Ccorr : ℝ := ∑ m ∈ Finset.range s,
    |EulerMaclaurin.saw (m + 2) 0| * (2 * A (m + 1))
  let Crem : ℝ := EulerMaclaurin.sawBound (s + 1) * A (s + 1)
  have hCrem : 0 ≤ Crem := by
    dsimp [Crem, A]
    exact mul_nonneg EulerMaclaurin.sawBound_nonneg (by positivity)
  change _ ≤ 1 + Ccorr + Crem
  calc
    _ ≤ 1 +
        (∑ m ∈ Finset.range s, |EulerMaclaurin.saw (m + 2) 0| *
          (2 * (A (m + 1) * (K * T / N) ^ (m + 1)))) +
        EulerMaclaurin.sawBound (s + 1) *
          (A (s + 1) * (K * T / N) ^ (s + 1)) * (b - a : ℕ) := hraw
    _ ≤ 1 + Ccorr + Crem := by
      have hcorr :
          (∑ m ∈ Finset.range s, |EulerMaclaurin.saw (m + 2) 0| *
            (2 * (A (m + 1) * (K * T / N) ^ (m + 1)))) ≤ Ccorr := by
        dsimp [Ccorr]
        apply Finset.sum_le_sum
        intro m hm
        have hA0 : 0 ≤ A (m + 1) := by dsimp [A]; positivity
        have hp : (K * T / N) ^ (m + 1) ≤ 1 :=
          pow_le_one₀ hrnonneg hratio
        gcongr
        simpa only [mul_one] using mul_le_mul_of_nonneg_left hp hA0
      have hlen : (b - a : ℝ) ≤ N := by
        linarith
      have hprod : (b - a : ℝ) * (K * T / N) ^ (s + 1) ≤ 1 :=
        (mul_le_mul_of_nonneg_right hlen (pow_nonneg hrnonneg _)).trans hremainder
      have hrem :
          EulerMaclaurin.sawBound (s + 1) *
              (A (s + 1) * (K * T / N) ^ (s + 1)) * (b - a : ℕ) ≤ Crem := by
        calc
          EulerMaclaurin.sawBound (s + 1) *
              (A (s + 1) * (K * T / N) ^ (s + 1)) * (b - a : ℕ) =
              Crem * ((b - a : ℝ) * (K * T / N) ^ (s + 1)) := by
            dsimp [Crem]
            rw [Nat.cast_sub hab.le]
            ring
          _ ≤ Crem := by
            simpa only [mul_one] using mul_le_mul_of_nonneg_left hprod hCrem
      linarith

/-- A uniform Euler–Maclaurin comparison between an oscillatory sum and its integral. The
constant depends only on the differentiation order. -/
theorem exists_norm_oscillatory_sum_sub_integral_le (s : ℕ) :
    ∃ C : ℝ, 1 ≤ C ∧
      ∀ {F : ℝ → ℝ} {T N K : ℝ} {a b : ℕ},
        a < b → 1 ≤ N → 1 ≤ T → 1 ≤ K →
        ContDiffOn ℝ ∞ F phaseInterval →
        HasPhaseDerivBound F (s + 1) K →
        N ≤ a → (b : ℝ) ≤ 2 * N →
        HasSmallEulerMaclaurinRemainder s T N K →
        ‖exponentialSumAt F T N a b -
          ∫ x in (a : ℝ)..b, oscillatory F T N x‖ ≤ C := by
  refine ⟨oscillatoryErrorConstant s, one_le_oscillatoryErrorConstant s, ?_⟩
  intro F T N K a b hab hN hT hK hF hderiv ha hb hscale
  exact norm_oscillatory_sum_sub_integral_le s hab hN hT hK hF hderiv
    ha hb hscale

private def oscillatorySumConstant (s : ℕ) (K c : ℝ) : ℝ :=
  oscillatoryErrorConstant s + (2 / c + K / c ^ 2)

private lemma one_le_oscillatorySumConstant
    (s : ℕ) {K c : ℝ} (hK : 1 ≤ K) (hc : 0 < c) :
    1 ≤ oscillatorySumConstant s K c := by
  dsimp only [oscillatorySumConstant]
  have hE := one_le_oscillatoryErrorConstant s
  have hD : 0 ≤ 2 / c + K / c ^ 2 := by positivity
  linarith

private lemma norm_oscillatory_sum_le
    (s : ℕ) {F : ℝ → ℝ} {T N K c : ℝ} {a b : ℕ}
    (hs : 1 ≤ s) (hab : a ≤ b) (hN : 1 ≤ N) (hT : 1 ≤ T) (hK : 1 ≤ K) (hc : 0 < c)
    (hF : ContDiffOn ℝ ∞ F phaseInterval)
    (hfirst : HasPhaseFirstDerivLowerBound F c)
    (hderiv : HasPhaseDerivBound F (s + 1) K)
    (ha : N ≤ a) (hb : (b : ℝ) ≤ 2 * N)
    (hscale : HasSmallEulerMaclaurinRemainder s T N K) :
    ‖exponentialSumAt F T N a b‖ ≤
      oscillatorySumConstant s K c * (1 + N / T) := by
  have hNpos : 0 < N := zero_lt_one.trans_le hN
  have hTpos : 0 < T := zero_lt_one.trans_le hT
  by_cases heq : a = b
  · subst b
    rw [exponentialSumAt_self, norm_oscillatory]
    have hconst : 1 ≤ oscillatorySumConstant s K c := by
      dsimp only [oscillatorySumConstant]
      have hE := one_le_oscillatoryErrorConstant s
      have hD : 0 ≤ 2 / c + K / c ^ 2 := by positivity
      linarith
    have hNT : 0 ≤ N / T := div_nonneg hNpos.le hTpos.le
    nlinarith
  have hablt : a < b := hab.lt_of_ne heq
  have herr := norm_oscillatory_sum_sub_integral_le s hablt hN hT hK hF hderiv
    ha hb hscale
  have hAB : (a : ℝ) / N ≤ (b : ℝ) / N :=
    div_le_div_of_nonneg_right (by exact_mod_cast hab) hNpos.le
  have hsub : Set.Icc ((a : ℝ) / N) ((b : ℝ) / N) ⊆ phaseInterval := by
    intro x hx
    constructor
    · exact ((le_div_iff₀ hNpos).2 (by simpa using ha)).trans hx.1
    · exact hx.2.trans ((div_le_iff₀ hNpos).2 (by simpa using hb))
  have hy := norm_integral_oscillatory_le hAB hsub
    (hF.of_le (ENat.natCast_le_of_coe_top_le_withTop le_rfl 2)) hTpos hc
    (fun x hx ↦ hfirst x (hsub hx))
    (fun x hx ↦ hderiv 2 (by norm_num) (by omega) x (hsub hx))
  have hcomp := intervalIntegral.integral_comp_div
    (fun y : ℝ ↦ (𝐞 (T * F y) : ℂ)) (a := (a : ℝ)) (b := (b : ℝ)) hNpos.ne'
  have hint :
      ‖∫ x in (a : ℝ)..b, oscillatory F T N x‖ ≤
        (2 / c + K / c ^ 2) * (N / T) := by
    rw [show (fun x ↦ oscillatory F T N x) =
      fun x ↦ (𝐞 (T * F (x / N)) : ℂ) from rfl,
      hcomp, norm_smul, Real.norm_eq_abs, abs_of_pos hNpos]
    calc
      N * ‖∫ x in (a : ℝ) / N..(b : ℝ) / N, (𝐞 (T * F x) : ℂ)‖ ≤
          N * ((2 / c + K / c ^ 2) / T) := by gcongr
      _ = (2 / c + K / c ^ 2) * (N / T) := by ring
  calc
    ‖exponentialSumAt F T N a b‖ ≤
        ‖exponentialSumAt F T N a b -
          ∫ x in (a : ℝ)..b, oscillatory F T N x‖ +
        ‖∫ x in (a : ℝ)..b, oscillatory F T N x‖ :=
      norm_le_norm_sub_add _ _
    _ ≤ oscillatoryErrorConstant s + (2 / c + K / c ^ 2) * (N / T) := by
      gcongr
    _ ≤ oscillatorySumConstant s K c * (1 + N / T) := by
      have hE := one_le_oscillatoryErrorConstant s
      have hD : 0 ≤ 2 / c + K / c ^ 2 := by positivity
      have hNT : 0 ≤ N / T := by positivity
      dsimp only [oscillatorySumConstant]
      nlinarith

/-- A first-derivative bound for an oscillatory sum, uniform at scales where the
Euler–Maclaurin remainder is bounded. -/
theorem exists_norm_oscillatory_sum_le
    (s : ℕ) (hs : 1 ≤ s) {K c : ℝ} (hK : 1 ≤ K) (hc : 0 < c) :
    ∃ C : ℝ, 1 ≤ C ∧
      ∀ {F : ℝ → ℝ} {T N : ℝ} {a b : ℕ},
        a ≤ b → 1 ≤ N → 1 ≤ T →
        ContDiffOn ℝ ∞ F phaseInterval →
        HasPhaseFirstDerivLowerBound F c →
        HasPhaseDerivBound F (s + 1) K →
        N ≤ a → (b : ℝ) ≤ 2 * N →
        HasSmallEulerMaclaurinRemainder s T N K →
        ‖exponentialSumAt F T N a b‖ ≤
          C * (1 + N / T) := by
  refine ⟨oscillatorySumConstant s K c,
    one_le_oscillatorySumConstant s hK hc, ?_⟩
  intro F T N a b hab hN hT hF hfirst hderiv ha hb hscale
  exact norm_oscillatory_sum_le s hs hab hN hT hK hc hF hfirst hderiv
    ha hb hscale


end Expdb

end
