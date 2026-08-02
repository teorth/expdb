module

public import Expdb.Basic.Definitions
public import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic

import Mathlib.Algebra.Order.Interval.Set.Group
import Mathlib.Analysis.Calculus.BumpFunction.InnerProduct
import Mathlib.Analysis.Distribution.SchwartzSpace.Fourier
import Mathlib.Analysis.PSeries

/-!
# L² integral estimate

This module formalizes the `L^2` integral estimate from Chapter 3 of the ANTEDB blueprint.
It proves both the blueprint's equality with a bounded error coefficient
and the corresponding absolute-error bound.
-/

@[expose] public section

open MeasureTheory Real Complex Filter Topology BigOperators
open scoped FourierTransform SchwartzMap ContDiff

noncomputable section

namespace Expdb

/-! ### A normalized smooth bump -/

private def rawBump : ContDiffBump (0 : ℝ) :=
  ⟨1 / 8, 1 / 4, by norm_num, by norm_num⟩

private def rawL2 : ℝ := ∫ x : ℝ, (rawBump x) ^ 2

private lemma rawL2_pos : 0 < rawL2 := by
  have hraw : Continuous (rawBump : ℝ → ℝ) :=
    (rawBump.contDiff (n := ⊤)).continuous
  have hsq : Continuous (fun x : ℝ => (rawBump x) ^ 2) := by
    convert (continuous_pow 2).comp hraw using 1
    ext x
    rfl
  apply integral_pos_of_integrable_nonneg_nonzero (x := 0)
  · exact hsq
  · apply hsq.integrable_of_hasCompactSupport
    apply HasCompactSupport.of_support_subset_isCompact rawBump.hasCompactSupport.isCompact
    intro x hx
    apply subset_tsupport rawBump
    simp only [Function.mem_support] at hx ⊢
    intro hzero
    apply hx
    simp [hzero]
  · intro x
    positivity
  · have hzero : (rawBump : ℝ → ℝ) 0 = 1 := by
      apply rawBump.one_of_mem_closedBall
      simp [rawBump, Metric.mem_closedBall]
    simp [hzero]

private def bump (x : ℝ) : ℝ := rawBump x / Real.sqrt rawL2

private lemma bump_smooth : ContDiff ℝ ∞ bump := by
  change ContDiff ℝ ∞ (fun x => rawBump x / Real.sqrt rawL2)
  exact (rawBump.contDiff (n := ⊤)).div_const (Real.sqrt rawL2)

private lemma bump_hasCompactSupport : HasCompactSupport bump := by
  apply HasCompactSupport.of_support_subset_isCompact rawBump.hasCompactSupport.isCompact
  intro x hx
  apply subset_tsupport rawBump
  simp only [Function.mem_support] at hx ⊢
  intro hzero
  apply hx
  simp [bump, hzero]

private lemma bump_supp (x : ℝ) (hx : bump x ≠ 0) : |x| ≤ 1 / 4 := by
  have hraw : (rawBump : ℝ → ℝ) x ≠ 0 := by
    intro hzero
    apply hx
    simp [bump, hzero]
  have hmem : x ∈ Metric.ball (0 : ℝ) rawBump.rOut := by
    rw [← rawBump.support_eq]
    exact hraw
  have : |x| < 1 / 4 := by
    simpa [rawBump, Metric.mem_ball, Real.dist_eq] using hmem
  exact this.le

private lemma bump_nonneg (x : ℝ) : 0 ≤ bump x :=
  div_nonneg (rawBump.nonneg' x) (Real.sqrt_nonneg rawL2)

private lemma bump_l2norm : ∫ x : ℝ, (bump x) ^ 2 = 1 := by
  rw [show (fun x : ℝ => (bump x) ^ 2) = fun x => (rawBump x) ^ 2 / rawL2 by
        funext x
        simp only [bump, div_pow]
        rw [Real.sq_sqrt rawL2_pos.le]]
  rw [integral_div]
  exact div_self (ne_of_gt rawL2_pos)

-- bump̂(u) = ∫_ℝ bump(x) e(-xu) dx
private def bumpFourier (u : ℝ) : ℂ :=
  ∫ x : ℝ, (bump x : ℂ) * 𝐞 (-(x * u))

-- bump is integrable
private lemma bump_integrable : Integrable bump :=
  bump_smooth.continuous.integrable_of_hasCompactSupport bump_hasCompactSupport

-- ∫ bump > 0  (from proof: "bump(t) ≥ 0 and ‖bump‖_{L²} = 1 so bump ≢ 0")
private lemma bump_integral_pos : 0 < ∫ x : ℝ, bump x := by
  have hne : ∃ x, bump x ≠ 0 := by
    by_contra h
    push Not at h
    have hsquare : ∫ x : ℝ, (bump x) ^ 2 = 0 := by
      calc
        ∫ x : ℝ, (bump x) ^ 2
            = ∫ x : ℝ, (0 : ℝ) ^ 2 := by
              congr 1
              ext x
              rw [h x]
        _ = 0 := by simp
    linarith [bump_l2norm, hsquare]
  obtain ⟨x, hx⟩ := hne
  exact integral_pos_of_integrable_nonneg_nonzero
    bump_smooth.continuous bump_integrable bump_nonneg hx

/-! ### Normalization and the exponential sum -/

private lemma normalize_coefficients {ι : Type*} [Fintype ι] (a : ι → ℂ) (M : ℝ)
    (hM : M = ∑ r, ‖a r‖ ^ 2)
    (hpos : 0 < M) :
    ∑ r, ‖(fun r => a r / (Real.sqrt M : ℂ)) r‖ ^ 2 = 1 := by
  simp_rw [norm_div]
  simp only [norm_real, Real.norm_eq_abs,
    abs_of_nonneg (Real.sqrt_nonneg M)]
  simp_rw [div_pow]
  rw [← Finset.sum_div]
  rw [← hM, Real.sq_sqrt hpos.le]
  exact div_self (ne_of_gt hpos)

private def expSum {ι : Type*} [Fintype ι] (a : ι → ℂ) (ξ : ι → ℝ) (t : ℝ) : ℂ :=
  ∑ r, a r * 𝐞 (ξ r * t)

private def expSumSq {ι : Type*} [Fintype ι] (a : ι → ℂ) (ξ : ι → ℝ) (t : ℝ) : ℝ :=
  ‖expSum a ξ t‖ ^ 2

private lemma expSumSq_nonneg {ι : Type*} [Fintype ι]
    (a : ι → ℂ) (ξ : ι → ℝ) (t : ℝ) :
    0 ≤ expSumSq a ξ t :=
  sq_nonneg _

private lemma expSumSq_continuous {ι : Type*} [Fintype ι]
    (a : ι → ℂ) (ξ : ι → ℝ) :
    Continuous (expSumSq a ξ) := by
  unfold expSumSq expSum
  fun_prop

/-! ### Plancherel and decay of the bump's Fourier transform -/

private lemma bump_complex_hasCompactSupport : HasCompactSupport (fun x : ℝ => (bump x : ℂ)) :=
  HasCompactSupport.of_support_subset_isCompact
    (isCompact_Icc (a := -1 / 4) (b := 1 / 4)) (by
      intro x hx
      simp only [Function.mem_support] at hx
      have h := bump_supp x (by exact_mod_cast hx)
      simp only [Set.mem_Icc, abs_le] at h ⊢
      simpa only [neg_div] using h)

private lemma bump_complex_smooth : ContDiff ℝ ∞ (fun x : ℝ => (bump x : ℂ)) := by
  change ContDiff ℝ ∞ (Complex.ofRealCLM ∘ bump)
  exact (Complex.ofRealCLM.contDiff (n := ∞)).comp bump_smooth

private def bumpSchwartz : 𝓢(ℝ, ℂ) :=
  bump_complex_hasCompactSupport.toSchwartzMap bump_complex_smooth

private lemma bumpFourier_eq_fourier : bumpFourier = 𝓕 (bumpSchwartz : ℝ → ℂ) := by
  funext u
  rw [Real.fourier_real_eq]
  simp only [bumpFourier, Real.fourierChar_apply, mul_neg, Complex.ofReal_neg,
    Complex.ofReal_mul, Complex.ofReal_ofNat, neg_mul, bumpSchwartz,
    HasCompactSupport.toSchwartzMap_toFun, Circle.smul_def, smul_eq_mul]
  apply integral_congr_ae
  filter_upwards with x
  ring

private lemma bumpFourier_continuous : Continuous bumpFourier := by
  rw [bumpFourier_eq_fourier]
  exact (𝓕 bumpSchwartz).continuous

private lemma bumpFourier_l2 : ∫ u : ℝ, ‖bumpFourier u‖ ^ 2 = 1 := by
  simp_rw [bumpFourier_eq_fourier]
  rw [← SchwartzMap.fourier_coe, SchwartzMap.integral_norm_sq_fourier]
  simpa [bumpSchwartz, Real.norm_eq_abs, abs_of_nonneg (bump_nonneg _)] using bump_l2norm

private lemma bumpFourier_sq_integrable :
    Integrable (fun u : ℝ => ‖bumpFourier u‖ ^ 2) := by
  by_contra h
  have hzero := bumpFourier_l2
  rw [integral_undef h] at hzero
  norm_num at hzero

private lemma bumpFourier_decay (K : ℕ) :
    ∃ C : ℝ, 0 < C ∧ ∀ u : ℝ, ‖bumpFourier u‖ ≤ C * (1 + |u|) ^ (-(K : ℝ)) := by
  let g : 𝓢(ℝ, ℂ) := 𝓕 bumpSchwartz
  have hfourier : ∀ u : ℝ, bumpFourier u = g u := by
    intro u
    exact congr_fun bumpFourier_eq_fourier u
  let c : ℝ := 2 ^ K * (Finset.Iic (K, 0)).sup
    (fun m => SchwartzMap.seminorm ℝ m.1 m.2) g
  refine ⟨|c| + 1, by positivity, ?_⟩
  intro u
  have hweight : (1 + |u|) ^ K * ‖bumpFourier u‖ ≤ c := by
    have hw := SchwartzMap.one_add_le_sup_seminorm_apply
      (𝕜 := ℝ) (m := (K, 0)) (k := K) (n := 0) le_rfl le_rfl g u
    rw [norm_iteratedFDeriv_zero, ← hfourier u] at hw
    simpa [c, Real.norm_eq_abs] using hw
  have hbase : 0 < 1 + |u| := by positivity
  have hpow : 0 < (1 + |u|) ^ K := by positivity
  have hrpow : (1 + |u|) ^ (-(K : ℝ)) = ((1 + |u|) ^ K)⁻¹ := by
    rw [← Real.rpow_natCast, Real.rpow_neg hbase.le]
  rw [hrpow, ← div_eq_mul_inv]
  apply (le_div_iff₀ hpow).2
  calc
    ‖bumpFourier u‖ * (1 + |u|) ^ K =
        (1 + |u|) ^ K * ‖bumpFourier u‖ := by ring
    _ ≤ c := hweight
    _ ≤ |c| + 1 := by linarith [le_abs_self c]

private lemma bumpFourier_lower_bound :
    ∃ c δ : ℝ, 0 < c ∧ 0 < δ ∧
    ∀ u : ℝ, |u| ≤ δ → c ≤ ‖bumpFourier u‖ ^ 2 := by
  have hcts : Continuous bumpFourier := bumpFourier_continuous
  have hpsi0_eq : (bumpFourier 0).re = ∫ x : ℝ, bump x := by
    simp only [bumpFourier, mul_zero, neg_zero]
    have hbumpc : Integrable (fun x : ℝ => (bump x : ℂ)) := bump_integrable.ofReal
    simpa using (integral_re hbumpc).symm
  have hpsi0 : 0 < (bumpFourier 0).re := by
    rw [hpsi0_eq]
    exact bump_integral_pos
  have hpos : 0 < ‖bumpFourier 0‖ := by
    rw [norm_pos_iff]
    intro hzero
    have : (bumpFourier 0).re = 0 := by rw [hzero]; rfl
    linarith
  set v₀ := ‖bumpFourier 0‖
  obtain ⟨δ, hδ, hball⟩ :=
    (Metric.continuousAt_iff.mp hcts.continuousAt) (v₀ / 2) (by linarith)
  refine ⟨(v₀ / 2) ^ 2, δ / 2, by positivity, by positivity, ?_⟩
  intro u hu
  have hdist : dist (bumpFourier u) (bumpFourier 0) < v₀ / 2 := by
    apply hball
    rw [Real.dist_eq]
    simp only [sub_zero]
    exact lt_of_le_of_lt hu (by linarith)
  have hbound : v₀ - ‖bumpFourier u‖ < v₀ / 2 := by
    calc
      v₀ - ‖bumpFourier u‖ = ‖bumpFourier 0‖ - ‖bumpFourier u‖ := rfl
      _ ≤ ‖bumpFourier 0 - bumpFourier u‖ := norm_sub_norm_le _ _
      _ = dist (bumpFourier u) (bumpFourier 0) := by
        rw [dist_eq_norm_sub, norm_sub_rev]
      _ < v₀ / 2 := hdist
  have hball' : v₀ / 2 < ‖bumpFourier u‖ := by linarith
  nlinarith [norm_nonneg (bumpFourier u)]

/-! ### Weighted Plancherel identity -/

-- Package translates of the concrete bump as Schwartz functions.
private def bumpShift (w : ℝ) : 𝓢(ℝ, ℂ) := by
  let f : ℝ → ℂ := fun x => (bump (x + w) : ℂ)
  have hcomp : HasCompactSupport f := by
    apply HasCompactSupport.of_support_subset_isCompact
      (isCompact_Icc (a := -w - 1 / 4) (b := -w + 1 / 4))
    intro x hx
    change (bump (x + w) : ℂ) ≠ 0 at hx
    have hx' : bump (x + w) ≠ 0 := by exact_mod_cast hx
    have h := abs_le.mp (bump_supp (x + w) hx')
    constructor <;> linarith
  have hsmooth : ContDiff ℝ ∞ f := by
    change ContDiff ℝ ∞ ((fun x : ℝ => (bump x : ℂ)) ∘ fun x => x + w)
    exact bump_complex_smooth.comp (contDiff_id.add contDiff_const)
  exact hcomp.toSchwartzMap hsmooth

private lemma bumpShift_apply (w x : ℝ) : bumpShift w x = bump (x + w) := rfl

private lemma fourier_bumpShift (w u : ℝ) :
    (𝓕 (bumpShift w : ℝ → ℂ)) u =
      𝐞 (w * u) * bumpFourier u := by
  have hcoe : (bumpShift w : ℝ → ℂ) = fun x : ℝ => (bump (x + w) : ℂ) := by
    ext x
    exact_mod_cast bumpShift_apply w x
  rw [hcoe, Real.fourier_real_eq]
  have hpsi : bumpFourier u =
      ∫ x : ℝ, 𝐞 (-(x * u)) * (bump x : ℂ) := by
    simp only [bumpFourier]
    apply integral_congr_ae
    filter_upwards with x
    ring
  rw [hpsi]
  have ht := congr_fun
    (Fourier.fourierIntegral_comp_add_right 𝐞 volume (fun x : ℝ => (bump x : ℂ)) w) u
  simpa [Fourier.fourierIntegral_def, Circle.smul_def, Real.fourierChar_apply] using ht

private lemma bumpShift_inner_eq_zero {v w : ℝ} (hvw : 1 ≤ |v - w|) :
    ∫ x : ℝ, inner ℂ (bumpShift v x) (bumpShift w x) = 0 := by
  apply integral_eq_zero_of_ae
  filter_upwards with x
  by_cases hxv : bump (x + v) = 0
  · simp [bumpShift_apply, hxv]
  by_cases hxw : bump (x + w) = 0
  · simp [bumpShift_apply, hxw]
  exfalso
  have hv := bump_supp (x + v) hxv
  have hw := bump_supp (x + w) hxw
  have hbound : |v - w| ≤ 1 / 2 := by
    calc
      |v - w| = |(x + v) - (x + w)| := by ring_nf
      _ ≤ |x + v| + |x + w| := abs_sub _ _
      _ ≤ 1 / 4 + 1 / 4 := by linarith
      _ = 1 / 2 := by norm_num
  linarith

private lemma bumpShift_inner_self (w : ℝ) :
    ∫ x : ℝ, inner ℂ (bumpShift w x) (bumpShift w x) = 1 := by
  calc
    ∫ x : ℝ, inner ℂ (bumpShift w x) (bumpShift w x) =
        ∫ x : ℝ, (((bump (x + w)) ^ 2 : ℝ) : ℂ) := by
          apply integral_congr_ae
          filter_upwards with x
          simp [bumpShift_apply, abs_of_nonneg (bump_nonneg _)]
    _ = ((∫ x : ℝ, (bump (x + w)) ^ 2 : ℝ) : ℂ) := integral_ofReal
    _ = ((∫ x : ℝ, (bump x) ^ 2 : ℝ) : ℂ) := by
      rw [integral_add_right_eq_self (fun x : ℝ => (bump x) ^ 2) w]
    _ = 1 := by rw [bump_l2norm]; norm_num

private lemma bumpShift_orthonormal {ι : Type*} [Finite ι]
    (ξ : ι → ℝ) (N : ℝ) (hN : 0 < N)
    (hsep : IsSeparatedFamily (1 / N) ξ) :
    Orthonormal ℂ (fun r => (bumpShift (N * ξ r)).toLp 2) := by
  classical
  letI := Fintype.ofFinite ι
  rw [orthonormal_iff_ite]
  intro r s
  rw [SchwartzMap.inner_toL2_toL2_eq
    (bumpShift (N * ξ r)) (bumpShift (N * ξ s)) volume]
  by_cases hrs : r = s
  · subst s
    rw [if_pos rfl, bumpShift_inner_self]
  · rw [if_neg hrs]
    apply bumpShift_inner_eq_zero
    rw [show N * ξ r - N * ξ s = N * (ξ r - ξ s) by ring, abs_mul, abs_of_pos hN]
    calc
      1 = N * (1 / N) := by field_simp
      _ ≤ N * |ξ r - ξ s| := mul_le_mul_of_nonneg_left (hsep hrs) hN.le

private theorem weighted_l2_identity {ι : Type*} [Fintype ι] (a : ι → ℂ) (ξ : ι → ℝ)
    (N : ℝ) (hN : 0 < N) (t₀ : ℝ)
    (hnorm : ∑ r, ‖a r‖ ^ 2 = 1)
    (hsep : IsSeparatedFamily (1 / N) ξ) :
    ∫ t : ℝ, expSumSq a ξ t * ‖bumpFourier ((t - t₀) / N)‖ ^ 2 = N := by
  classical
  let c : ι → ℂ := fun r => a r * 𝐞 (ξ r * t₀)
  let G : 𝓢(ℝ, ℂ) := ∑ r, c r • bumpShift (N * ξ r)
  have hfourier (u : ℝ) :
      (𝓕 G : 𝓢(ℝ, ℂ)) u = expSum a ξ (t₀ + N * u) * bumpFourier u := by
    change (SchwartzMap.fourierTransformCLM ℂ
      (∑ r, c r • bumpShift (N * ξ r))) u = _
    rw [map_sum]
    simp_rw [map_smul]
    have hsum_apply (s : Finset ι) (f : ι → 𝓢(ℝ, ℂ)) :
        (∑ r ∈ s, f r) u = ∑ r ∈ s, f r u := by
      induction s using Finset.induction_on with
      | empty => simp
      | @insert r s hrs ih =>
          simp [Finset.sum_insert, hrs, ih]
    rw [hsum_apply Finset.univ]
    simp_rw [smul_apply, smul_eq_mul]
    have hshift (r : ι) :
        SchwartzMap.fourierTransformCLM ℂ (bumpShift (N * ξ r)) u =
          𝐞 ((N * ξ r) * u) * bumpFourier u := by
      rw [SchwartzMap.fourierTransformCLM_apply]
      rw [SchwartzMap.fourier_coe]
      exact fourier_bumpShift (N * ξ r) u
    simp_rw [hshift]
    rw [show (∑ r, c r * (𝐞 ((N * ξ r) * u) * bumpFourier u)) =
        (∑ r, c r * 𝐞 ((N * ξ r) * u)) * bumpFourier u by
      rw [Finset.sum_mul]
      apply Finset.sum_congr rfl
      intro r _
      ring]
    congr 1
    apply Finset.sum_congr rfl
    intro r _
    simp only [c]
    calc
      a r * 𝐞 (ξ r * t₀) * 𝐞 ((N * ξ r) * u) =
          a r * (𝐞 (ξ r * t₀) * 𝐞 ((N * ξ r) * u)) := by ring
      _ = a r * 𝐞 (ξ r * (t₀ + N * u)) := by
        simp only [Real.fourierChar_apply]
        rw [← Complex.exp_add]
        congr 2
        push_cast
        ring
  have hG_toLp :
      G.toLp 2 = ∑ r, c r • (bumpShift (N * ξ r)).toLp 2 := by
    change SchwartzMap.toLpCLM ℂ ℂ 2 volume G = _
    simp [G]
  have hG_norm : ‖G.toLp 2‖ ^ 2 = 1 := by
    have horth := (bumpShift_orthonormal ξ N hN hsep).inner_sum c c Finset.univ
    have hc_norm : ∑ r, ‖c r‖ ^ 2 = 1 := by
      simpa only [c, norm_mul, Circle.norm_coe, mul_one] using hnorm
    have hsum_norm :
        ‖∑ r, c r • (bumpShift (N * ξ r)).toLp 2‖ ^ 2 = ∑ r, ‖c r‖ ^ 2 := by
      calc
        ‖∑ r, c r • (bumpShift (N * ξ r)).toLp 2‖ ^ 2 =
            (inner ℂ (∑ r, c r • (bumpShift (N * ξ r)).toLp 2)
              (∑ r, c r • (bumpShift (N * ξ r)).toLp 2)).re :=
                by
                  exact norm_sq_eq_re_inner (𝕜 := ℂ)
                    (∑ r, c r • (bumpShift (N * ξ r)).toLp 2)
        _ = (∑ r, (starRingEnd ℂ) (c r) * c r).re := congrArg Complex.re horth
        _ = ∑ r, ‖c r‖ ^ 2 := by
          rw [Complex.re_sum]
          apply Finset.sum_congr rfl
          intro r _
          rw [Complex.conj_mul']
          rw [← Complex.ofReal_pow, Complex.ofReal_re]
    rw [hG_toLp, hsum_norm, hc_norm]
  have hG_inner : inner ℂ (G.toLp 2) (G.toLp 2) = 1 := by
    rw [inner_self_eq_norm_sq_to_K]
    rw [← RCLike.ofReal_pow]
    exact congrArg (RCLike.ofReal : ℝ → ℂ) hG_norm
  have hG_l2 : ∫ x : ℝ, ‖G x‖ ^ 2 = 1 := by
    have hinner_integral : ∫ x : ℝ, inner ℂ (G x) (G x) = 1 := by
      rw [← SchwartzMap.inner_toL2_toL2_eq G G volume]
      exact hG_inner
    apply Complex.ofRealLI.injective
    simpa [← LinearIsometry.integral_comp_comm, inner_self_eq_norm_sq_to_K] using hinner_integral
  let H : ℝ → ℝ := fun u => ‖(𝓕 G : 𝓢(ℝ, ℂ)) u‖ ^ 2
  have hrewrite : (fun t : ℝ =>
      expSumSq a ξ t * ‖bumpFourier ((t - t₀) / N)‖ ^ 2) =
      fun t => H ((1 / N) * t + (-t₀ / N)) := by
    funext t
    rw [show (1 / N) * t + -t₀ / N = (t - t₀) / N by
      field_simp
      ring]
    have hu : t₀ + N * ((t - t₀) / N) = t := by
      field_simp
      ring
    simp only [H, hfourier, hu, expSumSq, norm_mul, mul_pow]
  rw [hrewrite]
  let K : ℝ → ℝ := fun x => H (x + (-t₀ / N))
  change ∫ t : ℝ, K ((1 / N) * t) = N
  rw [Measure.integral_comp_mul_left K (1 / N)]
  rw [show (∫ y : ℝ, K y) = ∫ y : ℝ, H y by
    exact integral_add_right_eq_self H (-t₀ / N)]
  rw [show (∫ y : ℝ, H y) = 1 by
    simpa only [H] using (SchwartzMap.integral_norm_sq_fourier G).trans hG_l2]
  simp [abs_of_pos hN]

/-! ### Local L² bound -/

private theorem local_l2_bound :
    ∃ C : ℝ, 0 < C ∧
    ∀ {ι : Type*} [Fintype ι] (a : ι → ℂ) (ξ : ι → ℝ) (N : ℝ),
    0 < N →
    (∑ r, ‖a r‖ ^ 2 = 1) →
    IsSeparatedFamily (1 / N) ξ →
    ∀ j₀ : ℝ,
    ∫ t in Set.Icc j₀ (j₀ + N), expSumSq a ξ t ≤ C * N := by
  obtain ⟨c, δ, hc, hδ, hlb⟩ := bumpFourier_lower_bound
  refine ⟨(1 + 1 / δ) / c, by positivity, ?_⟩
  intro ι _ a ξ N hN hnorm hsep j₀
  -- Enlarge the smoothing scale so that the whole interval lies in the
  -- neighbourhood on which `bumpFourier_lower_bound` applies.
  set M := N * (1 + 1 / δ)
  have hM : 0 < M := by
    dsimp [M]
    positivity
  have hNM : N ≤ M := by
    have hinvδ : 0 ≤ 1 / δ := by positivity
    dsimp [M]
    nlinarith
  have hsepM : IsSeparatedFamily (1 / M) ξ := by
    intro r s hrs
    apply le_trans _ (hsep hrs)
    apply (div_le_div_iff₀ hM hN).2
    simpa using hNM
  set t₀ := j₀ + N / 2
  have h2 := weighted_l2_identity a ξ M hM t₀ hnorm hsepM
  have hlb_J : ∀ t ∈ Set.Icc j₀ (j₀ + N),
      c ≤ ‖bumpFourier ((t - t₀) / M)‖ ^ 2 := by
    intro t ht
    apply hlb
    rw [abs_div, abs_of_pos hM]
    apply (div_le_iff₀ hM).2
    have habs : |t - t₀| ≤ N / 2 := by
      rw [abs_le]
      constructor <;> simp only [t₀] <;> linarith [ht.1, ht.2]
    apply habs.trans
    have hscale : δ * M = δ * N + N := by
      dsimp [M]
      field_simp
    rw [hscale]
    nlinarith
  have hFcont : Continuous (expSumSq a ξ) := expSumSq_continuous a ξ
  have hweighted : Integrable (fun t : ℝ =>
      expSumSq a ξ t * ‖bumpFourier ((t - t₀) / M)‖ ^ 2) := by
    by_contra h
    rw [integral_undef h] at h2
    linarith
  have hFub : c * ∫ t in Set.Icc j₀ (j₀ + N), expSumSq a ξ t ≤ M := by
    calc c * ∫ t in Set.Icc j₀ (j₀ + N), expSumSq a ξ t
        = ∫ t in Set.Icc j₀ (j₀ + N), c * expSumSq a ξ t :=
            (integral_const_mul _ _).symm
      _ ≤ ∫ t in Set.Icc j₀ (j₀ + N),
            expSumSq a ξ t * ‖bumpFourier ((t - t₀) / M)‖ ^ 2 := by
          apply setIntegral_mono_on
          · exact (continuous_const.mul hFcont).integrableOn_Icc
          · exact hweighted.integrableOn
          · exact measurableSet_Icc
          · intro t ht
            calc
              c * expSumSq a ξ t = expSumSq a ξ t * c := by ring
              _ ≤ expSumSq a ξ t * ‖bumpFourier ((t - t₀) / M)‖ ^ 2 :=
                mul_le_mul_of_nonneg_left (hlb_J t ht) (sq_nonneg _)
      _ ≤ ∫ t : ℝ, expSumSq a ξ t * ‖bumpFourier ((t - t₀) / M)‖ ^ 2 :=
          setIntegral_le_integral hweighted (ae_of_all _ fun t =>
            mul_nonneg (sq_nonneg _) (sq_nonneg _))
      _ = M := h2
  rw [show (1 + 1 / δ) / c * N = M / c by
    dsimp [M]
    ring]
  apply (le_div_iff₀ hc).2
  nlinarith

/-! ### Smoothing identity -/

-- E(t) = 1/N ∫_I |bump̂((t-t₀)/N)|² dt₀ - 1_I(t)
private def smoothingKernel (N : ℝ) (left right : ℝ) (t : ℝ) : ℝ :=
  (1 / N) * (∫ t₀ in Set.Icc left right, ‖bumpFourier ((t - t₀) / N)‖ ^ 2) -
  Set.indicator (Set.Icc left right) (fun _ => (1 : ℝ)) t

private lemma smoothingKernelAverage_eq_intervalIntegral (N : ℝ) (hN : 0 < N)
    (left right : ℝ) (hleft_right : left ≤ right) (t : ℝ) :
    (1 / N) * (∫ t₀ in Set.Icc left right, ‖bumpFourier ((t - t₀) / N)‖ ^ 2) =
    ∫ u in (t - right) / N..(t - left) / N, ‖bumpFourier u‖ ^ 2 := by
  rw [integral_Icc_eq_integral_Ioc, ← intervalIntegral.integral_of_le hleft_right]
  have hchange : (fun t₀ : ℝ => ‖bumpFourier ((t - t₀) / N)‖ ^ 2) =
      fun t₀ => ‖bumpFourier (t / N - t₀ / N)‖ ^ 2 := by
    funext t₀
    rw [sub_div]
  rw [hchange]
  simpa only [one_div, smul_eq_mul] using
    (intervalIntegral.inv_smul_integral_comp_sub_div
      (f := fun u : ℝ => ‖bumpFourier u‖ ^ 2) (a := left) (b := right) N (t / N)).trans (by
        congr 1 <;> field_simp)

private theorem smoothing_identity {ι : Type*} [Fintype ι] (a : ι → ℂ) (ξ : ι → ℝ)
    (N : ℝ) (hN : 0 < N)
    (hnorm : ∑ r, ‖a r‖ ^ 2 = 1)
    (hsep : IsSeparatedFamily (1 / N) ξ)
    (left right : ℝ) (T : ℝ) (hT : T = right - left) (hleft_right : left ≤ right) :
    ∫ t in Set.Icc left right, expSumSq a ξ t =
    T - ∫ t : ℝ, expSumSq a ξ t * smoothingKernel N left right t := by
  let F : ℝ → ℝ := expSumSq a ξ
  let K : ℝ → ℝ → ℝ := fun t₀ t => ‖bumpFourier ((t - t₀) / N)‖ ^ 2
  let H : ℝ × ℝ → ℝ := fun p => F p.2 * K p.1 p.2
  have hF_nonneg : ∀ t, 0 ≤ F t := expSumSq_nonneg a ξ
  have hK_nonneg : ∀ t₀ t, 0 ≤ K t₀ t := fun t₀ t => sq_nonneg _
  have hH_nonneg : ∀ p, 0 ≤ H p := fun p =>
    mul_nonneg (hF_nonneg p.2) (hK_nonneg p.1 p.2)
  have hF_cont : Continuous F := expSumSq_continuous a ξ
  have hbumpFourier_cont : Continuous bumpFourier := bumpFourier_continuous
  have hK_cont : Continuous (fun p : ℝ × ℝ => K p.1 p.2) := by
    exact ((hbumpFourier_cont.comp
      ((continuous_snd.sub continuous_fst).div_const N)).norm.pow 2)
  have hH_cont : Continuous H :=
    (hF_cont.comp continuous_snd).mul hK_cont
  -- From (3.1): ∫_ℝ F(t)|bump̂((t-t₀)/N)|² dt = N for each t₀.
  have h31 : ∀ t₀ : ℝ,
      ∫ t : ℝ, H (t₀, t) = N := by
    intro t₀
    simpa only [H, F, K] using weighted_l2_identity a ξ N hN t₀ hnorm hsep
  have hslice : ∀ t₀ : ℝ, Integrable (fun t => H (t₀, t)) := by
    intro t₀
    by_contra h
    have hzero := h31 t₀
    rw [integral_undef h] at hzero
    linarith
  -- The product-space integrability follows from the already evaluated slices.
  have hH_int : Integrable H ((volume.restrict (Set.Icc left right)).prod volume) := by
    apply (integrable_prod_iff hH_cont.aestronglyMeasurable).2
    constructor
    · exact ae_of_all _ hslice
    · have hmarginal : (fun t₀ : ℝ => ∫ t : ℝ, ‖H (t₀, t)‖) = fun _ => N := by
        funext t₀
        rw [show (fun t : ℝ => ‖H (t₀, t)‖) = fun t => H (t₀, t) by
          funext t
          exact Real.norm_of_nonneg (hH_nonneg (t₀, t))]
        exact h31 t₀
      rw [hmarginal]
      exact integrableOn_const measure_Icc_lt_top.ne
  -- Integrate the constant marginal over t₀ ∈ I.
  have hNT : ∫ _ in Set.Icc left right, N = N * T := by
    rw [setIntegral_const, Real.volume_real_Icc_of_le hleft_right, smul_eq_mul, hT]
    ring
  -- Fubini, with the interval restriction built into the first measure.
  have hswap :
      (∫ t₀ in Set.Icc left right, ∫ t : ℝ, H (t₀, t)) =
      ∫ t : ℝ, ∫ t₀ in Set.Icc left right, H (t₀, t) := by
    exact integral_integral_swap hH_int
  have hFubini : ∫ t : ℝ, F t * (∫ t₀ in Set.Icc left right, K t₀ t) = N * T := by
    calc
      ∫ t : ℝ, F t * (∫ t₀ in Set.Icc left right, K t₀ t) =
          ∫ t : ℝ, ∫ t₀ in Set.Icc left right, H (t₀, t) := by
            congr 1
            ext t
            change F t * (∫ t₀ in Set.Icc left right, K t₀ t) =
              ∫ t₀ in Set.Icc left right, F t * K t₀ t
            rw [integral_const_mul]
      _ = ∫ t₀ in Set.Icc left right, ∫ t : ℝ, H (t₀, t) := hswap.symm
      _ = ∫ _ in Set.Icc left right, N := by
        apply setIntegral_congr_fun measurableSet_Icc
        intro t₀ _
        exact h31 t₀
      _ = N * T := hNT
  have hFK_int : Integrable (fun t : ℝ => F t *
      (∫ t₀ in Set.Icc left right, K t₀ t)) := by
    have hmarginal := hH_int.integral_prod_right
    apply hmarginal.congr
    apply ae_of_all
    intro t
    change (∫ t₀ in Set.Icc left right, F t * K t₀ t) =
      F t * (∫ t₀ in Set.Icc left right, K t₀ t)
    rw [integral_const_mul]
  have hFI_int : Integrable ((Set.Icc left right).indicator F) :=
    hF_cont.integrableOn_Icc.integrable_indicator measurableSet_Icc
  have hexpand : (fun t : ℝ => expSumSq a ξ t * smoothingKernel N left right t) =
      fun t => (1 / N) * (F t * (∫ t₀ in Set.Icc left right, K t₀ t)) -
        (Set.Icc left right).indicator F t := by
    funext t
    simp only [smoothingKernel, F, K]
    by_cases ht : t ∈ Set.Icc left right
    · rw [Set.indicator_of_mem ht, Set.indicator_of_mem ht]
      ring
    · rw [Set.indicator_of_notMem ht, Set.indicator_of_notMem ht]
      ring
  rw [hexpand, integral_sub (hFK_int.const_mul _) hFI_int,
    integral_const_mul, hFubini, integral_indicator measurableSet_Icc]
  field_simp
  ring

/-! ### Decay of the smoothing error -/

private theorem smoothing_kernel_decay :
    ∃ C : ℝ, 0 < C ∧
    ∀ N : ℝ, 0 < N →
    ∀ left right : ℝ, left ≤ right →
    ∀ t : ℝ,
    |smoothingKernel N left right t| ≤
    C * (1 + min (|t - left| / N) (|t - right| / N)) ^ (-(10 : ℝ)) := by
  let f : ℝ → ℝ := fun u => ‖bumpFourier u‖ ^ 2
  have hf_nonneg : ∀ u, 0 ≤ f u := fun u => sq_nonneg _
  have hf_int : Integrable f := by
    simpa only [f] using bumpFourier_sq_integrable
  obtain ⟨C₀, hC₀, hdecay⟩ := bumpFourier_decay 6
  let w : ℝ → ℝ := fun u => (1 + |u|) ^ (-(2 : ℝ))
  have hw_nonneg : ∀ u, 0 ≤ w u := fun u => Real.rpow_nonneg (by positivity) _
  have hw_int : Integrable w := by
    simpa only [w, Real.norm_eq_abs] using
      (integrable_one_add_norm (E := ℝ) (μ := volume) (r := 2) (by norm_num))
  let A := ∫ u : ℝ, w u
  have hA : 0 ≤ A := integral_nonneg hw_nonneg
  let C := C₀ ^ 2 * (A + 1)
  have hC : 0 < C := by
    dsimp [C]
    positivity
  have htail : ∀ d : ℝ, 0 ≤ d →
      ∫ u in {u : ℝ | d ≤ |u|}, f u ≤ C * (1 + d) ^ (-(10 : ℝ)) := by
    intro d hd
    let q := (1 + d) ^ (-(10 : ℝ))
    have hq : 0 ≤ q := Real.rpow_nonneg (by positivity) _
    have hmajorant : Integrable (fun u => C₀ ^ 2 * q * w u) :=
      by simpa only [mul_assoc] using (hw_int.const_mul q).const_mul (C₀ ^ 2)
    calc
      ∫ u in {u : ℝ | d ≤ |u|}, f u ≤
          ∫ u in {u : ℝ | d ≤ |u|}, C₀ ^ 2 * q * w u := by
        apply setIntegral_mono_on hf_int.integrableOn hmajorant.integrableOn
          (measurableSet_le measurable_const (continuous_abs.measurable))
        intro u hu
        change d ≤ |u| at hu
        have hsquare : f u ≤ C₀ ^ 2 * (1 + |u|) ^ (-(12 : ℝ)) := by
          dsimp [f]
          have hd' := hdecay u
          have hp : 0 ≤ (1 + |u|) ^ (-(6 : ℝ)) := Real.rpow_nonneg (by positivity) _
          have hp_sq : (1 + |u|) ^ (-(12 : ℝ)) =
              ((1 + |u|) ^ (-(6 : ℝ))) ^ 2 := by
            rw [pow_two, ← Real.rpow_add (by positivity)]
            norm_num
          calc
            ‖bumpFourier u‖ ^ 2 ≤ (C₀ * (1 + |u|) ^ (-(6 : ℝ))) ^ 2 :=
              (sq_le_sq₀ (norm_nonneg _) (mul_nonneg hC₀.le hp)).2 hd'
            _ = C₀ ^ 2 * ((1 + |u|) ^ (-(6 : ℝ))) ^ 2 := by ring
            _ = C₀ ^ 2 * (1 + |u|) ^ (-(12 : ℝ)) := by rw [hp_sq]
        have hfactor : (1 + |u|) ^ (-(12 : ℝ)) =
            w u * (1 + |u|) ^ (-(10 : ℝ)) := by
          dsimp [w]
          rw [← Real.rpow_add (by positivity)]
          norm_num
        have hmono : (1 + |u|) ^ (-(10 : ℝ)) ≤ q := by
          dsimp [q]
          exact Real.rpow_le_rpow_of_nonpos (by positivity) (by linarith) (by norm_num)
        calc
          f u ≤ C₀ ^ 2 * (1 + |u|) ^ (-(12 : ℝ)) := hsquare
          _ = C₀ ^ 2 * (w u * (1 + |u|) ^ (-(10 : ℝ))) := by rw [hfactor]
          _ ≤ C₀ ^ 2 * (w u * q) := by gcongr
          _ = C₀ ^ 2 * q * w u := by ring
      _ ≤ ∫ u : ℝ, C₀ ^ 2 * q * w u :=
        setIntegral_le_integral hmajorant (ae_of_all _ fun u => by positivity)
      _ = C₀ ^ 2 * q * A := by simp only [integral_const_mul, A]
      _ ≤ C * q := by
        dsimp [C]
        nlinarith [mul_nonneg (sq_nonneg C₀) hq]
      _ = C * (1 + d) ^ (-(10 : ℝ)) := rfl
  refine ⟨C, hC, ?_⟩
  intro N hN left right hleft_right t
  let uLower := (t - right) / N
  let uUpper := (t - left) / N
  let d := min |uUpper| |uLower|
  have hLowerUpper : uLower ≤ uUpper := by
    dsimp [uLower, uUpper]
    exact div_le_div_of_nonneg_right (sub_le_sub_left hleft_right t) hN.le
  have hd : 0 ≤ d := le_min (abs_nonneg _) (abs_nonneg _)
  have hsubst : (1 / N) * (∫ t₀ in Set.Icc left right, f ((t - t₀) / N)) =
      ∫ u in Set.Icc uLower uUpper, f u := by
    simpa only [f, uLower, uUpper, intervalIntegral.integral_of_le hLowerUpper,
      integral_Icc_eq_integral_Ioc] using
      smoothingKernelAverage_eq_intervalIntegral N hN left right hleft_right t
  have hd_eq : d = min (|t - left| / N) (|t - right| / N) := by
    dsimp [d, uLower, uUpper]
    rw [abs_div, abs_div, abs_of_pos hN]
  by_cases ht : t ∈ Set.Icc left right
  · have huLower_nonpos : uLower ≤ 0 := by
      dsimp [uLower]
      exact div_nonpos_of_nonpos_of_nonneg (sub_nonpos.mpr ht.2) hN.le
    have huUpper_nonneg : 0 ≤ uUpper := by
      dsimp [uUpper]
      exact div_nonneg (sub_nonneg.mpr ht.1) hN.le
    have hsubset : (Set.Icc uLower uUpper)ᶜ ⊆ {u : ℝ | d ≤ |u|} := by
      intro u hu
      simp only [Set.mem_compl_iff, Set.mem_Icc, not_and_or, not_le] at hu
      simp only [Set.mem_setOf_eq]
      rcases hu with hu | hu
      · rw [abs_of_nonpos (hu.le.trans huLower_nonpos)]
        exact (min_le_right |uUpper| |uLower|).trans (by
          rw [abs_of_nonpos huLower_nonpos]
          linarith)
      · rw [abs_of_nonneg (huUpper_nonneg.trans hu.le)]
        exact (min_le_left |uUpper| |uLower|).trans (by
          rw [abs_of_nonneg huUpper_nonneg]
          linarith)
    have hcomp : (∫ u in Set.Icc uLower uUpper, f u) - 1 =
        -(∫ u in (Set.Icc uLower uUpper)ᶜ, f u) := by
      rw [setIntegral_compl measurableSet_Icc hf_int, bumpFourier_l2]
      ring
    rw [smoothingKernel, Set.indicator_of_mem ht]
    change |(1 / N) * (∫ t₀ in Set.Icc left right, f ((t - t₀) / N)) - 1| ≤ _
    rw [hsubst, hcomp, abs_neg, abs_of_nonneg
      (setIntegral_nonneg measurableSet_Icc.compl fun u _ => hf_nonneg u)]
    rw [← hd_eq]
    exact (setIntegral_mono_set hf_int.integrableOn (ae_of_all _ hf_nonneg)
      (ae_of_all _ hsubset)).trans (htail d hd)
  · have hzero : (0 : ℝ) ∉ Set.Icc uLower uUpper := by
      intro h0
      apply ht
      constructor
      · apply sub_nonneg.mp
        rcases div_nonneg_iff.mp h0.2 with h | h
        · exact h.1
        · exact (hN.not_ge h.2).elim
      · apply sub_nonpos.mp
        rcases div_nonpos_iff.mp h0.1 with h | h
        · exact (hN.not_ge h.2).elim
        · exact h.1
    have hsubset : Set.Icc uLower uUpper ⊆ {u : ℝ | d ≤ |u|} := by
      intro u hu
      rcases hu with ⟨hu1, hu2⟩
      simp only [Set.mem_setOf_eq]
      have hzero' : 0 < uLower ∨ uUpper < 0 := by
        simpa only [Set.mem_Icc, not_and_or, not_le] using hzero
      rcases hzero' with hLowerPos | hUpperNeg
      · rw [abs_of_nonneg (hLowerPos.le.trans hu1)]
        exact (min_le_right |uUpper| |uLower|).trans (by
          rw [abs_of_nonneg hLowerPos.le]
          linarith)
      · rw [abs_of_nonpos (hu2.trans hUpperNeg.le)]
        exact (min_le_left |uUpper| |uLower|).trans (by
          rw [abs_of_nonpos hUpperNeg.le]
          linarith)
    rw [smoothingKernel, Set.indicator_of_notMem ht]
    change |(1 / N) * (∫ t₀ in Set.Icc left right, f ((t - t₀) / N)) - 0| ≤ _
    rw [sub_zero, hsubst, abs_of_nonneg
      (setIntegral_nonneg measurableSet_Icc fun u _ => hf_nonneg u)]
    rw [← hd_eq]
    exact (setIntegral_mono_set hf_int.integrableOn (ae_of_all _ hf_nonneg)
      (ae_of_all _ hsubset)).trans (htail d hd)

private lemma smoothingKernelAverage_continuous (N : ℝ) (hN : 0 < N)
    (left right : ℝ) (hleft_right : left ≤ right) :
    Continuous (fun t : ℝ =>
      (1 / N) * (∫ t₀ in Set.Icc left right, ‖bumpFourier ((t - t₀) / N)‖ ^ 2)) := by
  let f : ℝ → ℝ := fun u => ‖bumpFourier u‖ ^ 2
  have hf_int : Integrable f := by
    simpa only [f] using bumpFourier_sq_integrable
  let P : ℝ → ℝ := fun x => ∫ u in 0..x, f u
  have hP_cont : Continuous P :=
    intervalIntegral.continuous_primitive (fun x y => hf_int.intervalIntegrable) 0
  have havg : (fun t : ℝ =>
      (1 / N) * (∫ t₀ in Set.Icc left right, ‖bumpFourier ((t - t₀) / N)‖ ^ 2)) =
      fun t => P ((t - left) / N) - P ((t - right) / N) := by
    funext t
    rw [smoothingKernelAverage_eq_intervalIntegral N hN left right hleft_right t]
    dsimp [P]
    exact (intervalIntegral.integral_interval_sub_left
      hf_int.intervalIntegrable hf_int.intervalIntegrable).symm
  rw [havg]
  exact (hP_cont.comp ((continuous_id.sub continuous_const).div_const N)).sub
    (hP_cont.comp ((continuous_id.sub continuous_const).div_const N))

private lemma smoothingKernel_aestronglyMeasurable (N : ℝ) (hN : 0 < N)
    (left right : ℝ) (hleft_right : left ≤ right) :
    AEStronglyMeasurable (smoothingKernel N left right) := by
  have hind_meas : AEStronglyMeasurable
      (Set.indicator (Set.Icc left right) (fun _ => (1 : ℝ))) :=
    (Measurable.indicator measurable_const measurableSet_Icc).aestronglyMeasurable
  exact (smoothingKernelAverage_continuous N hN left right hleft_right).aestronglyMeasurable.sub
    hind_meas

/-! ### Integrated smoothing-error bound -/

private theorem smoothing_error_bound :
    ∃ C : ℝ, 0 < C ∧
    ∀ {ι : Type*} [Fintype ι] (a : ι → ℂ) (ξ : ι → ℝ) (N : ℝ),
    0 < N →
    (∑ r, ‖a r‖ ^ 2 = 1) →
    IsSeparatedFamily (1 / N) ξ →
    ∀ left right : ℝ, left ≤ right →
    |∫ t : ℝ, expSumSq a ξ t * smoothingKernel N left right t| ≤ C * N := by
  obtain ⟨C₃, hC₃, hlocal_bound⟩ := local_l2_bound
  obtain ⟨C₅, hC₅, hkernel_bound⟩ := smoothing_kernel_decay
  let b : ℤ → ℝ := fun k => |(k : ℝ) + 1 / 2| ^ (-(10 : ℝ))
  have hb_nonneg : ∀ k, 0 ≤ b k :=
    fun k => Real.rpow_nonneg (abs_nonneg _) _
  have hb_summable : Summable b := by
    have h := (Real.summable_one_div_int_add_rpow (1 / 2) 10).2 (by norm_num)
    apply h.congr
    intro k
    dsimp [b]
    rw [Real.rpow_neg (abs_nonneg _)]
    rw [one_div]
  let B := ∑' k : ℤ, b k
  have hB : 0 ≤ B := tsum_nonneg hb_nonneg
  refine ⟨2 * C₅ * C₃ * (B + 1), by positivity, ?_⟩
  intro ι _ a ξ N hN hnorm hsep left right hleft_right
  let F : ℝ → ℝ := expSumSq a ξ
  have hF_nonneg : ∀ t, 0 ≤ F t := expSumSq_nonneg a ξ
  have hF_cont : Continuous F := expSumSq_continuous a ξ
  have hlocal := hlocal_bound a ξ N hN hnorm hsep
  have hE := hkernel_bound N hN left right hleft_right
  let weight : ℝ → ℝ → ℝ :=
    fun c t => (1 + |t - c| / N) ^ (-(10 : ℝ))
  have hweight_nonneg : ∀ c t, 0 ≤ weight c t :=
    fun c t => Real.rpow_nonneg (by positivity) _
  have hweight_cont : ∀ c, Continuous (weight c) := by
    intro c
    dsimp [weight]
    apply Continuous.rpow_const
    · fun_prop
    · intro t
      left
      positivity
  let W : ℝ → ℝ → ℝ := fun c t => F t * weight c t
  have hW_nonneg : ∀ c t, 0 ≤ W c t :=
    fun c t => mul_nonneg (hF_nonneg t) (hweight_nonneg c t)
  have hW_cont : ∀ c, Continuous (W c) :=
    fun c => hF_cont.mul (hweight_cont c)
  let cell : ℝ → ℤ → Set ℝ :=
    fun c k => Set.Ico (c + k • N) (c + (k + 1) • N)
  have hcell_meas : ∀ c k, MeasurableSet (cell c k) :=
    fun c k => measurableSet_Ico
  have hcell_pairwise : ∀ c, Pairwise (fun i j => Disjoint (cell c i) (cell c j)) := by
    intro c
    simpa only [cell] using Set.pairwise_disjoint_Ico_add_zsmul c N
  have hcell_cover : ∀ c, ⋃ k : ℤ, cell c k = Set.univ := by
    intro c
    simpa only [cell] using iUnion_Ico_add_zsmul hN c
  have hcell_weight : ∀ c k t, t ∈ cell c k → weight c t ≤ b k := by
    intro c k t ht
    have ht' : c + (k : ℝ) * N ≤ t ∧ t < c + ((k : ℝ) + 1) * N := by
      change c + k • N ≤ t ∧ t < c + (k + 1) • N at ht
      simpa only [zsmul_eq_mul, Int.cast_add, Int.cast_one] using ht
    let x := (t - c) / N
    have hx_left : (k : ℝ) ≤ x := by
      apply (le_div_iff₀ hN).2
      linarith [ht'.1]
    have hx_right : x < (k : ℝ) + 1 := by
      apply (div_lt_iff₀ hN).2
      linarith [ht'.2]
    have hmid_pos : 0 < |(k : ℝ) + 1 / 2| := by
      rw [abs_pos]
      intro hk
      have hk' : (2 : ℝ) * (k : ℝ) = -1 := by linarith
      have hk'' : (2 * k : ℤ) = -1 := by exact_mod_cast hk'
      omega
    have hmid : |((k : ℝ) + 1 / 2) - x| ≤ 1 / 2 := by
      rw [abs_le]
      constructor <;> linarith
    have hbase : |(k : ℝ) + 1 / 2| ≤ 1 + |x| := by
      calc
        |(k : ℝ) + 1 / 2| = |(((k : ℝ) + 1 / 2) - x) + x| := by ring_nf
        _ ≤ |((k : ℝ) + 1 / 2) - x| + |x| := abs_add_le _ _
        _ ≤ 1 / 2 + |x| := by linarith
        _ ≤ 1 + |x| := by linarith
    have hrpow := Real.rpow_le_rpow_of_nonpos hmid_pos hbase (by norm_num : (-(10 : ℝ)) ≤ 0)
    have hxabs : |t - c| / N = |x| := by
      dsimp [x]
      rw [abs_div, abs_of_pos hN]
    simpa only [weight, b, hxabs] using hrpow
  have hcell_bound : ∀ c k,
      ∫ t in cell c k, W c t ≤ b k * (C₃ * N) := by
    intro c k
    let cellLeft := c + k • N
    have hend : c + (k + 1) • N = cellLeft + N := by
      dsimp [cellLeft]
      rw [add_smul, one_smul]
      ring
    have hcell_eq : cell c k = Set.Ico cellLeft (cellLeft + N) := by
      dsimp [cell]
      rw [hend]
    have hWcell : IntegrableOn (W c) (cell c k) :=
      (hW_cont c).integrableOn_Icc.mono_set (by
        rw [hcell_eq]
        exact Set.Ico_subset_Icc_self)
    have hbFcell : IntegrableOn (fun t => b k * F t) (cell c k) :=
      ((continuous_const.mul hF_cont).integrableOn_Icc).mono_set (by
        rw [hcell_eq]
        exact Set.Ico_subset_Icc_self)
    have hFcell :
        ∫ t in cell c k, F t ≤ ∫ t in Set.Icc cellLeft (cellLeft + N), F t := by
      apply setIntegral_mono_set hF_cont.integrableOn_Icc
        (ae_of_all _ fun t => hF_nonneg t)
      apply ae_of_all
      rw [hcell_eq]
      exact Set.Ico_subset_Icc_self
    calc
      ∫ t in cell c k, W c t ≤ ∫ t in cell c k, b k * F t := by
        apply setIntegral_mono_on hWcell hbFcell (hcell_meas c k)
        intro t ht
        dsimp [W]
        calc
          F t * weight c t ≤ F t * b k :=
            mul_le_mul_of_nonneg_left (hcell_weight c k t ht) (hF_nonneg t)
          _ = b k * F t := by ring
      _ = b k * ∫ t in cell c k, F t := integral_const_mul _ _
      _ ≤ b k * ∫ t in Set.Icc cellLeft (cellLeft + N), F t :=
        mul_le_mul_of_nonneg_left hFcell (hb_nonneg k)
      _ ≤ b k * (C₃ * N) :=
        mul_le_mul_of_nonneg_left (hlocal cellLeft) (hb_nonneg k)
  have hW : ∀ c, Integrable (W c) ∧ ∫ t : ℝ, W c t ≤ B * (C₃ * N) := by
    intro c
    have hcell_int : ∀ k, IntegrableOn (W c) (cell c k) := by
      intro k
      exact (hW_cont c).integrableOn_Icc.mono_set Set.Ico_subset_Icc_self
    have hmajor_summable : Summable (fun k => b k * (C₃ * N)) :=
      hb_summable.mul_right (C₃ * N)
    have hsum_norm : Summable (fun k => ∫ t in cell c k, ‖W c t‖) := by
      refine Summable.of_nonneg_of_le
        (f := fun k => b k * (C₃ * N)) (g := fun k => ∫ t in cell c k, ‖W c t‖) ?_ ?_
        hmajor_summable
      · intro k
        exact setIntegral_nonneg (hcell_meas c k) fun t _ => norm_nonneg _
      · intro k
        rw [show (fun t => ‖W c t‖) = W c by
          funext t
          exact Real.norm_of_nonneg (hW_nonneg c t)]
        exact hcell_bound c k
    have hsum : Summable (fun k => ∫ t in cell c k, W c t) := by
      apply hsum_norm.congr
      intro k
      rw [show (fun t => ‖W c t‖) = W c by
        funext t
        exact Real.norm_of_nonneg (hW_nonneg c t)]
    have hW_int : Integrable (W c) := by
      have hu := integrableOn_iUnion_of_summable_integral_norm hcell_int hsum_norm
      rw [hcell_cover c] at hu
      simpa only [integrableOn_univ] using hu
    refine ⟨hW_int, ?_⟩
    have hpartition :
        ∫ t : ℝ, W c t = ∑' k : ℤ, ∫ t in cell c k, W c t := by
      rw [← setIntegral_univ, ← hcell_cover c]
      exact integral_iUnion (hcell_meas c) (hcell_pairwise c) hW_int.integrableOn
    rw [hpartition]
    calc
      ∑' k : ℤ, ∫ t in cell c k, W c t ≤ ∑' k : ℤ, b k * (C₃ * N) :=
        hsum.tsum_le_tsum (hcell_bound c) hmajor_summable
      _ = B * (C₃ * N) := by
        exact hb_summable.tsum_mul_right (C₃ * N)
  have hmin_split : ∀ t,
      (1 + min (|t - left| / N) (|t - right| / N)) ^ (-(10 : ℝ)) ≤
      weight left t + weight right t := by
    intro t
    rcases le_total (|t - left| / N) (|t - right| / N) with h | h
    · rw [min_eq_left h]
      dsimp [weight]
      exact le_add_of_nonneg_right (Real.rpow_nonneg (by positivity) _)
    · rw [min_eq_right h]
      dsimp [weight]
      exact le_add_of_nonneg_left (Real.rpow_nonneg (by positivity) _)
  have hmajor : ∀ t,
      F t * |smoothingKernel N left right t| ≤ C₅ * (W left t + W right t) := by
    intro t
    calc
      F t * |smoothingKernel N left right t| ≤
          F t * (C₅ * (1 + min (|t - left| / N) (|t - right| / N)) ^ (-(10 : ℝ))) :=
        mul_le_mul_of_nonneg_left (hE t) (hF_nonneg t)
      _ ≤ F t * (C₅ * (weight left t + weight right t)) :=
        mul_le_mul_of_nonneg_left
          (mul_le_mul_of_nonneg_left (hmin_split t) hC₅.le) (hF_nonneg t)
      _ = C₅ * (W left t + W right t) := by
        dsimp [W]
        ring
  have hWLeft := hW left
  have hWRight := hW right
  have hmajor_int : Integrable (fun t => C₅ * (W left t + W right t)) :=
    (hWLeft.1.add hWRight.1).const_mul C₅
  have hkernel_meas := smoothingKernel_aestronglyMeasurable N hN left right hleft_right
  have hkernel_abs_meas : AEStronglyMeasurable
      (fun t => |smoothingKernel N left right t|) := by
    simpa only [Real.norm_eq_abs] using hkernel_meas.norm
  have hFEabs_meas : AEStronglyMeasurable
      (fun t => F t * |smoothingKernel N left right t|) :=
    hF_cont.aestronglyMeasurable.mul hkernel_abs_meas
  have hFEabs : Integrable (fun t => F t * |smoothingKernel N left right t|) := by
    apply hmajor_int.mono' hFEabs_meas
    filter_upwards with t
    rw [Real.norm_eq_abs, abs_of_nonneg
      (mul_nonneg (hF_nonneg t) (abs_nonneg _))]
    exact hmajor t
  have hFE : Integrable (fun t => F t * smoothingKernel N left right t) := by
    apply hFEabs.mono' (hF_cont.aestronglyMeasurable.mul hkernel_meas)
    filter_upwards with t
    change |F t * smoothingKernel N left right t| ≤ F t * |smoothingKernel N left right t|
    rw [abs_mul, abs_of_nonneg (hF_nonneg t)]
  calc
    |∫ t : ℝ, expSumSq a ξ t * smoothingKernel N left right t| =
        |∫ t : ℝ, F t * smoothingKernel N left right t| := by rfl
    _ ≤ ∫ t : ℝ, |F t * smoothingKernel N left right t| :=
      abs_integral_le_integral_abs
    _ = ∫ t : ℝ, F t * |smoothingKernel N left right t| := by
      congr 1
      funext t
      rw [abs_mul, abs_of_nonneg (hF_nonneg t)]
    _ ≤ ∫ t : ℝ, C₅ * (W left t + W right t) :=
      integral_mono hFEabs hmajor_int hmajor
    _ = C₅ * ((∫ t : ℝ, W left t) + ∫ t : ℝ, W right t) := by
      rw [integral_const_mul, integral_add hWLeft.1 hWRight.1]
    _ ≤ C₅ * (B * (C₃ * N) + B * (C₃ * N)) := by
      gcongr
      · exact hWLeft.2
      · exact hWRight.2
    _ ≤ 2 * C₅ * C₃ * (B + 1) * N := by
      calc
        C₅ * (B * (C₃ * N) + B * (C₃ * N)) = (2 * C₅ * C₃ * N) * B := by ring
        _ ≤ (2 * C₅ * C₃ * N) * (B + 1) :=
          mul_le_mul_of_nonneg_left (le_add_of_nonneg_right zero_le_one) (by positivity)
        _ = 2 * C₅ * C₃ * (B + 1) * N := by ring

/-! ### Lemma 3.1 -/

/-- **Lemma 3.1 (L² integral estimate).** If `ξ` is a finite `1 / N`-separated
  family of real numbers, then over any interval of length `T`,
  `∫ |∑ r, a r * 𝐞 (ξ r * t)|² dt =
    (T + O(N)) * ∑ r, ‖a r‖²`.
  The conclusion expresses `O(N)` as `θ * N`, with `θ` bounded by a universal
  constant. -/
theorem l2_integral_estimate :
    ∃ C : ℝ, 0 < C ∧
    ∀ {ι : Type*} [Fintype ι] (a : ι → ℂ) (ξ : ι → ℝ) (N : ℝ),
    0 < N →
    IsSeparatedFamily (1 / N) ξ →
    ∀ left right T : ℝ,
    T = right - left →
    left ≤ right →
    ∃ θ : ℝ, |θ| ≤ C ∧
    ∫ t in Set.Icc left right,
      ‖∑ r, a r * 𝐞 (ξ r * t)‖ ^ 2 =
    (T + θ * N) * ∑ r, ‖a r‖ ^ 2 := by
  obtain ⟨C, hC, herror⟩ := smoothing_error_bound
  refine ⟨C, hC, ?_⟩
  intro ι _ a ξ N hN hsep left right T hT hleft_right
  change ∃ θ : ℝ, |θ| ≤ C ∧
    ∫ t in Set.Icc left right, expSumSq a ξ t =
    (T + θ * N) * ∑ r, ‖a r‖ ^ 2
  set M := ∑ r, ‖a r‖ ^ 2
  by_cases hM0 : M = 0
  · refine ⟨0, by simpa using hC.le, ?_⟩
    have hsum_zero : ∑ r, ‖a r‖ ^ 2 = 0 := by
      simpa only [M] using hM0
    have hterm_zero : ∀ r, ‖a r‖ ^ 2 = 0 := by
      intro r
      exact ((Finset.sum_eq_zero_iff_of_nonneg
        (s := Finset.univ) (f := fun i => ‖a i‖ ^ 2)
        (fun i _ => sq_nonneg _)).1 hsum_zero) r (Finset.mem_univ r)
    have hzero : ∀ r, a r = 0 := by
      intro r
      apply norm_eq_zero.mp
      nlinarith [hterm_zero r, norm_nonneg (a r)]
    simp [expSumSq, expSum, hzero, hM0]
  · have hMpos : 0 < M :=
      lt_of_le_of_ne (Finset.sum_nonneg fun r _ => sq_nonneg _) (Ne.symm hM0)
    set A : ι → ℂ := fun r => a r / (Real.sqrt M : ℂ)
    have hAnorm : ∑ r, ‖A r‖ ^ 2 = 1 := normalize_coefficients a M rfl hMpos
    have hsqrt_ne : (Real.sqrt M : ℂ) ≠ 0 :=
      Complex.ofReal_ne_zero.mpr (Real.sqrt_ne_zero'.mpr hMpos)
    have hsqrt_mul_A : ∀ r, (Real.sqrt M : ℂ) * A r = a r := by
      intro r
      dsimp [A]
      rw [mul_comm]
      exact div_mul_cancel₀ _ hsqrt_ne
    have hsum_scale : ∀ t,
        expSum a ξ t = (Real.sqrt M : ℂ) * expSum A ξ t := by
      intro t
      simp only [expSum]
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro r _
      rw [← mul_assoc, hsqrt_mul_A r]
    have hrescale : ∀ t, expSumSq a ξ t = M * expSumSq A ξ t := by
      intro t
      simp only [expSumSq]
      rw [hsum_scale t, norm_mul, Complex.norm_real,
        Real.norm_eq_abs, abs_of_nonneg (Real.sqrt_nonneg M), mul_pow,
        Real.sq_sqrt hMpos.le]
    have h4 := smoothing_identity A ξ N hN hAnorm hsep left right T hT hleft_right
    have h6 := herror A ξ N hN hAnorm hsep left right hleft_right
    set err := ∫ t : ℝ, expSumSq A ξ t * smoothingKernel N left right t
    have hA_eq : ∫ t in Set.Icc left right, expSumSq A ξ t = T - err := by
      simpa only [err] using h4
    have herr_bd : |err| ≤ C * N := by
      simpa only [err] using h6
    refine ⟨-err / N, ?_, ?_⟩
    · rw [abs_div, abs_neg, abs_of_pos hN, div_le_iff₀ hN]
      exact herr_bd
    · calc
        ∫ t in Set.Icc left right, expSumSq a ξ t =
            ∫ t in Set.Icc left right, M * expSumSq A ξ t :=
          setIntegral_congr_fun measurableSet_Icc (fun t _ => hrescale t)
        _ = M * ∫ t in Set.Icc left right, expSumSq A ξ t :=
          integral_const_mul _ _
        _ = M * (T - err) := by rw [hA_eq]
        _ = (T + (-err / N) * N) * M := by
          rw [div_mul_cancel₀ (-err) hN.ne']
          ring

/-- Absolute-error form of the $L^2$ integral estimate. -/
theorem l2_integral_estimate_error :
    ∃ C : ℝ, 0 < C ∧
    ∀ {ι : Type*} [Fintype ι] (a : ι → ℂ) (ξ : ι → ℝ) (N : ℝ),
    0 < N →
    IsSeparatedFamily (1 / N) ξ →
    ∀ left right T : ℝ,
    T = right - left →
    left ≤ right →
    |(∫ t in Set.Icc left right, ‖∑ r, a r * 𝐞 (ξ r * t)‖ ^ 2) -
        T * ∑ r, ‖a r‖ ^ 2| ≤
      C * N * ∑ r, ‖a r‖ ^ 2 := by
  obtain ⟨C, hC, hestimate⟩ := l2_integral_estimate
  refine ⟨C, hC, ?_⟩
  intro ι _ a ξ N hN hsep left right T hT hleft_right
  obtain ⟨θ, hθ, hidentity⟩ :=
    hestimate a ξ N hN hsep left right T hT hleft_right
  let S := ∑ r, ‖a r‖ ^ 2
  have hS : 0 ≤ S := Finset.sum_nonneg fun r _ => sq_nonneg ‖a r‖
  rw [hidentity]
  calc
    |(T + θ * N) * S - T * S| = |θ| * N * S := by
      rw [show (T + θ * N) * S - T * S = θ * N * S by ring,
        abs_mul, abs_mul, abs_of_pos hN, abs_of_nonneg hS]
    _ ≤ C * N * S := by gcongr

end Expdb

end
