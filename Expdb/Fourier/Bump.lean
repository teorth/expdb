module

public import Expdb.Basic.Definitions

public import Mathlib.Analysis.Calculus.BumpFunction.InnerProduct
public import Mathlib.Analysis.Distribution.SchwartzSpace.Fourier

/-!
# A normalized smooth bump and its Fourier transform

This module collects the concrete bump function used in the proof of the `L²` integral
estimate of the blueprint's Basic Fourier estimates chapter (`l2-chapter`), together with its
Fourier transform and the orthonormal family formed by its translates.

The bump `bump` is a fixed smooth, nonnegative, compactly supported function with
`‖bump‖_{L²} = 1` and support in `[-1/4, 1/4]`. Its translates by `1`-separated shifts are
therefore orthonormal, which is what the `L²` estimate uses.

These declarations are auxiliary to `Expdb.Fourier.L2Integral` and are collected in the
namespace `Expdb.L2Bump`.
-/

@[expose] public section

open MeasureTheory Real Complex Filter Topology BigOperators
open scoped FourierTransform SchwartzMap ContDiff

noncomputable section

namespace Expdb.L2Bump

/-! ## A normalized smooth bump -/

def rawBump : ContDiffBump (0 : ℝ) :=
  ⟨1 / 8, 1 / 4, by norm_num, by norm_num⟩

def rawL2 : ℝ := ∫ x : ℝ, (rawBump x) ^ 2

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

def bump (x : ℝ) : ℝ := rawBump x / Real.sqrt rawL2

lemma bump_smooth : ContDiff ℝ ∞ bump := by
  change ContDiff ℝ ∞ (fun x => rawBump x / Real.sqrt rawL2)
  exact (rawBump.contDiff (n := ⊤)).div_const (Real.sqrt rawL2)

lemma bump_hasCompactSupport : HasCompactSupport bump := by
  apply HasCompactSupport.of_support_subset_isCompact rawBump.hasCompactSupport.isCompact
  intro x hx
  apply subset_tsupport rawBump
  simp only [Function.mem_support] at hx ⊢
  intro hzero
  apply hx
  simp [bump, hzero]

lemma bump_supp (x : ℝ) (hx : bump x ≠ 0) : |x| ≤ 1 / 4 := by
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

lemma bump_nonneg (x : ℝ) : 0 ≤ bump x :=
  div_nonneg (rawBump.nonneg' x) (Real.sqrt_nonneg rawL2)

lemma bump_l2norm : ∫ x : ℝ, (bump x) ^ 2 = 1 := by
  rw [show (fun x : ℝ => (bump x) ^ 2) = fun x => (rawBump x) ^ 2 / rawL2 by
        funext x
        simp only [bump, div_pow]
        rw [Real.sq_sqrt rawL2_pos.le]]
  rw [integral_div]
  exact div_self (ne_of_gt rawL2_pos)

-- bump̂(u) = ∫_ℝ bump(x) e(-xu) dx
def bumpFourier (u : ℝ) : ℂ :=
  ∫ x : ℝ, (bump x : ℂ) * 𝐞 (-(x * u))

-- bump is integrable
lemma bump_integrable : Integrable bump :=
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

/-! ## Plancherel and decay of the bump's Fourier transform -/

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

lemma bumpFourier_continuous : Continuous bumpFourier := by
  rw [bumpFourier_eq_fourier]
  exact (𝓕 bumpSchwartz).continuous

lemma bumpFourier_l2 : ∫ u : ℝ, ‖bumpFourier u‖ ^ 2 = 1 := by
  simp_rw [bumpFourier_eq_fourier]
  rw [← SchwartzMap.fourier_coe, SchwartzMap.integral_norm_sq_fourier]
  simpa [bumpSchwartz, Real.norm_eq_abs, abs_of_nonneg (bump_nonneg _)] using bump_l2norm

lemma bumpFourier_sq_integrable :
    Integrable (fun u : ℝ => ‖bumpFourier u‖ ^ 2) := by
  by_contra h
  have hzero := bumpFourier_l2
  rw [integral_undef h] at hzero
  norm_num at hzero

lemma bumpFourier_decay (K : ℕ) :
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

lemma bumpFourier_lower_bound :
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

/-! ## Weighted Plancherel identity -/

-- Package translates of the concrete bump as Schwartz functions.
def bumpShift (w : ℝ) : 𝓢(ℝ, ℂ) := by
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
    exact
      ((Complex.ofRealCLM.contDiff (n := ∞)).comp bump_smooth).comp
        (contDiff_id.add contDiff_const)
  exact hcomp.toSchwartzMap hsmooth

lemma bumpShift_apply (w x : ℝ) : bumpShift w x = bump (x + w) := rfl

lemma fourier_bumpShift (w u : ℝ) :
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

lemma bumpShift_orthonormal {ι : Type*} [Finite ι]
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

end Expdb.L2Bump

end
