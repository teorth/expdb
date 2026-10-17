module

public import Expdb.Basic.DerivativeAsymptotics
public import Expdb.ExponentialSums.PhaseExtension
public import Expdb.ExponentialSums.PhaseFunctions
public import Expdb.Mathlib.StationaryPhase
public import Expdb.Mathlib.IteratedDeriv
public import Mathlib.Analysis.Calculus.BumpFunction.InnerProduct
public import Mathlib.Analysis.Calculus.ContDiff.Basic
import Mathlib.Analysis.Calculus.ContDiff.Operations
public import Mathlib.Analysis.Calculus.Deriv.Basic
import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.Calculus.Deriv.Support
public import Mathlib.Analysis.Calculus.IteratedDeriv.Lemmas
import Mathlib.Analysis.Calculus.MeanValue
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv
import Mathlib.Data.Nat.Find
public import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
import Mathlib.Topology.Algebra.Order.Field

/-!
# Regularization and the dual phase for the B-process

This module regularizes model phases and constructs the normalized dual model phases
used by the van der Corput B-process.
-/

@[expose] public section

namespace Expdb.BProcess

open StationaryPhase

open Filter Topology Set MeasureTheory
open scoped ContDiff

noncomputable section

section Primitives

/-! ## Primitives -/

private theorem contDiff_integral_from {h : ℝ → ℝ} (hh : ContDiff ℝ ∞ h) (x₀ : ℝ) :
    ContDiff ℝ ∞ (fun x ↦ ∫ t in x₀..x, h t) := by
  refine contDiff_infty_iff_deriv.2 ⟨fun x ↦
    (hh.continuous.integral_hasStrictDerivAt x₀ x).hasDerivAt.differentiableAt, ?_⟩
  rw [show deriv (fun x ↦ ∫ t in x₀..x, h t) = h from
    funext fun x ↦ Continuous.deriv_integral h hh.continuous x₀ x]
  exact hh

/-- The function with second derivative `h`, value `y₀` and slope `y₁` at `x₀`. -/
private def secondPrimitive (h : ℝ → ℝ) (x₀ y₀ y₁ : ℝ) (x : ℝ) : ℝ :=
  y₀ + ∫ t in x₀..x, (y₁ + ∫ s in x₀..t, h s)

private theorem deriv_secondPrimitive {h : ℝ → ℝ} (hh : Continuous h) (x₀ y₀ y₁ : ℝ) :
    deriv (secondPrimitive h x₀ y₀ y₁) = fun x ↦ y₁ + ∫ s in x₀..x, h s := by
  have hc : Continuous fun t ↦ y₁ + ∫ s in x₀..t, h s :=
    continuous_const.add (continuous_iff_continuousAt.2 fun t ↦
      (hh.integral_hasStrictDerivAt x₀ t).hasDerivAt.continuousAt)
  funext x
  have : HasDerivAt (secondPrimitive h x₀ y₀ y₁) (y₁ + ∫ s in x₀..x, h s) x :=
    (hc.integral_hasStrictDerivAt x₀ x).hasDerivAt.const_add y₀
  exact this.deriv

private theorem contDiff_secondPrimitive {h : ℝ → ℝ} (hh : ContDiff ℝ ∞ h) (x₀ y₀ y₁ : ℝ) :
    ContDiff ℝ ∞ (secondPrimitive h x₀ y₀ y₁) := by
  unfold secondPrimitive
  exact contDiff_const.add (contDiff_integral_from
    (contDiff_const.add (contDiff_integral_from hh x₀)) x₀)

private theorem iteratedDeriv_two_secondPrimitive {h : ℝ → ℝ} (hh : Continuous h) (x₀ y₀ y₁ : ℝ) :
    iteratedDeriv 2 (secondPrimitive h x₀ y₀ y₁) = h := by
  rw [iteratedDeriv_succ, iteratedDeriv_one, deriv_secondPrimitive hh]
  funext x
  have : HasDerivAt (fun x ↦ y₁ + ∫ s in x₀..x, h s) (h x) x :=
    (hh.integral_hasStrictDerivAt x₀ x).hasDerivAt.const_add y₁
  exact this.deriv

private theorem iteratedDeriv_add_two_secondPrimitive {h : ℝ → ℝ} (hh : Continuous h)
    (x₀ y₀ y₁ : ℝ) (k : ℕ) :
    iteratedDeriv (k + 2) (secondPrimitive h x₀ y₀ y₁) = iteratedDeriv k h := by
  rw [← iteratedDeriv_iteratedDeriv, iteratedDeriv_two_secondPrimitive hh]

/-! ## Plateau cutoffs -/

/-- A smooth cutoff equal to one on `[a, b]` and supported in `[a - η, b + η]`. -/
private def plateauBump {a b η : ℝ} (hab : a < b) (hη : 0 < η) : ContDiffBump ((a + b) / 2) :=
  ⟨(b - a) / 2, (b - a) / 2 + η, by linarith, by linarith⟩

private theorem plateauBump_eq_one {a b η : ℝ} (hab : a < b) (hη : 0 < η) {x : ℝ}
    (hx : x ∈ Icc a b) : plateauBump hab hη x = 1 := by
  apply ContDiffBump.one_of_mem_closedBall
  rw [Real.closedBall_eq_Icc]
  change x ∈ Icc ((a + b) / 2 - (b - a) / 2) ((a + b) / 2 + (b - a) / 2)
  constructor <;> linarith [hx.1, hx.2]

private theorem tsupport_plateauBump {a b η : ℝ} (hab : a < b) (hη : 0 < η) :
    tsupport (plateauBump hab hη : ℝ → ℝ) = Icc (a - η) (b + η) := by
  rw [ContDiffBump.tsupport_eq, Real.closedBall_eq_Icc]
  change Icc ((a + b) / 2 - ((b - a) / 2 + η)) ((a + b) / 2 + ((b - a) / 2 + η)) = _
  congr 1 <;> ring

/-! ## Products with cutoffs -/

/-- All derivatives of a cutoff times any function vanish outside the cutoff's support. -/
private theorem iteratedDeriv_mul_eq_zero_of_notMem_tsupport {χ g : ℝ → ℝ} {x : ℝ}
    (hx : x ∉ tsupport χ) (k : ℕ) : iteratedDeriv k (fun y ↦ χ y * g y) x = 0 := by
  have hzero : (fun y ↦ χ y * g y) =ᶠ[𝓝 x] fun _ ↦ 0 := by
    filter_upwards [(notMem_tsupport_iff_eventuallyEq.1 hx)] with y hy
    simp [show χ y = 0 from hy]
  rw [hzero.iteratedDeriv_eq]
  simp

/-! ## Integrals of functions supported in an interval -/

/-- An integral of a continuous function supported in `[a, b]` is bounded by the interval
length times a bound on the function. -/
private theorem norm_integral_le_of_vanishing_outside {f : ℝ → ℝ} (hf : Continuous f)
    {a b x₀ ε : ℝ} (hab : a ≤ b)
    (hzero : ∀ t, t ∉ Icc a b → f t = 0) (hbound : ∀ t ∈ Icc a b, ‖f t‖ ≤ ε) (x : ℝ) :
    ‖∫ t in x₀..x, f t‖ ≤ ε * (b - a) := by
  have hsupp : Function.support f ⊆ Icc a b := by
    intro y hy
    by_contra h
    exact hy (hzero y h)
  have hint : Integrable f := hf.integrable_of_hasCompactSupport
    (HasCompactSupport.of_support_subset_isCompact isCompact_Icc hsupp)
  calc
    ‖∫ t in x₀..x, f t‖ ≤ ∫ t in uIoc x₀ x, ‖f t‖ :=
      intervalIntegral.norm_integral_le_integral_norm_uIoc
    _ ≤ ∫ t, ‖f t‖ := setIntegral_le_integral hint.norm
      (Eventually.of_forall fun _ ↦ norm_nonneg _)
    _ ≤ ε * (b - a) := integral_norm_le_of_support_Icc hab hzero hbound

end Primitives

section Regularization

/-- Global regularizations `H i` of the phases `G i`, with limit `H₀`. -/
structure RegularizedPhases (σ a b : ℝ) (G H : ℕ → ℝ → ℝ) (H₀ : ℝ → ℝ) : Prop where
  smooth : ∀ i, ContDiff ℝ ∞ (H i)
  smooth_lim : ContDiff ℝ ∞ H₀
  eqOn : ∀ᶠ i in atTop, EqOn (H i) (G i) (Icc a b)
  curv : ∃ c C : ℝ, 0 < c ∧
    (∀ i x, -C ≤ iteratedDeriv 2 (H i) x ∧ iteratedDeriv 2 (H i) x ≤ -c) ∧
    ∀ x, -C ≤ iteratedDeriv 2 H₀ x ∧ iteratedDeriv 2 H₀ x ≤ -c
  conv : ∀ k, 1 ≤ k → ∀ ε > 0, ∀ᶠ i in atTop, ∀ x,
    ‖iteratedDeriv k (H i) x - iteratedDeriv k H₀ x‖ ≤ ε
  bdd : ∀ k, 2 ≤ k → ∃ K, ∀ x, ‖iteratedDeriv k H₀ x‖ ≤ K
  model : EqOn (deriv H₀) (modelPhase σ) (Icc a b)

/-- The value of `secondPrimitive` agrees with any `C²` function having the same second
derivative on an interval and the same initial data. -/
private theorem secondPrimitive_eqOn {h F : ℝ → ℝ} (hF : ContDiff ℝ ∞ F)
    {a b x₀ : ℝ} (hx₀ : x₀ ∈ Icc a b) (heq : EqOn h (iteratedDeriv 2 F) (Icc a b)) :
    EqOn (secondPrimitive h x₀ (F x₀) (deriv F x₀)) F (Icc a b) := by
  have hF1 : ContDiff ℝ ∞ (deriv F) := hF.deriv' (n := ∞)
  have hF2 : ContDiff ℝ ∞ (deriv (deriv F)) := hF1.deriv' (n := ∞)
  have hsub : ∀ x ∈ Icc a b, uIcc x₀ x ⊆ Icc a b := fun x hx ↦
    uIcc_subset_Icc hx₀ hx
  have hderiv : ∀ x ∈ Icc a b, deriv F x₀ + ∫ s in x₀..x, h s = deriv F x := by
    intro x hx
    have hcongr : ∫ s in x₀..x, h s = ∫ s in x₀..x, deriv (deriv F) s := by
      apply intervalIntegral.integral_congr
      intro s hs
      rw [heq (hsub x hx hs), iteratedDeriv_succ, iteratedDeriv_one]
    rw [hcongr, intervalIntegral.integral_deriv_eq_sub
      (fun s _ ↦ (hF1.differentiable (by simp)) s) (hF2.continuous.intervalIntegrable _ _)]
    ring
  intro x hx
  unfold secondPrimitive
  have hcongr : ∫ t in x₀..x, (deriv F x₀ + ∫ s in x₀..t, h s) =
      ∫ t in x₀..x, deriv F t := by
    apply intervalIntegral.integral_congr
    intro t ht
    exact hderiv t (hsub x hx ht)
  rw [hcongr, intervalIntegral.integral_deriv_eq_sub
    (fun s _ ↦ (hF.differentiable (by simp)) s) (hF1.continuous.intervalIntegrable _ _)]
  ring

private def regSecond (χ : ℝ → ℝ) (c : ℝ) (g : ℝ → ℝ) (x : ℝ) : ℝ :=
  χ x * g x - c * (1 - χ x)

private theorem regSecond_bounds {χ g : ℝ → ℝ} {c C x : ℝ} (hcC : c ≤ C)
    (h0 : 0 ≤ χ x) (h1 : χ x ≤ 1) (hg : χ x ≠ 0 → -C ≤ g x ∧ g x ≤ -c) :
    -C ≤ regSecond χ c g x ∧ regSecond χ c g x ≤ -c := by
  unfold regSecond
  by_cases hχ : χ x = 0
  · rw [hχ]
    constructor <;> linarith
  · obtain ⟨hl, hu⟩ := hg hχ
    constructor <;> nlinarith

/-- Phases whose derivatives converge to those of the reference phase on a neighbourhood
of `[a, b]` admit global concave regularizations. -/
theorem exists_regularizedPhases {σ : ℝ} (hσ : 0 < σ) {a b : ℝ} (ha : 0 < a) (hab : a < b)
    {G : ℕ → ℝ → ℝ} (hG : ∀ i, ContDiff ℝ ∞ (G i))
    (hconv : ∀ p : ℕ, ∀ ε > 0, ∀ᶠ i in atTop, ∀ x ∈ Icc (a / 2) (b + 1),
      ‖iteratedDeriv (p + 1) (G i) x - iteratedDeriv p (modelPhase σ) x‖ ≤ ε) :
    ∃ H H₀, RegularizedPhases σ a b G H H₀ := by
  classical
  -- The cutoff.
  set η : ℝ := min (a / 4) (1 / 2) with hη_def
  have hη : 0 < η := lt_min (by linarith) (by norm_num)
  let χb := plateauBump hab hη
  let χ : ℝ → ℝ := χb
  have hχs : ContDiff ℝ ∞ χ := χb.contDiff (n := ⊤)
  have hχ0 : ∀ x, 0 ≤ χ x := fun x ↦ χb.nonneg
  have hχ1 : ∀ x, χ x ≤ 1 := fun x ↦ χb.le_one
  have hχone : ∀ x ∈ Icc a b, χ x = 1 := fun x hx ↦ plateauBump_eq_one hab hη hx
  have htsupp : tsupport χ ⊆ Icc (a / 2) (b + 1) := by
    rw [tsupport_plateauBump hab hη]
    apply Icc_subset_Icc
    · have : η ≤ a / 4 := min_le_left _ _
      linarith
    · have : η ≤ 1 / 2 := min_le_right _ _
      linarith
  have hχzero : ∀ x, x ∉ Icc (a / 2) (b + 1) → χ x = 0 := by
    intro x hx
    by_contra h
    exact hx (htsupp (subset_tsupport _ h))
  have hnotsupp : ∀ x, x ∉ Icc (a / 2) (b + 1) → x ∉ tsupport χ :=
    fun x hx h ↦ hx (htsupp h)
  have hχcs : HasCompactSupport χ := χb.hasCompactSupport
  -- The reference curvature.
  let φ : ℝ → ℝ := modelPhase σ
  let g₀ : ℝ → ℝ := deriv φ
  have hφsmooth : ContDiffOn ℝ ∞ φ (Ioi 0) := fun x hx ↦
    (Real.contDiffAt_rpow_const_of_ne (ne_of_gt hx)).contDiffWithinAt
  have hg₀smooth : ContDiffOn ℝ ∞ g₀ (Ioi 0) :=
    hφsmooth.deriv_of_isOpen isOpen_Ioi (by simp)
  have hg₀eq : ∀ x, g₀ x = -σ * x ^ (-σ - 1) := by
    intro x
    exact Real.deriv_rpow_const x (-σ)
  set c₁ : ℝ := σ * (b + 1) ^ (-σ - 1) / 2 with hc₁_def
  have hb1 : 0 < b + 1 := by linarith
  have hc₁ : 0 < c₁ := by positivity
  set C₁ : ℝ := σ * (a / 2) ^ (-σ - 1) + c₁ with hC₁_def
  have hg₀bounds : ∀ x ∈ Icc (a / 2) (b + 1),
      -σ * (a / 2) ^ (-σ - 1) ≤ g₀ x ∧ g₀ x ≤ -2 * c₁ := by
    intro x hx
    have hx0 : 0 < x := lt_of_lt_of_le (by linarith) hx.1
    have h1 : x ^ (-σ - 1) ≤ (a / 2) ^ (-σ - 1) :=
      Real.rpow_le_rpow_of_nonpos (by linarith) hx.1 (by linarith)
    have h2 : (b + 1) ^ (-σ - 1) ≤ x ^ (-σ - 1) :=
      Real.rpow_le_rpow_of_nonpos hx0 hx.2 (by linarith)
    rw [hg₀eq]
    constructor
    · nlinarith
    · rw [hc₁_def]
      nlinarith
  have hcC : c₁ ≤ C₁ := by
    rw [hC₁_def]
    have : 0 ≤ σ * (a / 2) ^ (-σ - 1) := by positivity
    linarith
  have hr₀bounds : ∀ x, -C₁ ≤ regSecond χ c₁ g₀ x ∧ regSecond χ c₁ g₀ x ≤ -c₁ := by
    intro x
    apply regSecond_bounds hcC (hχ0 x) (hχ1 x)
    intro hx
    have hxI : x ∈ Icc (a / 2) (b + 1) := by
      by_contra h
      exact hx (hχzero x h)
    have h2 := hg₀bounds x hxI
    rw [hC₁_def]
    constructor <;> linarith [h2.1, h2.2]
  -- Good indices.
  let good : ℕ → Prop := fun i ↦ ∀ x ∈ Icc (a / 2) (b + 1),
    ‖iteratedDeriv 2 (G i) x - g₀ x‖ ≤ c₁
  have hgood : ∀ᶠ i in atTop, good i := by
    filter_upwards [hconv 1 c₁ hc₁] with i hi x hx
    simpa [g₀, φ, iteratedDeriv_one] using hi x hx
  -- The regularized second derivatives.
  let r : ℕ → ℝ → ℝ := fun i ↦ regSecond χ c₁ (iteratedDeriv 2 (G i))
  let r₀ : ℝ → ℝ := regSecond χ c₁ g₀
  have hG2 : ∀ i, ContDiff ℝ ∞ (iteratedDeriv 2 (G i)) := fun i ↦ by
    rw [iteratedDeriv_succ, iteratedDeriv_one]
    exact (((hG i).deriv' (n := ∞))).deriv' (n := ∞)
  have hχg₀ : ContDiff ℝ ∞ (fun x ↦ χ x * g₀ x) :=
    SmoothExtension.contDiff_mul_of_tsupport isOpen_Ioi hχs hg₀smooth (fun x hx ↦
      show (0 : ℝ) < x from lt_of_lt_of_le (by linarith) (htsupp hx).1)
  have hrs : ∀ i, ContDiff ℝ ∞ (r i) := fun i ↦
    (hχs.mul (hG2 i)).sub (contDiff_const.mul (contDiff_const.sub hχs))
  have hr₀s : ContDiff ℝ ∞ r₀ := hχg₀.sub (contDiff_const.mul (contDiff_const.sub hχs))
  -- The regularized phases.
  let H₀ : ℝ → ℝ := secondPrimitive r₀ a 0 (φ a)
  let H : ℕ → ℝ → ℝ := fun i ↦
    if good i then secondPrimitive (r i) a (G i a) (deriv (G i) a) else H₀
  have hH₀s : ContDiff ℝ ∞ H₀ := contDiff_secondPrimitive hr₀s _ _ _
  have hHs : ∀ i, ContDiff ℝ ∞ (H i) := by
    intro i
    by_cases hi : good i
    · simp only [H, if_pos hi]
      exact contDiff_secondPrimitive (hrs i) _ _ _
    · simp only [H, if_neg hi]
      exact hH₀s
  have hH₀2 : iteratedDeriv 2 H₀ = r₀ := iteratedDeriv_two_secondPrimitive hr₀s.continuous _ _ _
  have hHgood : ∀ i, good i →
      H i = secondPrimitive (r i) a (G i a) (deriv (G i) a) := fun i hi ↦ by
    simp only [H, if_pos hi]
  have hHbad : ∀ i, ¬ good i → H i = H₀ := fun i hi ↦ by simp only [H, if_neg hi]
  have ha_mem : a ∈ Icc (a / 2) (b + 1) := ⟨by linarith, by linarith⟩
  refine ⟨H, H₀, ⟨hHs, hH₀s, ?_, ?_, ?_, ?_, ?_⟩⟩
  · -- Agreement on `[a, b]`.
    filter_upwards [hgood] with i hi
    rw [hHgood i hi]
    apply secondPrimitive_eqOn (hG i) ⟨le_rfl, hab.le⟩
    intro x hx
    simp only [r, regSecond, hχone x hx]
    ring
  · -- Curvature bounds.
    refine ⟨c₁, C₁, hc₁, fun i x ↦ ?_, fun x ↦ ?_⟩
    · by_cases hi : good i
      · rw [hHgood i hi, iteratedDeriv_two_secondPrimitive (hrs i).continuous]
        apply regSecond_bounds hcC (hχ0 x) (hχ1 x)
        intro hx
        have hxI : x ∈ Icc (a / 2) (b + 1) := by
          by_contra h
          exact hx (hχzero x h)
        have h1 := hi x hxI
        have h2 := hg₀bounds x hxI
        rw [Real.norm_eq_abs, abs_le] at h1
        rw [hC₁_def]
        constructor <;> linarith [h1.1, h1.2, h2.1, h2.2]
      · rw [hHbad i hi, hH₀2]
        exact hr₀bounds x
    · rw [hH₀2]
      exact hr₀bounds x
  · -- Convergence of the derivatives.
    intro k hk ε hε
    rcases Nat.lt_or_ge k 2 with hk2 | hk2
    · -- First derivatives.
      obtain rfl : k = 1 := by omega
      set δ : ℝ := ε / (b + 3) with hδ_def
      have hδ : 0 < δ := div_pos hε (by linarith)
      filter_upwards [hgood, hconv 0 δ hδ, hconv 1 δ hδ] with i hi h0 h1
      intro x
      rw [hHgood i hi, iteratedDeriv_one, iteratedDeriv_one,
        deriv_secondPrimitive (hrs i).continuous, deriv_secondPrimitive hr₀s.continuous]
      have hint : (∫ s in a..x, r i s) - ∫ s in a..x, r₀ s = ∫ s in a..x, (r i s - r₀ s) :=
        (intervalIntegral.integral_sub ((hrs i).continuous.intervalIntegrable _ _)
          (hr₀s.continuous.intervalIntegrable _ _)).symm
      have hdiff : ∀ s, r i s - r₀ s = χ s * (iteratedDeriv 2 (G i) s - g₀ s) := by
        intro s
        simp only [r, r₀, regSecond]
        ring
      have hbound := norm_integral_le_of_vanishing_outside (x₀ := a)
        (f := fun s ↦ r i s - r₀ s)
        ((hrs i).continuous.sub hr₀s.continuous) (by linarith : a / 2 ≤ b + 1)
        (fun s hs ↦ by rw [hdiff, hχzero s hs, zero_mul])
        (fun s hs ↦ by
          rw [hdiff, norm_mul, Real.norm_of_nonneg (hχ0 s)]
          have := h1 s hs
          simp only [iteratedDeriv_one] at this
          calc χ s * ‖iteratedDeriv 2 (G i) s - g₀ s‖ ≤ 1 * δ :=
                mul_le_mul (hχ1 s) this (norm_nonneg _) zero_le_one
            _ = δ := one_mul δ) x
      have hfirst := h0 a ha_mem
      simp only [zero_add, iteratedDeriv_one, iteratedDeriv_zero] at hfirst
      change ‖deriv (G i) a + (∫ s in a..x, r i s) - (φ a + ∫ s in a..x, r₀ s)‖ ≤ ε
      rw [show deriv (G i) a + (∫ s in a..x, r i s) - (φ a + ∫ s in a..x, r₀ s) =
        (deriv (G i) a - φ a) + ((∫ s in a..x, r i s) - ∫ s in a..x, r₀ s) by ring, hint]
      have hδε : δ * (b + 3) = ε := by
        rw [hδ_def, div_mul_cancel₀ _ (by linarith : b + 3 ≠ 0)]
      calc ‖(deriv (G i) a - φ a) + ∫ s in a..x, (r i s - r₀ s)‖
          ≤ ‖deriv (G i) a - φ a‖ + ‖∫ s in a..x, (r i s - r₀ s)‖ := norm_add_le _ _
        _ ≤ δ + δ * (b + 1 - a / 2) := add_le_add hfirst hbound
        _ ≤ ε := by nlinarith
    · -- Higher derivatives.
      obtain ⟨j, rfl⟩ : ∃ j, k = j + 2 := ⟨k - 2, by omega⟩
      have hKχ : ∀ m, ∃ K, ∀ x, ‖iteratedDeriv m χ x‖ ≤ K := fun m ↦
        exists_bound_iteratedDeriv_of_hasCompactSupport hχs hχcs m
      choose Kχ hKχ using hKχ
      set L : ℝ := ∑ m ∈ Finset.range (j + 1), (j.choose m : ℝ) * Kχ m with hL_def
      have hKχ0 : ∀ m, 0 ≤ Kχ m := fun m ↦ (norm_nonneg _).trans (hKχ m 0)
      have hL0 : 0 ≤ L := Finset.sum_nonneg fun m _ ↦ mul_nonneg (Nat.cast_nonneg _) (hKχ0 m)
      set δ : ℝ := ε / (L + 1) with hδ_def
      have hδ : 0 < δ := div_pos hε (by linarith)
      have hall : ∀ᶠ i in atTop, ∀ m ∈ Finset.range (j + 1), ∀ x ∈ Icc (a / 2) (b + 1),
          ‖iteratedDeriv (j - m + 1 + 1) (G i) x - iteratedDeriv (j - m + 1) φ x‖ ≤ δ :=
        (Finset.range (j + 1)).eventually_all.2 fun m _ ↦ hconv (j - m + 1) δ hδ
      filter_upwards [hgood, hall] with i hi hconvi
      intro x
      rw [hHgood i hi, iteratedDeriv_add_two_secondPrimitive (hrs i).continuous,
        iteratedDeriv_add_two_secondPrimitive hr₀s.continuous]
      -- Remove the common constant part.
      have hsplit : ∀ (g : ℝ → ℝ), ContDiff ℝ ∞ (fun y ↦ χ y * g y) →
          iteratedDeriv j (regSecond χ c₁ g) x =
            iteratedDeriv j (fun y ↦ χ y * g y) x -
              iteratedDeriv j (fun y ↦ c₁ * (1 - χ y)) x := by
        intro g hg
        exact iteratedDeriv_fun_sub
          (hg.contDiffAt.of_le (by exact_mod_cast le_top))
          ((contDiff_const.mul (contDiff_const.sub hχs)).contDiffAt.of_le
            (by exact_mod_cast le_top))
      have hr_eq : r i = regSecond χ c₁ (iteratedDeriv 2 (G i)) := rfl
      rw [hr_eq, hsplit _ (hχs.mul (hG2 i)), hsplit _ hχg₀, sub_sub_sub_cancel_right]
      by_cases hx : x ∈ Icc (a / 2) (b + 1)
      · have hx0 : 0 < x := lt_of_lt_of_le (by linarith) hx.1
        have hg₀at : ContDiffAt ℝ ∞ g₀ x := hg₀smooth.contDiffAt (isOpen_Ioi.mem_nhds hx0)
        have hbound := norm_iteratedDeriv_mul_sub_le (n := j) (x := x)
          (hχs.contDiffAt.of_le (by exact_mod_cast le_top))
          ((hG2 i).contDiffAt.of_le (by exact_mod_cast le_top))
          (hg₀at.of_le (by exact_mod_cast le_top))
          (A := Kχ) (ε := δ) (fun m _ ↦ hKχ m x) (fun m hm ↦ by
            rw [iteratedDeriv_iteratedDeriv, show iteratedDeriv m g₀ =
              iteratedDeriv (m + 1) φ from (by rw [iteratedDeriv_succ'])]
            simpa only [show j - (j - m) = m by omega, Nat.add_assoc] using
              hconvi (j - m) (Finset.mem_range.2 (by omega)) x hx)
        calc _ ≤ L * δ := hbound
          _ ≤ ε := by
            rw [hδ_def, mul_div_assoc', div_le_iff₀ (by linarith)]
            nlinarith
      · rw [iteratedDeriv_mul_eq_zero_of_notMem_tsupport (hnotsupp x hx),
          iteratedDeriv_mul_eq_zero_of_notMem_tsupport (hnotsupp x hx), sub_zero, norm_zero]
        exact hε.le
  · -- Bounded derivatives of the limit.
    intro k hk
    obtain ⟨j, rfl⟩ : ∃ j, k = j + 2 := ⟨k - 2, by omega⟩
    have hcs : HasCompactSupport (fun y ↦ χ y * (g₀ y + c₁)) :=
      hχcs.mul_right
    have hsm : ContDiff ℝ ∞ (fun y ↦ χ y * (g₀ y + c₁)) := by
      have : (fun y ↦ χ y * (g₀ y + c₁)) = fun y ↦ χ y * g₀ y + c₁ * χ y := by
        funext y
        ring
      rw [this]
      exact hχg₀.add (contDiff_const.mul hχs)
    obtain ⟨K, hK⟩ := exists_bound_iteratedDeriv_of_hasCompactSupport hsm hcs j
    refine ⟨K + |c₁|, fun x ↦ ?_⟩
    rw [iteratedDeriv_add_two_secondPrimitive hr₀s.continuous]
    have hr₀eq : r₀ = fun y ↦ χ y * (g₀ y + c₁) - c₁ := by
      funext y
      simp only [r₀, regSecond]
      ring
    rw [hr₀eq]
    have := iteratedDeriv_sub (n := j) (x := x) (f := fun y ↦ χ y * (g₀ y + c₁))
      (g := fun _ ↦ c₁) (hsm.contDiffAt.of_le (by exact_mod_cast le_top))
      contDiffAt_const
    change iteratedDeriv j ((fun y ↦ χ y * (g₀ y + c₁)) - fun _ ↦ c₁) x = _ at this
    rw [show (fun y ↦ χ y * (g₀ y + c₁) - c₁) =
      (fun y ↦ χ y * (g₀ y + c₁)) - fun _ ↦ c₁ from rfl, this]
    calc ‖iteratedDeriv j (fun y ↦ χ y * (g₀ y + c₁)) x - iteratedDeriv j (fun _ ↦ c₁) x‖
        ≤ ‖iteratedDeriv j (fun y ↦ χ y * (g₀ y + c₁)) x‖ +
          ‖iteratedDeriv j (fun _ ↦ c₁) x‖ := norm_sub_le _ _
      _ ≤ K + |c₁| := by
          apply add_le_add (hK x)
          rcases j with _ | j
          · simp
          · rw [iteratedDeriv_const]
            simp
  · -- The limit has the reference derivative on `[a, b]`.
    intro x hx
    rw [deriv_secondPrimitive hr₀s.continuous]
    have hsub : uIcc a x ⊆ Icc a b := uIcc_subset_Icc ⟨le_rfl, hab.le⟩ hx
    have hcongr : ∫ s in a..x, r₀ s = ∫ s in a..x, deriv φ s := by
      apply intervalIntegral.integral_congr
      intro s hs
      simp only [r₀, regSecond, hχone s (hsub hs), g₀]
      ring
    have hpos : ∀ s ∈ uIcc a x, 0 < s := fun s hs ↦ lt_of_lt_of_le ha (hsub hs).1
    change φ a + (∫ s in a..x, r₀ s) = modelPhase σ x
    rw [hcongr, intervalIntegral.integral_deriv_eq_sub
      (fun s hs ↦ ((hφsmooth.contDiffAt (isOpen_Ioi.mem_nhds (hpos s hs))).differentiableAt
        (by simp)))
      (ContinuousOn.intervalIntegrable (hg₀smooth.continuousOn.mono fun s hs ↦ hpos s hs))]
    simp [φ]

end Regularization

section Extension

/-- A model phase admits global extensions converging to the reference phase in every
derivative of positive order, uniformly on `[l, u]`. -/
theorem exists_global_extensions {F : VariableFunction (VariableObject.fixed ℝ) ℝ} {σ : ℝ}
    (hF : IsModelPhaseFunctionWith F σ) {l u : ℝ} (hl : 0 < l) (hl1 : l < 1) (hu : 2 < u) :
    ∃ G : ℕ → ℝ → ℝ, (∀ i, ContDiff ℝ ∞ (G i)) ∧
      (∀ᶠ i in atTop, EqOn (G i) (F i) phaseInterval) ∧
      ∀ p : ℕ, ∀ ε > 0, ∀ᶠ i in atTop, ∀ x ∈ Icc l u,
        ‖iteratedDeriv (p + 1) (G i) x - iteratedDeriv p (modelPhase σ) x‖ ≤ ε := by
  classical
  choose C hC hext using fun P : ℕ ↦ exists_global_approximateModelPhase_extension P
  let δ : ℕ → ℝ := fun P ↦ 1 / ((P + 1) * C P)
  have hCpos : ∀ P, 0 < C P := fun P ↦ lt_of_lt_of_le one_pos (hC P)
  have hδ : ∀ P, 0 < δ P := fun P ↦ by
    have := hCpos P
    positivity
  have hCδ : ∀ P, C P * δ P = 1 / (P + 1) := fun P ↦ by
    have := (hCpos P).ne'
    simp only [δ]
    field_simp
  let good : ℕ → ℕ → Prop := fun P i ↦
    ∀ Q ≤ P, IsApproximateModelPhaseFunction (F i) σ Q (δ Q)
  have hgood : ∀ P, ∀ᶠ i in atTop, good P i := by
    intro P
    have := (Finset.range (P + 1)).eventually_all.2 fun Q _ ↦
      hF.eventually_isApproximate Q (hδ Q)
    filter_upwards [this] with i hi Q hQ
    exact hi Q (Finset.mem_range.2 (by omega))
  let Pidx : ℕ → ℕ := fun i ↦ Nat.findGreatest (fun P ↦ good P i) i
  let G : ℕ → ℝ → ℝ := fun i ↦
    if h : good (Pidx i) i then
      Classical.choose (hext (Pidx i) (δ (Pidx i)) (hδ _) σ (F i) (h (Pidx i) le_rfl) l u
        hl hl1 hu)
    else fun _ ↦ 0
  have hGspec : ∀ i (h : good (Pidx i) i),
      ContDiff ℝ ∞ (G i) ∧ EqOn (G i) (F i) phaseInterval ∧
        ∀ p ≤ Pidx i, ∀ x ∈ Icc l u,
          ‖iteratedDeriv (p + 1) (G i) x - iteratedDeriv p (modelPhase σ) x‖ ≤
            C (Pidx i) * δ (Pidx i) := by
    intro i h
    simp only [G, dif_pos h]
    exact Classical.choose_spec (hext (Pidx i) (δ (Pidx i)) (hδ _) σ (F i)
      (h (Pidx i) le_rfl) l u hl hl1 hu)
  have hgoodP : ∀ᶠ i in atTop, good (Pidx i) i := by
    filter_upwards [hgood 0] with i hi
    exact Nat.findGreatest_spec (P := fun P ↦ good P i) (Nat.zero_le i) hi
  have hPlarge : ∀ P, ∀ᶠ i in atTop, P ≤ Pidx i := by
    intro P
    filter_upwards [hgood P, eventually_ge_atTop P] with i hi hPi
    exact Nat.le_findGreatest (P := fun P ↦ good P i) hPi hi
  refine ⟨G, fun i ↦ ?_, ?_, ?_⟩
  · by_cases h : good (Pidx i) i
    · exact (hGspec i h).1
    · simp only [G, dif_neg h]
      exact contDiff_const
  · filter_upwards [hgoodP] with i hi
    exact (hGspec i hi).2.1
  · intro p ε hε
    obtain ⟨M, hM⟩ := exists_nat_gt (1 / ε)
    filter_upwards [hgoodP, hPlarge (max p M)] with i hi hPi x hx
    refine ((hGspec i hi).2.2 p (le_of_max_le_left hPi) x hx).trans ?_
    rw [hCδ]
    have hM' : (M : ℝ) ≤ Pidx i := by exact_mod_cast le_of_max_le_right hPi
    rw [div_le_iff₀ (by positivity)]
    rw [div_lt_iff₀ hε] at hM
    nlinarith

end Extension

section DualPhase

/-! ## Normalized dual phases -/

/-- The normalized dual phase on the block `[κ, 2κ]`. -/
def dualPhase (H : ℝ → ℝ) (κ S : ℝ) (u : ℝ) : ℝ := -legendre H (κ * u) / S

private theorem contDiff_dualPhase {H : ℝ → ℝ} (hH : ContDiff ℝ ∞ H) {c : ℝ} (hc : 0 < c)
    (hcurv : ∀ x, iteratedDeriv 2 H x ≤ -c) (κ S : ℝ) : ContDiff ℝ ∞ (dualPhase H κ S) :=
  ((contDiff_legendre hH hc hcurv).comp (contDiff_const.mul contDiff_id)).neg.div_const S

private theorem deriv_dualPhase {H : ℝ → ℝ} (hH : ContDiff ℝ ∞ H) {c : ℝ} (hc : 0 < c)
    (hcurv : ∀ x, iteratedDeriv 2 H x ≤ -c) (κ S : ℝ) :
    deriv (dualPhase H κ S) = fun u ↦ κ / S * phaseInverse H (κ * u) := by
  funext u
  have h1 : HasDerivAt (fun u ↦ κ * u) κ u := by
    simpa using (hasDerivAt_id u).const_mul κ
  have h : HasDerivAt (dualPhase H κ S) (-(-phaseInverse H (κ * u) * κ) / S) u :=
    ((hasDerivAt_legendre hH hc hcurv (κ * u)).comp u h1).neg.div_const S
  rw [h.deriv]
  ring

/-! ## Regularized phases -/

namespace RegularizedPhases

variable {σ a b : ℝ} {G H : ℕ → ℝ → ℝ} {H₀ : ℝ → ℝ}

/-- The limit stationary point is an explicit power on the reference range. -/
private theorem phaseInverse_lim_eq (hR : RegularizedPhases σ a b G H H₀) (hσ : 0 < σ)
    (ha : 0 < a) (hab : a < b) {ν : ℝ} (hν : ν ∈ Icc (b ^ (-σ)) (a ^ (-σ))) :
    phaseInverse H₀ ν = ν ^ (-1 / σ) := by
  obtain ⟨c, C, hc, -, hcurv₀⟩ := hR.curv
  obtain ⟨-, hinv, -⟩ := phaseInverse_spec hR.smooth_lim hc (fun x ↦ (hcurv₀ x).2)
  have hb : 0 < b := ha.trans hab
  have hν0 : 0 < ν := lt_of_lt_of_le (Real.rpow_pos_of_pos hb _) hν.1
  set x := ν ^ (-1 / σ) with hx_def
  have hxa : a ≤ x := by
    have := Real.rpow_le_rpow_of_nonpos hν0 hν.2 (show -1 / σ ≤ 0 by
      have : 0 < 1 / σ := by positivity
      rw [neg_div]; linarith)
    rwa [← Real.rpow_mul ha.le, show -σ * (-1 / σ) = 1 by field_simp, Real.rpow_one] at this
  have hxb : x ≤ b := by
    have := Real.rpow_le_rpow_of_nonpos (Real.rpow_pos_of_pos hb _) hν.1 (show -1 / σ ≤ 0 by
      have : 0 < 1 / σ := by positivity
      rw [neg_div]; linarith)
    rwa [← Real.rpow_mul hb.le, show -σ * (-1 / σ) = 1 by field_simp, Real.rpow_one] at this
  have hderiv : deriv H₀ x = ν := by
    rw [hR.model ⟨hxa, hxb⟩, modelPhase, hx_def, ← Real.rpow_mul hν0.le,
      show -1 / σ * -σ = 1 by field_simp, Real.rpow_one]
  rw [← hderiv, hinv]

/-- The stationary points converge, uniformly on the line. -/
private theorem tendsto_phaseInverse_sub (hR : RegularizedPhases σ a b G H H₀) (ν : ℕ → ℝ) :
    Tendsto (fun i ↦ phaseInverse (H i) (ν i) - phaseInverse H₀ (ν i)) atTop (𝓝 0) := by
  obtain ⟨c, C, hc, hcurv, hcurv₀⟩ := hR.curv
  have hspec := fun i ↦ phaseInverse_spec (hR.smooth i) hc (fun x ↦ (hcurv i x).2)
  have hspec₀ := phaseInverse_spec hR.smooth_lim hc (fun x ↦ (hcurv₀ x).2)
  have hd₀ : ∀ x, deriv (deriv H₀) x ≤ -c := fun x ↦ by
    simpa [iteratedDeriv_succ, iteratedDeriv_one] using (hcurv₀ x).2
  rw [Metric.tendsto_nhds]
  intro ε hε
  filter_upwards [hR.conv 1 le_rfl (ε * c / 2) (by positivity)] with i hi
  rw [dist_zero_right]
  set y := phaseInverse (H i) (ν i)
  have hy : deriv (H i) y = ν i := (hspec i).2.2 _
  have hy₀ : phaseInverse H₀ (deriv H₀ y) = y := hspec₀.2.1 y
  have hlip := abs_inverse_sub_le_of_neg
    ((hR.smooth_lim.deriv' (n := ∞)).differentiable (by simp)) hc hd₀ hspec₀.2.2
    (deriv H₀ y) (ν i)
  rw [hy₀] at hlip
  have hclose := hi y
  simp only [iteratedDeriv_one] at hclose
  rw [Real.norm_eq_abs]
  calc |y - phaseInverse H₀ (ν i)| ≤ |deriv H₀ y - ν i| / c := hlip
    _ = ‖deriv (H i) y - deriv H₀ y‖ / c := by
        rw [← hy, Real.norm_eq_abs, abs_sub_comm]
    _ ≤ (ε * c / 2) / c := div_le_div_of_nonneg_right hclose hc.le
    _ < ε := by
        rw [show ε * c / 2 / c = ε / 2 by field_simp]
        linarith

/-- Derivatives of the stationary-point maps converge and are bounded at bounded points. -/
private theorem derivativesClose_phaseInverse (hR : RegularizedPhases σ a b G H H₀) {ν : ℕ → ℝ}
    (hν : IsBoundedUnder (· ≤ ·) atTop (fun i ↦ ‖ν i‖)) (n : ℕ) :
    DerivativesClose (fun i ↦ phaseInverse (H i)) (fun _ ↦ phaseInverse H₀) ν ν n := by
  obtain ⟨c, C, hc, hcurv, hcurv₀⟩ := hR.curv
  have hspec := fun i ↦ phaseInverse_spec (hR.smooth i) hc (fun x ↦ (hcurv i x).2)
  have hspec₀ := phaseInverse_spec hR.smooth_lim hc (fun x ↦ (hcurv₀ x).2)
  have hd₀ : ∀ x, deriv (deriv H₀) x ≤ -c := fun x ↦ by
    simpa [iteratedDeriv_succ, iteratedDeriv_one] using (hcurv₀ x).2
  -- Values.
  have h0 : ((fun i ↦ phaseInverse (H i) (ν i)) =o (fun i ↦ phaseInverse H₀ (ν i)) ∧
      VariableObject.IsBounded (fun i ↦ phaseInverse H₀ (ν i))) := by
    refine ⟨(isEqUpToInfinitesimal_iff_tendsto_sub _ _).2 (hR.tendsto_phaseInverse_sub ν), ?_⟩
    obtain ⟨B, hB⟩ := hν
    refine ⟨|phaseInverse H₀ 0| + B / c, ?_⟩
    filter_upwards [hB] with i hi
    have hlip := abs_inverse_sub_le_of_neg
      ((hR.smooth_lim.deriv' (n := ∞)).differentiable (by simp)) hc hd₀ hspec₀.2.2
      (ν i) 0
    rw [Real.norm_eq_abs]
    calc |phaseInverse H₀ (ν i)| ≤ |phaseInverse H₀ (ν i) - phaseInverse H₀ 0| +
          |phaseInverse H₀ 0| := by
          have := abs_add_le (phaseInverse H₀ (ν i) - phaseInverse H₀ 0) (phaseInverse H₀ 0)
          simpa using this
      _ ≤ |ν i - 0| / c + |phaseInverse H₀ 0| := add_le_add hlip le_rfl
      _ ≤ B / c + |phaseInverse H₀ 0| := by
          apply add_le_add _ le_rfl
          rw [sub_zero]
          exact div_le_div_of_nonneg_right (by simpa [Real.norm_eq_abs] using hi) hc.le
      _ = _ := by ring
  -- Derivative bounds and convergence of the curvature at the stationary points.
  have hcurvature : DerivativesClose (fun i ↦ deriv (deriv (H i))) (fun _ ↦ deriv (deriv H₀))
      (fun i ↦ phaseInverse (H i) (ν i)) (fun i ↦ phaseInverse H₀ (ν i)) n := by
    have h2 : ∀ (F : ℝ → ℝ), deriv (deriv F) = iteratedDeriv 2 F := fun F ↦ by
      rw [iteratedDeriv_succ, iteratedDeriv_one]
    simp only [h2]
    apply DerivativesClose.of_uniform
    · rw [← h2]
      exact ((hR.smooth_lim.deriv' (n := ∞))).deriv' (n := ∞)
    · intro k ε hε
      filter_upwards [hR.conv (k + 2) (by omega) ε hε] with i hi x
      simpa only [iteratedDeriv_iteratedDeriv] using hi x
    · intro k
      obtain ⟨K, hK⟩ := hR.bdd (k + 2) (by omega)
      refine ⟨K, fun x ↦ ?_⟩
      simpa only [iteratedDeriv_iteratedDeriv] using hK x
    · exact hR.tendsto_phaseInverse_sub ν
  have hinv := DerivativesClose.inverse (X := fun i ↦ phaseInverse (H i))
    (X' := fun _ ↦ phaseInverse H₀) (Y := fun i ↦ deriv (H i)) (Y' := fun _ ↦ deriv H₀)
    (x := ν) (x' := ν) (n := n)
    (fun i ↦ (hspec i).1) (fun _ ↦ hspec₀.1)
    (fun i ↦ (hR.smooth i).deriv' (n := ∞))
    (fun _ ↦ hR.smooth_lim.deriv' (n := ∞))
    (fun i y ↦ (hspec i).2.2 y) (fun _ y ↦ hspec₀.2.2 y) hcurvature h0.1 h0.2 hc
    (Eventually.of_forall fun i ↦ by
      rw [Real.norm_eq_abs]
      have := hd₀ (phaseInverse H₀ (ν i))
      rw [abs_of_neg (by linarith)]
      linarith)
  exact hinv.mono (Nat.le_succ n)

/-- On every block inside the reference range, the normalized dual phases form a model
phase of exponent `1 / σ`. -/
theorem isModelPhaseFunctionWith_dualPhase (hR : RegularizedPhases σ a b G H H₀)
    (hσ : 0 < σ) (ha : 0 < a) (hab : a < b) {κ : ℝ} (hκ : 0 < κ) (hlo : b ^ (-σ) < κ)
    (hhi : 2 * κ < a ^ (-σ)) :
    IsModelPhaseFunctionWith (fun i ↦ dualPhase (H i) κ (κ ^ (1 - 1 / σ))) (1 / σ) := by
  obtain ⟨c, C, hc, hcurv, hcurv₀⟩ := hR.curv
  set S := κ ^ (1 - 1 / σ) with hS_def
  have hS : 0 < S := Real.rpow_pos_of_pos hκ _
  have hsmooth : ∀ i, ContDiff ℝ ∞ (dualPhase (H i) κ S) := fun i ↦
    contDiff_dualPhase (hR.smooth i) hc (fun x ↦ (hcurv i x).2) κ S
  have hψs : ∀ i, ContDiff ℝ ∞ (phaseInverse (H i)) := fun i ↦
    (phaseInverse_spec (hR.smooth i) hc (fun x ↦ (hcurv i x).2)).1
  have hψ₀s : ContDiff ℝ ∞ (phaseInverse H₀) :=
    (phaseInverse_spec hR.smooth_lim hc (fun x ↦ (hcurv₀ x).2)).1
  -- Derivatives of the dual phase.
  have hder : ∀ (ψ : ℝ → ℝ), ContDiff ℝ ∞ ψ → ∀ (p : ℕ) (u : ℝ),
      iteratedDeriv p (fun u ↦ κ / S * ψ (κ * u)) u =
        κ / S * (κ ^ p * iteratedDeriv p ψ (κ * u)) := by
    intro ψ hψ p u
    rw [iteratedDeriv_const_mul (f := fun u ↦ ψ (κ * u)) _ ((hψ.comp (contDiff_const.mul
      contDiff_id)).contDiffAt.of_le
      (by exact_mod_cast le_top))]
    rw [iteratedDeriv_comp_const_mul (hψ.of_le (by exact_mod_cast le_top)) κ]
  -- The reference phase is the limit dual derivative near `[1, 2]`.
  set U : Set ℝ := Ioo (b ^ (-σ) / κ) (a ^ (-σ) / κ) with hU_def
  have hUopen : IsOpen U := isOpen_Ioo
  have hIU : phaseInterval ⊆ U := by
    intro u hu
    constructor
    · rw [div_lt_iff₀ hκ]
      nlinarith [hu.1]
    · rw [lt_div_iff₀ hκ]
      nlinarith [hu.2]
  have hrefeq : ∀ u ∈ U, modelPhase (1 / σ) u = κ / S * phaseInverse H₀ (κ * u) := by
    intro u hu
    have hκu : κ * u ∈ Icc (b ^ (-σ)) (a ^ (-σ)) := by
      constructor
      · have := hu.1
        rw [div_lt_iff₀ hκ] at this
        linarith
      · have := hu.2
        rw [lt_div_iff₀ hκ] at this
        linarith
    have hu0 : 0 < u := by
      have hb0 : 0 < b ^ (-σ) / κ := div_pos (Real.rpow_pos_of_pos (ha.trans hab) _) hκ
      linarith [hu.1]
    rw [hR.phaseInverse_lim_eq hσ ha hab hκu, Real.mul_rpow hκ.le hu0.le, modelPhase, hS_def]
    rw [show -(1 / σ) = -1 / σ by ring]
    rw [show κ / κ ^ (1 - 1 / σ) * (κ ^ (-1 / σ) * u ^ (-1 / σ)) =
      (κ ^ (1 : ℝ) / κ ^ (1 - 1 / σ) * κ ^ (-1 / σ)) * u ^ (-1 / σ) by rw [Real.rpow_one]; ring]
    rw [div_mul_eq_mul_div, ← Real.rpow_add hκ, ← Real.rpow_sub hκ,
      show (1 : ℝ) + -1 / σ - (1 - 1 / σ) = 0 by ring, Real.rpow_zero, one_mul]
  refine ⟨fun i ↦ (hsmooth i).contDiffOn, ?_⟩
  intro p u
  rw [VariableObject.IsInfinitesimal]
  set x : ℕ → ℝ := fun i ↦ κ * (u i : ℝ) with hx_def
  have hxb : IsBoundedUnder (· ≤ ·) atTop (fun i ↦ ‖x i‖) := by
    refine ⟨κ * 2, ?_⟩
    simp only [eventually_map]
    exact Eventually.of_forall fun i ↦ by
      rw [Real.norm_eq_abs, abs_of_nonneg (by
        have := (u i).2.1
        simp only [hx_def]
        positivity)]
      simp only [hx_def]
      nlinarith [(u i).2.2]
  have hderivatives := (hR.derivativesClose_phaseInverse hxb p) p le_rfl
  have hlim := ((isEqUpToInfinitesimal_iff_tendsto_sub _ _).1 hderivatives.1).const_mul
    (κ / S * κ ^ p)
  rw [mul_zero] at hlim
  have hlim' := hlim.norm
  rw [norm_zero] at hlim'
  refine hlim'.congr fun i ↦ ?_
  congr 1
  have hu : (u i : ℝ) ∈ phaseInterval := (u i).2
  simp only [modelPhaseError_apply, modelPhaseErrorAt]
  rw [iteratedDerivWithin_eq_iteratedDeriv uniqueDiffOn_phaseInterval
      ((hsmooth i).contDiffAt.of_le (by exact_mod_cast le_top)) hu,
    iteratedDerivWithin_eq_iteratedDeriv (f := modelPhase (1 / σ)) uniqueDiffOn_phaseInterval
      (show ContDiffAt ℝ p (modelPhase (1 / σ)) (u i : ℝ) from
        Real.contDiffAt_rpow_const_of_ne (ne_of_gt (lt_of_lt_of_le one_pos hu.1))) hu]
  rw [iteratedDeriv_succ', deriv_dualPhase (hR.smooth i) hc (fun x ↦ (hcurv i x).2),
    hder _ (hψs i)]
  have hevent : modelPhase (1 / σ) =ᶠ[𝓝 (u i : ℝ)]
      fun u ↦ κ / S * phaseInverse H₀ (κ * u) := by
    filter_upwards [hUopen.mem_nhds (hIU hu)] with v hv
    exact hrefeq v hv
  rw [hevent.iteratedDeriv_eq, hder _ hψ₀s]
  simp only [hx_def]
  ring

end RegularizedPhases

end DualPhase

end

end Expdb.BProcess
