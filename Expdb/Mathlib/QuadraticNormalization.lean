module

public import Expdb.Mathlib.InverseFunction

public import Mathlib.Analysis.Calculus.IteratedDeriv.Lemmas
public import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic
import Mathlib.Analysis.Calculus.TaylorIntegral
import Mathlib.Analysis.Calculus.ParametricIntegral
import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.Tactic.Ring

/-!
# The quadratic factor of a phase near a stationary point

Mathlib's integral Taylor theorem factors a phase minus its tangent as a
quadratic polynomial times an averaged second derivative. This representation
includes the stationary point itself, where dividing by `(x - v)²` would be
singular. Curvature bounds give a smooth quadratic coordinate and derivative estimates.
-/

@[expose] public section
namespace Expdb

open MeasureTheory
open scoped ContDiff Topology
noncomputable section

/-- Twice the averaged second derivative along the segment from `v` to `x`. -/
def phaseQuadraticFactor (F : ℝ → ℝ) (v x : ℝ) : ℝ :=
  2 * ∫ t : ℝ in 0..1, (1 - t) * iteratedDeriv 2 F (v + t * (x - v))

/-- The phase minus its tangent has an exact quadratic factor, including at
`x = v`. Only twice continuous differentiability along the segment is needed. -/
theorem phase_sub_tangent_eq_quadraticFactor (F : ℝ → ℝ) (v x : ℝ)
    (hF : ∀ t ∈ Set.Icc (0 : ℝ) 1, ContDiffAt ℝ 2 F (v + t * (x - v))) :
    F x - F v - deriv F v * (x - v) =
      (x - v) ^ 2 * phaseQuadraticFactor F v x / 2 := by
  have h := map_add_eq_sum_add_integral_iteratedFDeriv
    (n := 1) (x := v) (y := x - v) hF
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add,
    iteratedFDeriv_apply_eq_iteratedDeriv_mul_prod, Finset.prod_const,
    Finset.card_univ, Fintype.card_fin, pow_zero, pow_one, Nat.factorial_zero,
    Nat.factorial_one, Nat.cast_one, inv_one, iteratedDeriv_zero,
    iteratedDeriv_one, smul_eq_mul, one_mul, Nat.reduceAdd] at h
  rw [show v + (x - v) = x by ring] at h
  have he : (∫ t : ℝ in 0..1,
      (1 - t) * ((x - v) ^ 2 * iteratedDeriv 2 F (v + t * (x - v)))) =
      (x - v) ^ 2 * ∫ t : ℝ in 0..1,
        (1 - t) * iteratedDeriv 2 F (v + t * (x - v)) := by
    rw [← intervalIntegral.integral_const_mul]
    congr 1
    funext t
    ring
  rw [he] at h
  dsimp [phaseQuadraticFactor]
  linarith

private theorem integral_one_sub : (∫ t : ℝ in 0..1, 1 - t) = 1 / 2 := by
  rw [intervalIntegral.integral_sub (f := fun _ : ℝ => (1 : ℝ))
    (g := fun t : ℝ => t) intervalIntegrable_const
    (continuous_id.intervalIntegrable _ _), intervalIntegral.integral_const, integral_id]
  norm_num

/-- At the stationary point the factor is the second derivative itself. -/
@[simp] theorem phaseQuadraticFactor_self (F : ℝ → ℝ) (v : ℝ) :
    phaseQuadraticFactor F v v = iteratedDeriv 2 F v := by
  simp only [phaseQuadraticFactor, sub_self, mul_zero, add_zero]
  rw [intervalIntegral.integral_mul_const]
  rw [integral_one_sub]
  ring

/-- Bounds on curvature along the segment also bound the quadratic factor.
In particular, positive curvature keeps its square root away from zero. -/
theorem phaseQuadraticFactor_bounds (F : ℝ → ℝ) (v x : ℝ) {c C : ℝ}
    (hcontinuous : ContinuousOn (fun t : ℝ => iteratedDeriv 2 F (v + t * (x - v)))
      (Set.Icc 0 1))
    (hlower : ∀ t ∈ Set.Icc (0 : ℝ) 1, c ≤ iteratedDeriv 2 F (v + t * (x - v)))
    (hupper : ∀ t ∈ Set.Icc (0 : ℝ) 1, iteratedDeriv 2 F (v + t * (x - v)) ≤ C) :
    c ≤ phaseQuadraticFactor F v x ∧ phaseQuadraticFactor F v x ≤ C := by
  have hi : IntervalIntegrable
      (fun t : ℝ => (1 - t) * iteratedDeriv 2 F (v + t * (x - v))) volume 0 1 := by
    have hc : ContinuousOn (fun t : ℝ => (1 - t) *
        iteratedDeriv 2 F (v + t * (x - v))) (Set.Icc 0 1) :=
      (continuous_const.sub continuous_id).continuousOn.mul hcontinuous
    apply ContinuousOn.intervalIntegrable
    simpa only [Set.uIcc_of_le (show (0 : ℝ) ≤ 1 by norm_num)] using hc
  have hconst (a : ℝ) : (∫ t : ℝ in 0..1, (1 - t) * a) = a / 2 := by
    rw [intervalIntegral.integral_mul_const]
    rw [integral_one_sub]
    ring
  have hl := intervalIntegral.integral_mono_on (show (0 : ℝ) ≤ 1 by norm_num)
    ((by fun_prop : Continuous (fun t : ℝ => (1 - t) * c)).intervalIntegrable 0 1) hi
    (fun t ht => mul_le_mul_of_nonneg_left (hlower t ht) (sub_nonneg.2 ht.2))
  have hu := intervalIntegral.integral_mono_on (show (0 : ℝ) ≤ 1 by norm_num) hi
    ((by fun_prop : Continuous (fun t : ℝ => (1 - t) * C)).intervalIntegrable 0 1)
    (fun t ht => mul_le_mul_of_nonneg_left (hupper t ht) (sub_nonneg.2 ht.2))
  rw [hconst] at hl hu
  dsimp [phaseQuadraticFactor]
  constructor <;> linarith

private theorem hasDerivAt_average {f f' : ℝ → ℝ → ℝ} (x : ℝ)
    (hf : Continuous (Function.uncurry f))
    (hf' : Continuous (Function.uncurry f'))
    (hd : ∀ x t, HasDerivAt (fun y => f y t) (f' x t) x) :
    HasDerivAt (fun y => ∫ t : ℝ in 0..1, f y t)
      (∫ t : ℝ in 0..1, f' x t) x := by
  obtain ⟨B, hB⟩ := ((isCompact_Icc : IsCompact (Set.Icc (x - 1) (x + 1))).prod
    (isCompact_Icc : IsCompact (Set.Icc (0 : ℝ) 1))).bddAbove_image
    (hf'.norm.continuousOn : ContinuousOn
      (fun p : ℝ × ℝ => ‖f' p.1 p.2‖)
      (Set.Icc (x - 1) (x + 1) ×ˢ Set.Icc 0 1))
  have hm (y : ℝ) : Continuous (f y) :=
    hf.comp (continuous_const.prodMk continuous_id)
  have hm' : Continuous (f' x) :=
    hf'.comp (continuous_const.prodMk continuous_id)
  have hb : ∀ᵐ t ∂volume.restrict (Set.Ioc (0 : ℝ) 1),
      ∀ y ∈ Set.Ioo (x - 1) (x + 1), ‖f' y t‖ ≤ B := by
    filter_upwards [ae_restrict_mem measurableSet_Ioc] with t ht
    intro y hy
    have hp : (y, t) ∈ Set.Icc (x - 1) (x + 1) ×ˢ Set.Icc (0 : ℝ) 1 :=
      ⟨⟨hy.1.le, hy.2.le⟩, ⟨ht.1.le, ht.2⟩⟩
    exact hB (Set.mem_image_of_mem (fun p : ℝ × ℝ => ‖f' p.1 p.2‖) hp)
  have hd' : ∀ᵐ t ∂volume.restrict (Set.Ioc (0 : ℝ) 1),
      ∀ y ∈ Set.Ioo (x - 1) (x + 1), HasDerivAt (fun z => f z t) (f' y t) y :=
    Filter.Eventually.of_forall (fun t y _ => hd y t)
  have hi := hasDerivAt_integral_of_dominated_loc_of_deriv_le
    (μ := volume.restrict (Set.Ioc (0 : ℝ) 1))
    (bound := fun _ => B)
    (Ioo_mem_nhds (by linarith : x - 1 < x) (by linarith : x < x + 1))
    (Filter.Eventually.of_forall (fun y => (hm y).aestronglyMeasurable))
    ((hm x).intervalIntegrable 0 1).1 hm'.aestronglyMeasurable hb
    (intervalIntegrable_const : IntervalIntegrable (fun _ : ℝ => B) volume 0 1).1 hd'
  simpa only [intervalIntegral.integral_of_le (by norm_num : (0 : ℝ) ≤ 1)] using hi.2

/-- The averaged Taylor factor has a derivative without a singular quotient
at the stationary point. -/
theorem hasDerivAt_phaseQuadraticFactor (F : ℝ → ℝ) (v x : ℝ)
    (hF : ContDiff ℝ 3 F) :
    HasDerivAt (phaseQuadraticFactor F v)
      (2 * ∫ t : ℝ in 0..1,
        (1 - t) * t * iteratedDeriv 3 F (v + t * (x - v))) x := by
  have h₂ := hF.continuous_iteratedDeriv 2 (by norm_num)
  have h₃ := hF.continuous_iteratedDeriv 3 le_rfl
  have hd₂ := hF.differentiable_iteratedDeriv 2 (by norm_num)
  have hd (y t : ℝ) : HasDerivAt
      (fun z : ℝ => (1 - t) * iteratedDeriv 2 F (v + t * (z - v)))
      ((1 - t) * t * iteratedDeriv 3 F (v + t * (y - v))) y := by
    have he : HasDerivAt (fun z : ℝ => v + t * (z - v)) t y := by
      convert! (((hasDerivAt_id y).sub_const v).const_mul t).const_add v using 1
      ring
    have hh := ((hd₂ _).hasDerivAt.comp y he).const_mul (1 - t)
    rw [← iteratedDeriv_succ] at hh
    convert! hh using 1
    ring
  exact (hasDerivAt_average x (by fun_prop) (by fun_prop) hd).const_mul 2

private def quadraticFactorJet (F : ℝ → ℝ) (v : ℝ) (j : ℕ) (x : ℝ) : ℝ :=
  2 * ∫ t : ℝ in 0..1,
    (1 - t) * t ^ j * iteratedDeriv (j + 2) F (v + t * (x - v))

private theorem hasDerivAt_quadraticFactorJet (F : ℝ → ℝ) (v x : ℝ) (j : ℕ)
    (hF : ContDiff ℝ ∞ F) :
    HasDerivAt (quadraticFactorJet F v j) (quadraticFactorJet F v (j + 1) x) x := by
  have h₂ := hF.continuous_iteratedDeriv (j + 2)
    (ENat.natCast_le_of_coe_top_le_withTop le_rfl _)
  have h₃ := hF.continuous_iteratedDeriv (j + 3)
    (ENat.natCast_le_of_coe_top_le_withTop le_rfl _)
  have hd₂ := (hF.of_le (ENat.natCast_le_of_coe_top_le_withTop le_rfl (j + 2 +
    1))).differentiable_iteratedDeriv'
    (j + 2)
  have hd (y t : ℝ) : HasDerivAt
      (fun z : ℝ => (1 - t) * t ^ j * iteratedDeriv (j + 2) F (v + t * (z - v)))
      ((1 - t) * t ^ (j + 1) * iteratedDeriv (j + 3) F (v + t * (y - v))) y := by
    have he : HasDerivAt (fun z : ℝ => v + t * (z - v)) t y := by
      convert! (((hasDerivAt_id y).sub_const v).const_mul t).const_add v using 1
      ring
    have hh := ((hd₂ _).hasDerivAt.comp y he).const_mul ((1 - t) * t ^ j)
    rw [← iteratedDeriv_succ] at hh
    convert! hh using 1
    simp only [pow_succ, Nat.add_assoc]
    ring
  simpa only [quadraticFactorJet, show j + 1 + 2 = j + 3 by omega] using!
    (hasDerivAt_average x (by fun_prop) (by fun_prop) hd).const_mul 2

private theorem iteratedDeriv_quadraticFactorJet (F : ℝ → ℝ) (v : ℝ)
    (hF : ContDiff ℝ ∞ F) (j : ℕ) :
    iteratedDeriv j (phaseQuadraticFactor F v) = quadraticFactorJet F v j := by
  induction j with
  | zero =>
    funext x
    simp [quadraticFactorJet, phaseQuadraticFactor]
  | succ j ih =>
    rw [iteratedDeriv_succ, ih]
    funext x
    exact (hasDerivAt_quadraticFactorJet F v x j hF).deriv

/-- Derivatives of the quadratic factor are averages of higher phase derivatives.
Only the displayed derivative order enters the integrand. -/
theorem iteratedDeriv_phaseQuadraticFactor (F : ℝ → ℝ) (v x : ℝ) (j : ℕ)
    (hF : ContDiff ℝ ∞ F) :
    iteratedDeriv j (phaseQuadraticFactor F v) x =
      2 * ∫ t : ℝ in 0..1,
        (1 - t) * t ^ j * iteratedDeriv (j + 2) F (v + t * (x - v)) := by
  rw [iteratedDeriv_quadraticFactorJet F v hF j]
  rfl

/-- The integral factor stays smooth across the stationary point. -/
theorem contDiff_phaseQuadraticFactor (F : ℝ → ℝ) (v : ℝ)
    (hF : ContDiff ℝ ∞ F) : ContDiff ℝ ∞ (phaseQuadraticFactor F v) := by
  apply contDiff_iff_iteratedDeriv.mpr
  constructor
  · intro j _
    rw [iteratedDeriv_quadraticFactorJet F v hF j]
    exact continuous_iff_continuousAt.mpr
      (fun x => (hasDerivAt_quadraticFactorJet F v x j hF).continuousAt)
  · intro j _
    rw [iteratedDeriv_quadraticFactorJet F v hF j]
    exact fun x => (hasDerivAt_quadraticFactorJet F v x j hF).differentiableAt

/-- The increasing quadratic coordinate for a phase of negative curvature. -/
def phaseQuadraticCoordinate (F : ℝ → ℝ) (v x : ℝ) : ℝ :=
  (x - v) * Real.sqrt (-phaseQuadraticFactor F v x)

@[simp] theorem phaseQuadraticCoordinate_self (F : ℝ → ℝ) (v : ℝ) :
    phaseQuadraticCoordinate F v v = 0 := by
  simp [phaseQuadraticCoordinate]

/-- The coordinate puts the phase minus its tangent in exact negative quadratic
form. At a stationary point this is precisely the stationary-phase normal form. -/
theorem phase_sub_tangent_eq_neg_coordinate_sq (F : ℝ → ℝ) (v x : ℝ)
    (hF : ∀ t ∈ Set.Icc (0 : ℝ) 1, ContDiffAt ℝ 2 F (v + t * (x - v)))
    (hq : phaseQuadraticFactor F v x ≤ 0) :
    F x - F v - deriv F v * (x - v) =
      -(phaseQuadraticCoordinate F v x) ^ 2 / 2 := by
  rw [phase_sub_tangent_eq_quadraticFactor F v x hF]
  simp only [phaseQuadraticCoordinate, mul_pow,
    Real.sq_sqrt (neg_nonneg.mpr hq)]
  ring

/-- Nonzero curvature makes the coordinate smooth, including at its centre. -/
theorem contDiffOn_phaseQuadraticCoordinate (F : ℝ → ℝ) (v : ℝ) {s : Set ℝ}
    (hF : ContDiff ℝ ∞ F)
    (hq : ∀ x ∈ s, phaseQuadraticFactor F v x < 0) :
    ContDiffOn ℝ ∞ (phaseQuadraticCoordinate F v) s := by
  exact (contDiff_id.sub contDiff_const).contDiffOn.mul
    ((contDiff_phaseQuadraticFactor F v hF).neg.contDiffOn.sqrt
      (fun x hx => ne_of_gt (neg_pos.mpr (hq x hx))))

/-- Differentiating the coordinate does not divide by `x - v`. -/
theorem hasDerivAt_phaseQuadraticCoordinate (F : ℝ → ℝ) (v x : ℝ)
    (hF : ContDiff ℝ 3 F) (hq : phaseQuadraticFactor F v x < 0) :
    HasDerivAt (phaseQuadraticCoordinate F v)
      (Real.sqrt (-phaseQuadraticFactor F v x) -
        (x - v) * (2 * ∫ t : ℝ in 0..1,
          (1 - t) * t * iteratedDeriv 3 F (v + t * (x - v))) /
          (2 * Real.sqrt (-phaseQuadraticFactor F v x))) x := by
  have hs := (hasDerivAt_phaseQuadraticFactor F v x hF).neg.sqrt
    (ne_of_gt (neg_pos.mpr hq))
  convert! ((hasDerivAt_id x).sub_const v).mul hs using 1
  simp only [Pi.neg_apply, id_eq]
  ring

/-- The derivative at the stationary point is the square root of its curvature. -/
theorem hasDerivAt_phaseQuadraticCoordinate_self (F : ℝ → ℝ) (v : ℝ)
    (hF : ContDiff ℝ 3 F) (hcurvature : iteratedDeriv 2 F v < 0) :
    HasDerivAt (phaseQuadraticCoordinate F v) (Real.sqrt (-iteratedDeriv 2 F v)) v := by
  simpa only [phaseQuadraticFactor_self, sub_self, zero_mul, zero_div, sub_zero] using
    hasDerivAt_phaseQuadraticCoordinate F v v hF
      (by simpa only [phaseQuadraticFactor_self] using hcurvature)

/-- Differentiating the normal form identifies the coordinate's Jacobian. -/
theorem phaseQuadraticCoordinate_mul_deriv (F : ℝ → ℝ) (v x : ℝ)
    (hF : ContDiff ℝ 3 F) (hq : phaseQuadraticFactor F v x < 0) :
    phaseQuadraticCoordinate F v x * deriv (phaseQuadraticCoordinate F v) x =
      -(deriv F x - deriv F v) := by
  have hd := hasDerivAt_phaseQuadraticCoordinate F v x hF hq
  have he : (fun y => -(phaseQuadraticCoordinate F v y) ^ 2 / 2) =ᶠ[𝓝 x]
      (fun y => F y - F v - deriv F v * (y - v)) := by
    have hn := (hasDerivAt_phaseQuadraticFactor F v x hF).continuousAt.eventually_lt_const hq
    filter_upwards [hn] with y hy
    exact (phase_sub_tangent_eq_neg_coordinate_sq F v y
      (fun _ _ => (hF.of_le (by norm_num : (2 : ℕ∞ω) ≤ 3)).contDiffAt) hy.le).symm
  have ht := ((hF.differentiable (by norm_num) x).hasDerivAt).sub_const (F v)
  have ht' := ht.sub (((hasDerivAt_id x).sub_const v).const_mul (deriv F v))
  have hn := ((hd.pow 2).neg.div_const 2).congr_of_eventuallyEq he.symm
  have hu := hn.unique ht'
  rw [hd.deriv] at ⊢
  simp only [mul_one] at hu
  nlinarith

/-- Uniform negative curvature gives a uniform positive lower bound on the
Jacobian, including at the stationary point. -/
theorem phaseQuadraticCoordinate_deriv_lower (F : ℝ → ℝ) (v : ℝ) {s : Set ℝ}
    (hs : Convex ℝ s) (hv : v ∈ s) {c C : ℝ} (hc : 0 < c)
    (hF : ContDiff ℝ 3 F)
    (hlower : ∀ x ∈ s, -C ≤ iteratedDeriv 2 F x)
    (hupper : ∀ x ∈ s, iteratedDeriv 2 F x ≤ -c)
    (x : ℝ) (hx : x ∈ s) :
    c / Real.sqrt C ≤ deriv (phaseQuadraticCoordinate F v) x := by
  have hseg (t : ℝ) (ht : t ∈ Set.Icc (0 : ℝ) 1) : v + t * (x - v) ∈ s :=
    hs.add_smul_sub_mem hv hx ht
  have hcont := hF.continuous_iteratedDeriv 2 (by norm_num : (2 : ℕ∞ω) ≤ 3)
  have hcontseg : ContinuousOn (fun t : ℝ => iteratedDeriv 2 F (v + t * (x - v)))
      (Set.Icc 0 1) := (hcont.comp (by fun_prop)).continuousOn
  obtain ⟨hqlo, hqhi⟩ := phaseQuadraticFactor_bounds F v x hcontseg
    (fun t ht => hlower _ (hseg t ht)) (fun t ht => hupper _ (hseg t ht))
  have hq : phaseQuadraticFactor F v x < 0 := lt_of_le_of_lt hqhi (neg_neg_of_pos hc)
  have hC : 0 < C := by linarith [hlower v hv, hupper v hv]
  have hsqrtC : 0 < Real.sqrt C := Real.sqrt_pos.mpr hC
  have hsqc := Real.sq_sqrt hc.le
  have hsqq := Real.sq_sqrt (show 0 ≤ -phaseQuadraticFactor F v x by linarith)
  have hj := phaseQuadraticCoordinate_mul_deriv F v x hF hq
  have hdx := hF.differentiable_iteratedDeriv 1 (by norm_num : (1 : ℕ∞ω) < 3)
  have hdy : Continuous (deriv F) := by simpa only [iteratedDeriv_one] using hdx.continuous
  have hddy : Differentiable ℝ (deriv F) := by simpa only [iteratedDeriv_one] using hdx
  have hm (a : ℝ) (ha : a ∈ s) (b : ℝ) (hb : b ∈ s) (hab : a ≤ b) :
      deriv F b - deriv F a ≤ -c * (b - a) := by
    apply hs.image_sub_le_mul_sub_of_deriv_le hdy.continuousOn hddy.differentiableOn
      (fun y hy => ?_) a ha b hb hab
    simpa only [iteratedDeriv_succ' (n := 1), iteratedDeriv_one] using
      hupper y (interior_subset hy)
  have hroot := Real.sqrt_le_sqrt (show -phaseQuadraticFactor F v x ≤ C by linarith)
  have hrootpos : 0 < Real.sqrt (-phaseQuadraticFactor F v x) :=
    Real.sqrt_pos.mpr (neg_pos.mpr hq)
  by_cases he : x = v
  · subst x
    rw [(hasDerivAt_phaseQuadraticCoordinate_self F v hF
      (lt_of_le_of_lt (hupper v hv) (neg_neg_of_pos hc))).deriv]
    apply (div_le_iff₀ hsqrtC).mpr
    have hs := Real.sqrt_le_sqrt (show c ≤ -iteratedDeriv 2 F v by linarith [hupper v hv])
    have hsC := Real.sqrt_le_sqrt (show c ≤ C by linarith [hlower v hv, hupper v hv])
    nlinarith [Real.sqrt_nonneg c, Real.sqrt_nonneg (-iteratedDeriv 2 F v)]
  · have hxv : x < v ∨ v < x := lt_or_gt_of_ne he
    rcases hxv with hlt | hgt
    · have hg := hm x hx v hv hlt.le
      have hy : phaseQuadraticCoordinate F v x < 0 := by
        exact mul_neg_of_neg_of_pos (sub_neg.mpr hlt) hrootpos
      have hp : 0 < deriv (phaseQuadraticCoordinate F v) x := by
        nlinarith
      apply (div_le_iff₀ hsqrtC).mpr
      dsimp [phaseQuadraticCoordinate] at hj
      nlinarith [mul_le_mul_of_nonneg_right hroot hp.le]
    · have hg := hm v hv x hx hgt.le
      have hy : 0 < phaseQuadraticCoordinate F v x := by
        exact mul_pos (sub_pos.mpr hgt) hrootpos
      have hp : 0 < deriv (phaseQuadraticCoordinate F v) x := by
        nlinarith
      apply (div_le_iff₀ hsqrtC).mpr
      dsimp [phaseQuadraticCoordinate] at hj
      nlinarith [mul_le_mul_of_nonneg_right hroot hp.le]

/-- Uniform bounds on phase derivatives give the same bounds on the derivatives
of its averaged quadratic factor; the constant does not depend on the centre. -/
theorem norm_iteratedDeriv_phaseQuadraticFactor_le (F : ℝ → ℝ) (v x : ℝ) (j : ℕ)
    (hF : ContDiff ℝ ∞ F) {C : ℝ}
    (hbound : ∀ t ∈ Set.Icc (0 : ℝ) 1,
      ‖iteratedDeriv (j + 2) F (v + t * (x - v))‖ ≤ C) :
    ‖iteratedDeriv j (phaseQuadraticFactor F v) x‖ ≤ C := by
  have hC : 0 ≤ C := (norm_nonneg _).trans (hbound 0 (by norm_num))
  rw [iteratedDeriv_phaseQuadraticFactor F v x j hF, norm_mul, Real.norm_of_nonneg
    (by norm_num : (0 : ℝ) ≤ 2)]
  have hb : ∀ᵐ t : ℝ ∂volume, t ∈ Set.Ioc (0 : ℝ) 1 →
      ‖(1 - t) * t ^ j * iteratedDeriv (j + 2) F (v + t * (x - v))‖ ≤ (1 - t) * C := by
    apply Filter.Eventually.of_forall
    intro t ht
    have ht' : t ∈ Set.Icc (0 : ℝ) 1 := ⟨ht.1.le, ht.2⟩
    rw [norm_mul, norm_mul, Real.norm_of_nonneg (sub_nonneg.mpr ht.2),
      Real.norm_of_nonneg (pow_nonneg ht.1.le _), mul_assoc]
    apply mul_le_mul_of_nonneg_left _ (sub_nonneg.mpr ht.2)
    exact (mul_le_mul_of_nonneg_left (hbound t ht') (pow_nonneg ht.1.le _)).trans
      (by simpa using mul_le_mul_of_nonneg_right (pow_le_one₀ ht.1.le ht.2) hC)
  have hi := intervalIntegral.norm_integral_le_of_norm_le
    (by norm_num : (0 : ℝ) ≤ 1) hb
    ((by fun_prop : Continuous (fun t : ℝ => (1 - t) * C)).intervalIntegrable 0 1)
  rw [intervalIntegral.integral_mul_const, integral_one_sub] at hi
  linarith

/-- Global curvature bounds bound the quadratic factor. -/
theorem phaseQuadraticFactor_mem {H : ℝ → ℝ} (hH : ContDiff ℝ ∞ H) {c C : ℝ}
    (hcurv : ∀ x, -C ≤ iteratedDeriv 2 H x ∧ iteratedDeriv 2 H x ≤ -c) (v x : ℝ) :
    -C ≤ phaseQuadraticFactor H v x ∧ phaseQuadraticFactor H v x ≤ -c :=
  phaseQuadraticFactor_bounds H v x
    ((hH.continuous_iteratedDeriv 2 (ENat.natCast_le_of_coe_top_le_withTop le_rfl 2)).comp
      (by fun_prop)).continuousOn
    (fun _ _ ↦ (hcurv _).1) (fun _ _ ↦ (hcurv _).2)

/-- The quadratic coordinate of a globally concave phase is smooth. -/
theorem contDiff_phaseQuadraticCoordinate {H : ℝ → ℝ} (hH : ContDiff ℝ ∞ H) {c C : ℝ}
    (hc : 0 < c) (hcurv : ∀ x, -C ≤ iteratedDeriv 2 H x ∧ iteratedDeriv 2 H x ≤ -c)
    (v : ℝ) : ContDiff ℝ ∞ (phaseQuadraticCoordinate H v) :=
  contDiffOn_univ.1 (contDiffOn_phaseQuadraticCoordinate H v hH
    (fun x _ ↦ (phaseQuadraticFactor_mem hH hcurv v x).2.trans_lt (by linarith)))

/-- The quadratic coordinate of a globally concave phase has derivative bounded below. -/
theorem le_deriv_phaseQuadraticCoordinate {H : ℝ → ℝ} (hH : ContDiff ℝ ∞ H) {c C : ℝ}
    (hc : 0 < c) (hcurv : ∀ x, -C ≤ iteratedDeriv 2 H x ∧ iteratedDeriv 2 H x ≤ -c)
    (v x : ℝ) : c / Real.sqrt C ≤ deriv (phaseQuadraticCoordinate H v) x :=
  phaseQuadraticCoordinate_deriv_lower H v convex_univ (Set.mem_univ v) hc
    (hH.of_le (ENat.natCast_le_of_coe_top_le_withTop le_rfl 3)) (fun x _ ↦ (hcurv x).1) (fun x
      _ ↦ (hcurv x).2)
    x (Set.mem_univ x)

/-! ## The inverse quadratic coordinate -/

/-- The inverse of the quadratic normal-form coordinate centred at `v`. -/
def quadInverse (H : ℝ → ℝ) (v : ℝ) : ℝ → ℝ :=
  Function.invFun (phaseQuadraticCoordinate H v)

theorem quadInverse_spec {H : ℝ → ℝ} (hH : ContDiff ℝ ∞ H) {c C : ℝ} (hc : 0 < c)
    (hcurv : ∀ x, -C ≤ iteratedDeriv 2 H x ∧ iteratedDeriv 2 H x ≤ -c) (v : ℝ) :
    ContDiff ℝ ∞ (quadInverse H v) ∧
      (∀ y, phaseQuadraticCoordinate H v (quadInverse H v y) = y) ∧
      (∀ x, quadInverse H v (phaseQuadraticCoordinate H v x) = x) := by
  have hC : 0 < C := by linarith [(hcurv 0).1, (hcurv 0).2]
  obtain ⟨g, hg, hgY, hYg⟩ := exists_smooth_inverse_of_le_deriv
    (contDiff_phaseQuadraticCoordinate hH hc hcurv v) (div_pos hc (Real.sqrt_pos.2 hC))
    (le_deriv_phaseQuadraticCoordinate hH hc hcurv v)
  have heq : quadInverse H v = g :=
    Function.invFun_eq_of_injective_of_rightInverse
      (show Function.LeftInverse g _ from hgY).injective hYg
  rw [heq]
  exact ⟨hg, hYg, hgY⟩

/-- The quadratic coordinate grows at most linearly. -/
theorem abs_phaseQuadraticCoordinate_le {H : ℝ → ℝ} (hH : ContDiff ℝ ∞ H) {c C : ℝ}
    (hcurv : ∀ x, -C ≤ iteratedDeriv 2 H x ∧ iteratedDeriv 2 H x ≤ -c) (v x : ℝ) :
    |phaseQuadraticCoordinate H v x| ≤ |x - v| * Real.sqrt C := by
  have hq := phaseQuadraticFactor_mem hH hcurv v x
  rw [phaseQuadraticCoordinate, abs_mul, abs_of_nonneg (Real.sqrt_nonneg _)]
  exact mul_le_mul_of_nonneg_left (Real.sqrt_le_sqrt (by linarith)) (abs_nonneg _)

end

end Expdb
