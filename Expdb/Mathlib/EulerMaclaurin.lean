module

public import Mathlib.Analysis.Calculus.IteratedDeriv.Lemmas
public import Mathlib.NumberTheory.ZetaValues

import Mathlib.MeasureTheory.Integral.DivergenceTheorem
import Mathlib.Tactic.Abel
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Module
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Order

/-!
# Arbitrary-order Euler--Maclaurin summation

This file is adapted from
`Interval/EulerMaclaurin/Bernoulli.lean` and
`Interval/EulerMaclaurin/EulerMaclaurin.lean` in
[`girving/interval`](https://github.com/girving/interval) at commit
`f4a3231fd735cdf6d5265864d9a9207744a7aca2`, licensed under Apache-2.0.

The source has been modified, namespaced for Expdb, and updated for newer version of Mathlib.

-/

namespace Expdb.EulerMaclaurin

open Filter Set
open MeasureTheory (volume)
open scoped Real Topology
open intervalIntegral

noncomputable section

@[expose] public section

/-! ## Periodic Bernoulli functions -/

variable {s : ℕ} {a : ℤ}

private def presaw (s : ℕ) (a : ℤ) (x : ℝ) : ℝ :=
  (s.factorial : ℝ)⁻¹ • bernoulliFun s (x - a)

private lemma contDiff_presaw : ContDiff ℝ ⊤ (presaw s a) := by
  exact ((contDiff_bernoulliFun _).comp (contDiff_id.sub contDiff_const)).const_smul _

private lemma hasDerivAt_presaw {x : ℝ} :
    HasDerivAt (presaw (s + 1) a) (presaw s a x) x := by
  have e : presaw (s + 1) a = (fun x ↦ presaw (s + 1) a x) := rfl
  simp only [presaw, Nat.factorial_succ, mul_comm _ s.factorial, Nat.cast_mul,
    mul_inv, ← smul_smul, e]
  apply HasDerivAt.fun_const_smul
  have s0 : ((s + 1 : ℕ) : ℝ) ≠ 0 := by
    simp only [Nat.cast_ne_zero]
    omega
  have sc : s = s + 1 - 1 := by
    omega
  rw [← inv_smul_smul₀ s0 (x := bernoulliFun s (x - ↑a))]
  nth_rw 5 [sc]
  apply HasDerivAt.fun_const_smul
  rw [Nat.cast_smul_eq_nsmul]
  rw [← mul_one ((s + 1) • bernoulliFun (s + 1 - 1) (x - ↑a))]
  simp only [nsmul_eq_mul]
  exact HasDerivAt.comp _ (hasDerivAt_bernoulliFun (s + 1) _) (h := fun x : ℝ ↦ x - a)
    ((hasDerivAt_id' x).sub_const _)

@[simp] private lemma presaw_zero {x : ℝ} : presaw 0 a x = 1 := by
  simp only [presaw, Nat.factorial_zero, Nat.cast_one, inv_one, bernoulliFun_zero, smul_eq_mul,
    mul_one]

/-- The scaled periodic Bernoulli function of order `s`. -/
def saw (s : ℕ) (x : ℝ) : ℝ :=
  (s.factorial : ℝ)⁻¹ * periodizedBernoulli s (x : AddCircle (1 : ℝ))

private lemma saw_eq_bernoulliFun_fract (s : ℕ) (x : ℝ) :
    saw s x = (s.factorial : ℝ)⁻¹ * bernoulliFun s (Int.fract x) := by
  rw [saw, ← AddCircle.coe_fract x, periodizedBernoulli, AddCircle.liftIco_coe_apply]
  constructor
  · exact Int.fract_nonneg x
  · simpa using Int.fract_lt_one x

private lemma saw_eqOn {a : ℤ} : EqOn (saw s) (presaw s a) (Ico a (a + 1)) := by
  intro x m
  simp only [mem_Ico, ← Int.floor_eq_iff] at m
  simp only [saw_eq_bernoulliFun_fract, presaw, smul_eq_mul, Int.fract, m]

@[simp] private lemma presaw_coe_same {a : ℤ} : presaw s a a = saw s 0 := by
  rw [← saw_eqOn (a := a)]
  all_goals simp [saw_eq_bernoulliFun_fract]

@[simp] private lemma presaw_coe_succ {a : ℤ} :
    presaw (s + 2) a (a + 1) = saw (s + 2) 0 := by
  simp [presaw, bernoulliFun_endpoints_eq_of_ne_one, saw_eq_bernoulliFun_fract]

@[simp] private lemma presaw_one_coe_succ {a : ℤ} : presaw 1 a (a + 1) = 1 / 2 := by
  simp only [presaw, Nat.factorial_one, Nat.cast_one, inv_one, add_sub_cancel_left,
    bernoulliFun_one, one_div, smul_eq_mul, one_mul]
  norm_num

private lemma saw_one {x : ℝ} : saw 1 x = Int.fract x - 1 / 2 := by
  simp only [saw_eq_bernoulliFun_fract, Nat.factorial_one, Nat.cast_one, inv_one,
    bernoulliFun_one, one_div, one_mul]

@[simp] private lemma saw_one_zero : saw 1 0 = -2⁻¹ := by
  simp only [saw_one, Int.fract_zero, one_div, zero_sub]

/-! ## Non-explicit bounds -/

private lemma exists_max_bernoulli (s : ℕ) :
    ∃ x ∈ Icc 0 1, IsMaxOn (fun x ↦ |bernoulliFun s x|) (Icc 0 1) x := by
  apply isCompact_Icc.exists_isMaxOn
  · exact nonempty_Icc.mpr zero_le_one
  · exact (continuous_bernoulliFun s).abs.continuousOn

private def bernoulliBound (s : ℕ) : ℝ :=
  |bernoulliFun s (exists_max_bernoulli s).choose|

private lemma abs_bernoulliFun_le (s : ℕ) (x : ℝ) (m : x ∈ Icc 0 1) :
    |bernoulliFun s x| ≤ bernoulliBound s := by
  simp only [bernoulliBound]
  obtain ⟨_, hmax⟩ := (exists_max_bernoulli s).choose_spec
  exact hmax m

private lemma bddAbove_range_abs_saw (s : ℕ) :
    BddAbove (range fun x ↦ |saw s x|) := by
  refine ⟨(s.factorial : ℝ)⁻¹ * bernoulliBound s, ?_⟩
  rintro _ ⟨x, rfl⟩
  have sp : 0 < (s.factorial : ℝ)⁻¹ := inv_pos.mpr (Nat.cast_pos.mpr (Nat.factorial_pos _))
  simp only [saw_eq_bernoulliFun_fract, abs_mul, abs_of_pos sp]
  refine mul_le_mul_of_nonneg_left ?_ sp.le
  exact abs_bernoulliFun_le _ _ (unitInterval.fract_mem x)

/-- A uniform bound for the absolute value of `saw s`. -/
def sawBound (s : ℕ) : ℝ :=
  sSup (range fun x ↦ |saw s x|)

private lemma abs_saw_le (s : ℕ) (x : ℝ) : |saw s x| ≤ sawBound s := by
  exact le_csSup (bddAbove_range_abs_saw s) (mem_range_self x)

/-- The uniform bound for `saw` is nonnegative. -/
@[simp] lemma sawBound_nonneg {s : ℕ} : 0 ≤ sawBound s := by
  exact (abs_nonneg (saw s 0)).trans (abs_saw_le s 0)

/-! ## Euler--Maclaurin on one unit interval -/

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
variable {f : ℝ → E} {t : Set ℝ} {b c : ℤ} {n : ℕ}

private lemma intervalIntegrable_presaw_smul
    (fc : ContDiffOn ℝ n f t) {a b : ℝ} (ab : a ≤ b)
    (u : UniqueDiffOn ℝ t) (abt : Icc a b ⊆ t) {s : ℕ} :
    IntervalIntegrable (fun x ↦ presaw s c x • iteratedDerivWithin n f t x) volume a b := by
  refine (contDiff_presaw.continuous.continuousOn.smul ?_).intervalIntegrable
  simp only [uIcc_of_le ab]
  exact (fc.continuousOn_iteratedDerivWithin le_rfl u).mono abt

private lemma integral_saw_eq_integral_presaw :
    ∫ x in a..a + 1, saw s x • iteratedDerivWithin s f t x =
      ∫ x in a..a + 1, presaw s a x • iteratedDerivWithin s f t x := by
  apply intervalIntegral.integral_congr_ae
  have e : uIoc (a : ℝ) (a + 1) =ᵐ[volume] Ioo a (a + 1) := by
    rw [uIoc_of_le]
    · exact MeasureTheory.Ioo_ae_eq_Ioc.symm
    · simp only [(lt_add_one _).le]
  simp only [← MeasureTheory.ae_restrict_iff' measurableSet_uIoc,
    MeasureTheory.ae_restrict_congr_set e,
    MeasureTheory.ae_restrict_iff' measurableSet_Ioo]
  filter_upwards
  intro x m
  rw [saw_eqOn (a := a)]
  exact Ioo_subset_Ico_self m

private lemma presaw_smul_iteratedDeriv_by_parts [CompleteSpace E]
    (fc : ContDiffOn ℝ (s + 1) f t) (u : UniqueDiffOn ℝ t)
    (abt : Icc (a : ℝ) (a + 1) ⊆ t) :
    ∫ x in a..a + 1, presaw s c x • iteratedDerivWithin s f t x =
      presaw (s + 1) c (a + 1) • iteratedDerivWithin s f t (a + 1) -
      presaw (s + 1) c a • iteratedDerivWithin s f t a -
      ∫ x in a..a + 1, presaw (s + 1) c x • iteratedDerivWithin (s + 1) f t x := by
  have i0 := intervalIntegrable_presaw_smul (s := s) (c := c) fc.of_succ (by linarith) u abt
  have i1 : IntervalIntegrable (fun x ↦ presaw (s + 1) c x •
      iteratedDerivWithin (s + 1) f t x) volume (a : ℝ) (a + 1) :=
    intervalIntegrable_presaw_smul (s := s + 1) (c := c) fc (by linarith) u abt
  rw [eq_sub_iff_add_eq, ← intervalIntegral.integral_add i0 i1]
  set g := fun x ↦ presaw (s + 1) c x • iteratedDerivWithin s f t x
  set g' := fun x ↦ presaw s c x • iteratedDerivWithin s f t x +
    presaw (s + 1) c x • iteratedDerivWithin (s + 1) f t x
  have df : ∀ x ∈ Ioo (a : ℝ) (a + 1) \ ∅, HasDerivAt g (g' x) x := by
    intro x m
    simp only [sdiff_empty] at m
    simp only [g, g', add_comm (presaw s c _ • _) _]
    apply HasDerivAt.smul
    · exact hasDerivAt_presaw
    · have mt := abt (Ioo_subset_Icc_self m)
      apply HasDerivWithinAt.hasDerivAt
      · rw [iteratedDerivWithin_succ, hasDerivWithinAt_derivWithin_iff]
        apply (fc.contDiffWithinAt mt).differentiableWithinAt_iteratedDerivWithin
        · simp only [← Nat.cast_add_one, Nat.cast_lt, Nat.lt_add_one]
        · simp only [mt, insert_eq_of_mem, u]
      · exact Filter.monotone_mem (subset_trans Ioo_subset_Icc_self abt) (isOpen_Ioo.mem_nhds m)
  refine Eq.trans (MeasureTheory.integral_eq_of_hasDerivAt_off_countable_of_le
    (f := g) (f' := g') (Hd := df) (by linarith) countable_empty ?_ ?_) ?_
  · apply contDiff_presaw.continuous.continuousOn.smul
    exact (fc.continuousOn_iteratedDerivWithin le_self_add u).mono abt
  · exact i0.add i1
  · simp only [g]

private lemma presaw_iterated_by_parts [CompleteSpace E] (fc : ContDiffOn ℝ (s + 1) f t)
    (u : UniqueDiffOn ℝ t) (abt : Icc (a : ℝ) (a + 1) ⊆ t) :
    ∫ x in a..a + 1, f x = (2⁻¹ : ℝ) • (f a + f (a + 1)) -
      ∑ m ∈ Finset.range s, (-1 : ℝ) ^ m • saw (m + 2) 0 •
        (iteratedDerivWithin (m + 1) f t (a + 1) - iteratedDerivWithin (m + 1) f t a) -
      (-1 : ℝ) ^ s •
        ∫ x in a..a + 1, presaw (s + 1) a x • iteratedDerivWithin (s + 1) f t x := by
  induction s
  · simpa only [zero_add, presaw_one_coe_succ, presaw_coe_same, saw_one_zero, neg_smul,
      sub_neg_eq_add, one_div, iteratedDerivWithin_zero, presaw_zero, smul_add,
      Finset.range_zero, Finset.sum_empty, sub_zero, ← smul_add, Int.reduceNeg, pow_zero,
      one_smul, Int.cast_add, Int.cast_one, add_comm _ (f a)] using
      presaw_smul_iteratedDeriv_by_parts fc u abt (c := a)
  · rename_i s h
    refine (h fc.of_succ).trans ?_
    clear h
    have z : (-1 : ℝ) ^ s ≠ 0 := by simp
    simp only [sub_sub, add_left_cancel_iff, sub_right_inj, Finset.sum_range_succ, add_assoc,
      ← smul_smul, pow_succ, ← smul_add, neg_smul, ← sub_eq_add_neg, one_smul,
      smul_right_inj z]
    refine (presaw_smul_iteratedDeriv_by_parts fc u abt (c := a)).trans ?_
    simp only [(by omega : s + 1 + 1 = s + 2), presaw_coe_same, Nat.reduceAdd,
      presaw_coe_succ, ← smul_sub]

/-! ## The full formula -/

/-- The trapezoidal rule for integer step size on `[a, a + n]`. -/
private def trapezoid_sum (f : ℝ → E) (a : ℤ) : ℕ → E
  | 0 => 0
  | n + 1 => trapezoid_sum f a n + (2⁻¹ : ℝ) • (f (a + n) + f (a + n + 1))

@[simp] private lemma trapezoid_sum_zero : trapezoid_sum f a 0 = 0 := rfl

private lemma trapezoid_sum_succ :
    trapezoid_sum f a (n + 1) =
      trapezoid_sum f a n + (2⁻¹ : ℝ) • (f (a + n) + f (a + n + 1)) := rfl

/-- The trapezoidal sum is the ordinary inclusive sum with half of each endpoint removed. -/
private lemma trapezoid_sum_eq_sum_Icc_sub (f : ℝ → E) (a : ℤ) (n : ℕ) :
    trapezoid_sum f a n =
      ∑ k ∈ Finset.Icc a (a + n), f k - (2⁻¹ : ℝ) • (f a + f (a + n)) := by
  induction n with
  | zero =>
      rw [trapezoid_sum_zero]
      simp only [Nat.cast_zero, add_zero, Finset.Icc_self,
        Finset.sum_singleton]
      module
  | succ n ih =>
      simp only [Nat.cast_add, Nat.cast_one] at ih ⊢
      rw [trapezoid_sum_succ, ih]
      have hinsert :
          Finset.Icc a (a + ((n : ℤ) + 1)) =
            insert (a + ((n : ℤ) + 1)) (Finset.Icc a (a + n)) := by
        have hle : a ≤ Order.succ (a + (n : ℤ)) := by
          rw [Order.succ_eq_add_one]
          omega
        have h := Finset.insert_Icc_right_eq_Icc_succ hle
        rw [Order.succ_eq_add_one] at h
        simpa only [add_assoc] using h.symm
      rw [hinsert, Finset.sum_insert]
      · push_cast
        simp only [add_assoc]
        module
      · simp

private lemma intervalIntegrable_saw_smul
    (fc : ContDiffOn ℝ s f t) (u : UniqueDiffOn ℝ t)
    (abt : Icc (a : ℝ) (a + n) ⊆ t) :
    IntervalIntegrable (fun x ↦ saw s x • iteratedDerivWithin s f t x) volume a (a + n) := by
  induction n
  · simp
  · rename_i n h
    simp only [Nat.cast_add, Nat.cast_one, ← add_assoc]
    refine (h (subset_trans (Icc_subset_Icc (by rfl) (by simp)) abt)).trans ?_
    have ab1 : Icc ((a + n : ℤ) : ℝ) ((a + n : ℤ) + 1) ⊆ t := by
      refine subset_trans (Icc_subset_Icc ?_ ?_) abt
      all_goals simp [← add_assoc]
    have i := intervalIntegrable_presaw_smul (s := s) (c := a + n) fc (by simp) u ab1
    simp only [Int.cast_add, Int.cast_natCast, intervalIntegrable_iff,
      le_add_iff_nonneg_right, zero_le_one, uIoc_of_le] at i ⊢
    apply i.congr
    simp only [MeasureTheory.ae_restrict_iff' measurableSet_Ioc, Filter.EventuallyEq]
    filter_upwards [MeasureTheory.Measure.ae_ne volume (a + n + 1 : ℝ)]
    intro x ne m
    rw [saw_eqOn (a := a + n)]
    simp only [mem_Ioc, ne_eq, Int.cast_add, Int.cast_natCast, mem_Ico] at m ne ⊢
    exact ⟨m.1.le, Ne.lt_of_le ne m.2⟩

/-- The arbitrary-order Euler--Maclaurin formula. -/
private lemma trapezoid_sum_eq_integral_add [CompleteSpace E]
    (fc : ContDiffOn ℝ (s + 1) f t) (u : UniqueDiffOn ℝ t)
    (abt : Icc (a : ℝ) (a + n) ⊆ t) :
    trapezoid_sum f a n = (∫ x in a..a + n, f x) +
      ∑ m ∈ Finset.range s, (-1 : ℝ) ^ m • saw (m + 2) 0 •
        (iteratedDerivWithin (m + 1) f t (a + n) - iteratedDerivWithin (m + 1) f t a) +
      (-1 : ℝ) ^ s •
        ∫ x in a..a + n, saw (s + 1) x • iteratedDerivWithin (s + 1) f t x := by
  induction n
  · simp only [trapezoid_sum_zero, CharP.cast_eq_zero, add_zero, integral_same, sub_self,
      smul_zero, Finset.sum_const_zero]
  · rename_i n h
    generalize hb : a + n = b
    have hb' : (a + n : ℝ) = b := by simp [← hb]
    have ab0 : Icc (a : ℝ) (a + n) ⊆ t :=
      subset_trans (Icc_subset_Icc (by rfl) (by simp)) abt
    have ab1 : Icc (a + n : ℝ) (a + n + 1) ⊆ t :=
      subset_trans (Icc_subset_Icc (by simp) (by simp [add_assoc])) abt
    simp only [hb', Nat.cast_add, Nat.cast_one, ← add_assoc] at ab0 ab1 h ⊢
    simp only [trapezoid_sum_succ]
    rw [← intervalIntegral.integral_add_adjacent_intervals (b := a + n),
      ← intervalIntegral.integral_add_adjacent_intervals (a := a) (b := b) (c := b + 1),
      h ab0]
    · simp only [presaw_iterated_by_parts (a := b) fc u ab1,
        ← integral_saw_eq_integral_presaw, smul_sub, smul_add, Finset.sum_sub_distrib, hb']
      abel
    · rw [← hb']
      exact intervalIntegrable_saw_smul (by simpa) u (by simpa only [hb'])
    · rw [← Nat.cast_one (R := ℝ)]
      exact intervalIntegrable_saw_smul (by simpa) u (by simpa only [Nat.cast_one])
    · exact (fc.continuousOn.mono (by simpa only [hb'])).intervalIntegrable_of_Icc (by linarith)
    · exact (fc.continuousOn.mono (by simpa only [hb'])).intervalIntegrable_of_Icc (by linarith)

/-- The arbitrary-order Euler--Maclaurin formula for a sum over consecutive integers. -/
theorem sum_Icc_eq_integral_add [CompleteSpace E]
    (fc : ContDiffOn ℝ (s + 1) f t) (u : UniqueDiffOn ℝ t)
    (abt : Icc (a : ℝ) (a + n) ⊆ t) :
    ∑ k ∈ Finset.Icc a (a + n), f k =
      (∫ x in a..a + n, f x) +
      (2⁻¹ : ℝ) • (f a + f (a + n)) +
      ∑ m ∈ Finset.range s, (-1 : ℝ) ^ m • saw (m + 2) 0 •
        (iteratedDerivWithin (m + 1) f t (a + n) - iteratedDerivWithin (m + 1) f t a) +
      (-1 : ℝ) ^ s •
        ∫ x in a..a + n, saw (s + 1) x • iteratedDerivWithin (s + 1) f t x := by
  have h := trapezoid_sum_eq_integral_add fc u abt
  rw [trapezoid_sum_eq_sum_Icc_sub, sub_eq_iff_eq_add] at h
  rw [h]
  module

private lemma norm_integral_saw_smul_iteratedDerivWithin_le
    {C : ℝ}
    (hderiv : ∀ x ∈ Icc (a : ℝ) (a + n),
      ‖iteratedDerivWithin (s + 1) f t x‖ ≤ C) :
    ‖∫ x in a..a + n, saw (s + 1) x • iteratedDerivWithin (s + 1) f t x‖ ≤
      sawBound (s + 1) * C * n := by
  refine (intervalIntegral.norm_integral_le_of_norm_le_const
    (C := sawBound (s + 1) * C) ?_).trans ?_
  · intro x hx
    have hxIcc : x ∈ Icc (a : ℝ) (a + n) := by
      rw [uIoc_of_le (by simp), mem_Ioc] at hx
      exact ⟨hx.1.le, hx.2⟩
    rw [norm_smul, Real.norm_eq_abs]
    exact mul_le_mul (abs_saw_le _ _) (hderiv x hxIcc) (norm_nonneg _)
      sawBound_nonneg
  · simp only [add_sub_cancel_left]
    rw [abs_of_nonneg (Nat.cast_nonneg n)]

private lemma norm_sum_Icc_int_sub_integral_le [CompleteSpace E]
    (fc : ContDiffOn ℝ (s + 1) f t) (u : UniqueDiffOn ℝ t)
    (abt : Icc (a : ℝ) (a + n) ⊆ t)
    {C₀ : ℝ} (ha : ‖f a‖ ≤ C₀) (hb : ‖f (a + n)‖ ≤ C₀)
    (C : ℕ → ℝ)
    (hderiv : ∀ k, 1 ≤ k → k ≤ s + 1 →
      ∀ x ∈ Icc (a : ℝ) (a + n), ‖iteratedDerivWithin k f t x‖ ≤ C k) :
    ‖(∑ k ∈ Finset.Icc a (a + n), f k) - ∫ x in a..a + n, f x‖ ≤ C₀ +
        (∑ m ∈ Finset.range s, |saw (m + 2) 0| * (2 * C (m + 1))) +
        sawBound (s + 1) * C (s + 1) * n := by
  let endpoint : E := (2⁻¹ : ℝ) • (f a + f (a + n))
  let correction : E := ∑ m ∈ Finset.range s, (-1 : ℝ) ^ m • saw (m + 2) 0 •
    (iteratedDerivWithin (m + 1) f t (a + n) - iteratedDerivWithin (m + 1) f t a)
  let remainder : E := (-1 : ℝ) ^ s • ∫ x in a..a + n,
    saw (s + 1) x • iteratedDerivWithin (s + 1) f t x
  have hformula : (∑ k ∈ Finset.Icc a (a + n), f k) = (∫ x in a..a + n, f x) +
      endpoint + correction + remainder := by
    simpa only [endpoint, correction, remainder] using sum_Icc_eq_integral_add fc u abt
  have hendpoint : ‖endpoint‖ ≤ C₀ := by
    calc
      ‖endpoint‖ ≤ |(2⁻¹ : ℝ)| * (‖f a‖ + ‖f (a + n)‖) := by
        simpa only [endpoint, norm_smul, Real.norm_eq_abs] using
          mul_le_mul_of_nonneg_left (norm_add_le (f a) (f (a + n))) (abs_nonneg _)
      _ ≤ |(2⁻¹ : ℝ)| * (C₀ + C₀) := by gcongr
      _ = C₀ := by ring_nf
  have hcorrection : ‖correction‖ ≤
      ∑ m ∈ Finset.range s, |saw (m + 2) 0| * (2 * C (m + 1)) := by
    calc
      ‖correction‖ ≤ ∑ m ∈ Finset.range s,
          ‖(-1 : ℝ) ^ m • saw (m + 2) 0 •
            (iteratedDerivWithin (m + 1) f t (a + n) -
              iteratedDerivWithin (m + 1) f t a)‖ := norm_sum_le _ _
      _ ≤ ∑ m ∈ Finset.range s, |saw (m + 2) 0| * (2 * C (m + 1)) := by
        apply Finset.sum_le_sum
        intro m hm
        have hmS : m < s := Finset.mem_range.mp hm
        rw [norm_smul, norm_smul, Real.norm_eq_abs, Real.norm_eq_abs,
          abs_pow, abs_neg, abs_one, one_pow, one_mul]
        calc
          |saw (m + 2) 0| * ‖_ - _‖ ≤
              |saw (m + 2) 0| * (C (m + 1) + C (m + 1)) := by
            gcongr
            exact (norm_sub_le _ _).trans (add_le_add
              (hderiv (m + 1) (by omega) (by omega) (a + n) (by simp))
              (hderiv (m + 1) (by omega) (by omega) a (by simp)))
          _ = |saw (m + 2) 0| * (2 * C (m + 1)) := by ring
  have hremainder : ‖remainder‖ ≤ sawBound (s + 1) * C (s + 1) * n := by
    rw [show ‖remainder‖ =
        ‖∫ x in a..a + n, saw (s + 1) x • iteratedDerivWithin (s + 1) f t x‖ by
      dsimp only [remainder]
      rw [norm_smul, Real.norm_eq_abs, abs_pow, abs_neg, abs_one, one_pow, one_mul]]
    exact norm_integral_saw_smul_iteratedDerivWithin_le
      (fun x hx ↦ hderiv (s + 1) (by omega) le_rfl x hx)
  rw [hformula]
  rw [show ((∫ x in a..a + n, f x) + endpoint + correction + remainder) -
      (∫ x in a..a + n, f x) = endpoint + correction + remainder by abel]
  exact (norm_add_le _ _).trans (add_le_add (norm_add_le _ _) le_rfl) |>.trans <|
    add_le_add (add_le_add hendpoint hcorrection) hremainder

private lemma map_Icc_natCast_int (a b : ℕ) :
    (Finset.Icc a b).map
      ⟨(fun n : ℕ ↦ (n : ℤ)), fun _ _ h ↦ Int.ofNat_inj.mp h⟩ =
      Finset.Icc (a : ℤ) (b : ℤ) := by
  ext z
  simp only [Finset.mem_map, Finset.mem_Icc]
  constructor
  · rintro ⟨n, hn, rfl⟩
    change (a : ℤ) ≤ (n : ℤ) ∧ (n : ℤ) ≤ (b : ℤ)
    exact_mod_cast hn
  · intro hz
    have hz0 : 0 ≤ z := le_trans (Int.natCast_nonneg a) hz.1
    have hzt : (z.toNat : ℤ) = z := Int.toNat_of_nonneg hz0
    refine ⟨z.toNat, ?_, hzt⟩
    constructor
    · exact_mod_cast (show (a : ℤ) ≤ (z.toNat : ℤ) by simpa [hzt] using hz.1)
    · exact_mod_cast (show (z.toNat : ℤ) ≤ (b : ℤ) by simpa [hzt] using hz.2)

omit [NormedSpace ℝ E] in
private lemma sum_Icc_int_eq_sum_Icc_nat
    (f : ℝ → E) (a b : ℕ) :
    ∑ z ∈ Finset.Icc (a : ℤ) (b : ℤ), f z =
      ∑ n ∈ Finset.Icc a b, f n := by
  rw [← map_Icc_natCast_int]
  simp

/-- A convenient norm form of Euler--Maclaurin for natural endpoints. Uniform bounds on the
function and its derivatives control the difference between the inclusive sum and its integral;
the right-hand side records the endpoint, correction-term, and remainder contributions. -/
theorem norm_sum_Icc_nat_sub_integral_le [CompleteSpace E]
    {a b s : ℕ} (hab : a ≤ b)
    (fc : ContDiffOn ℝ (s + 1) f t) (u : UniqueDiffOn ℝ t)
    (abt : Icc (a : ℝ) (b : ℝ) ⊆ t)
    {C₀ : ℝ} (hendpoint : ∀ x ∈ Icc (a : ℝ) (b : ℝ), ‖f x‖ ≤ C₀)
    (C : ℕ → ℝ)
    (hderiv : ∀ k, 1 ≤ k → k ≤ s + 1 →
      ∀ x ∈ Icc (a : ℝ) (b : ℝ), ‖iteratedDerivWithin k f t x‖ ≤ C k) :
    ‖(∑ n ∈ Finset.Icc a b, f n) - ∫ x in (a : ℝ)..b, f x‖ ≤
      C₀ +
        (∑ m ∈ Finset.range s, |saw (m + 2) 0| * (2 * C (m + 1))) +
        sawBound (s + 1) * C (s + 1) * (b - a : ℕ) := by
  have hendNat : a + (b - a) = b := Nat.add_sub_of_le hab
  have hendInt : (a : ℤ) + (b - a : ℕ) = (b : ℤ) := by exact_mod_cast hendNat
  have hendReal : ((a : ℤ) : ℝ) + (b - a : ℕ) = ((b : ℤ) : ℝ) := by
    exact_mod_cast hendNat
  have hendRealNat : (a : ℝ) + (b - a : ℕ) = (b : ℝ) := by
    exact_mod_cast hendNat
  have hraw := norm_sum_Icc_int_sub_integral_le
    (s := s) (a := (a : ℤ)) (n := b - a) (f := f) (t := t)
    fc u (by simpa only [Int.cast_natCast, hendRealNat] using abt)
    (hendpoint a ⟨le_rfl, by exact_mod_cast hab⟩)
    (by simpa only [Int.cast_natCast, hendRealNat] using
      hendpoint b ⟨by exact_mod_cast hab, le_rfl⟩)
    C (fun k hk hks x hx ↦ hderiv k hk hks x (by
      simpa only [Int.cast_natCast, hendRealNat] using hx))
  rw [hendInt, hendReal, sum_Icc_int_eq_sum_Icc_nat] at hraw
  simpa only [Int.cast_natCast] using hraw

end
end
end Expdb.EulerMaclaurin
