module

public import Expdb.ExponentialSums.Oscillatory
public import Mathlib.Algebra.Order.Floor.Div
public import Mathlib.Analysis.Fourier.ZMod

import Mathlib.Data.Nat.Cast.Order.Field

/-!
# Finite scale transfer for exponential sums

Concrete finite Fourier and residue-class identities used to transfer exponential-sum estimates
between scales.
-/

@[expose] public section

open scoped Expdb FourierTransform

noncomputable section

namespace Expdb

/-! ## Finite Fourier orthogonality -/

private lemma sum_stdAddChar_mul_eq_ite (q m : ℕ) [NeZero q] :
    (∑ h : ZMod q, ZMod.stdAddChar (h * (m : ZMod q))) =
      if q ∣ m then (q : ℂ) else 0 := by
  rw [AddChar.sum_mulShift (m : ZMod q) (ZMod.isPrimitive_stdAddChar q)]
  simp only [ZMod.natCast_eq_zero_iff m q, ZMod.card]
  split_ifs <;> norm_num

private lemma stdAddChar_mul_nat_eq_fourierChar (q m : ℕ) [NeZero q] (h : ZMod q) :
    ZMod.stdAddChar (h * (m : ZMod q)) =
      (𝐞 ((h.val : ℝ) * (m : ℝ) / q) : ℂ) := by
  rw [← h.natCast_zmod_val, ← Nat.cast_mul, ZMod.stdAddChar_apply,
    ZMod.toCircle_natCast, Real.fourierChar_apply]
  push_cast
  rw [ZMod.val_natCast_of_lt h.val_lt]
  congr 1
  ring_nf

/-! ## Reindexing exponential sums -/

/-- Finite Fourier inversion expresses a sum at scale `N` as an average of sums dilated by `q`. -/
lemma exponentialSumAt_dilate (F : ℝ → ℝ) (T N : ℝ) (a b q : ℕ) [NeZero q]
    (hq : 0 < q) (hT : T ≠ 0) (hN : N ≠ 0) :
    exponentialSumAt F T N a b =
      (q : ℂ)⁻¹ * ∑ h : ZMod q,
        exponentialSumAt (fun u ↦ F u + (h.val : ℝ) * N / T * u)
          T (q * N) (q * a) (q * b) := by
  rw [exponentialSumAt]
  simp_rw [exponentialSumAt, oscillatory]
  rw [Finset.sum_comm]
  classical
  symm
  calc
    (q : ℂ)⁻¹ * ∑ m ∈ Finset.Icc (q * a) (q * b), ∑ h : ZMod q,
        (𝐞 (T * (F ((m : ℝ) / (q * N)) +
          (h.val : ℝ) * N / T * ((m : ℝ) / (q * N)))) : ℂ) =
        ∑ m ∈ Finset.Icc (q * a) (q * b),
          (q : ℂ)⁻¹ * ((𝐞 (T * F ((m : ℝ) / (q * N))) : ℂ) *
            ∑ h : ZMod q, ZMod.stdAddChar (h * (m : ZMod q))) := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro m hm
      apply congrArg (fun z : ℂ ↦ (q : ℂ)⁻¹ * z)
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro h _
      have harg :
          T * (F ((m : ℝ) / (q * N)) +
            (h.val : ℝ) * N / T * ((m : ℝ) / (q * N))) =
            T * F ((m : ℝ) / (q * N)) + (h.val : ℝ) * m / q := by
        field_simp [show (q : ℝ) ≠ 0 by exact_mod_cast hq.ne', hT, hN]
      rw [harg, AddChar.map_add_eq_mul, stdAddChar_mul_nat_eq_fourierChar]
      exact Circle.coe_mul _ _
    _ = ∑ m ∈ Finset.Icc (q * a) (q * b),
          if q ∣ m then (𝐞 (T * F ((m : ℝ) / (q * N))) : ℂ) else 0 := by
      apply Finset.sum_congr rfl
      intro m hm
      rw [sum_stdAddChar_mul_eq_ite]
      split_ifs with hdiv
      · field_simp [show (q : ℂ) ≠ 0 by exact_mod_cast hq.ne']
      · simp
    _ = ∑ n ∈ Finset.Icc a b, (𝐞 (T * F ((n : ℝ) / N)) : ℂ) := by
      symm
      rw [← Finset.sum_filter]
      apply Finset.sum_bij (fun n _ ↦ q * n)
      · intro n hn
        simp only [Finset.mem_Icc] at hn
        simp only [Finset.mem_filter, Finset.mem_Icc]
        exact ⟨⟨Nat.mul_le_mul_left q hn.1, Nat.mul_le_mul_left q hn.2⟩,
          dvd_mul_right q n⟩
      · intro n₁ hn₁ n₂ hn₂ heq
        exact Nat.eq_of_mul_eq_mul_left hq heq
      · intro m hm
        simp only [Finset.mem_filter, Finset.mem_Icc] at hm
        obtain ⟨n, rfl⟩ := hm.2
        refine ⟨n, ?_, rfl⟩
        simp only [Finset.mem_Icc]
        constructor
        · exact Nat.le_of_mul_le_mul_left hm.1.1 hq
        · exact Nat.le_of_mul_le_mul_left hm.1.2 hq
      · intro n hn
        push_cast
        congr 3
        field_simp [show (q : ℝ) ≠ 0 by exact_mod_cast hq.ne']

/-- Triangle-inequality form of `exponentialSumAt_dilate`. -/
lemma norm_exponentialSumAt_dilate_le (F : ℝ → ℝ) (T N : ℝ) (a b q : ℕ)
    [NeZero q]
    (hq : 0 < q) (hT : T ≠ 0) (hN : N ≠ 0) :
    ‖exponentialSumAt F T N a b‖ ≤
      ∑ h : ZMod q,
        ‖exponentialSumAt (fun u ↦ F u + (h.val : ℝ) * N / T * u)
          T (q * N) (q * a) (q * b)‖ / q := by
  rw [exponentialSumAt_dilate F T N a b q hq hT hN, norm_mul, norm_inv,
    Complex.norm_natCast]
  calc
    (q : ℝ)⁻¹ * ‖∑ h : ZMod q,
        exponentialSumAt (fun u ↦ F u + (h.val : ℝ) * N / T * u)
          T (q * N) (q * a) (q * b)‖ ≤
        (q : ℝ)⁻¹ * ∑ h : ZMod q,
          ‖exponentialSumAt (fun u ↦ F u + (h.val : ℝ) * N / T * u)
            T (q * N) (q * a) (q * b)‖ :=
      mul_le_mul_of_nonneg_left (norm_sum_le _ _) (by positivity)
    _ = _ := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro h _
      simp only [div_eq_mul_inv]
      ring

/-- The scale induced by restricting a sum to the residue class `r mod q`. -/
def residueScale (N : ℝ) (q r : ℕ) : ℝ := (N - r) / q

/-- The affine phase induced by restricting a sum to the residue class `r mod q`. -/
def residuePhase (F : ℝ → ℝ) (N : ℝ) (r : ℕ) : ℝ → ℝ :=
  fun u ↦ F ((1 - (r : ℝ) / N) * u + (r : ℝ) / N)

/-- Evaluation of the residue-class phase at its rescaled integer arguments. -/
lemma residuePhase_identity (F : ℝ → ℝ) (N : ℝ) (q r m : ℕ)
    (hN : N ≠ 0) (hrN : (r : ℝ) < N) :
    residuePhase F N r ((m : ℝ) / residueScale N q r) =
      F (((q * m + r : ℕ) : ℝ) / N) := by
  simp only [residuePhase, residueScale]
  push_cast
  field_simp [sub_ne_zero.mpr (ne_of_gt hrN)]

/-- Reindex one residue class as an exponential sum with `residuePhase` and `residueScale`. -/
lemma exponentialSumAt_residue (F : ℝ → ℝ) (T N : ℝ) (a b q r : ℕ)
    (hq : 0 < q) (hr : r < q) (hra : r ≤ a) (hab : a ≤ b)
    (hN : N ≠ 0) (hrN : (r : ℝ) < N) :
    (∑ n ∈ Finset.Icc a b with n % q = r, oscillatory F T N n) =
      exponentialSumAt (residuePhase F N r) T (residueScale N q r)
        ((a - r) ⌈/⌉ q) ((b - r) / q) := by
  rw [exponentialSumAt]
  classical
  symm
  apply Finset.sum_bij (fun m _ ↦ q * m + r)
  · intro m hm
    simp only [Finset.mem_Icc] at hm
    simp only [Finset.mem_filter, Finset.mem_Icc]
    constructor
    · constructor
      · have ha : a - r ≤ q * m :=
          (ceilDiv_le_iff_le_mul hq).mp hm.1
        omega
      · have hb : q * m ≤ b - r := by
          simpa [Nat.mul_comm] using (Nat.le_div_iff_mul_le hq).mp hm.2
        omega
    · simp [Nat.add_mod, Nat.mod_eq_of_lt hr]
  · intro m₁ hm₁ m₂ hm₂ heq
    exact Nat.eq_of_mul_eq_mul_left hq (Nat.add_right_cancel heq)
  · intro n hn
    simp only [Finset.mem_filter, Finset.mem_Icc] at hn
    have hnrepr : n = q * (n / q) + r := by
      have hnmod := hn.2
      have hdecomp := Nat.mod_add_div n q
      omega
    refine ⟨n / q, ?_, ?_⟩
    · simp only [Finset.mem_Icc]
      constructor
      · apply (ceilDiv_le_iff_le_mul hq).mpr
        omega
      · apply (Nat.le_div_iff_mul_le hq).mpr
        rw [Nat.mul_comm]
        omega
    · exact hnrepr.symm
  · intro m hm
    simp only [oscillatory]
    rw [residuePhase_identity F N q r m hN hrN]

private lemma residue_upper_endpoint_le (N : ℝ) (b q r : ℕ)
    (hq : 0 < q) (hr : r < q) (hb : (b : ℝ) ≤ 2 * N) :
    (b - r) / q ≤ ⌊2 * residueScale N q r⌋₊ + 1 := by
  by_cases hrb : r ≤ b
  · have hdiv : (((b - r) / q : ℕ) : ℝ) < (2 * residueScale N q r) + 1 := by
      have hnatdiv : (((b - r) / q : ℕ) : ℝ) ≤ ((b - r : ℕ) : ℝ) / q :=
        Nat.cast_div_le
      refine lt_of_le_of_lt hnatdiv ?_
      simp only [residueScale]
      rw [Nat.cast_sub hrb]
      have hqR : (0 : ℝ) < q := by exact_mod_cast hq
      have hrR : (r : ℝ) < q := by exact_mod_cast hr
      rw [div_lt_iff₀ hqR]
      field_simp [show (q : ℝ) ≠ 0 by positivity]
      nlinarith
    have hfloor := Nat.lt_floor_add_one (2 * residueScale N q r)
    have hcast : (((b - r) / q : ℕ) : ℝ) <
        ((⌊2 * residueScale N q r⌋₊ + 2 : ℕ) : ℝ) := by
      push_cast
      linarith
    have hnat : (b - r) / q < ⌊2 * residueScale N q r⌋₊ + 2 := by
      exact_mod_cast hcast
    omega
  · have hbr : b ≤ r := by omega
    simp [Nat.sub_eq_zero_of_le hbr]

/-- Split an exponential sum into residue classes, trimming each upper endpoint at a cost of at
most one term. -/
lemma norm_exponentialSumAt_le_sum_residues
    (F : ℝ → ℝ) (T N : ℝ) (a b q : ℕ)
    (hq : 0 < q) (hqN : (q : ℝ) ≤ N) (hN : N ≠ 0)
    (ha : N ≤ (a : ℝ)) (hb : (b : ℝ) ≤ 2 * N) :
    ‖exponentialSumAt F T N a b‖ ≤
      ∑ r ∈ Finset.range q,
        (‖exponentialSumAt (residuePhase F N r) T (residueScale N q r)
          ((a - r) ⌈/⌉ q)
          (min ((b - r) / q) ⌊2 * residueScale N q r⌋₊)‖ + 1) := by
  by_cases hab : a ≤ b
  swap
  · rw [exponentialSumAt, Finset.Icc_eq_empty hab]
    simp only [Finset.sum_empty, norm_zero]
    apply Finset.sum_nonneg
    intro r _
    exact add_nonneg (norm_nonneg (exponentialSumAt
      (residuePhase F N r) T (residueScale N q r)
      ((a - r) ⌈/⌉ q) (min ((b - r) / q) ⌊2 * residueScale N q r⌋₊))) zero_le_one
  rw [exponentialSumAt]
  classical
  have hpartition :
      ∑ n ∈ Finset.Icc a b, oscillatory F T N n =
        ∑ r ∈ Finset.range q,
          ∑ n ∈ Finset.Icc a b with n % q = r, oscillatory F T N n := by
    rw [Finset.sum_fiberwise_of_maps_to]
    intro n hn
    exact Finset.mem_range.mpr (Nat.mod_lt n hq)
  rw [hpartition]
  calc
    ‖∑ r ∈ Finset.range q, ∑ n ∈ Finset.Icc a b with n % q = r,
        oscillatory F T N n‖ ≤
        ∑ r ∈ Finset.range q,
          ‖∑ n ∈ Finset.Icc a b with n % q = r, oscillatory F T N n‖ :=
      norm_sum_le _ _
    _ ≤ _ := by
      apply Finset.sum_le_sum
      intro r hr
      have hrq := Finset.mem_range.mp hr
      have hrN : (r : ℝ) < N := lt_of_lt_of_le (by exact_mod_cast hrq) hqN
      have hra : r ≤ a := by
        have : (r : ℝ) < a := hrN.trans_le ha
        exact_mod_cast this.le
      rw [exponentialSumAt_residue F T N a b q r hq hrq hra hab hN hrN]
      let A := (a - r) ⌈/⌉ q
      let B := (b - r) / q
      let B' := min B ⌊2 * residueScale N q r⌋₊
      have hBB' : B ≤ B' + 1 := by
        dsimp [B, B']
        have hend := residue_upper_endpoint_le N b q r hq hrq hb
        rw [min_def]
        split_ifs <;> omega
      rw [exponentialSumAt]
      change ‖∑ n ∈ Finset.Icc A B,
          oscillatory (residuePhase F N r) T (residueScale N q r) n‖ ≤
        ‖∑ n ∈ Finset.Icc A B',
          oscillatory (residuePhase F N r) T (residueScale N q r) n‖ + 1
      by_cases hAB0 : A ≤ B
      swap
      · rw [Finset.Icc_eq_empty hAB0]
        simp only [Finset.sum_empty, norm_zero]
        positivity
      by_cases hle : B ≤ B'
      · have : B = B' := le_antisymm hle (min_le_left _ _)
        rw [this]
        exact le_add_of_nonneg_right (by norm_num)
      · have hsucc : B = B' + 1 := by omega
        by_cases hAB : A ≤ B'
        · rw [hsucc, Finset.sum_Icc_succ_top (hAB.trans (Nat.le_succ _))]
          calc
            _ ≤ ‖∑ n ∈ Finset.Icc A B',
                  oscillatory (residuePhase F N r) T (residueScale N q r) n‖ +
                ‖oscillatory (residuePhase F N r) T (residueScale N q r)
                  ((B' + 1 : ℕ) : ℝ)‖ :=
              norm_add_le _ _
            _ = _ := by
              rw [norm_oscillatory]
        · have hAeq : A = B' + 1 := by omega
          rw [hAeq, hsucc]
          simp [norm_oscillatory]

end Expdb
