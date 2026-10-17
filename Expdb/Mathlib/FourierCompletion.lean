module

public import Mathlib.Analysis.Fourier.ZMod
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Bounds
import Mathlib.NumberTheory.Harmonic.Bounds
import Mathlib.Algebra.Ring.GeomSum
import Mathlib.Order.Interval.Finset.Nat

/-!
# Finite Fourier completion of intervals

An interval indicator has Fourier coefficients with logarithmically bounded total mass.
Consequently a partial sum is controlled by a full sum with an additive-character twist.
-/

@[expose] public section
open scoped FourierTransform
noncomputable section
namespace Expdb

private theorem norm_fourierChar_sub_one (x : ℝ) :
    ‖(𝐞 x : ℂ) - 1‖ = 2 * |Real.sin (Real.pi * x)| := by
  rw [Real.fourierChar_apply]
  have h := Complex.norm_exp_I_mul_ofReal_sub_one (2 * Real.pi * x)
  simpa only [Complex.ofReal_mul, Complex.ofReal_ofNat, mul_comm Complex.I,
    show 2 * Real.pi * x / 2 = Real.pi * x by ring, Real.norm_eq_abs,
    abs_mul, abs_of_pos (by norm_num : (0 : ℝ) < 2)] using h

private theorem fourierChar_sub_one_lower {x : ℝ} (hx : 0 ≤ x) (hx1 : x ≤ 1) :
    4 * min x (1 - x) ≤ ‖(𝐞 x : ℂ) - 1‖ := by
  rw [norm_fourierChar_sub_one]
  have hs := Real.sin_nonneg_of_nonneg_of_le_pi
    (mul_nonneg Real.pi_pos.le hx) (by nlinarith [Real.pi_pos])
  rw [abs_of_nonneg hs]
  by_cases hh : x ≤ 1 / 2
  · rw [min_eq_left (by linarith)]
    have h := Real.mul_le_sin (x := Real.pi * x)
      (by positivity) (by nlinarith [Real.pi_pos])
    have he : 2 / Real.pi * (Real.pi * x) = 2 * x := by field_simp
    rw [he] at h
    linarith
  · rw [min_eq_right (by linarith)]
    have h := Real.mul_le_sin (x := Real.pi * (1 - x))
      (by positivity) (by nlinarith [Real.pi_pos])
    have he : 2 / Real.pi * (Real.pi * (1 - x)) = 2 * (1 - x) := by field_simp
    rw [he, show Real.pi * (1 - x) = Real.pi - Real.pi * x by ring, Real.sin_pi_sub] at h
    linarith

private theorem norm_geometric_interval_le (z : ℂ) (hz : ‖z‖ = 1) (a b : ℕ) :
    ‖∑ n ∈ Finset.Icc a b, z ^ n‖ * ‖z - 1‖ ≤ 2 := by
  by_cases hab : a ≤ b
  · rw [← norm_mul, ← Finset.Ico_add_one_right_eq_Icc,
      geom_sum_Ico_mul z (by omega)]
    calc
      ‖z ^ (b + 1) - z ^ a‖ ≤ ‖z ^ (b + 1)‖ + ‖z ^ a‖ := norm_sub_le _ _
      _ = 2 := by simp only [norm_pow, hz, one_pow]; norm_num
  · simp [Finset.Icc_eq_empty hab]

private def intervalCoefficient (q a b j : ℕ) : ℂ :=
  (q : ℂ)⁻¹ * ∑ n ∈ Finset.Icc a b, (𝐞 ((j : ℝ) * n / q) : ℂ)

private theorem coefficient_zero_le {q a b : ℕ} (hq : 0 < q) (hb : b < q) :
    ‖intervalCoefficient q a b 0‖ ≤ 1 := by
  simp only [intervalCoefficient, Nat.cast_zero, zero_mul, zero_div,
    Finset.sum_const, nsmul_eq_mul,
    norm_mul, norm_inv, Complex.norm_natCast]
  apply (inv_mul_le_iff₀ (by exact_mod_cast hq : (0 : ℝ) < q)).2
  have hc : (Finset.Icc a b).card ≤ q := by rw [Nat.card_Icc]; omega
  simpa using (show ((Finset.Icc a b).card : ℝ) ≤ q by exact_mod_cast hc)

private theorem coefficient_nonzero_le {q a b j : ℕ} (hj : 0 < j) (hjq : j < q) :
    ‖intervalCoefficient q a b j‖ ≤ 1 / (j : ℝ) + 1 / ((q - j : ℕ) : ℝ) := by
  have hq : (0 : ℝ) < q := by exact_mod_cast (hj.trans hjq)
  have hjr : (0 : ℝ) < j := by exact_mod_cast hj
  have hdiff : (0 : ℝ) < (q : ℝ) - (j : ℝ) := by
    have : (j : ℝ) < q := by exact_mod_cast hjq
    linarith
  let x : ℝ := (j : ℝ) / q
  have hx : 0 < x := div_pos hjr hq
  have hx1 : x < 1 := (div_lt_one hq).2 (by exact_mod_cast hjq)
  have he (n : ℕ) : (𝐞 ((j : ℝ) * n / q) : ℂ) = (𝐞 x : ℂ) ^ n := by
    have hh := AddChar.map_nsmul_eq_pow Real.fourierChar n x
    have hn : n • x = (j : ℝ) * n / q := by simp only [nsmul_eq_mul, x]; ring
    rw [hn] at hh
    exact congrArg (fun z : Circle ↦ (z : ℂ)) hh
  have hgeom := norm_geometric_interval_le (𝐞 x : ℂ) (Circle.norm_coe _) a b
  have hlower := fourierChar_sub_one_lower hx.le hx1.le
  have hm : 0 < min x (1 - x) := lt_min hx (by linarith)
  have hbnd : ‖∑ n ∈ Finset.Icc a b, (𝐞 x : ℂ) ^ n‖ ≤ 1 / (2 * min x (1 - x)) := by
    apply (le_div_iff₀ (by positivity)).2
    nlinarith [norm_nonneg (∑ n ∈ Finset.Icc a b, (𝐞 x : ℂ) ^ n)]
  simp only [intervalCoefficient, he, norm_mul, norm_inv, Complex.norm_natCast]
  have hd : ((q - j : ℕ) : ℝ) = (q : ℝ) - j := Nat.cast_sub hjq.le
  rw [hd]
  calc
    (q : ℝ)⁻¹ * ‖∑ n ∈ Finset.Icc a b, (𝐞 x : ℂ) ^ n‖ ≤
        (q : ℝ)⁻¹ * (1 / (2 * min x (1 - x))) :=
      mul_le_mul_of_nonneg_left hbnd (inv_nonneg.mpr hq.le)
    _ ≤ 1 / (j : ℝ) + 1 / ((q : ℝ) - j) := by
      by_cases hh : x ≤ 1 - x
      · rw [min_eq_left hh]
        dsimp [x]
        field_simp
        nlinarith
      · rw [min_eq_right (le_of_not_ge hh)]
        dsimp [x]
        field_simp
        nlinarith

private theorem coefficient_mass_le {q a b : ℕ} (hq : 0 < q) (hb : b < q) :
    ∑ j ∈ Finset.range q, ‖intervalCoefficient q a b j‖ ≤ 3 + 2 * Real.log q := by
  have hsplit : Finset.range q = insert 0 (Finset.Ico 1 q) := by
    ext j
    simp only [Finset.mem_range, Finset.mem_insert, Finset.mem_Ico]
    omega
  rw [hsplit, Finset.sum_insert (by simp)]
  have hsum := Finset.sum_le_sum (fun (j : ℕ) (hj : j ∈ Finset.Ico 1 q) ↦ coefficient_nonzero_le
    (by have := (Finset.mem_Ico.mp hj).1; omega) (Finset.mem_Ico.mp hj).2 (a := a) (b := b))
  have hreverse : ∑ j ∈ Finset.Ico 1 q, 1 / ((q - j : ℕ) : ℝ) =
      ∑ j ∈ Finset.Ico 1 q, 1 / (j : ℝ) := by
    simpa only [Nat.add_sub_cancel, Nat.add_sub_cancel_left] using
      (Finset.sum_Ico_reflect (fun j ↦ 1 / (j : ℝ)) 1 (n := q) (m := q) (by omega))
  rw [Finset.sum_add_distrib, hreverse] at hsum
  have hharm : ∑ j ∈ Finset.Ico 1 q, 1 / (j : ℝ) ≤ 1 + Real.log q := by
    have h := harmonic_le_one_add_log q
    rw [harmonic_eq_sum_Icc, Rat.cast_sum] at h
    simp only [Rat.cast_inv, Rat.cast_natCast, one_div] at h ⊢
    apply le_trans ?_ h
    apply Finset.sum_le_sum_of_subset_of_nonneg
    · intro j hj
      simp only [Finset.mem_Ico, Finset.mem_Icc] at hj ⊢
      omega
    · intro j _ _
      positivity
  linarith [coefficient_zero_le (a := a) hq hb]

private theorem stdAddChar_mul_nat {q : ℕ} [NeZero q] (h : ZMod q) (n : ℕ) :
    ZMod.stdAddChar (h * (n : ZMod q)) = (𝐞 ((h.val : ℝ) * n / q) : ℂ) := by
  rw [← h.natCast_zmod_val, ← Nat.cast_mul, ZMod.stdAddChar_apply,
    ZMod.toCircle_natCast, Real.fourierChar_apply]
  push_cast
  rw [ZMod.val_natCast_of_lt h.val_lt]
  congr 1
  ring

private theorem sum_character_eq {q m n : ℕ} [NeZero q] (hm : m < q) (hn : n < q) :
    ∑ j ∈ Finset.range q,
      (𝐞 ((j : ℝ) * m / q) : ℂ) * (𝐞 (-((j : ℝ) * n / q)) : ℂ) =
      if m = n then (q : ℂ) else 0 := by
  have he (h : ZMod q) : ZMod.stdAddChar (h * ((m : ZMod q) - (n : ZMod q))) =
      (𝐞 ((h.val : ℝ) * m / q) : ℂ) * (𝐞 (-((h.val : ℝ) * n / q)) : ℂ) := by
    rw [mul_sub, AddChar.map_sub_eq_div, stdAddChar_mul_nat, stdAddChar_mul_nat]
    rw [← Circle.coe_div, ← AddChar.map_sub_eq_div, ← Circle.coe_mul,
      ← AddChar.map_add_eq_mul]
    congr 2
  have hsum := AddChar.sum_mulShift ((m : ZMod q) - (n : ZMod q)) (ZMod.isPrimitive_stdAddChar q)
  have hzero : (m : ZMod q) - (n : ZMod q) = 0 ↔ m = n := by
    rw [sub_eq_zero]
    constructor
    · intro h
      have hv := congrArg ZMod.val h
      simpa [ZMod.val_natCast_of_lt hm, ZMod.val_natCast_of_lt hn] using hv
    · rintro rfl
      rfl
  simp only [hzero, ZMod.card, Nat.cast_ite, Nat.cast_zero] at hsum
  rw [← hsum]
  rw [← Fin.sum_univ_eq_sum_range]
  symm
  apply Fintype.sum_equiv (ZMod.finEquiv q).toEquiv.symm
  intro h
  rw [he]
  cases q with
  | zero => exact (NeZero.ne 0 rfl).elim
  | succ q => rfl

private theorem interval_completion_identity {q A B a b : ℕ} [NeZero q]
    (hAa : A ≤ a) (hbB : b ≤ B) (hBq : B < q) (f : ℕ → ℂ) :
    ∑ n ∈ Finset.Icc a b, f n =
      ∑ j ∈ Finset.range q, intervalCoefficient q a b j *
        ∑ n ∈ Finset.Icc A B, f n * (𝐞 (-((j : ℝ) * n / q)) : ℂ) := by
  classical
  simp only [intervalCoefficient]
  simp_rw [Finset.mul_sum, Finset.sum_mul]
  rw [Finset.sum_comm]
  simp_rw [Finset.sum_comm (s := Finset.range q) (t := Finset.Icc a b)]
  have hterm (n : ℕ) (hn : n ∈ Finset.Icc A B) :
      ∑ m ∈ Finset.Icc a b, ∑ j ∈ Finset.range q,
        ((q : ℂ)⁻¹ * (𝐞 ((j : ℝ) * m / q) : ℂ)) *
          (f n * (𝐞 (-((j : ℝ) * n / q)) : ℂ)) =
      if n ∈ Finset.Icc a b then f n else 0 := by
    have hnq : n < q := (Finset.mem_Icc.mp hn).2.trans_lt hBq
    have ht (m : ℕ) (hm : m ∈ Finset.Icc a b) :
        ∑ j ∈ Finset.range q,
          ((q : ℂ)⁻¹ * (𝐞 ((j : ℝ) * m / q) : ℂ)) *
            (f n * (𝐞 (-((j : ℝ) * n / q)) : ℂ)) =
          if m = n then f n else 0 := by
      have hmq : m < q := (Finset.mem_Icc.mp hm).2.trans hbB |>.trans_lt hBq
      calc
        _ = (q : ℂ)⁻¹ * f n * (∑ j ∈ Finset.range q,
            (𝐞 ((j : ℝ) * m / q) : ℂ) * (𝐞 (-((j : ℝ) * n / q)) : ℂ)) := by
          rw [Finset.mul_sum]
          apply Finset.sum_congr rfl
          intro j _
          ring
        _ = _ := by
          rw [sum_character_eq hmq hnq]
          have hqne : (q : ℂ) ≠ 0 := by exact_mod_cast NeZero.ne q
          split_ifs
          · rw [mul_assoc, mul_comm (f n) (q : ℂ), ← mul_assoc,
              inv_mul_cancel₀ hqne, one_mul]
          · simp
    simp_rw [Finset.sum_congr rfl ht]
    simp
  rw [Finset.sum_congr rfl hterm]
  rw [← Finset.sum_filter]
  congr 1
  ext n
  simp only [Finset.mem_filter, Finset.mem_Icc]
  omega

/-- Fourier completion bounds a partial sum by a logarithmic factor times one twisted
full sum. No regularity assumptions on the summands are needed. -/
theorem exists_norm_sum_Icc_le_log_mul_twist {q A B a b : ℕ}
    (hq : 0 < q) (hAa : A ≤ a) (hbB : b ≤ B) (hBq : B < q) (f : ℕ → ℂ) :
    ∃ r : ℝ, 0 ≤ r ∧ r ≤ 1 ∧
      ‖∑ n ∈ Finset.Icc a b, f n‖ ≤ (3 + 2 * Real.log q) *
        ‖∑ n ∈ Finset.Icc A B, f n * (𝐞 (-(r * n)) : ℂ)‖ := by
  letI : NeZero q := ⟨hq.ne'⟩
  let S (j : ℕ) : ℝ := ‖∑ n ∈ Finset.Icc A B, f n * (𝐞 (-((j : ℝ) * n / q)) : ℂ)‖
  obtain ⟨j, hj, hmax⟩ := Finset.exists_max_image (Finset.range q) S (by simp; omega)
  refine ⟨(j : ℝ) / q, (div_nonneg (Nat.cast_nonneg _) (Nat.cast_nonneg _)), ?_, ?_⟩
  · apply (div_le_one (by exact_mod_cast hq : (0 : ℝ) < q)).2
    exact_mod_cast (Finset.mem_range.mp hj).le
  have hS : 0 ≤ S j := norm_nonneg _
  have he : ∑ n ∈ Finset.Icc A B, f n * (𝐞 (-((j : ℝ) / q * n)) : ℂ) =
      ∑ n ∈ Finset.Icc A B, f n * (𝐞 (-((j : ℝ) * n / q)) : ℂ) := by
    apply Finset.sum_congr rfl
    intro n _
    congr 3
    ring
  rw [he]
  calc
    _ = ‖∑ k ∈ Finset.range q, intervalCoefficient q a b k *
        ∑ n ∈ Finset.Icc A B, f n * (𝐞 (-((k : ℝ) * n / q)) : ℂ)‖ := by
      rw [interval_completion_identity hAa hbB hBq]
    _ ≤ ∑ k ∈ Finset.range q, ‖intervalCoefficient q a b k‖ * S k := by
      exact (norm_sum_le _ _).trans_eq (by simp [S])
    _ ≤ ∑ k ∈ Finset.range q, ‖intervalCoefficient q a b k‖ * S j := by
      apply Finset.sum_le_sum
      intro k hk
      exact mul_le_mul_of_nonneg_left (hmax k hk) (norm_nonneg _)
    _ = (∑ k ∈ Finset.range q, ‖intervalCoefficient q a b k‖) * S j :=
      (Finset.sum_mul _ _ _).symm
    _ ≤ (3 + 2 * Real.log q) * S j :=
      mul_le_mul_of_nonneg_right (coefficient_mass_le hq (hbB.trans_lt hBq)) hS

end Expdb
