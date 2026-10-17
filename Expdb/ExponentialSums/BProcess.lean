module

import Expdb.ExponentialSums.DualPhase
public import Expdb.ExponentialSums.ExponentSumGrowth
public import Expdb.ExponentialSums.ExponentialSumTransformation
public import Expdb.ExponentialSums.FixedExponentialSum
import Expdb.ExponentialSums.TrivialBounds
import Expdb.Mathlib.PartialSummation
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Analysis.Calculus.Deriv.Mul
import Mathlib.Analysis.Calculus.IteratedDeriv.Lemmas
public import Mathlib.Analysis.Complex.Basic

/-!
# The van der Corput B-process for growth exponents

This module proves the growth-bound transfer underlying the blueprint's reflection lemma.
-/

@[expose] public section

namespace Expdb.BProcess

open StationaryPhase

open Filter Topology Set MeasureTheory Complex PartialSummation
open scoped ContDiff FourierTransform NNReal

noncomputable section

section DualBlocks

/-- A sum over `[P 0, P J)` splits into the blocks `[P j, P (j + 1))`. -/
private theorem sum_Ico_eq_sum_blocks {E : Type*} [AddCommMonoid E] (f : ℕ → E) {P : ℕ → ℕ}
    (hP : Monotone P) (J : ℕ) :
    ∑ m ∈ Finset.Ico (P 0) (P J), f m =
      ∑ j ∈ Finset.range J, ∑ m ∈ Finset.Ico (P j) (P (j + 1)), f m := by
  induction J with
  | zero => simp
  | succ J ih =>
    rw [Finset.sum_range_succ, ← ih, Finset.sum_Ico_consecutive _ (hP (Nat.zero_le J))
      (hP (Nat.le_succ J))]

/-- Partial sums of the dual exponentials are conjugates of exponential sums with the
normalized dual phase. -/
private theorem norm_sum_Ico_fourierChar_legendre_eq (H : ℝ → ℝ) {T N κ S : ℝ} (hT : 0 < T)
    (hκ : 0 < κ) (hS : 0 < S) {P k : ℕ} (hP : 1 ≤ P) :
    ‖∑ m ∈ Finset.Ico P k, (𝐞 (T * legendre H (N * m / T)) : ℂ)‖ =
      ‖exponentialSumAt (dualPhase H κ S) (T * S) (κ * (T / N)) P (k - 1)‖ := by
  have hIco : Finset.Ico P k = Finset.Icc P (k - 1) := by
    ext m; simp only [Finset.mem_Ico, Finset.mem_Icc]; omega
  rw [exponentialSumAt, hIco, ← Complex.norm_conj, map_sum]
  congr 1
  refine Finset.sum_congr rfl fun m _ ↦ ?_
  rw [oscillatory]
  have hconj (x : ℝ) : (𝐞 (-x) : ℂ) = (starRingEnd ℂ) (𝐞 x : ℂ) := by
    rw [AddChar.map_neg_eq_inv, Circle.coe_inv_eq_conj]
  rw [← hconj]
  congr 2
  simp only [dualPhase]
  field_simp

section

variable {H : ℝ → ℝ} {c C : ℝ}

/-- The main amplitude is bounded. -/
private theorem abs_mainAmplitude_le (hc : 0 < c)
    (hcurv : ∀ x, -C ≤ iteratedDeriv 2 H x ∧ iteratedDeriv 2 H x ≤ -c) {W : ℝ}
    (hW : ∀ x, |phaseCutoff x| ≤ W) (ν : ℝ) :
    |mainAmplitude H ν| ≤ W / Real.sqrt c := by
  have hs : Real.sqrt c ≤ Real.sqrt (-iteratedDeriv 2 H (phaseInverse H ν)) :=
    Real.sqrt_le_sqrt (by linarith [(hcurv (phaseInverse H ν)).2])
  have hsc : 0 < Real.sqrt c := Real.sqrt_pos.2 hc
  rw [mainAmplitude, abs_div, abs_of_nonneg (Real.sqrt_nonneg _)]
  have hW0 : 0 ≤ W := (abs_nonneg _).trans (hW 0)
  calc |phaseCutoff (phaseInverse H ν)| / Real.sqrt (-iteratedDeriv 2 H (phaseInverse H ν))
      ≤ W / Real.sqrt (-iteratedDeriv 2 H (phaseInverse H ν)) :=
        div_le_div_of_nonneg_right (hW _) (Real.sqrt_nonneg _)
    _ ≤ W / Real.sqrt c := div_le_div_of_nonneg_left hW0 hsc hs

/-- The main amplitude is Lipschitz. -/
private theorem abs_mainAmplitude_sub_le (hH : ContDiff ℝ ∞ H) (hc : 0 < c)
    (hcurv : ∀ x, -C ≤ iteratedDeriv 2 H x ∧ iteratedDeriv 2 H x ≤ -c) {W K3 : ℝ}
    (hW : ∀ x, |phaseCutoff x| ≤ W) (hW1 : ∀ x, |deriv (phaseCutoff : ℝ → ℝ) x| ≤ W)
    (h3 : ∀ x, |iteratedDeriv 3 H x| ≤ K3) (ν ν' : ℝ) :
    |mainAmplitude H ν - mainAmplitude H ν'| ≤
      (W / Real.sqrt c + W * (K3 / (2 * c * Real.sqrt c))) / c * |ν - ν'| := by
  obtain ⟨-, -, hfg⟩ := phaseInverse_spec hH hc (fun x ↦ (hcurv x).2)
  have hH1 : ContDiff ℝ ∞ (deriv H) := hH.deriv' (n := ∞)
  have hH2 : ContDiff ℝ ∞ (iteratedDeriv 2 H) := by
    rw [iteratedDeriv_succ, iteratedDeriv_one]
    exact hH1.deriv' (n := ∞)
  have hd2 : ∀ x, deriv (deriv H) x ≤ -c := fun x ↦ by
    simpa [iteratedDeriv_succ, iteratedDeriv_one] using (hcurv x).2
  set x := phaseInverse H ν
  set x' := phaseInverse H ν'
  have hxx : |x - x'| ≤ |ν - ν'| / c :=
    abs_inverse_sub_le_of_neg (hH1.differentiable (by simp)) hc hd2 hfg ν ν'
  have hW0 : 0 ≤ W := (abs_nonneg _).trans (hW 0)
  have hK3 : 0 ≤ K3 := (abs_nonneg _).trans (h3 0)
  have hsc : 0 < Real.sqrt c := Real.sqrt_pos.2 hc
  -- Lipschitz bound for the cutoff
  have hwl : |phaseCutoff x - phaseCutoff x'| ≤ W * |x - x'| := by
    have := Convex.norm_image_sub_le_of_norm_deriv_le (f := (phaseCutoff : ℝ → ℝ)) (s := univ)
      (x := x') (y := x) (C := W)
      (fun z _ ↦ ((phaseCutoff.contDiff (n := 1)).differentiable (by simp)) z)
      (fun z _ ↦ by simpa [Real.norm_eq_abs] using hW1 z) convex_univ (mem_univ _) (mem_univ _)
    simpa [Real.norm_eq_abs] using this
  -- Lipschitz bound for the curvature
  have hsl : |iteratedDeriv 2 H x - iteratedDeriv 2 H x'| ≤ K3 * |x - x'| := by
    have := Convex.norm_image_sub_le_of_norm_deriv_le (f := iteratedDeriv 2 H) (s := univ)
      (x := x') (y := x) (C := K3)
      (fun z _ ↦ (hH2.differentiable (by simp)) z)
      (fun z _ ↦ by
        rw [← iteratedDeriv_succ]
        simpa [Real.norm_eq_abs] using h3 z) convex_univ (mem_univ _) (mem_univ _)
    simpa [Real.norm_eq_abs] using this
  set s := -iteratedDeriv 2 H x
  set s' := -iteratedDeriv 2 H x'
  have hs : c ≤ s := by linarith [(hcurv x).2]
  have hs' : c ≤ s' := by linarith [(hcurv x').2]
  have hrs : Real.sqrt c ≤ Real.sqrt s := Real.sqrt_le_sqrt hs
  have hrs' : Real.sqrt c ≤ Real.sqrt s' := Real.sqrt_le_sqrt hs'
  have hss : |s - s'| ≤ K3 * |x - x'| := by
    simp only [s, s']
    rw [show -iteratedDeriv 2 H x - -iteratedDeriv 2 H x' =
      -(iteratedDeriv 2 H x - iteratedDeriv 2 H x') by ring, abs_neg]
    exact hsl
  have hg : |1 / Real.sqrt s - 1 / Real.sqrt s'| ≤ K3 / (2 * c * Real.sqrt c) * |x - x'| := by
    calc _ ≤ |s - s'| / (2 * c * Real.sqrt c) := abs_inv_sqrt_sub_le hc hs hs'
      _ ≤ K3 * |x - x'| / (2 * c * Real.sqrt c) :=
        div_le_div_of_nonneg_right hss (by positivity)
      _ = _ := by ring
  have hgb : |1 / Real.sqrt s| ≤ 1 / Real.sqrt c := by
    rw [abs_of_pos (one_div_pos.2 (lt_of_lt_of_le hsc hrs))]
    exact div_le_div_of_nonneg_left zero_le_one hsc hrs
  have hwb : |phaseCutoff x'| ≤ W := hW x'
  have hA : mainAmplitude H ν - mainAmplitude H ν' =
      (phaseCutoff x - phaseCutoff x') * (1 / Real.sqrt s) +
        phaseCutoff x' * (1 / Real.sqrt s - 1 / Real.sqrt s') := by
    simp only [mainAmplitude, s, s', x, x']
    ring
  rw [hA]
  refine (abs_add_le _ _).trans ?_
  rw [abs_mul, abs_mul]
  have e1 : |phaseCutoff x - phaseCutoff x'| * |1 / Real.sqrt s| ≤
      W * |x - x'| * (1 / Real.sqrt c) :=
    mul_le_mul hwl hgb (abs_nonneg _) (by positivity)
  have e2 : |phaseCutoff x'| * |1 / Real.sqrt s - 1 / Real.sqrt s'| ≤
      W * (K3 / (2 * c * Real.sqrt c) * |x - x'|) :=
    mul_le_mul hwb hg (abs_nonneg _) hW0
  have hxx' : |x - x'| * c ≤ |ν - ν'| := by rwa [le_div_iff₀ hc] at hxx
  calc _ ≤ W * |x - x'| * (1 / Real.sqrt c) + W * (K3 / (2 * c * Real.sqrt c) * |x - x'|) :=
        add_le_add e1 e2
    _ = (W / Real.sqrt c + W * (K3 / (2 * c * Real.sqrt c))) * |x - x'| := by ring
    _ ≤ (W / Real.sqrt c + W * (K3 / (2 * c * Real.sqrt c))) * (|ν - ν'| / c) :=
        mul_le_mul_of_nonneg_left hxx (by positivity)
    _ = _ := by ring

/-- A block of main terms is bounded by the partial sums of the dual exponentials. -/
private theorem norm_sum_mainTerm_block_le {T N : ℝ} (hT : 0 < T) (hN : 0 < N) {Amax LA : ℝ}
    (hA : ∀ ν, |mainAmplitude H ν| ≤ Amax)
    (hLip : ∀ ν ν', |mainAmplitude H ν - mainAmplitude H ν'| ≤ LA * |ν - ν'|)
    {P Q : ℕ} (hPQ : P ≤ Q) {D : ℝ}
    (hD : ∀ k ∈ Finset.Icc P Q,
      ‖∑ m ∈ Finset.Ico P k, (𝐞 (T * legendre H (N * m / T)) : ℂ)‖ ≤ D) :
    ‖∑ m ∈ Finset.Ico P Q, mainTerm H T (N * m / T)‖ ≤
      (Amax + LA * ((Q : ℝ) - P) / (T / N)) * D / Real.sqrt T := by
  have hM : 0 < T / N := div_pos hT hN
  have hν : ∀ m : ℕ, N * (m : ℝ) / T = (m : ℝ) / (T / N) := fun m ↦ by field_simp
  have hsum : ∑ m ∈ Finset.Ico P Q, mainTerm H T (N * m / T) =
      (Complex.exp (-(Real.pi : ℂ) / 4 * I) / (Real.sqrt T : ℂ)) *
        ∑ m ∈ Finset.Ico P Q, ((mainAmplitude H ((m : ℝ) / (T / N)) : ℝ) : ℂ) *
          (𝐞 (T * legendre H (N * m / T)) : ℂ) := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun m _ ↦ ?_
    rw [mainTerm, hν m]
    ring
  rw [hsum, norm_mul, norm_div, Complex.norm_exp]
  have hre : (-(Real.pi : ℂ) / 4 * I).re = 0 := by simp
  rw [hre, Real.exp_zero, Complex.norm_real, Real.norm_of_nonneg (Real.sqrt_nonneg _)]
  have hb := norm_sum_Ico_lipschitz_mul_le (mainAmplitude H) hM hA hLip
    (fun m ↦ (𝐞 (T * legendre H (N * m / T)) : ℂ)) hPQ hD
  rw [one_div_mul_eq_div]
  exact div_le_div_of_nonneg_right hb (Real.sqrt_nonneg _)

end

end DualBlocks

section Assembly

/-- Block ratios: `κ 0 = 3 ^ (-σ)`, `κ J = 2 ^ σ`, and consecutive ratios at most `2`. -/
private theorem exists_block_constants {σ : ℝ} (hσ : 0 < σ) :
    ∃ (J : ℕ) (κ : ℕ → ℝ), κ 0 = (3 : ℝ) ^ (-σ) ∧ κ J = (2 : ℝ) ^ σ ∧ StrictMono κ ∧
      (∀ j, 0 < κ j) ∧ ∀ j, κ (j + 1) ≤ 2 * κ j := by
  set J : ℕ := ⌈3 * σ⌉₊ + 1 with hJ
  have hJpos : (0 : ℝ) < J := by positivity
  have hJ3 : 3 * σ ≤ J := by
    have := Nat.le_ceil (3 * σ)
    simp only [hJ, Nat.cast_add, Nat.cast_one]
    linarith
  set ρ : ℝ := (6 : ℝ) ^ (σ / J) with hρ
  have hρ1 : 1 < ρ := Real.one_lt_rpow (by norm_num) (div_pos hσ hJpos)
  have hρ2 : ρ ≤ 2 := by
    calc ρ ≤ (8 : ℝ) ^ (σ / J) :=
          Real.rpow_le_rpow (by norm_num) (by norm_num) (div_pos hσ hJpos).le
      _ = (2 : ℝ) ^ (3 * σ / J) := by
          rw [show (8 : ℝ) = 2 ^ (3 : ℝ) by norm_num, ← Real.rpow_mul (by norm_num)]
          ring_nf
      _ ≤ (2 : ℝ) ^ (1 : ℝ) :=
          Real.rpow_le_rpow_of_exponent_le (by norm_num) ((div_le_one hJpos).2 hJ3)
      _ = 2 := Real.rpow_one 2
  refine ⟨J, fun j ↦ (3 : ℝ) ^ (-σ) * ρ ^ j, by simp, ?_, ?_, ?_, ?_⟩
  · simp only
    rw [hρ, ← Real.rpow_natCast, ← Real.rpow_mul (by norm_num),
      div_mul_cancel₀ σ hJpos.ne', Real.rpow_neg (by norm_num),
      show (6 : ℝ) = 3 * 2 by norm_num, Real.mul_rpow (by norm_num) (by norm_num)]
    field_simp
  · intro j k hjk
    exact mul_lt_mul_of_pos_left (pow_lt_pow_right₀ hρ1 hjk) (by positivity)
  · intro j; positivity
  · intro j
    simp only [pow_succ]
    have : 0 < (3 : ℝ) ^ (-σ) * ρ ^ j := by positivity
    nlinarith

/-- Separation of frequencies outside `[κ₀, κ_J)` from the values of `H'`. -/
private theorem sep_of_deriv_mem {H : ℝ → ℝ} {κ₀ κJ g : ℝ} (hκ₀ : 0 < κ₀) (hκ : κ₀ < κJ)
    (hg : 0 < g) (hg1 : g ≤ 1) {s : Set ℝ}
    (hmem : ∀ x ∈ s, κ₀ + g ≤ deriv H x ∧ deriv H x ≤ κJ - g)
    {ν : ℝ} (hν : ν < κ₀ ∨ κJ ≤ ν) (x : ℝ) (hx : x ∈ s) :
    g / (1 + κJ) * (1 + |ν|) ≤ |deriv H x - ν| := by
  obtain ⟨h1, h2⟩ := hmem x hx
  have hκJ : 0 < 1 + κJ := by linarith
  rw [div_mul_eq_mul_div, div_le_iff₀ hκJ]
  rcases hν with hν | hν
  · rw [abs_of_pos (by linarith : 0 < deriv H x - ν)]
    rcases le_or_gt 0 ν with h0 | h0
    · rw [abs_of_nonneg h0]
      nlinarith
    · rw [abs_of_neg h0]
      nlinarith
  · rw [abs_of_neg (by linarith : deriv H x - ν < 0), abs_of_pos (by linarith : 0 < ν)]
    nlinarith

/-- Stationary points of a uniformly concave phase grow at most linearly. -/
private theorem abs_phaseInverse_le {H : ℝ → ℝ} (hH : ContDiff ℝ ∞ H) {c : ℝ} (hc : 0 < c)
    (hcurv : ∀ x, iteratedDeriv 2 H x ≤ -c) (ν : ℝ) :
    |phaseInverse H ν| ≤ 1 + (|ν| + |deriv H 1|) / c := by
  obtain ⟨-, hgf, hfg⟩ := phaseInverse_spec hH hc hcurv
  have hd2 : ∀ x, deriv (deriv H) x ≤ -c := fun x ↦ by
    simpa [iteratedDeriv_succ, iteratedDeriv_one] using hcurv x
  have h := abs_inverse_sub_le_of_neg ((hH.deriv' (n := ∞)).differentiable (by simp))
    hc hd2 hfg ν (deriv H 1)
  rw [hgf] at h
  have h' : |ν - deriv H 1| ≤ |ν| + |deriv H 1| := abs_sub _ _
  calc |phaseInverse H ν| = |(phaseInverse H ν - 1) + 1| := by ring_nf
    _ ≤ |phaseInverse H ν - 1| + 1 := by simpa using abs_add_le (phaseInverse H ν - 1) 1
    _ ≤ (|ν| + |deriv H 1|) / c + 1 := by
        linarith [div_le_div_of_nonneg_right h' hc.le]
    _ = _ := by ring

/-- The fixed-parameter B-process inequality. -/
private theorem norm_smooth_sum_le_blocks {H : ℝ → ℝ} (hH : ContDiff ℝ ∞ H) {c C : ℝ} (hc : 0 < c)
    (hcurv : ∀ x, -C ≤ iteratedDeriv 2 H x ∧ iteratedDeriv 2 H x ≤ -c)
    {T N : ℝ} (hN1 : 1 ≤ N) (hNT : N ≤ T) {Bc : ℕ} (hBc : 9 / 4 * N ≤ Bc)
    {J : ℕ} {κ : ℕ → ℝ} (hκmono : Monotone κ) (hκpos : ∀ j, 0 < κ j)
    {M0 Bamp δ K3 W : ℝ} (hM0 : 0 ≤ M0) (hδ : 0 < δ)
    (hstat : ∀ ν ∈ Icc (κ 0) (κ J), |phaseInverse H ν| ≤ M0)
    (hamp : ∀ v, |v| ≤ M0 → ∀ y, ‖iteratedDeriv 2 (stationaryAmplitude H v) y‖ ≤ Bamp ∧
      ‖iteratedDeriv 4 (stationaryAmplitude H v) y‖ ≤ Bamp)
    (hsep : ∀ ν, (ν < κ 0 ∨ κ J ≤ ν) → ∀ x ∈ Icc (5 / 8 : ℝ) (19 / 8),
      δ * (1 + |ν|) ≤ |deriv H x - ν|)
    (h3 : ∀ x, |iteratedDeriv 3 H x| ≤ K3)
    (hW : ∀ k ≤ 2, ∀ x, |iteratedDeriv k (phaseCutoff : ℝ → ℝ) x| ≤ W)
    {D : ℕ → ℝ}
    (hD : ∀ j < J, ∀ k ∈ Finset.Icc ⌈κ j * (T / N)⌉₊ ⌈κ (j + 1) * (T / N)⌉₊,
      ‖∑ m ∈ Finset.Ico ⌈κ j * (T / N)⌉₊ k, (𝐞 (T * legendre H (N * m / T)) : ℂ)‖ ≤ D j) :
    ‖∑ n ∈ Finset.Icc 0 Bc, (phaseCutoff ((n : ℝ) / N) : ℂ) * 𝐞 (T * H ((n : ℝ) / N))‖ ≤
      N / Real.sqrt T * ((W / Real.sqrt c + (W / Real.sqrt c +
          W * (K3 / (2 * c * Real.sqrt c))) / c * (κ J + 1)) * ∑ j ∈ Finset.range J, D j) +
        ((κ J + 1) * (4 * Real.pi ^ 2 * ((9 / 4 + M0) * Real.sqrt C) * Bamp) +
          (7 / 4 * W * (1 + 3 * C / δ + (K3 / δ + 3 * C ^ 2 / δ ^ 2))) *
            (∑' m : ℤ, 1 / (1 + |(m : ℝ)|) ^ 2) / δ ^ 2) := by
  have hN : 0 < N := zero_lt_one.trans_le hN1
  have hT : 0 < T := hN.trans_le hNT
  have hT1 : 1 ≤ T := hN1.trans hNT
  set M : ℝ := T / N with hMdef
  have hM : 0 < M := div_pos hT hN
  have hM1 : 1 ≤ M := (one_le_div hN).2 hNT
  set P : ℕ → ℕ := fun j ↦ ⌈κ j * M⌉₊ with hPdef
  have hPmono : Monotone P := fun j k hjk ↦
    Nat.ceil_mono (mul_le_mul_of_nonneg_right (hκmono hjk) hM.le)
  have hC : 0 ≤ C := by linarith [(hcurv 0).1, (hcurv 0).2]
  have hW0 : 0 ≤ W := (abs_nonneg _).trans (hW 0 (by norm_num) 0)
  have hBamp : 0 ≤ Bamp := (norm_nonneg _).trans (hamp 0 (by simpa using hM0) 0).1
  have hS0 : 0 ≤ ∑' m : ℤ, 1 / (1 + |(m : ℝ)|) ^ 2 := tsum_nonneg fun m ↦ by positivity
  have hsc : 0 < Real.sqrt c := Real.sqrt_pos.2 hc
  have hνeq : ∀ m : ℝ, N * m / T = m / M := fun m ↦ by rw [hMdef]; field_simp
  -- stationary frequencies
  have hstat' : ∀ m ∈ Finset.Ico (P 0) (P J), |phaseInverse H (N * m / T)| ≤ M0 := by
    intro m hm
    obtain ⟨h1, h2⟩ := Finset.mem_Ico.1 hm
    apply hstat
    rw [hνeq]
    constructor
    · rw [le_div_iff₀ hM]
      exact (Nat.le_ceil _).trans (by exact_mod_cast h1)
    · rw [div_le_iff₀ hM]
      exact (Nat.lt_ceil.1 h2).le
  -- nonstationary frequencies
  have hsep' : ∀ m : ℤ, (m < P 0 ∨ ((P J : ℕ) : ℤ) ≤ m) → ∀ x ∈ Icc (5 / 8 : ℝ) (19 / 8),
      δ * (1 + |N * m / T|) ≤ |deriv H x - N * m / T| := by
    intro m hm x hx
    apply hsep _ _ x hx
    rw [hνeq]
    rcases hm with hm | hm
    · left
      rw [div_lt_iff₀ hM]
      rcases lt_or_ge m 0 with hneg | hnn
      · have : (m : ℝ) < 0 := by exact_mod_cast hneg
        nlinarith [hκpos 0]
      · have hm' : m.toNat < P 0 := by omega
        have := Nat.lt_ceil.1 hm'
        have hcast : ((m.toNat : ℕ) : ℝ) = (m : ℝ) := by
          rw [← Int.cast_natCast, Int.toNat_of_nonneg hnn]
        rw [hcast] at this
        exact this
    · right
      rw [le_div_iff₀ hM]
      exact (Nat.le_ceil _).trans (by exact_mod_cast hm)
  have hK3 : 0 ≤ K3 := (abs_nonneg _).trans (h3 0)
  have hdec := norm_smooth_sum_sub_main_le hH hc hcurv hT hN hNT hBc hδ hstat' hamp hsep'
    (fun x _ ↦ h3 x) hW
  -- main terms
  have hblocks : ∑ m ∈ Finset.Ico (P 0) (P J), mainTerm H T (N * m / T) =
      ∑ j ∈ Finset.range J, ∑ m ∈ Finset.Ico (P j) (P (j + 1)), mainTerm H T (N * m / T) :=
    sum_Ico_eq_sum_blocks _ hPmono J
  set Amax : ℝ := W / Real.sqrt c with hAmax
  set LA : ℝ := (W / Real.sqrt c + W * (K3 / (2 * c * Real.sqrt c))) / c with hLA
  have hLA0 : 0 ≤ LA := by positivity
  have hAb : ∀ ν, |mainAmplitude H ν| ≤ Amax := abs_mainAmplitude_le hc hcurv
    (fun x ↦ by simpa using hW 0 (by norm_num) x)
  have hLip : ∀ ν ν', |mainAmplitude H ν - mainAmplitude H ν'| ≤ LA * |ν - ν'| :=
    abs_mainAmplitude_sub_le hH hc hcurv (fun x ↦ by simpa using hW 0 (by norm_num) x)
      (fun x ↦ by simpa [iteratedDeriv_one] using hW 1 (by norm_num) x) h3
  have hblock : ∀ j ∈ Finset.range J,
      ‖∑ m ∈ Finset.Ico (P j) (P (j + 1)), mainTerm H T (N * m / T)‖ ≤
        (Amax + LA * (κ J + 1)) * D j / Real.sqrt T := by
    intro j hj
    have hjJ := Finset.mem_range.1 hj
    have hDj := hD j hjJ
    have hD0 : 0 ≤ D j := by
      have := hDj (P j) (Finset.mem_Icc.2 ⟨le_rfl, hPmono (Nat.le_succ j)⟩)
      simp only [P, Finset.Ico_self, Finset.sum_empty, norm_zero] at this
      exact this
    have h1 := norm_sum_mainTerm_block_le (H := H) hT hN hAb hLip (hPmono (Nat.le_succ j)) hDj
    refine h1.trans (div_le_div_of_nonneg_right ?_ (Real.sqrt_nonneg _))
    apply mul_le_mul_of_nonneg_right _ hD0
    apply add_le_add le_rfl
    rw [mul_div_assoc]
    apply mul_le_mul_of_nonneg_left _ hLA0
    rw [div_le_iff₀ hM]
    have hQ : ((P (j + 1) : ℕ) : ℝ) < κ (j + 1) * M + 1 :=
      Nat.ceil_lt_add_one (mul_nonneg (hκpos _).le hM.le)
    have hP : κ j * M ≤ ((P j : ℕ) : ℝ) := Nat.le_ceil _
    have hκj : κ (j + 1) ≤ κ J := hκmono (by omega)
    have : 0 ≤ κ j * M := mul_nonneg (hκpos _).le hM.le
    nlinarith
  have hmain : ‖(N : ℂ) * ∑ m ∈ Finset.Ico (P 0) (P J), mainTerm H T (N * m / T)‖ ≤
      N / Real.sqrt T * ((Amax + LA * (κ J + 1)) * ∑ j ∈ Finset.range J, D j) := by
    rw [norm_mul, Complex.norm_real, Real.norm_of_nonneg hN.le, hblocks]
    calc N * ‖∑ j ∈ Finset.range J, ∑ m ∈ Finset.Ico (P j) (P (j + 1)),
          mainTerm H T (N * m / T)‖
        ≤ N * ∑ j ∈ Finset.range J, (Amax + LA * (κ J + 1)) * D j / Real.sqrt T :=
          mul_le_mul_of_nonneg_left ((norm_sum_le _ _).trans (Finset.sum_le_sum hblock)) hN.le
      _ = _ := by rw [← Finset.sum_div, ← Finset.mul_sum]; ring
  -- errors
  set R : ℝ := 4 * Real.pi ^ 2 * ((9 / 4 + M0) * Real.sqrt C) * Bamp with hR
  have hR0 : 0 ≤ R := by positivity
  have herr1 : N * (((P J - P 0 : ℕ) : ℝ) * (R / (T * Real.sqrt T))) ≤ (κ J + 1) * R := by
    have hcount : ((P J - P 0 : ℕ) : ℝ) ≤ κ J * M + 1 := by
      have h1 : ((P J - P 0 : ℕ) : ℝ) ≤ ((P J : ℕ) : ℝ) := by exact_mod_cast Nat.sub_le _ _
      have h2 : ((P J : ℕ) : ℝ) < κ J * M + 1 :=
        Nat.ceil_lt_add_one (mul_nonneg (hκpos _).le hM.le)
      linarith
    have hsT : 1 ≤ Real.sqrt T := Real.one_le_sqrt.2 hT1
    have hsTpos : 0 < Real.sqrt T := zero_lt_one.trans_le hsT
    calc N * (((P J - P 0 : ℕ) : ℝ) * (R / (T * Real.sqrt T)))
        ≤ N * ((κ J * M + 1) * (R / (T * Real.sqrt T))) := by
          gcongr
      _ = (κ J * T + N) * R / (T * Real.sqrt T) := by
          rw [hMdef]; field_simp
      _ ≤ (κ J * T + T) * R / (T * Real.sqrt T) := by gcongr
      _ = (κ J + 1) * R / Real.sqrt T := by field_simp
      _ ≤ (κ J + 1) * R := div_le_self (mul_nonneg (by linarith [hκpos J]) hR0) hsT
  set Kn : ℝ := 7 / 4 * W * (1 + 3 * C / δ + (K3 / δ + 3 * C ^ 2 / δ ^ 2)) with hKn
  have hKn0 : 0 ≤ Kn := by positivity
  have herr2 : Kn * (∑' m : ℤ, 1 / (1 + |(m : ℝ)|) ^ 2) / (δ ^ 2 * N) ≤
      Kn * (∑' m : ℤ, 1 / (1 + |(m : ℝ)|) ^ 2) / δ ^ 2 := by
    apply div_le_div_of_nonneg_left (by positivity) (by positivity)
    nlinarith [sq_nonneg δ, pow_pos hδ 2]
  have htri := norm_le_insert' (∑ n ∈ Finset.Icc 0 Bc,
      (phaseCutoff ((n : ℝ) / N) : ℂ) * 𝐞 (T * H ((n : ℝ) / N)))
    ((N : ℂ) * ∑ m ∈ Finset.Ico (P 0) (P J), mainTerm H T (N * m / T))
  linarith

end Assembly

section Setup

/-- Uniform constants for the fixed-parameter B-process inequality. -/
private theorem RegularizedPhases.eventually_setup {σ ar : ℝ} {G H : ℕ → ℝ → ℝ} {H₀ : ℝ → ℝ}
    (hR : RegularizedPhases σ ar 4 G H H₀) (hσ : 0 < σ) (har0 : 0 < ar) (har : ar ≤ 1 / 2)
    {J : ℕ} {κ : ℕ → ℝ} (hκ0 : κ 0 = (3 : ℝ) ^ (-σ)) (hκJ : κ J = (2 : ℝ) ^ σ)
    (hκpos : ∀ j, 0 < κ j) :
    ∃ c C M0 Bamp δ K3 : ℝ, 0 < c ∧ 0 ≤ M0 ∧ 0 < δ ∧
      (∀ i x, -C ≤ iteratedDeriv 2 (H i) x ∧ iteratedDeriv 2 (H i) x ≤ -c) ∧
      ∀ᶠ i in atTop,
        (∀ ν ∈ Icc (κ 0) (κ J), |phaseInverse (H i) ν| ≤ M0) ∧
        (∀ v, |v| ≤ M0 → ∀ y, ‖iteratedDeriv 2 (stationaryAmplitude (H i) v) y‖ ≤ Bamp ∧
          ‖iteratedDeriv 4 (stationaryAmplitude (H i) v) y‖ ≤ Bamp) ∧
        (∀ ν, (ν < κ 0 ∨ κ J ≤ ν) → ∀ x ∈ Icc (5 / 8 : ℝ) (19 / 8),
          δ * (1 + |ν|) ≤ |deriv (H i) x - ν|) ∧
        (∀ x, |iteratedDeriv 3 (H i) x| ≤ K3) := by
  obtain ⟨c, C, hc, hcurv, -⟩ := hR.curv
  have hbdd : ∀ k, 2 ≤ k → ∃ K, ∀ᶠ i in atTop, ∀ x, ‖iteratedDeriv k (H i) x‖ ≤ K := by
    intro k hk
    obtain ⟨K, hK⟩ := hR.bdd k hk
    refine ⟨K + 1, (hR.conv k (by omega) 1 one_pos).mono fun i hi x ↦ ?_⟩
    have h1 := hi x
    have h2 := hK x
    calc ‖iteratedDeriv k (H i) x‖
        = ‖(iteratedDeriv k (H i) x - iteratedDeriv k H₀ x) + iteratedDeriv k H₀ x‖ := by
          ring_nf
      _ ≤ ‖iteratedDeriv k (H i) x - iteratedDeriv k H₀ x‖ + ‖iteratedDeriv k H₀ x‖ :=
          norm_add_le _ _
      _ ≤ K + 1 := by linarith
  obtain ⟨K3, hK3⟩ := hbdd 3 (by norm_num)
  -- the stationary points
  have hJ : (0 : ℝ) < κ J := hκpos J
  set M0 : ℝ := 1 + (κ J + 2) / c with hM0
  have hM0pos : 0 ≤ M0 := by positivity
  have h1mem : (1 : ℝ) ∈ Icc ar 4 := ⟨by linarith, by norm_num⟩
  have hd1 : ∀ᶠ i in atTop, |deriv (H i) 1| ≤ 2 := by
    filter_upwards [hR.conv 1 le_rfl 1 one_pos] with i hi
    have h := hi 1
    simp only [iteratedDeriv_one] at h
    have hm : deriv H₀ 1 = 1 := by
      rw [hR.model h1mem]
      simp [modelPhase]
    rw [hm, Real.norm_eq_abs] at h
    have := abs_sub_abs_le_abs_sub (deriv (H i) 1) 1
    rw [abs_one] at this
    linarith
  -- the amplitudes
  have hw0 : ∀ x, x ∉ Icc (3 / 4 : ℝ) (9 / 4) → phaseCutoff x = 0 :=
    fun x hx ↦ phaseCutoff_eq_zero fun h ↦ hx ⟨h.1.le, h.2.le⟩
  obtain ⟨B2, hB2⟩ := exists_amplitude_bound hR.smooth hc hcurv hbdd phaseCutoff.contDiff
    (p := 3 / 4) (q := 9 / 4) (by norm_num) hw0 (M := M0) 2
  obtain ⟨B4, hB4⟩ := exists_amplitude_bound hR.smooth hc hcurv hbdd phaseCutoff.contDiff
    (p := 3 / 4) (q := 9 / 4) (by norm_num) hw0 (M := M0) 4
  -- the separation
  have hA : (3 : ℝ) ^ (-σ) < (19 / 8 : ℝ) ^ (-σ) :=
    Real.rpow_lt_rpow_of_neg (by norm_num) (by norm_num) (by linarith)
  have h58 : (5 / 8 : ℝ) ^ (-σ) = (8 / 5 : ℝ) ^ σ := by
    rw [Real.rpow_neg (by norm_num), ← Real.inv_rpow (by norm_num)]
    norm_num
  have hB : (8 / 5 : ℝ) ^ σ < (2 : ℝ) ^ σ := Real.rpow_lt_rpow (by norm_num) (by norm_num) hσ
  set g : ℝ := min 1 (min ((19 / 8 : ℝ) ^ (-σ) - (3 : ℝ) ^ (-σ))
    ((2 : ℝ) ^ σ - (8 / 5 : ℝ) ^ σ)) / 2 with hg
  have hgpos : 0 < g := by
    simp only [hg]
    apply div_pos _ two_pos
    exact lt_min one_pos (lt_min (by linarith) (by linarith))
  have hg1 : g ≤ 1 := by
    have := min_le_left (1 : ℝ) (min ((19 / 8 : ℝ) ^ (-σ) - (3 : ℝ) ^ (-σ))
      ((2 : ℝ) ^ σ - (8 / 5 : ℝ) ^ σ))
    simp only [hg]; linarith
  have hgA : 2 * g ≤ (19 / 8 : ℝ) ^ (-σ) - (3 : ℝ) ^ (-σ) := by
    have := (min_le_right (1 : ℝ) _).trans (min_le_left ((19 / 8 : ℝ) ^ (-σ) - (3 : ℝ) ^ (-σ))
      ((2 : ℝ) ^ σ - (8 / 5 : ℝ) ^ σ))
    simp only [hg]; linarith
  have hgB : 2 * g ≤ (2 : ℝ) ^ σ - (8 / 5 : ℝ) ^ σ := by
    have := (min_le_right (1 : ℝ) _).trans (min_le_right ((19 / 8 : ℝ) ^ (-σ) - (3 : ℝ) ^ (-σ))
      ((2 : ℝ) ^ σ - (8 / 5 : ℝ) ^ σ))
    simp only [hg]; linarith
  have hκlt : κ 0 < κ J := by
    rw [hκ0, hκJ]
    have : (3 : ℝ) ^ (-σ) < 1 := Real.rpow_lt_one_of_one_lt_of_neg (by norm_num) (by linarith)
    have : (1 : ℝ) < 2 ^ σ := Real.one_lt_rpow (by norm_num) hσ
    linarith
  have hmemd : ∀ᶠ i in atTop, ∀ x ∈ Icc (5 / 8 : ℝ) (19 / 8),
      κ 0 + g ≤ deriv (H i) x ∧ deriv (H i) x ≤ κ J - g := by
    filter_upwards [hR.conv 1 le_rfl g hgpos] with i hi x hx
    have h := hi x
    simp only [iteratedDeriv_one, Real.norm_eq_abs] at h
    have hxm : x ∈ Icc ar 4 := ⟨by linarith [hx.1], by linarith [hx.2]⟩
    rw [hR.model hxm] at h
    have hx0 : 0 < x := by linarith [hx.1]
    have hlo : (19 / 8 : ℝ) ^ (-σ) ≤ modelPhase σ x :=
      Real.rpow_le_rpow_of_nonpos hx0 hx.2 (by linarith)
    have hhi : modelPhase σ x ≤ (8 / 5 : ℝ) ^ σ := by
      rw [← h58]; exact Real.rpow_le_rpow_of_nonpos (by norm_num) hx.1 (by linarith)
    rw [abs_le] at h
    rw [hκ0, hκJ]
    constructor <;> linarith [h.1, h.2]
  refine ⟨c, C, M0, max B2 B4, g / (1 + κ J), K3, hc, hM0pos, by positivity, hcurv, ?_⟩
  filter_upwards [hd1, hB2, hB4, hmemd, hK3] with i hi1 hiB2 hiB4 himemd hiK3
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro ν hν
    have hνabs : |ν| ≤ κ J := by
      rw [abs_of_pos (lt_of_lt_of_le (hκpos 0) hν.1)]; exact hν.2
    refine (abs_phaseInverse_le (hR.smooth i) hc (fun x ↦ (hcurv i x).2) ν).trans ?_
    simp only [hM0]
    gcongr
  · intro v hv y
    exact ⟨(hiB2 v hv y).trans (le_max_left _ _), (hiB4 v hv y).trans (le_max_right _ _)⟩
  · intro ν hν x hx
    exact sep_of_deriv_mem (hκpos 0) hκlt hgpos hg1 himemd hν x hx
  · intro x
    simpa [Real.norm_eq_abs] using hiK3 x

end Setup

section PowerBounds

/-- Power bounds for the dual exponential sums of a B-process block. -/
private theorem isPowerBounded_dual {α : ℝ≥0} (hα1 : α < 1) {B σ' κ S : ℝ}
    (hB : IsExponentSumBound (1 - α) B) (hσ' : 0 < σ') (hκ : 0 < κ) (hS : 0 < S)
    {Hd : ℕ → ℝ → ℝ} (hmodel : IsModelPhaseFunctionWith Hd σ')
    {N T : VariableObject ℝ} (hT : ∀ i, 1 ≤ T i)
    (hTunb : T.IsUnbounded) (hNT : IsPowerAsymptotic N T α)
    {a' b' : ℕ → ℕ}
    (hab : ∀ᶠ i in atTop, IsDyadicSubinterval (κ * (T i / N i)) (a' i) (b' i)) :
    IsPowerBounded
      (fun i ↦ exponentialSumAt (Hd i) (T i * S) (κ * (T i / N i)) (a' i) (b' i)) T B := by
  set T' : VariableObject ℝ := fun i ↦ T i * S with hT'
  set N' : VariableObject ℝ := fun i ↦ κ * (T i / N i) with hN'
  have hTtop : Tendsto T atTop atTop :=
    (VariableObject.isUnbounded_iff_tendsto_atTop (fun i ↦ zero_le_one.trans (hT i))).1 hTunb
  have hT'top : Tendsto T' atTop atTop := hTtop.atTop_mul_const hS
  have hα1' : (α : ℝ) < 1 := by exact_mod_cast hα1
  have hTN : IsPowerAsymptotic (T / N) T (1 - α) := (isPowerAsymptotic_self T).div hNT hT
  have hN'T : IsPowerAsymptotic N' T (1 - α) := by
    have := (isPowerAsymptotic_zero_of_tendsto_pos hκ tendsto_const_nhds hT hTunb).mul hTN hT
    rw [zero_add] at this
    exact this
  have hT'T : IsPowerAsymptotic T' T 1 := by
    have := (isPowerAsymptotic_self T).mul (isPowerAsymptotic_zero_of_tendsto_pos hS
      tendsto_const_nhds hT hTunb) hT
    rw [add_zero] at this
    exact this
  have hN'T' : IsPowerAsymptotic N' T' ((1 - α : ℝ≥0) : ℝ) := by
    have := hN'T.change_base hT'T one_pos hT
    rw [div_one] at this
    rwa [NNReal.coe_sub hα1.le, NNReal.coe_one]
  have hN'top : Tendsto N' atTop atTop :=
    hN'T.tendsto_atTop_of_pos (by linarith) hT hTunb
  have hsetup : ∀ᶠ i in atTop, 1 ≤ N' i ∧ 1 ≤ T' i ∧
      IsDyadicSubinterval (N' i) (a' i) (b' i) := by
    filter_upwards [hN'top.eventually (eventually_ge_atTop 1),
      hT'top.eventually (eventually_ge_atTop 1), hab] with i h1 h2 h3
    exact ⟨h1, h2, h3⟩
  have h := IsExponentSumBound.isPowerBounded_of_eventually_setup hB hσ' hT'top hN'T' hmodel hsetup
  have h' := h.change_base hT'T hT
  rw [one_mul] at h'
  exact h'

/-- A logarithm of a quantity of polynomial size has power exponent zero. -/
private theorem isPowerBounded_three_add_two_log {X T : VariableObject ℝ} {C p : ℝ} (hC : 1 ≤ C)
    (hp : 0 ≤ p) (hT : ∀ i, 1 ≤ T i) (hTunb : T.IsUnbounded)
    (hX1 : ∀ᶠ i in atTop, 1 ≤ X i) (hXle : ∀ᶠ i in atTop, X i ≤ C * T i ^ p) :
    IsPowerBounded (fun i ↦ 3 + 2 * Real.log (X i)) T 0 := by
  have hlog := isPowerBounded_log_of_le_const_rpow hC hp hT hTunb hX1 hXle
  have hconst := isPowerBounded_const T 3
  simpa only [Pi.add_def, max_self] using hconst.add (hlog.const_mul 2) hT

/-- Twisting by a linear function only changes the first derivative. -/
private theorem iteratedDeriv_succ_sub_linear {g : ℝ → ℝ} (hg : ContDiff ℝ ∞ g) (c : ℝ) (p : ℕ)
    (x : ℝ) :
    iteratedDeriv (p + 1) (fun y ↦ g y - c * y) x =
      iteratedDeriv (p + 1) g x - if p = 0 then c else 0 := by
  rw [iteratedDeriv_fun_sub
    (hg.of_le (ENat.natCast_le_of_coe_top_le_withTop le_rfl (p + 1))).contDiffAt (by fun_prop),
    iteratedDeriv_const_mul_field]
  change _ - c * iteratedDeriv (p + 1) id x = _
  simp only [iteratedDeriv_id]
  split_ifs <;> simp_all

end PowerBounds

end

end Expdb.BProcess

namespace Expdb

open Filter Topology Set Complex
open scoped NNReal ContDiff FourierTransform
open BProcess StationaryPhase

noncomputable section

/-- **The B-process transfer.** An admissible exponential-sum exponent `B` at scale `1 - α`
gives the admissible exponent `α - 1 / 2 + B` at scale `α`, for `0 < α < 1`. -/
theorem isExponentSumBound_bProcess {α : ℝ≥0} (hα0 : 0 < α) (hα1 : α < 1) {B : ℝ}
    (hB : IsExponentSumBound (1 - α) B) :
    IsExponentSumBound α ((α : ℝ) - 1 / 2 + B) := by
  intro N T F a b hN hT hTunb hNT hF hab
  obtain ⟨σ, hσ, hFσ⟩ := hF.exists_with
  have hTtop : Tendsto T atTop atTop :=
    (VariableObject.isUnbounded_iff_tendsto_atTop (fun i ↦ zero_le_one.trans (hT i))).1 hTunb
  have hα0' : (0 : ℝ) < α := by exact_mod_cast hα0
  have hα1' : (α : ℝ) < 1 := by exact_mod_cast hα1
  have hNpos : ∀ i, 0 < N i := fun i ↦ zero_lt_one.trans_le (hN i)
  have hTpos : ∀ i, 0 < T i := fun i ↦ zero_lt_one.trans_le (hT i)
  have hBlow : (1 - (α : ℝ)) / 2 ≤ B := by
    have h := hB.fullInterval.half_le (tsub_pos_of_lt hα1) tsub_le_self
    rwa [NNReal.coe_sub hα1.le, NNReal.coe_one] at h
  have hexp : 0 ≤ (α : ℝ) - 1 / 2 + B := by linarith
  -- the length of the dual sums
  have hTN : IsPowerAsymptotic (T / N) T (1 - α) := (isPowerAsymptotic_self T).div hNT hT
  have hMtop : Tendsto (fun i ↦ T i / N i) atTop atTop :=
    hTN.tendsto_atTop_of_pos (by linarith) hT hTunb
  have hNleT : ∀ᶠ i in atTop, N i ≤ T i := by
    filter_upwards [hMtop.eventually (eventually_ge_atTop 1)] with i hi
    rwa [le_div_iff₀ (hNpos i), one_mul] at hi
  have hNT0 : Tendsto (fun i ↦ N i / T i) atTop (𝓝 0) := by
    have := tendsto_inv_atTop_zero.comp hMtop
    refine this.congr fun i ↦ ?_
    simp
  -- block constants
  obtain ⟨J, κ, hκ0, hκJ, hκmono, hκpos, hκ2⟩ := exists_block_constants hσ
  set ar : ℝ := min (1 / 2) ((2 * 2 ^ σ) ^ (-1 / σ)) with har_def
  have har0 : 0 < ar := lt_min (by norm_num) (by positivity)
  have har : ar ≤ 1 / 2 := min_le_left _ _
  have harσ : 2 * (2 : ℝ) ^ σ ≤ ar ^ (-σ) := by
    have hy : 0 < (2 * (2 : ℝ) ^ σ) := by positivity
    have h1 : ((2 * (2 : ℝ) ^ σ) ^ (-1 / σ)) ^ (-σ) = 2 * (2 : ℝ) ^ σ := by
      rw [← Real.rpow_mul hy.le, show -1 / σ * -σ = 1 by field_simp, Real.rpow_one]
    rw [← h1]
    exact Real.rpow_le_rpow_of_nonpos har0 (min_le_right _ _) (by linarith)
  -- phase extension
  obtain ⟨G, hGs, hGeq, hGconv⟩ := exists_global_extensions hFσ (l := ar / 2) (u := 5)
    (by positivity) (by linarith) (by norm_num)
  -- finite Fourier completion
  set Bc : ℕ → ℕ := fun i ↦ ⌈9 / 4 * N i⌉₊ with hBc
  have hcomp : ∀ i, ∃ r : ℝ, 0 ≤ r ∧ r ≤ 1 ∧
      ‖exponentialSumAt (G i) (T i) (N i) (a i) (b i)‖ ≤
        (3 + 2 * Real.log ((Bc i + 1 : ℕ) : ℝ)) *
          ‖∑ n ∈ Finset.Icc 0 (Bc i), (phaseCutoff ((n : ℝ) / N i) : ℂ) *
            (𝐞 (T i * G i ((n : ℝ) / N i) - r * n) : ℂ)‖ := by
    intro i
    refine exists_norm_exponentialSumAt_le_smooth_twist (G i) (T i) (hNpos i) (q := Bc i + 1)
      (A := 0) (B := Bc i) (by omega) (Nat.zero_le _) ?_ (by omega) (hab i).1 (hab i).2
    have h1 : (b i : ℝ) ≤ Bc i := by
      calc (b i : ℝ) ≤ 2 * N i := (hab i).2
        _ ≤ 9 / 4 * N i := by linarith [hNpos i]
        _ ≤ Bc i := Nat.le_ceil _
    exact_mod_cast h1
  choose r hr0 hr1 hrb using hcomp
  -- the twisted phases
  set cT : ℕ → ℝ := fun i ↦ r i * N i / T i with hcT
  have hcT0 : Tendsto cT atTop (𝓝 0) := by
    refine squeeze_zero (fun i ↦ ?_) (fun i ↦ ?_) hNT0
    · exact div_nonneg (mul_nonneg (hr0 i) (hNpos i).le) (hTpos i).le
    · simp only [hcT]
      exact div_le_div_of_nonneg_right
        (by nlinarith [hr1 i, hNpos i]) (hTpos i).le
  set Gt : ℕ → ℝ → ℝ := fun i x ↦ G i x - cT i * x with hGt
  have hGts : ∀ i, ContDiff ℝ ∞ (Gt i) := fun i ↦ (hGs i).sub (contDiff_const.mul contDiff_id)
  have hGtconv : ∀ p : ℕ, ∀ ε > 0, ∀ᶠ i in atTop, ∀ x ∈ Icc (ar / 2) (4 + 1),
      ‖iteratedDeriv (p + 1) (Gt i) x - iteratedDeriv p (modelPhase σ) x‖ ≤ ε := by
    intro p ε hε
    filter_upwards [hGconv p (ε / 2) (half_pos hε),
      (hcT0.norm.eventually
        (ge_mem_nhds (show ‖(0 : ℝ)‖ < ε / 2 by simpa using half_pos hε)))] with i hi hci x hx
    have hx' : x ∈ Icc (ar / 2) 5 := ⟨hx.1, by linarith [hx.2]⟩
    simp only [hGt]
    rw [iteratedDeriv_succ_sub_linear (hGs i)]
    have h1 := hi x hx'
    have h2 : ‖(if p = 0 then cT i else 0)‖ ≤ ε / 2 := by
      split_ifs
      · simpa using hci
      · simp; linarith
    calc ‖iteratedDeriv (p + 1) (G i) x - (if p = 0 then cT i else 0) -
          iteratedDeriv p (modelPhase σ) x‖
        = ‖(iteratedDeriv (p + 1) (G i) x - iteratedDeriv p (modelPhase σ) x) -
            (if p = 0 then cT i else 0)‖ := by ring_nf
      _ ≤ ‖iteratedDeriv (p + 1) (G i) x - iteratedDeriv p (modelPhase σ) x‖ +
            ‖(if p = 0 then cT i else 0)‖ := norm_sub_le _ _
      _ ≤ ε := by linarith
  -- regularization
  obtain ⟨H, H₀, hR⟩ := exists_regularizedPhases hσ har0 (by linarith) hGts hGtconv
  obtain ⟨c, C, M0, Bamp, δ, K3, hc, hM0, hδ, hcurv, hsetup⟩ :=
    hR.eventually_setup hσ har0 har hκ0 hκJ hκpos
  obtain ⟨W, hW⟩ := exists_phaseCutoff_bound 2
  -- the dual phases
  set S : ℕ → ℝ := fun j ↦ κ j ^ (1 - 1 / σ) with hS
  have hmodel : ∀ j < J, IsModelPhaseFunctionWith
      (fun i ↦ dualPhase (H i) (κ j) (S j)) (1 / σ) := by
    intro j hj
    apply hR.isModelPhaseFunctionWith_dualPhase hσ har0 (by linarith) (hκpos j)
    · have h1 : (4 : ℝ) ^ (-σ) < (3 : ℝ) ^ (-σ) :=
        Real.rpow_lt_rpow_of_neg (by norm_num) (by norm_num) (by linarith)
      have h2 : κ 0 ≤ κ j := hκmono.monotone (Nat.zero_le j)
      rw [hκ0] at h2
      linarith
    · have h1 : κ j < κ J := hκmono hj
      rw [hκJ] at h1
      linarith
  -- maximizing partial sums
  set Pf : ℕ → ℕ → ℕ := fun j i ↦ ⌈κ j * (T i / N i)⌉₊ with hPf
  have hkex : ∀ j i, ∃ k ∈ Finset.Icc (Pf j i) (Pf (j + 1) i),
      ∀ k' ∈ Finset.Icc (Pf j i) (Pf (j + 1) i),
        ‖∑ m ∈ Finset.Ico (Pf j i) k', (𝐞 (T i * legendre (H i) (N i * m / T i)) : ℂ)‖ ≤
          ‖∑ m ∈ Finset.Ico (Pf j i) k, (𝐞 (T i * legendre (H i) (N i * m / T i)) : ℂ)‖ := by
    intro j i
    apply Finset.exists_max_image
    exact ⟨Pf j i, Finset.mem_Icc.2 ⟨le_rfl, Nat.ceil_mono
      (mul_le_mul_of_nonneg_right (hκmono.monotone (Nat.le_succ j)) (div_pos (hTpos i) (hNpos
        i)).le)⟩⟩
  choose kk hkk hkkmax using hkex
  set D : ℕ → ℕ → ℝ := fun j i ↦ ‖exponentialSumAt (dualPhase (H i) (κ j) (S j))
    (T i * S j) (κ j * (T i / N i)) (Pf j i) (kk j i - 1)‖ with hD
  have hDeq : ∀ j i, ‖∑ m ∈ Finset.Ico (Pf j i) (kk j i),
      (𝐞 (T i * legendre (H i) (N i * m / T i)) : ℂ)‖ = D j i := by
    intro j i
    have hP : 1 ≤ Pf j i := Nat.one_le_iff_ne_zero.2 (Nat.pos_iff_ne_zero.1
      (Nat.ceil_pos.2 (mul_pos (hκpos j) (div_pos (hTpos i) (hNpos i)))))
    exact norm_sum_Ico_fourierChar_legendre_eq (H i) (hTpos i) (hκpos j)
      (Real.rpow_pos_of_pos (hκpos j) _) hP
  have hDpow : ∀ j ∈ Finset.range J, IsPowerBounded (D j) T B := by
    intro j hj
    have hjJ := Finset.mem_range.1 hj
    refine (isPowerBounded_dual hα1 hB (one_div_pos.2 hσ) (hκpos j)
      (Real.rpow_pos_of_pos (hκpos j) _) (hmodel j hjJ) hT hTunb hNT
      (a' := Pf j) (b' := fun i ↦ kk j i - 1) (Eventually.of_forall fun i ↦ ?_)).norm
    have hM : 0 < T i / N i := div_pos (hTpos i) (hNpos i)
    refine ⟨Nat.le_ceil _, ?_⟩
    have hk := (Finset.mem_Icc.1 (hkk j i)).2
    have hQ : ((Pf (j + 1) i : ℕ) : ℝ) < κ (j + 1) * (T i / N i) + 1 :=
      Nat.ceil_lt_add_one (mul_nonneg (hκpos _).le hM.le)
    have hκ2' : κ (j + 1) * (T i / N i) ≤ 2 * (κ j * (T i / N i)) := by
      nlinarith [hκ2 j]
    have hk1 : ((kk j i - 1 : ℕ) : ℝ) ≤ κ (j + 1) * (T i / N i) := by
      rcases Nat.eq_zero_or_pos (kk j i) with h0 | hpos
      · rw [h0]; simp only [Nat.zero_sub, Nat.cast_zero]
        exact mul_nonneg (hκpos _).le hM.le
      · rw [Nat.cast_sub hpos, Nat.cast_one]
        have : ((kk j i : ℕ) : ℝ) ≤ ((Pf (j + 1) i : ℕ) : ℝ) := by exact_mod_cast hk
        linarith
    linarith
  -- the main inequality
  set KA : ℝ := W / Real.sqrt c + (W / Real.sqrt c + W * (K3 / (2 * c * Real.sqrt c))) / c *
    (κ J + 1) with hKA
  set KE : ℝ := (κ J + 1) * (4 * Real.pi ^ 2 * ((9 / 4 + M0) * Real.sqrt C) * Bamp) +
    (7 / 4 * W * (1 + 3 * C / δ + (K3 / δ + 3 * C ^ 2 / δ ^ 2))) *
      (∑' m : ℤ, 1 / (1 + |(m : ℝ)|) ^ 2) / δ ^ 2 with hKE
  have hmain : ∀ᶠ i in atTop, ‖exponentialSum F T N a b i‖ ≤
      (3 + 2 * Real.log ((Bc i + 1 : ℕ) : ℝ)) *
        (N i / Real.sqrt (T i) * (KA * ∑ j ∈ Finset.range J, D j i) + KE) := by
    filter_upwards [hGeq, hR.eqOn, hsetup, hNleT] with i hGi hHi hsi hNTi
    obtain ⟨hstat, hamp, hsep, h3⟩ := hsi
    have h1 : exponentialSum F T N a b i = exponentialSumAt (G i) (T i) (N i) (a i) (b i) := by
      simp only [exponentialSum_apply, exponentialSumAt]
      refine Finset.sum_congr rfl fun n hn ↦ ?_
      have hn' := Finset.mem_Icc.1 hn
      have hx : (n : ℝ) / N i ∈ phaseInterval := by
        constructor
        · rw [le_div_iff₀ (hNpos i), one_mul]
          exact (hab i).1.trans (by exact_mod_cast hn'.1)
        · rw [div_le_iff₀ (hNpos i)]
          have : (n : ℝ) ≤ b i := by exact_mod_cast hn'.2
          linarith [(hab i).2]
      simp only [oscillatory, hGi hx]
    have h2 : ∑ n ∈ Finset.Icc 0 (Bc i), (phaseCutoff ((n : ℝ) / N i) : ℂ) *
          (𝐞 (T i * G i ((n : ℝ) / N i) - r i * n) : ℂ) =
        ∑ n ∈ Finset.Icc 0 (Bc i), (phaseCutoff ((n : ℝ) / N i) : ℂ) *
          𝐞 (T i * H i ((n : ℝ) / N i)) := by
      refine Finset.sum_congr rfl fun n _ ↦ ?_
      by_cases hw : phaseCutoff ((n : ℝ) / N i) = 0
      · simp [hw]
      · have hmem : (n : ℝ) / N i ∈ Ioo (3 / 4 : ℝ) (9 / 4) := by
          by_contra h
          exact hw (phaseCutoff_eq_zero h)
        have hx : (n : ℝ) / N i ∈ Icc ar 4 := ⟨by linarith [hmem.1], by linarith [hmem.2]⟩
        rw [hHi hx]
        have : T i * Gt i ((n : ℝ) / N i) = T i * G i ((n : ℝ) / N i) - r i * n := by
          simp only [hGt, hcT]
          have := hNpos i
          have := hTpos i
          field_simp
        rw [this]
    have h3' := norm_smooth_sum_le_blocks (hR.smooth i) hc (hcurv i) (hN i) hNTi
      (Nat.le_ceil _) hκmono.monotone hκpos hM0 hδ hstat hamp hsep h3 hW
      (D := fun j ↦ D j i) (fun j _ k hk ↦ (hkkmax j i k hk).trans (le_of_eq (hDeq j i)))
    have hlog : 0 ≤ 3 + 2 * Real.log ((Bc i + 1 : ℕ) : ℝ) := by
      have : 0 ≤ Real.log ((Bc i + 1 : ℕ) : ℝ) :=
        Real.log_nonneg (by exact_mod_cast Nat.le_add_left 1 (Bc i))
      linarith
    rw [h1]
    refine (hrb i).trans ?_
    rw [h2]
    exact mul_le_mul_of_nonneg_left h3' hlog
  -- power bounds
  have hX1 : IsPowerBounded (fun i ↦ 3 + 2 * Real.log ((Bc i + 1 : ℕ) : ℝ)) T 0 := by
    apply isPowerBounded_three_add_two_log (C := 5) (p := 1) (by norm_num) zero_le_one hT hTunb
    · exact Eventually.of_forall fun i ↦ by exact_mod_cast Nat.le_add_left 1 (Bc i)
    · filter_upwards [hNleT] with i hi
      rw [Real.rpow_one]
      have h1 : ((Bc i : ℕ) : ℝ) < 9 / 4 * N i + 1 :=
        Nat.ceil_lt_add_one (by linarith [hNpos i])
      push_cast
      linarith [hT i]
  have hsqrt : IsPowerAsymptotic (fun i ↦ Real.sqrt (T i)) T (1 / 2) :=
    ⟨VariableObject.fixed (1 / 2), IsEqUpToInfinitesimal.refl _,
      Eventually.of_forall fun i ↦ Real.sqrt_eq_rpow (T i)⟩
  have hNs : IsPowerBounded (fun i ↦ N i / Real.sqrt (T i)) T ((α : ℝ) - 1 / 2) :=
    (hNT.div hsqrt hT).isPowerBounded
  have hsum := IsPowerBounded.const_mul (isPowerBounded_finset_sum (Finset.range J) hT hDpow) KA
  have hprod := hNs.mul hsum hT
  have hconst := (isPowerBounded_const T KE).mono (Eventually.of_forall hT) hexp
  have hX2 := hprod.add hconst hT
  rw [max_self] at hX2
  have hX := hX1.mul hX2 hT
  rw [zero_add] at hX
  exact hX.of_norm_le hmain

end

end Expdb
