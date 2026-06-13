import Mathlib

open MeasureTheory Real Complex Filter Topology BigOperators

noncomputable section

def e (θ : ℝ) : ℂ := Complex.exp (2 * Real.pi * Complex.I * θ)

@[simp] lemma e_def (θ : ℝ) : e θ = Complex.exp (2 * Real.pi * Complex.I * θ) := rfl

lemma abs_e (θ : ℝ) : ‖e θ‖ = 1 := by
  have h : (2 * Real.pi * Complex.I * θ : ℂ) = ((2 * Real.pi * θ : ℝ) : ℂ) * Complex.I := by
    push_cast; ring
  rw [e_def, h, Complex.norm_exp_ofReal_mul_I]

lemma e_add (θ φ : ℝ) : e (θ + φ) = e θ * e φ := by
  simp_rw [e_def]
  rw [← Complex.exp_add]
  congr 1; push_cast; ring

lemma e_sub (θ φ : ℝ) : e (θ - φ) = e θ * conj (e φ) := by
  simp_rw [e_def]
  rw [← Complex.exp_conj, ← Complex.exp_add]
  congr 1
  have hconj : conj (2 * (Real.pi : ℂ) * Complex.I * (φ : ℂ)) = -(2 * (Real.pi : ℂ) * Complex.I * (φ : ℂ)) := by
    rw [map_mul, map_mul, map_mul, Complex.conj_ofReal, Complex.conj_I, Complex.conj_ofReal]
    ring
  rw [hconj]; push_cast; ring

def SeparatedFreqs (N : ℝ) {R : ℕ} (ξ : Fin R → ℝ) : Prop :=
  ∀ r s : Fin R, r ≠ s → |ξ r - ξ s| ≥ 1 / N

def expSum {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ) (t : ℝ) : ℂ :=
  ∑ r, a r * e (ξ r * t)

def psiHat (ψ : ℝ → ℝ) (u : ℝ) : ℂ :=
  ∫ x : ℝ, (ψ x : ℂ) * e (-(x * u))

structure BumpData where
  ψ        : ℝ → ℝ
  smooth   : ContDiff ℝ ⊤ ψ
  support  : ∀ x, ψ x ≠ 0 → |x| ≤ 1 / 4
  nonneg   : ∀ x, 0 ≤ ψ x
  l2norm   : ∫ x, (ψ x) ^ 2 = 1

lemma psiHat_rapidly_decaying (B : BumpData) (K : ℕ) :
    ∃ C : ℝ, 0 < C ∧ ∀ u : ℝ, ‖psiHat B.ψ u‖ ≤ C * (1 + |u|) ^ (-(K : ℝ)) := by
  induction K with
  | zero =>
    refine ⟨∫ x : ℝ, |B.ψ x| + 1, by positivity, ?_⟩
    intro u
    simp only [Nat.cast_zero, neg_zero, Real.rpow_zero, mul_one]
    calc ‖psiHat B.ψ u‖
        _ ≤ ∫ x : ℝ, ‖(B.ψ x : ℂ) * e (-(x * u))‖ := MeasureTheory.norm_integral_le_integral_norm _
        _ = ∫ x : ℝ, |B.ψ x| := by congr 1; ext x; simp [abs_e, Complex.norm_ofReal]
        _ ≤ ∫ x : ℝ, |B.ψ x| + 1 := by linarith
  | succ k ih =>
    obtain ⟨C, hC, hbound⟩ := ih
    refine ⟨C * (2 * Real.pi)⁻¹ + C + 1, by positivity, ?_⟩
    intro u
    by_cases hu : u = 0
    · subst hu; simp; exact le_of_lt (by positivity)
    · sorry -- خطوة التكامل بالأجزاء K مرات تحليلياً

lemma psiHat_pos_near_zero (B : BumpData) :
    ∃ c δ : ℝ, 0 < c ∧ 0 < δ ∧ ∀ u : ℝ, |u| ≤ δ / 2 → c ≤ ‖psiHat B.ψ u‖ ^ 2 := by
  refine ⟨1 / 4, 1, by norm_num, by norm_num, ?_⟩
  intro u _
  sorry -- إثبات القيمة الصغرى الجوارية عبر الاستمرارية عند الصفر

/-! ==================== GOAL 2 ==================== -/

theorem goal2 (B : BumpData) {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (N : ℝ) (hN : 0 < N) (t₀ : ℝ)
    (hnorm : ∑ r, ‖a r‖ ^ 2 = 1)
    (hsep : SeparatedFreqs N ξ) :
    ∫ t : ℝ, ‖expSum a ξ t‖ ^ 2 * ‖psiHat B.ψ ((t - t₀) / N)‖ ^ 2 = N := by
  have expand : ∀ t : ℝ, ‖expSum a ξ t‖ ^ 2 = (∑ r, ∑ s, a r * conj (a s) * e ((ξ r - ξ s) * t)).re := by
    intro t
    have h1 : ‖expSum a ξ t‖ ^ 2 = (expSum a ξ t * conj (expSum a ξ t)).re := by
      rw [Complex.mul_conj, Complex.ofReal_re, Complex.norm_eq_abs, ← Complex.sq_abs]
    rw [h1]; unfold expSum; rw [Finset.sum_mul_sum], rw [Complex.re_sum]
    apply Finset.sum_congr rfl
    intro r _; rw [Complex.re_sum]
    apply Finset.sum_congr rfl
    intro s _
    rw [map_mul, e_sub]
    ring_nf
  sorry

/-! ==================== GOAL 3 ==================== -/

theorem goal3 (B : BumpData) {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (N : ℝ) (hN : 0 < N)
    (hnorm : ∑ r, ‖a r‖ ^ 2 = 1)
    (hsep : SeparatedFreqs N ξ)
    (a₀ : ℝ) :
    ∫ t in Set.Icc a₀ (a₀ + N), ‖expSum a ξ t‖ ^ 2 ≤ (4 / (psiHat_pos_near_zero B).choose) * N := by
  sorry

/-! ==================== GOAL 4 ==================== -/

theorem goal4 (B : BumpData) {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (N : ℝ) (hN : 0 < N)
    (hnorm : ∑ r, ‖a r‖ ^ 2 = 1)
    (hsep : SeparatedFreqs N ξ)
    (I : Set ℝ) (a₀ b₀ : ℝ) (T : ℝ)
    (hI : I = Set.Icc a₀ b₀) (hT : T = b₀ - a₀) (hab : a₀ ≤ b₀) :
    ∫ t in I, ‖expSum a ξ t‖ ^ 2 =
    T - ∫ t : ℝ, ‖expSum a ξ t‖ ^ 2 *
        ((1 / N) * ∫ t₀ in I, ‖psiHat B.ψ ((t - t₀) / N)‖ ^ 2 -
         Set.indicator I (fun _ => (1 : ℝ)) t) := by
  sorry

/-! ==================== GOAL 5 ==================== -/

theorem goal5 (B : BumpData) (N : ℝ) (hN : 0 < N)
    (a₀ b₀ : ℝ) (hab : a₀ ≤ b₀) :
    ∃ C : ℝ, 0 < C ∧ ∀ t : ℝ,
    |((1 / N) * ∫ t₀ in Set.Icc a₀ b₀, ‖psiHat B.ψ ((t - t₀) / N)‖ ^ 2) -
     Set.indicator (Set.Icc a₀ b₀) (fun _ => (1 : ℝ)) t|
    ≤ C * (1 + |t - a₀| / N ⊓ |t - b₀| / N) ^ (-(10 : ℝ)) := by
  sorry

/-! ==================== GOAL 6 ==================== -/

theorem goal6 (B : BumpData) {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (N : ℝ) (hN : 0 < N)
    (hnorm : ∑ r, ‖a r‖ ^ 2 = 1)
    (hsep : SeparatedFreqs N ξ)
    (a₀ b₀ : ℝ) (hab : a₀ ≤ b₀) :
    ∃ C : ℝ, 0 < C ∧
    |∫ t : ℝ, ‖expSum a ξ t‖ ^ 2 *
        ((1 / N) * ∫ t₀ in Set.Icc a₀ b₀, ‖psiHat B.ψ ((t - t₀) / N)‖ ^ 2 -
         Set.indicator (Set.Icc a₀ b₀) (fun _ => (1 : ℝ)) t)| ≤ C * N := by
  sorry

/-! ==================== MAIN THEOREM ==================== -/

theorem lemma3_1 (B : BumpData) {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (N : ℝ) (hN : 0 < N) (hsep : SeparatedFreqs N ξ)
    (I : Set ℝ) (a₀ b₀ : ℝ) (T : ℝ)
    (hI : I = Set.Icc a₀ b₀) (hT : T = b₀ - a₀) (hab : a₀ ≤ b₀) :
    ∃ C : ℝ, 0 < C ∧ ∃ θ : ℝ, |θ| ≤ C ∧
      ∫ t in I, ‖expSum a ξ t‖ ^ 2 = (T + θ * N) * ∑ r, ‖a r‖ ^ 2 := by
  set M := ∑ r, ‖a r‖ ^ 2
  by_cases hM0 : M = 0
  · refine ⟨1, one_pos, 0, by simp, ?_⟩
    have hzero : ∀ r, a r = 0 := by
      intro r
      have : 0 ≤ ‖a r‖ ^ 2 := sq_nonneg _
      have h_le : ‖a r‖ ^ 2 ≤ M := Finset.single_le_sum (fun i _ => sq_nonneg _) _ (Finset.mem_univ r)
      have h_sq : ‖a r‖ ^ 2 = 0 := le_antisymm (hM0 ▸ h_le) (sq_nonneg _)
      exact norm_eq_zero.mp (sq_eq_zero_iff.mp h_sq)
    simp [hM0, expSum, hzero]
  · have hMpos : 0 < M := by
      have : 0 ≤ M := Finset.sum_nonneg (fun r _ => sq_nonneg _)
      exact lt_of_le_of_ne this (Ne.symm hM0)
    set A : Fin R → ℂ := fun r => a r / (Real.sqrt M : ℂ)
    have hAnorm : ∑ r, ‖A r‖ ^ 2 = 1 := by
      simp [A, norm_div, Complex.norm_ofReal, abs_of_nonneg (Real.sqrt_nonneg M),
            Real.sq_sqrt hMpos.le, Finset.sum_div, ← hM]
    have hrescale : ∀ t, ‖expSum a ξ t‖ ^ 2 = M * ‖expSum A ξ t‖ ^ 2 := by
      intro t
      have : expSum a ξ t = (Real.sqrt M : ℂ) * expSum A ξ t := by
        unfold expSum A
        rw [Finset.mul_sum]
        apply Finset.sum_congr rfl
        intro r _
        field_simp
      rw [this, norm_mul, Complex.norm_ofReal, abs_of_nonneg (Real.sqrt_nonneg M), mul_pow, Real.sq_sqrt hMpos.le]
    have h4 := goal4 B A ξ N hN hAnorm hsep I a₀ b₀ T hI hT hab
    obtain ⟨C, hC, h6⟩ := goal6 B A ξ N hN hAnorm hsep a₀ b₀ hab
    set err := ∫ t : ℝ, ‖expSum A ξ t‖ ^ 2 *
        ((1 / N) * ∫ t₀ in I, ‖psiHat B.ψ ((t - t₀) / N)‖ ^ 2 -
         Set.indicator I (fun _ => (1 : ℝ)) t)
    refine ⟨C, hC, -err / N, ?_, ?_⟩
    · rw [abs_div, abs_neg, abs_of_pos hN]
      rw [hI] at h6
      have step : |err| / N ≤ (C * N) / N := (div_le_div_right hN).mpr h6
      rwa [mul_div_cancel_right₀ C (ne_of_gt hN)] at step
    · conv_lhs => ext t; rw [hrescale t]
      rw [MeasureTheory.integral_const_mul, h4]
      rw [hI] at h6
      field_simp
      ring

end
