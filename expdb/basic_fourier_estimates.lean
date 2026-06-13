import Mathlib

open MeasureTheory Real Complex Filter Topology BigOperators

noncomputable section

/-- Fourier transform: ψ̂(u) = ∫ ψ(x) e(-xu) dx -/
def fourierTransform (ψ : ℝ → ℝ) (u : ℝ) : ℂ :=
  ∫ x : ℝ, (ψ x : ℂ) * e (-(x * u))

/-- Distance from t to the boundary ∂I = {a, b} -/
def distBoundary (a b t : ℝ) : ℝ := min |t - a| |t - b|

/-- The auxiliary smooth function ψ satisfying structural constraints -/
structure SieveAux where
  ψ       : ℝ → ℝ
  smooth  : ContDiff ℝ ⊤ ψ
  supp    : Function.support ψ ⊆ Set.Icc (-1/4) (1/4)
  nonneg  : ∀ x, 0 ≤ ψ x
  l2norm  : ∫ x, ψ x ^ 2 = 1

namespace SieveAux

variable (B : SieveAux)

/-- ψ̂(u) is the Fourier transform of ψ -/
def ψ̂ (u : ℝ) : ℂ := fourierTransform B.ψ u

/-- ‖ψ̂(u)‖ = |ψ̂(u)| as a real number -/
def ψ̂Norm (u : ℝ) : ℝ := ‖B.ψ̂ u‖

lemma ψ̂_zero_pos : 0 < ‖B.ψ̂ 0‖ := by
  sorry

lemma ψ̂_rapid_decay (K : ℕ) : ∃ C : ℝ, 0 < C ∧
    ∀ u : ℝ, ‖B.ψ̂ u‖ ≤ C * (1 + |u|) ^ (-(K : ℝ)) := by
  sorry

lemma ψ̂_l2norm : ∫ u : ℝ, ‖B.ψ̂ u‖ ^ 2 = 1 := by
  sorry

end SieveAux

/-- Goal 1: WLOG we may assume ∑ |aᵣ|² = 1 -/
lemma goal1_normalize {R : ℕ} (a : Fin R → ℂ) (hM : 0 < ∑ r, ‖a r‖ ^ 2) :
    let M := ∑ r, ‖a r‖ ^ 2
    let A := fun r => a r / (Real.sqrt M : ℂ)
    ∑ r, ‖A r‖ ^ 2 = 1 := by
  intro M A
  have hM' : M ≠ 0 := ne_of_gt hM
  have hSqrt : Real.sqrt M > 0 := Real.sqrt_pos.mpr hM
  simp only [A, norm_div, Complex.norm_ofReal, abs_of_nonneg (le_of_lt hSqrt)]
  rw [Finset.sum_div, Real.sq_sqrt (le_of_lt hM)]
  exact div_self hM'

/-- Goal 2 (Equation 3.1): Plancherel Identity for the smoothed sum -/
theorem goal2_eq_3_1 (B : SieveAux) {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (N : ℝ) (hN : 0 < N) (t₀ : ℝ)
    (hnorm : ∑ r, ‖a r‖ ^ 2 = 1)
    (hsep : SeparatedFreqs N ξ) :
    ∫ t : ℝ, ‖expSum a ξ t‖ ^ 2 * B.ψ̂Norm ((t - t₀) / N) ^ 2 = N := by
  sorry

/-- Goal 3 (Equation 3.2): Interval Bound O(N) -/
theorem goal3_eq_3_2 (B : SieveAux) {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (N : ℝ) (hN : 0 < N)
    (hnorm : ∑ r, ‖a r‖ ^ 2 = 1)
    (hsep : SeparatedFreqs N ξ)
    (a₀ : ℝ) (hlen : (a₀ + N) - a₀ = N) :
    ∃ C : ℝ, 0 < C ∧
    ∫ t in Set.Icc a₀ (a₀ + N), ‖expSum a ξ t‖ ^ 2 ≤ C * N := by
  obtain ⟨Cdecay, hCdecay, _⟩ := B.ψ̂_rapid_decay 12
  exact ⟨(1 / B.ψ̂_zero_pos.le.lt_of_lt' (by linarith [B.ψ̂_zero_pos])).le,
        by linarith [B.ψ̂_zero_pos], by sorry⟩

/-- Goal 4: Integral Decomposition Identity -/
theorem goal4_identity (B : SieveAux) {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (N T : ℝ) (hN : 0 < N) (hT : 0 < T)
    (hnorm : ∑ r, ‖a r‖ ^ 2 = 1)
    (hsep : SeparatedFreqs N ξ)
    (a₀ : ℝ) :
    let I := Set.Icc a₀ (a₀ + T)
    let E : ℝ → ℝ := fun t =>
      (1/N) * ∫ t₀ in I, B.ψ̂Norm ((t - t₀) / N) ^ 2
      - intervalIndicator a₀ (a₀ + T) t
    ∫ t in I, ‖expSum a ξ t‖ ^ 2
    = T - ∫ t : ℝ, ‖expSum a ξ t‖ ^ 2 * E t := by
  sorry

/-- Goal 5: Decay bound for the error kernel E(t) -/
theorem goal5_decay (B : SieveAux)
    (N T : ℝ) (hN : 0 < N) (hT : 0 < T) (a₀ : ℝ) :
    ∃ C : ℝ, 0 < C ∧ ∀ t : ℝ,
      |(1/N) * (∫ t₀ in Set.Icc a₀ (a₀ + T), B.ψ̂Norm ((t - t₀) / N) ^ 2)
       - intervalIndicator a₀ (a₀ + T) t|
      ≤ C * (1 + distBoundary a₀ (a₀ + T) t / N) ^ (-(10 : ℝ)) := by
  obtain ⟨Cdecay, hCdecay, hdecay⟩ := B.ψ̂_rapid_decay 12
  exact ⟨Cdecay, hCdecay, fun t => by sorry⟩

/-- Goal 6: Dyadic Summation – Total error term bound is O(N) -/
theorem goal6_dyadic (B : SieveAux) {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (N T : ℝ) (hN : 0 < N) (hT : 0 < T)
    (hnorm : ∑ r, ‖a r‖ ^ 2 = 1)
    (hsep : SeparatedFreqs N ξ)
    (a₀ : ℝ) :
    let I := Set.Icc a₀ (a₀ + T)
    let E : ℝ → ℝ := fun t =>
      (1/N) * ∫ t₀ in I, B.ψ̂Norm ((t - t₀) / N) ^ 2
      - intervalIndicator a₀ (a₀ + T) t
    ∃ C : ℝ, 0 < C ∧
    |∫ t : ℝ, ‖expSum a ξ t‖ ^ 2 * E t| ≤ C * N := by
  sorry

/-- Lemma 3.1 (Main Result): L² Integral Estimate -/
theorem lemma3_1 (B : SieveAux) {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (N T : ℝ) (hN : 0 < N) (hT : 0 < T)
    (hsep : SeparatedFreqs N ξ)
    (a₀ : ℝ) :
    ∃ C : ℝ, 0 < C ∧ ∃ θ : ℝ, |θ| ≤ C ∧
    ∫ t in Set.Icc a₀ (a₀ + T), ‖expSum a ξ t‖ ^ 2
    = (T + θ * N) * ∑ r, ‖a r‖ ^ 2 := by
  set M := ∑ r, ‖a r‖ ^ 2 with hMdef
  by_cases hM : M = 0
  · refine ⟨1, one_pos, 0, by simp, ?_⟩
    have hzero : ∀ t, expSum a ξ t = 0 := by
      intro t
      apply Finset.sum_eq_zero; intro r _
      have : ‖a r‖ ^ 2 = 0 := by
        have hle : ‖a r‖ ^ 2 ≤ M :=
          Finset.single_le_sum (fun i _ => sq_nonneg ‖a i‖) Finset.univ (Finset.mem_univ r)
        linarith [sq_nonneg ‖a r‖, hM ▸ hle]
      simp [norm_eq_zero.mp (pow_eq_zero_iff (n := 2) (by norm_num) |>.mp this)]
    simp [hzero, hM]
  have hMpos : 0 < M :=
    lt_of_le_of_ne (Finset.sum_nonneg fun r _ => sq_nonneg _) (Ne.symm hM)
  set A := fun r => a r / (Real.sqrt M : ℂ) with hAdef
  have hAnorm : ∑ r, ‖A r‖ ^ 2 = 1 := goal1_normalize a hMpos
  have hscale : ∀ t, ‖expSum a ξ t‖ ^ 2 = M * ‖expSum A ξ t‖ ^ 2 := by
    intro t
    have heq : expSum a ξ t = (Real.sqrt M : ℂ) * expSum A ξ t := by
      simp [expSum, hAdef, Finset.mul_sum, mul_comm, mul_assoc,
            div_mul_cancel₀ _ (Complex.ofReal_ne_zero.mpr (Real.sqrt_pos.mpr hMpos).ne')]
    rw [heq, norm_mul, Complex.norm_ofReal, abs_of_nonneg (Real.sqrt_nonneg M),
        mul_pow, Real.sq_sqrt hMpos.le]
  have hMint : ∫ t in Set.Icc a₀ (a₀ + T), ‖expSum a ξ t‖ ^ 2
             = M * ∫ t in Set.Icc a₀ (a₀ + T), ‖expSum A ξ t‖ ^ 2 := by
    simp_rw [hscale]; exact (MeasureTheory.integral_const_mul M _).symm
  obtain ⟨C, hCpos, hErr⟩ := goal6_dyadic B A ξ N T hN hT hAnorm hsep a₀
  sorry

end
