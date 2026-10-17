module

public import Expdb.ExponentialSums.FixedExponentialSum
public import Mathlib.Analysis.Calculus.BumpFunction.InnerProduct
import Expdb.Mathlib.FourierCompletion

/-!
# Smooth completion for the B-process

A fixed cutoff is one on `[1, 2]` and supported inside `(3/4, 9/4)`.
Finite Fourier completion then converts an arbitrary interval sum into a
smoothly weighted sum on a larger interval, with a logarithmic loss and a
linear twist. This uses the completion machinery already developed for
blueprint Lemma `interval-set`.
-/

@[expose] public section

namespace Expdb

open scoped FourierTransform
noncomputable section

/-- The fixed smooth cutoff used to complete an exponential sum. -/
def phaseCutoff : ContDiffBump (3 / 2 : ℝ) :=
  ⟨1 / 2, 3 / 4, by norm_num, by norm_num⟩

/-- The cutoff preserves every summand in the original phase interval. -/
lemma phaseCutoff_eq_one {x : ℝ} (hx : x ∈ phaseInterval) : phaseCutoff x = 1 := by
  apply phaseCutoff.one_of_mem_closedBall
  rw [Metric.mem_closedBall, Real.dist_eq, abs_le]
  change -(1 / 2 : ℝ) ≤ x - 3 / 2 ∧ x - 3 / 2 ≤ 1 / 2
  constructor <;> linarith [hx.1, hx.2]

/-- The cutoff vanishes outside the larger positive interval used for extension. -/
lemma support_phaseCutoff :
    Function.support (phaseCutoff : ℝ → ℝ) = Set.Ioo (3 / 4 : ℝ) (9 / 4) := by
  rw [phaseCutoff.support_eq]
  ext x
  rw [Metric.mem_ball, Real.dist_eq, abs_lt]
  change (-((3 / 4 : ℝ)) < x - 3 / 2 ∧ x - 3 / 2 < 3 / 4) ↔ _
  constructor <;> rintro ⟨h₁, h₂⟩ <;> constructor <;> linarith

/-- The cutoff vanishes outside its support. -/
lemma phaseCutoff_eq_zero {x : ℝ} (hx : x ∉ Set.Ioo (3 / 4 : ℝ) (9 / 4)) :
    phaseCutoff x = 0 := by
  rw [← support_phaseCutoff] at hx
  exact Function.notMem_support.1 hx

/-- Completing a sharp sum to the fixed smooth cutoff costs only a logarithm.
The outer interval should contain the cutoff's support at scale `N`. -/
theorem exists_norm_exponentialSumAt_le_smooth_twist
    (F : ℝ → ℝ) (T : ℝ) {N : ℝ} (hN : 0 < N)
    {q A B a b : ℕ} (hq : 0 < q) (hAa : A ≤ a) (hbB : b ≤ B) (hBq : B < q)
    (hNa : N ≤ a) (hbN : (b : ℝ) ≤ 2 * N) :
    ∃ r : ℝ, 0 ≤ r ∧ r ≤ 1 ∧
      ‖exponentialSumAt F T N a b‖ ≤ (3 + 2 * Real.log q) *
        ‖∑ n ∈ Finset.Icc A B, (phaseCutoff ((n : ℝ) / N) : ℂ) *
          (𝐞 (T * F ((n : ℝ) / N) - r * n) : ℂ)‖ := by
  let f : ℕ → ℂ := fun n => (phaseCutoff ((n : ℝ) / N) : ℂ) *
    (𝐞 (T * F ((n : ℝ) / N)) : ℂ)
  have heq : (∑ n ∈ Finset.Icc a b, f n) = exponentialSumAt F T N a b := by
    apply Finset.sum_congr rfl
    intro n hn
    have hn' := Finset.mem_Icc.mp hn
    have hx : (n : ℝ) / N ∈ phaseInterval := by
      constructor
      · apply (le_div_iff₀ hN).2
        simpa only [one_mul] using hNa.trans (by exact_mod_cast hn'.1)
      · exact (div_le_iff₀ hN).2 (by
          simpa only [mul_comm] using (show (n : ℝ) ≤ b by exact_mod_cast hn'.2).trans hbN)
    simp only [f, phaseCutoff_eq_one hx, Complex.ofReal_one, one_mul, oscillatory]
  obtain ⟨r, hr, hr', hcomplete⟩ :=
    exists_norm_sum_Icc_le_log_mul_twist hq hAa hbB hBq f
  refine ⟨r, hr, hr', ?_⟩
  rw [heq] at hcomplete
  convert hcomplete using 2
  congr 1
  apply Finset.sum_congr rfl
  intro n _
  dsimp [f]
  rw [mul_assoc, ← Circle.coe_mul, ← AddChar.map_add_eq_mul, sub_eq_add_neg]

end
end Expdb
