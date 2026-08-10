module

public import Expdb.ExponentialSums.ExponentSumGrowth

import Expdb.Basic.AutomaticUniformity

/-!
# Exponential-sum bounds over compact scale ranges

The definition of `exponentSumGrowthExponent α` applies when the power-scale exponent `α` is
fixed. The first theorem below is an operational upper-semicontinuity analogue: the scale
exponent may converge to `α` while the model phase family (and hence its model order) remains
fixed. In applications, a dyadic or pigeonholed scale may instead vary over a fixed compact set;
the second theorem obtains the required uniformity by passing to a convergent subsequence.

This operational compactness statement does not assert upper semicontinuity of
`exponentSumGrowthExponent` and requires no uniformity over model orders.
-/

@[expose] public section

open Filter Topology
open scoped Expdb NNReal

noncomputable section

namespace Expdb

/-- Failure of a big-O estimate produces a subsequence on which the ratio exceeds every
successive natural-number threshold. -/
private lemma extract_not_isBigO_subseq
    {E₁ E₂ : Type*} [SeminormedAddCommGroup E₁] [SeminormedAddCommGroup E₂]
    {X : ℕ → E₁} {Y : ℕ → E₂} (h : ¬ X =O[atTop] Y) :
    ∃ φ : ℕ → ℕ, StrictMono φ ∧
      ∀ n : ℕ, (n : ℝ) * ‖Y (φ n)‖ < ‖X (φ n)‖ := by
  rw [Asymptotics.isBigO_iff] at h
  push Not at h
  have step : ∀ (n prev : ℕ), ∃ i > prev,
      (n : ℝ) * ‖Y i‖ < ‖X i‖ := fun n prev => by
    obtain ⟨i, hi, hbad⟩ :=
      (Filter.frequently_atTop.mp (h (n : ℝ))) (prev + 1)
    exact ⟨i, by omega, hbad⟩
  let φ : ℕ → ℕ :=
    fun n => n.rec (step 0 0).choose (fun n prev => (step (n + 1) prev).choose)
  exact ⟨φ, strictMono_nat_of_lt_succ fun n => (step (n + 1) (φ n)).choose_spec.1,
    fun n => n.casesOn (step 0 0).choose_spec.2
      (fun n => (step (n + 1) (φ n)).choose_spec.2)⟩

/-- Restricting a model phase family to a subsequence preserves the model phase property. -/
private lemma IsModelPhaseFunction.comp_tendsto_atTop
    {F : VariableFunction (VariableObject.fixed ℝ) ℝ} (hF : IsModelPhaseFunction F)
    {φ : ℕ → ℕ} (hφ : Tendsto φ atTop atTop) :
    IsModelPhaseFunction (fun i ↦ F (φ i)) := by
  obtain ⟨hphase, σ, hσ, herror⟩ := hF
  refine ⟨fun i ↦ hphase (φ i), σ, hσ, fun p ↦ ?_⟩
  apply (VariableFunction.isPointwiseInfinitesimal_iff_forall_pos_uniform
    (VariableObject.fixed phaseInterval)
    (fun _ ↦ ⟨1, by simp [phaseInterval]⟩)
    (modelPhaseError (fun i ↦ F (φ i)) σ p)).2
  intro ε hε
  have horiginal :=
    (VariableFunction.isPointwiseInfinitesimal_iff_forall_pos_uniform
      (VariableObject.fixed phaseInterval)
      (fun _ ↦ ⟨1, by simp [phaseInterval]⟩)
      (modelPhaseError F σ p)).1 (herror p) ε hε
  filter_upwards [hφ horiginal] with i hi
  intro u
  simpa [modelPhaseError] using hi u

/-- Fixed-phase scale-limit bound, serving as the operational upper-semicontinuity analogue for
`exponentSumGrowthExponent`: if the varying scale exponent tends to `α`, then a single model
phase family satisfies every exponent bound valid at `α`.

The phase family is fixed before the scale limit is taken, so its model order cannot vary along
the sequence. -/
theorem isPowerBounded_exponentialSum_of_scale_tendsto
    {α : ℝ≥0} {B : ℝ} (hB : exponentSumGrowthExponent α ≤ B)
    (q : VariableObject ℝ≥0)
    (N T : VariableObject ℝ)
    (F : VariableFunction (VariableObject.fixed ℝ) ℝ)
    (a b : VariableObject ℕ)
    (hq : Tendsto q atTop (𝓝 α))
    (hN : ∀ i, 1 ≤ N i)
    (hT : ∀ i, 1 ≤ T i)
    (hTunbounded : T.IsUnbounded)
    (hNT : ∀ᶠ i in atTop, N i = T i ^ (q i : ℝ))
    (hF : IsModelPhaseFunction F)
    (hab : ∀ i, N i ≤ (a i : ℝ) ∧ (b i : ℝ) ≤ 2 * N i) :
    IsPowerBounded (exponentialSum F T N a b) T B := by
  have hq' : Tendsto (fun i ↦ (q i : ℝ)) atTop (𝓝 (α : ℝ)) :=
    (NNReal.continuous_coe.tendsto α).comp hq
  have hNT' : IsPowerAsymptotic N T (α : ℝ) := by
    refine ⟨fun i ↦ (q i : ℝ), ?_, hNT⟩
    rw [IsEqUpToInfinitesimal, VariableObject.IsInfinitesimal]
    have hsub : Tendsto (fun i ↦ (q i : ℝ) - (α : ℝ)) atTop (𝓝 0) := by
      simpa using hq'.sub
        (tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ (α : ℝ)) atTop (𝓝 (α : ℝ)))
    simpa [Function.comp_def, VariableObject.fixed] using
      (continuous_norm.tendsto 0).comp hsub
  exact ((isExponentSumBound_exponentSumGrowthExponent α).mono hB)
    N T F a b hN hT hTunbounded hNT' hF hab

/-- If a varying power-scale exponent stays in a compact set on which `B` bounds `β`, then a
single fixed model phase family satisfies the exponential-sum bound with exponent `B`.

This is the compact-scale uniformity used by dyadic decompositions: any putative bad sequence
has a subsequence on which the scale exponent converges to a fixed `α`, where the defining bound
for `exponentSumGrowthExponent α` applies. -/
theorem isPowerBounded_exponentialSum_of_compact_scale
    {K : Set ℝ≥0} (hK : IsCompact K) {B : ℝ}
    (hB : ∀ α ∈ K, exponentSumGrowthExponent α ≤ B)
    (q : VariableObject ℝ≥0)
    (N T : VariableObject ℝ)
    (F : VariableFunction (VariableObject.fixed ℝ) ℝ)
    (a b : VariableObject ℕ)
    (hqK : ∀ i, q i ∈ K)
    (hN : ∀ i, 1 ≤ N i)
    (hT : ∀ i, 1 ≤ T i)
    (hTunbounded : T.IsUnbounded)
    (hNT : ∀ᶠ i in atTop, N i = T i ^ (q i : ℝ))
    (hF : IsModelPhaseFunction F)
    (hab : ∀ i, N i ≤ (a i : ℝ) ∧ (b i : ℝ) ≤ 2 * N i) :
    IsPowerBounded (exponentialSum F T N a b) T B := by
  apply (isPowerBounded_iff_forall_pos
    (exponentialSum F T N a b) T B hT hTunbounded).2
  intro ε hε
  by_contra hbigO
  obtain ⟨φ, hφ, hbad⟩ := extract_not_isBigO_subseq hbigO
  obtain ⟨α, hαK, ψ, hψ, hq⟩ := hK.tendsto_subseq (fun n ↦ hqK (φ n))
  let ρ : ℕ → ℕ := φ ∘ ψ
  have hρ : StrictMono ρ := hφ.comp hψ
  have hρtop : Tendsto ρ atTop atTop := hρ.tendsto_atTop
  let N' : VariableObject ℝ := N ∘ ρ
  let T' : VariableObject ℝ := T ∘ ρ
  let F' : VariableFunction (VariableObject.fixed ℝ) ℝ := fun i ↦ F (ρ i)
  let a' : VariableObject ℕ := a ∘ ρ
  let b' : VariableObject ℕ := b ∘ ρ
  have hT'unbounded : T'.IsUnbounded := hTunbounded.comp hρtop
  have hbound : IsPowerBounded (exponentialSum F' T' N' a' b') T' B :=
    isPowerBounded_exponentialSum_of_scale_tendsto (hB α hαK)
      (q ∘ ρ) N' T' F' a' b'
      (by simpa [ρ, Function.comp_def] using hq)
      (fun i ↦ hN (ρ i)) (fun i ↦ hT (ρ i)) hT'unbounded
      (by simpa [N', T', ρ, Function.comp_def] using hρtop hNT)
      (hF.comp_tendsto_atTop hρtop) (fun i ↦ hab (ρ i))
  have hO := (isPowerBounded_iff_forall_pos
    (exponentialSum F' T' N' a' b') T' B
    (fun i ↦ hT (ρ i)) hT'unbounded).1 hbound ε hε
  rw [Asymptotics.isBigO_iff] at hO
  obtain ⟨C, hC⟩ := hO
  have hψtop : Tendsto ψ atTop atTop := hψ.tendsto_atTop
  have hClarge : ∀ᶠ i in atTop, C < (ψ i : ℝ) := by
    exact (tendsto_natCast_atTop_atTop.comp hψtop).eventually (eventually_gt_atTop C)
  obtain ⟨i, hiC, hibound⟩ := (Filter.Eventually.and hClarge hC).exists
  have hibad := hbad (ψ i)
  have hpowpos : 0 < ‖T (φ (ψ i)) ^ (B + ε)‖ := by
    rw [Real.norm_eq_abs, abs_pos]
    exact (Real.rpow_pos_of_pos (zero_lt_one.trans_le (hT (φ (ψ i)))) _).ne'
  have hrewriteX :
      exponentialSum F' T' N' a' b' i = exponentialSum F T N a b (ρ i) := by
    rfl
  have hrewriteY : T' i ^ (B + ε) = T (ρ i) ^ (B + ε) := by rfl
  rw [hrewriteX, hrewriteY] at hibound
  rw [show ρ i = φ (ψ i) by rfl] at hibound
  exact (not_lt_of_ge hibound) ((mul_lt_mul_of_pos_right hiC hpowpos).trans hibad)

end Expdb
