module

public import Expdb.ExponentialSums.ExponentSumGrowthBounds
import Mathlib.Analysis.Convex.Combination

/-!
# Exponent pairs

This module formalizes the blueprint's Exponent pairs chapter (`exponent-pairs-chapter`).

It defines exponent pairs and the van der Corput transforms, proves duality with the exponential
sum growth exponent, and shows that the set of exponent pairs is closed and convex.
-/

@[expose] public section

open Filter Topology
open scoped Expdb FourierTransform NNReal

noncomputable section

namespace Expdb

/-! ## Definition -/

/-- The analytic exponent-pair bound with the blueprint's subsequence convention made
explicit: every subsequence has a further subsequence on which the variable exponents
converge to the fixed parameters `k` and `ℓ`. All data use the same composed index map.
`isExponentPairBound_iff_forall_pos` recovers the epsilon estimates on arbitrary scales. -/
def IsExponentPairBound (k l : ℝ) : Prop :=
  ∀ (N T : VariableObject ℝ)
    (F : VariableFunction (VariableObject.fixed ℝ) ℝ)
    (a b : VariableObject ℕ),
    (∀ i, 1 ≤ N i) →
    (∀ i, N i ≤ T i) →
    IsModelPhaseFunction F →
    (∀ i, IsDyadicSubinterval (N i) (a i) (b i)) →
    ∀ φ : ℕ → ℕ, StrictMono φ →
      ∃ ψ : ℕ → ℕ, StrictMono ψ ∧
        ∃ k' l' : VariableObject ℝ,
          k' =o VariableObject.fixed k ∧ l' =o VariableObject.fixed l ∧
          Asymptotics.IsBigO atTop
            (fun j ↦ exponentialSum F T N a b (φ (ψ j)))
            (fun j ↦ (T (φ (ψ j)) / N (φ (ψ j))) ^ k' j * N (φ (ψ j)) ^ l' j)

/-- The blueprint's Definition `exp-pair-def`: `(k, ℓ)` is an exponent pair if it lies in the
triangle `{0 ≤ k ≤ 1/2 ≤ ℓ ≤ 1, k + ℓ ≤ 1}` and satisfies the exponential sum bound
`IsExponentPairBound`. -/
def IsExponentPair (k l : ℝ) : Prop :=
  (0 ≤ k ∧ k ≤ 1 / 2 ∧ 1 / 2 ≤ l ∧ l ≤ 1 ∧ k + l ≤ 1) ∧ IsExponentPairBound k l

/-! ## The triangle constraints -/

/-- The first coordinate of an exponent pair is nonnegative. -/
theorem IsExponentPair.k_nonneg {k l : ℝ} (h : IsExponentPair k l) : 0 ≤ k := h.1.1

/-- The first coordinate of an exponent pair is at most `1 / 2`. -/
theorem IsExponentPair.k_le_half {k l : ℝ} (h : IsExponentPair k l) : k ≤ 1 / 2 := h.1.2.1

/-- The second coordinate of an exponent pair is at least `1 / 2`. -/
theorem IsExponentPair.half_le_l {k l : ℝ} (h : IsExponentPair k l) : 1 / 2 ≤ l := h.1.2.2.1

/-- The second coordinate of an exponent pair is at most `1`. -/
theorem IsExponentPair.l_le_one {k l : ℝ} (h : IsExponentPair k l) : l ≤ 1 := h.1.2.2.2.1

/-- The coordinates of an exponent pair sum to at most `1`. -/
theorem IsExponentPair.add_le_one {k l : ℝ} (h : IsExponentPair k l) : k + l ≤ 1 :=
  h.1.2.2.2.2

/-- The first coordinate of an exponent pair does not exceed the second. -/
theorem IsExponentPair.k_le_l {k l : ℝ} (h : IsExponentPair k l) : k ≤ l :=
  h.k_le_half.trans h.half_le_l

/-- An exponent pair satisfies the exponential sum bound. -/
theorem IsExponentPair.bound {k l : ℝ} (h : IsExponentPair k l) : IsExponentPairBound k l :=
  h.2

/-! ## Subsequences -/

/-- Being a model phase function is inherited by subsequences. -/
theorem IsModelPhaseFunction.comp_strictMono
    {F : VariableFunction (VariableObject.fixed ℝ) ℝ}
    (hF : IsModelPhaseFunction F) {φ : ℕ → ℕ} (hφ : StrictMono φ) :
    IsModelPhaseFunction (fun i ↦ F (φ i)) := by
  classical
  obtain ⟨hphase, σ, hσ, herror⟩ := hF
  refine ⟨fun i ↦ hphase (φ i), σ, hσ, ?_⟩
  intro p y
  set x : ∀ i, (VariableObject.fixed phaseInterval) i := fun j ↦
    if h : ∃ i, φ i = j then y h.choose
      else ⟨1, by simp [phaseInterval]⟩ with hxdef
  have hx : ∀ i, x (φ i) = y i := by
    intro i
    have h : ∃ i', φ i' = φ i := ⟨i, rfl⟩
    simp only [hxdef, dif_pos h]
    congr 1
    exact hφ.injective h.choose_spec
  have hlimit := (herror p x).comp hφ.tendsto_atTop
  refine hlimit.congr fun i ↦ ?_
  simp only [Function.comp_apply, hx i]
  rfl

/-- An exponential sum over an interval inside `[N, 2N]` has at most `3N` unit-sized terms. -/
theorem norm_exponentialSum_le_three_mul
    (F : VariableFunction (VariableObject.fixed ℝ) ℝ)
    (T N : VariableObject ℝ) (a b : VariableObject ℕ)
    (hN : ∀ i, 1 ≤ N i)
    (hab : ∀ i, IsDyadicSubinterval (N i) (a i) (b i))
    (i : ℕ) :
    ‖exponentialSum F T N a b i‖ ≤ 3 * N i := by
  have h := norm_exponentialSumAt_le_add_one (F i) (T i) (N i) (a i) (b i)
  have h2 := (hab i).2
  have h3 := hN i
  simpa [exponentialSum] using h.trans (by linarith)

private theorem exponentPairMajorant_shift (N T k l δ : ℝ) (hN : 0 < N) (hT : 0 < T) :
    (T / N) ^ (k + δ) * N ^ (l + δ) =
      ((T / N) ^ k * N ^ l) * T ^ δ := by
  rw [Real.rpow_add (div_pos hT hN), Real.rpow_add hN]
  calc
    (T / N) ^ k * (T / N) ^ δ * (N ^ l * N ^ δ) =
        ((T / N) ^ k * N ^ l) * ((T / N) ^ δ * N ^ δ) := by ring
    _ = ((T / N) ^ k * N ^ l) * T ^ δ := by
      rw [← Real.mul_rpow (div_pos hT hN).le hN.le, div_mul_cancel₀ T hN.ne']

/-- Vanishing exponent losses imply every fixed positive loss, without a growth
assumption on the heights. -/
private theorem exponentPairBound_isBigO
    (S : VariableObject ℂ) (N T : VariableObject ℝ) (k l : ℝ)
    (hN : ∀ i, 1 ≤ N i) (hNT : ∀ i, N i ≤ T i)
    (hS : ∃ k' l' : VariableObject ℝ,
      k' =o VariableObject.fixed k ∧ l' =o VariableObject.fixed l ∧
      Asymptotics.IsBigO atTop S (fun i ↦ (T i / N i) ^ k' i * N i ^ l' i)) :
    ∀ ε : ℝ, 0 < ε → Asymptotics.IsBigO atTop S
      (fun i ↦ (T i / N i) ^ (k + ε) * N i ^ (l + ε)) := by
  have hNpos : ∀ i, 0 < N i := fun i ↦ zero_lt_one.trans_le (hN i)
  have hT : ∀ i, 1 ≤ T i := fun i ↦ (hN i).trans (hNT i)
  have hTpos : ∀ i, 0 < T i := fun i ↦ zero_lt_one.trans_le (hT i)
  have hdiv : ∀ i, 1 ≤ T i / N i := fun i ↦ (le_div_iff₀ (hNpos i)).mpr (by simpa using hNT i)
  obtain ⟨k', l', hk', hl', hbound⟩ := hS
  intro ε hε
  have hk := (isEqUpToInfinitesimal_iff_forall_pos k' (VariableObject.fixed k)).mp hk' ε hε
  have hl := (isEqUpToInfinitesimal_iff_forall_pos l' (VariableObject.fixed l)).mp hl' ε hε
  refine hbound.trans (Asymptotics.IsBigO.of_bound 1 ?_)
  filter_upwards [hk, hl] with i hki hli
  simp only [Real.norm_eq_abs, abs_of_nonneg (mul_nonneg
    (Real.rpow_nonneg (div_pos (hTpos i) (hNpos i)).le _) (Real.rpow_nonneg (hNpos i).le _)),
    one_mul]
  rw [Real.norm_eq_abs] at hki hli
  apply mul_le_mul
    (Real.rpow_le_rpow_of_exponent_le (hdiv i) (by
      dsimp [VariableObject.fixed] at hki; linarith [le_abs_self (k' i - k)]))
    (Real.rpow_le_rpow_of_exponent_le (hN i) (by
      dsimp [VariableObject.fixed] at hli; linarith [le_abs_self (l' i - l)]))
    (Real.rpow_nonneg (hNpos i).le _) (Real.rpow_nonneg (div_pos (hTpos i) (hNpos i)).le _)

/-- On divergent heights, vanishing exponent losses are equivalent to every fixed
positive exponent loss. The two variable exponents need not have the same error. -/
private theorem exponentPairBound_iff_forall_pos
    (S : VariableObject ℂ) (N T : VariableObject ℝ) (k l : ℝ)
    (hN : ∀ i, 1 ≤ N i) (hNT : ∀ i, N i ≤ T i) (hTunbounded : T.IsUnbounded) :
    (∃ k' l' : VariableObject ℝ,
      k' =o VariableObject.fixed k ∧ l' =o VariableObject.fixed l ∧
      Asymptotics.IsBigO atTop S (fun i ↦ (T i / N i) ^ k' i * N i ^ l' i)) ↔
    ∀ ε : ℝ, 0 < ε → Asymptotics.IsBigO atTop S
      (fun i ↦ (T i / N i) ^ (k + ε) * N i ^ (l + ε)) := by
  have hNpos : ∀ i, 0 < N i := fun i ↦ zero_lt_one.trans_le (hN i)
  have hT : ∀ i, 1 ≤ T i := fun i ↦ (hN i).trans (hNT i)
  have hTpos : ∀ i, 0 < T i := fun i ↦ zero_lt_one.trans_le (hT i)
  constructor
  · intro hS
    exact exponentPairBound_isBigO S N T k l hN hNT hS
  · intro hS
    let B : VariableObject ℝ := fun i ↦ (T i / N i) ^ k * N i ^ l
    have hBpos : ∀ i, 0 < B i := fun i ↦ mul_pos
      (Real.rpow_pos_of_pos (div_pos (hTpos i) (hNpos i)) _) (Real.rpow_pos_of_pos (hNpos i) _)
    let X : VariableObject ℝ := fun i ↦ ‖S i‖ / B i
    have hX : IsPowerBounded X T 0 := by
      apply (isPowerBounded_iff_forall_pos X T 0 hT hTunbounded).mpr
      intro ε hε
      obtain ⟨C, hCpos, hC⟩ := (hS ε hε).exists_pos
      refine Asymptotics.IsBigO.of_bound C ?_
      filter_upwards [hC.bound] with i hi
      have hshift := exponentPairMajorant_shift (N i) (T i) k l ε (hNpos i) (hTpos i)
      simp only [Real.norm_eq_abs, abs_of_nonneg (mul_nonneg
        (Real.rpow_nonneg (div_pos (hTpos i) (hNpos i)).le _)
        (Real.rpow_nonneg (hNpos i).le _))] at hi
      simp only [X, Real.norm_eq_abs, abs_of_nonneg (div_nonneg (norm_nonneg _) (hBpos i).le),
        zero_add, abs_of_nonneg (Real.rpow_nonneg (hTpos i).le _)]
      apply (div_le_iff₀ (hBpos i)).mpr
      rw [hshift] at hi
      dsimp [B]
      nlinarith [hi]
    obtain ⟨δ, hδ, hbound⟩ := hX
    refine ⟨fun i ↦ k + δ i, fun i ↦ l + δ i, ?_, ?_, ?_⟩
    · rw [isEqUpToInfinitesimal_iff_tendsto_sub] at hδ ⊢
      simpa [VariableObject.fixed] using hδ
    · rw [isEqUpToInfinitesimal_iff_tendsto_sub] at hδ ⊢
      simpa [VariableObject.fixed] using hδ
    · obtain ⟨C, hCpos, hC⟩ := hbound.exists_pos
      refine Asymptotics.IsBigO.of_bound C ?_
      filter_upwards [hC.bound] with i hi
      simp only [X, Real.norm_eq_abs, abs_of_nonneg (div_nonneg (norm_nonneg _) (hBpos i).le),
        abs_of_nonneg (Real.rpow_nonneg (hTpos i).le _)] at hi
      simp only [Real.norm_eq_abs, abs_of_nonneg (mul_nonneg
        (Real.rpow_nonneg (div_pos (hTpos i) (hNpos i)).le _) (Real.rpow_nonneg (hNpos i).le _))]
      rw [exponentPairMajorant_shift (N i) (T i) k l (δ i) (hNpos i) (hTpos i)]
      have := (div_le_iff₀ (hBpos i)).mp hi
      dsimp [B] at this
      nlinarith [this]

/-- The subsequence formulation of the vanishing-exponent bound is equivalent to epsilon
bounds on all scale sequences, including bounded or oscillating heights. -/
theorem isExponentPairBound_iff_forall_pos {k l : ℝ} :
    IsExponentPairBound k l ↔
    ∀ (N T : VariableObject ℝ)
      (F : VariableFunction (VariableObject.fixed ℝ) ℝ) (a b : VariableObject ℕ),
      (∀ i, 1 ≤ N i) → (∀ i, N i ≤ T i) → IsModelPhaseFunction F →
      (∀ i, IsDyadicSubinterval (N i) (a i) (b i)) →
      ∀ ε : ℝ, 0 < ε → Asymptotics.IsBigO atTop (exponentialSum F T N a b)
        (fun i ↦ (T i / N i) ^ (k + ε) * N i ^ (l + ε)) := by
  constructor
  · intro h N T F a b hN hNT hF hab ε hε
    by_contra hO
    have hfreq : ∀ n : ℕ, ∃ᶠ i in atTop,
        (n : ℝ) * ‖(T i / N i) ^ (k + ε) * N i ^ (l + ε)‖ <
          ‖exponentialSum F T N a b i‖ := by
      intro n
      rw [Asymptotics.isBigO_iff] at hO
      push Not at hO
      simpa using hO n
    obtain ⟨φ, hφ, hφprop⟩ := extraction_forall_of_frequently hfreq
    obtain ⟨ψ, hψ, k', l', hk', hl', hbound⟩ := h N T F a b hN hNT hF hab φ hφ
    have hsubε := exponentPairBound_isBigO
      (fun j ↦ exponentialSum F T N a b (φ (ψ j)))
      (N ∘ φ ∘ ψ) (T ∘ φ ∘ ψ) k l (fun j ↦ hN _) (fun j ↦ hNT _)
      ⟨k', l', hk', hl', hbound⟩ ε hε
    obtain ⟨C, hCpos, hC⟩ := hsubε.exists_pos
    have hcontra : ∀ᶠ j : ℕ in atTop, (ψ j : ℝ) < C := by
      filter_upwards [hC.bound] with j hj
      have hf := hφprop (ψ j)
      have hNpos : 0 < N (φ (ψ j)) := zero_lt_one.trans_le (hN _)
      have hTpos : 0 < T (φ (ψ j)) := hNpos.trans_le (hNT _)
      have hRpos : 0 < ‖(T (φ (ψ j)) / N (φ (ψ j))) ^ (k + ε) *
          N (φ (ψ j)) ^ (l + ε)‖ := by
        exact norm_pos_iff.mpr (mul_pos
          (Real.rpow_pos_of_pos (div_pos hTpos hNpos) _)
          (Real.rpow_pos_of_pos hNpos _)).ne'
      change ‖exponentialSum F T N a b (φ (ψ j))‖ ≤
        C * ‖(T (φ (ψ j)) / N (φ (ψ j))) ^ (k + ε) *
          N (φ (ψ j)) ^ (l + ε)‖ at hj
      exact lt_of_mul_lt_mul_right (hf.trans_le hj) hRpos.le
    have hlarge : ∀ᶠ j : ℕ in atTop, C ≤ (ψ j : ℝ) :=
      (tendsto_natCast_atTop_atTop.comp hψ.tendsto_atTop).eventually (eventually_ge_atTop C)
    obtain ⟨j, hj, hjC⟩ := (hcontra.and hlarge).exists
    exact (not_lt_of_ge hjC) hj
  · intro h N T F a b hN hNT hF hab φ hφ
    have hNpos : ∀ i, 0 < N i := fun i ↦ zero_lt_one.trans_le (hN i)
    have hT : ∀ i, 1 ≤ T i := fun i ↦ (hN i).trans (hNT i)
    have hTpos : ∀ i, 0 < T i := fun i ↦ zero_lt_one.trans_le (hT i)
    let B : VariableObject ℝ := fun i ↦ (T i / N i) ^ k * N i ^ l
    have hBpos : ∀ i, 0 < B i := fun i ↦ mul_pos
      (Real.rpow_pos_of_pos (div_pos (hTpos i) (hNpos i)) _) (Real.rpow_pos_of_pos (hNpos i) _)
    let X : VariableObject ℝ := fun i ↦ ‖exponentialSum F T N a b i‖ / B i
    have hXbound : ∀ i, X i ≤ 3 * T i ^ (1 + |k| + |l|) := by
      intro i
      have hdiv1 : 1 ≤ T i / N i := (le_div_iff₀ (hNpos i)).mpr (by simpa using hNT i)
      have hdivT : T i / N i ≤ T i := by
        apply (div_le_iff₀ (hNpos i)).mpr
        nlinarith [hN i, hTpos i]
      have hkbound : (T i / N i) ^ (-k) ≤ T i ^ |k| :=
        (Real.rpow_le_rpow_of_exponent_le hdiv1 (neg_le_abs k)).trans
          (Real.rpow_le_rpow (div_pos (hTpos i) (hNpos i)).le hdivT (abs_nonneg k))
      have hlbound : N i ^ (-l) ≤ T i ^ |l| :=
        (Real.rpow_le_rpow_of_exponent_le (hN i) (neg_le_abs l)).trans
          (Real.rpow_le_rpow (hNpos i).le (hNT i) (abs_nonneg l))
      have hSbound : ‖exponentialSum F T N a b i‖ ≤ 3 * T i :=
        (norm_exponentialSum_le_three_mul F T N a b hN hab i).trans (by linarith [hNT i])
      have hinv : (B i)⁻¹ = (T i / N i) ^ (-k) * N i ^ (-l) := by
        dsimp [B]
        rw [mul_inv_rev, Real.rpow_neg (div_pos (hTpos i) (hNpos i)).le,
          Real.rpow_neg (hNpos i).le]
        ring
      calc
        X i = ‖exponentialSum F T N a b i‖ *
            ((T i / N i) ^ (-k) * N i ^ (-l)) := by
          dsimp [X]
          rw [div_eq_mul_inv, hinv]
        _ ≤ (3 * T i) * (T i ^ |k| * T i ^ |l|) := by
          exact mul_le_mul hSbound
            (mul_le_mul hkbound hlbound (Real.rpow_nonneg (hNpos i).le _)
              (Real.rpow_nonneg (hTpos i).le _))
            (mul_nonneg (Real.rpow_nonneg (div_pos (hTpos i) (hNpos i)).le _)
              (Real.rpow_nonneg (hNpos i).le _)) (mul_nonneg (by norm_num) (hTpos i).le)
        _ = 3 * T i ^ (1 + |k| + |l|) := by
          rw [Real.rpow_add (hTpos i), Real.rpow_add (hTpos i), Real.rpow_one]
          ring
    by_cases htop : Tendsto (T ∘ φ) atTop atTop
    · refine ⟨id, strictMono_id, ?_⟩
      have hTunbounded : VariableObject.IsUnbounded (T ∘ φ) :=
        (VariableObject.isUnbounded_iff_tendsto_atTop
          (fun j ↦ (hTpos (φ j)).le)).mpr htop
      apply (exponentPairBound_iff_forall_pos
        (fun j ↦ exponentialSum F T N a b (φ j)) (N ∘ φ) (T ∘ φ) k l
        (fun j ↦ hN _) (fun j ↦ hNT _) hTunbounded).mpr
      intro ε hε
      exact (h N T F a b hN hNT hF hab ε hε).comp_tendsto hφ.tendsto_atTop
    · have hfreq : ∃ M : ℝ, ∃ᶠ j in atTop, T (φ j) < M := by
        simpa only [Filter.tendsto_atTop, not_forall, Filter.not_eventually, not_le,
          Function.comp_apply]
          using htop
      obtain ⟨M, hM⟩ := hfreq
      obtain ⟨ψ, hψ, hψbound⟩ := extraction_of_frequently_atTop hM
      refine ⟨ψ, hψ, VariableObject.fixed k, VariableObject.fixed l,
        IsEqUpToInfinitesimal.refl _, IsEqUpToInfinitesimal.refl _, ?_⟩
      let Q : ℝ := max M 1
      have hMQ : M ≤ Q := le_max_left _ _
      refine Asymptotics.IsBigO.of_bound (3 * Q ^ (1 + |k| + |l|)) ?_
      filter_upwards [] with j
      have hsmall : T (φ (ψ j)) ≤ Q := (hψbound j).le.trans hMQ
      have hb := Real.rpow_le_rpow (hTpos (φ (ψ j))).le hsmall
        (show 0 ≤ 1 + |k| + |l| by positivity)
      have hx : X (φ (ψ j)) ≤ 3 * Q ^ (1 + |k| + |l|) :=
        (hXbound _).trans (mul_le_mul_of_nonneg_left hb (by norm_num))
      have hS := (div_le_iff₀ (hBpos (φ (ψ j)))).mp hx
      simpa only [B, VariableObject.fixed, Real.norm_eq_abs,
        abs_of_pos (hBpos (φ (ψ j)))] using hS

/-! ## Consequences of the exponential sum bound -/

/-- The exponential sum bound may be applied when `N ≤ T` only holds eventually. -/
theorem IsExponentPairBound.isBigO_of_eventually_le
    {k l : ℝ} (h : IsExponentPairBound k l) {N T : VariableObject ℝ}
    {F : VariableFunction (VariableObject.fixed ℝ) ℝ} {a b : VariableObject ℕ}
    (hN : ∀ i, 1 ≤ N i) (hT : ∀ i, 1 ≤ T i) (hNT : ∀ᶠ i in atTop, N i ≤ T i)
    (hF : IsModelPhaseFunction F)
    (hab : ∀ i, IsDyadicSubinterval (N i) (a i) (b i))
    {ε : ℝ} (hε : 0 < ε) :
    Asymptotics.IsBigO atTop (exponentialSum F T N a b)
      (fun i ↦ (T i / N i) ^ (k + ε) * N i ^ (l + ε)) := by
  classical
  set N' : VariableObject ℝ := fun i ↦ if N i ≤ T i then N i else T i with hN'def
  set a' : VariableObject ℕ := fun i ↦ if N i ≤ T i then a i else ⌈T i⌉₊ with ha'def
  set b' : VariableObject ℕ := fun i ↦ if N i ≤ T i then b i else 0 with hb'def
  have hN'1 : ∀ i, 1 ≤ N' i := by
    intro i
    by_cases hi : N i ≤ T i <;> simp [hN'def, hi, hN i, hT i]
  have hN'T : ∀ i, N' i ≤ T i := by
    intro i
    by_cases hi : N i ≤ T i <;> simp [hN'def, hi]
  have hab' : ∀ i, IsDyadicSubinterval (N' i) (a' i) (b' i) := by
    intro i
    by_cases hi : N i ≤ T i
    · simpa [hN'def, ha'def, hb'def, hi] using hab i
    · refine ⟨?_, ?_⟩
      · simpa [hN'def, ha'def, hi] using Nat.le_ceil (T i)
      · have := hT i
        simp only [hN'def, hb'def, hi, if_false, Nat.cast_zero]
        linarith
  have hmain := isExponentPairBound_iff_forall_pos.mp h N' T F a' b' hN'1 hN'T hF hab' ε hε
  refine hmain.congr' ?_ ?_ <;> filter_upwards [hNT] with i hi
  · simp [exponentialSum, hN'def, ha'def, hb'def, hi]
  · simp [hN'def, hi]

/-! ## Duality with the exponential sum growth exponent -/

/-- Substitute the scale relation into the exponent-pair majorant. -/
private theorem exponentPairMajorant_eq_rpow {N T : ℝ} (hN : 0 < N) (hT : 0 < T)
    (k l a : ℝ) (hNa : N = T ^ a) :
    (T / N) ^ k * N ^ l = T ^ (k + (l - k) * a) := by
  rw [Real.div_rpow hT.le hN.le, div_mul_eq_mul_div, mul_div_assoc,
    ← Real.rpow_sub hN, hNa, ← Real.rpow_mul hT.le, ← Real.rpow_add hT]
  congr 1
  ring

/-- Rewrite a power bound at a fixed scale as an exponent-pair estimate, keeping the
vanishing errors as variable exponents rather than introducing epsilon losses. -/
private theorem exponentPairBound_of_isPowerBounded
    {S : VariableObject ℂ} {N T : VariableObject ℝ} {α k l : ℝ}
    (hN : ∀ i, 1 ≤ N i) (hT : ∀ i, 1 ≤ T i)
    (hNT : IsPowerAsymptotic N T α)
    (hS : IsPowerBounded S T (k + (l - k) * α)) :
    ∃ k' l' : VariableObject ℝ,
      k' =o VariableObject.fixed k ∧ l' =o VariableObject.fixed l ∧
      Asymptotics.IsBigO atTop S (fun i ↦ (T i / N i) ^ k' i * N i ^ l' i) := by
  obtain ⟨a, ha, hNa⟩ := hNT
  obtain ⟨e, he, hbound⟩ := hS
  let δ : VariableObject ℝ := fun i ↦ e i - (k + (l - k) * a i)
  have hδ : δ =o VariableObject.fixed 0 := by
    have h := he.sub ((IsEqUpToInfinitesimal.refl (VariableObject.fixed k)).add
      (ha.const_smul (l - k)))
    convert h using 1 <;> ext i <;> simp [δ, VariableObject.fixed, smul_eq_mul]
  refine ⟨fun i ↦ k + δ i, fun i ↦ l + δ i, ?_, ?_, ?_⟩
  · convert (IsEqUpToInfinitesimal.refl (VariableObject.fixed k)).add hδ using 1 <;>
      ext i <;> simp [VariableObject.fixed]
  · convert (IsEqUpToInfinitesimal.refl (VariableObject.fixed l)).add hδ using 1 <;>
      ext i <;> simp [VariableObject.fixed]
  · apply hbound.congr' Filter.EventuallyEq.rfl
    filter_upwards [hNa] with i hi
    rw [exponentPairMajorant_eq_rpow (zero_lt_one.trans_le (hN i))
      (zero_lt_one.trans_le (hT i)) (k + δ i) (l + δ i) (a i) hi]
    congr 1
    dsimp [δ]
    ring

/-- Half of the blueprint's Lemma `beta-duality`: an exponent pair gives the linear bound
`β(α) ≤ k + (ℓ - k) α`, here in the range `α < 1`. -/
theorem IsExponentPairBound.isExponentSumBound
    {k l : ℝ} (h : IsExponentPairBound k l)
    {α : ℝ≥0} (hα : (α : ℝ) < 1) :
    IsExponentSumBound α (k + (l - k) * α) := by
  intro N T F a b hN hT hTunbounded hNT hF hab
  -- For `α < 1`, the scale relation implies `N ≤ T` after discarding finitely many indices.
  have hbetween := hNT.eventually_between (Eventually.of_forall hT)
    (show 0 < (1 - (α : ℝ)) / 2 by linarith)
  have hNTle : ∀ᶠ i in atTop, N i ≤ T i := by
    filter_upwards [hbetween] with i hi
    calc N i ≤ T i ^ ((α : ℝ) + (1 - α) / 2) := hi.2
      _ ≤ T i ^ (1 : ℝ) := Real.rpow_le_rpow_of_exponent_le (hT i) (by linarith)
      _ = T i := Real.rpow_one _
  -- Fill the exceptional indices with admissible data; the eventual sum is unchanged.
  let N' : VariableObject ℝ := fun i ↦ if N i ≤ T i then N i else T i
  let a' : VariableObject ℕ := fun i ↦ if N i ≤ T i then a i else ⌈T i⌉₊
  let b' : VariableObject ℕ := fun i ↦ if N i ≤ T i then b i else 0
  have hN' : ∀ i, 1 ≤ N' i := by
    intro i
    by_cases hi : N i ≤ T i <;> simp [N', hi, hN i, hT i]
  have hN'T : ∀ i, N' i ≤ T i := by
    intro i
    by_cases hi : N i ≤ T i <;> simp [N', hi]
  have hab' : ∀ i, IsDyadicSubinterval (N' i) (a' i) (b' i) := by
    intro i
    by_cases hi : N i ≤ T i
    · simpa [N', a', b', hi] using hab i
    · simp only [N', a', b', if_neg hi]
      exact ⟨Nat.le_ceil _, by simp; linarith [hT i]⟩
  apply isPowerBounded_of_subsequences hT hTunbounded
  intro φ hφ
  obtain ⟨ψ, hψ, k', l', hk', hl', hbound⟩ := h N' T F a' b' hN' hN'T hF hab' φ hφ
  obtain ⟨c, hc, hNc⟩ := hNT.comp_tendsto (hφ.comp hψ).tendsto_atTop
  -- Substitute `N = T^(α + o(1))` into `(T/N)^(k + o(1)) N^(ℓ + o(1))`.
  let e : VariableObject ℝ := fun j ↦ k' j + (l' j - k' j) * c j
  have he := hk'.add ((hl'.sub hk').mul hc
    ⟨‖l - k‖, Eventually.of_forall fun _ ↦ le_rfl⟩
    ⟨‖(α : ℝ)‖, Eventually.of_forall fun _ ↦ le_rfl⟩)
  refine ⟨ψ, hψ, e, ?_, ?_⟩
  · convert he using 1 <;> ext j <;> rfl
  · have hgood := (hφ.comp hψ).tendsto_atTop.eventually hNTle
    refine hbound.congr' ?_ ?_
    · filter_upwards [hgood] with j hj
      simp only [exponentialSum, Function.comp_apply, N', a', b']
      dsimp [Function.comp_def] at hj
      simp only [if_pos hj]
    · filter_upwards [hgood, hNc] with j hj hjN
      dsimp [Function.comp_def] at hj hjN ⊢
      simp only [N', if_pos hj]
      exact exponentPairMajorant_eq_rpow (zero_lt_one.trans_le (hN _))
        (zero_lt_one.trans_le (hT _)) (k' j) (l' j) (c j) hjN

/-- An exponent pair gives the linear bound `β(α) ≤ k + (ℓ - k) α` on `[0, 1]`. -/
theorem IsExponentPair.exponentSumGrowthExponent_le
    {k l : ℝ} (h : IsExponentPair k l) {α : ℝ≥0} (hα : α ≤ 1) :
    exponentSumGrowthExponent α ≤ k + (l - k) * α := by
  -- The endpoint values of `β` and the triangle constraints settle `α = 0, 1`.
  by_cases hzero : α = 0
  · subst α
    simpa using h.k_nonneg
  by_cases hone : α = 1
  · subst α
    simpa using h.half_le_l
  exact exponentSumGrowthExponent_le_iff.mpr
    (h.bound.isExponentSumBound (by exact_mod_cast lt_of_le_of_ne hα hone))

/-- The converse half of the blueprint's Lemma `beta-duality`: the linear bound
`β(α) ≤ k + (ℓ - k) α` on `[0, 1]` implies the exponential sum bound. -/
theorem isExponentPairBound_of_forall_exponentSumGrowthExponent_le
    {k l : ℝ} (hk : 0 ≤ k) (hkl : k ≤ l)
    (hβ : ∀ α : ℝ≥0, α ≤ 1 → exponentSumGrowthExponent α ≤ k + (l - k) * α) :
    IsExponentPairBound k l := by
  intro N T F a b hN hNT hF hab φ hφ
  have hT1 : ∀ i, 1 ≤ T i := fun i ↦ (hN i).trans (hNT i)
  by_cases hNtop : Tendsto (N ∘ φ) atTop atTop
  · have hTtop : Tendsto (fun j ↦ T (φ j)) atTop atTop :=
      tendsto_atTop_mono (fun j ↦ hNT (φ j)) hNtop
    -- On growing scales, pass to a further subsequence with a fixed scale exponent.
    set x : ℕ → ℝ := fun j ↦ max 0 (min 1 (Real.logb (T (φ j)) (N (φ j)))) with hxdef
    have hxmem : ∀ j, x j ∈ Set.Icc (0 : ℝ) 1 := by
      intro j
      refine ⟨le_max_left _ _, max_le zero_le_one (min_le_left _ _)⟩
    obtain ⟨α0, hα0mem, ψ, hψ, hlimit⟩ := isCompact_Icc.tendsto_subseq hxmem
    set φ2 : ℕ → ℕ := fun j ↦ φ (ψ j) with hφ2def
    have hφ2 : StrictMono φ2 := hφ.comp hψ
    set N₂ : VariableObject ℝ := fun j ↦ N (φ2 j) with hN₂def
    set T₂ : VariableObject ℝ := fun j ↦ T (φ2 j) with hT₂def
    set F₂ : VariableFunction (VariableObject.fixed ℝ) ℝ := fun j ↦ F (φ2 j) with hF₂def
    set a₂ : VariableObject ℕ := fun j ↦ a (φ2 j) with ha₂def
    set b₂ : VariableObject ℕ := fun j ↦ b (φ2 j) with hb₂def
    have hN₂1 : ∀ j, 1 ≤ N₂ j := fun j ↦ hN _
    have hT₂1 : ∀ j, 1 ≤ T₂ j := fun j ↦ hT1 _
    have hT₂top : Tendsto T₂ atTop atTop := hTtop.comp hψ.tendsto_atTop
    have hT₂unbounded : T₂.IsUnbounded :=
      (VariableObject.isUnbounded_iff_tendsto_atTop
        (fun j ↦ zero_le_one.trans (hT₂1 j))).2 hT₂top
    have hT₂gt : ∀ᶠ j in atTop, 1 < T₂ j := hT₂top.eventually (eventually_gt_atTop 1)
    have hlogb : Tendsto (fun j ↦ Real.logb (T₂ j) (N₂ j)) atTop (𝓝 α0) := by
      refine hlimit.congr' ?_
      filter_upwards [hT₂gt] with j hj
      have hlogT : 0 < Real.log (T₂ j) := Real.log_pos hj
      have hlogN : 0 ≤ Real.log (N₂ j) := Real.log_nonneg (hN₂1 j)
      have hle : Real.log (N₂ j) ≤ Real.log (T₂ j) :=
        Real.log_le_log (lt_of_lt_of_le zero_lt_one (hN₂1 j)) (hNT _)
      have h0 : 0 ≤ Real.logb (T₂ j) (N₂ j) := by
        rw [Real.logb]; positivity
      have h1 : Real.logb (T₂ j) (N₂ j) ≤ 1 := by
        rw [Real.logb, div_le_one hlogT]; exact hle
      simp only [Function.comp_apply, hxdef]
      change max 0 (min 1 (Real.logb (T₂ j) (N₂ j))) = Real.logb (T₂ j) (N₂ j)
      rw [min_eq_right h1, max_eq_right h0]
    have hpa : IsPowerAsymptotic N₂ T₂ α0 :=
      isPowerAsymptotic_of_logb_tendsto hT₂gt
        (Filter.Eventually.of_forall fun j ↦ lt_of_lt_of_le zero_lt_one (hN₂1 j)) hlogb
    set αnn : ℝ≥0 := ⟨α0, hα0mem.1⟩ with hαnndef
    have hαnncoe : (αnn : ℝ) = α0 := rfl
    have hαnnle : αnn ≤ 1 := by
      rw [← NNReal.coe_le_coe, hαnncoe, NNReal.coe_one]
      exact hα0mem.2
    have hsumbound : IsExponentSumBound αnn (k + (l - k) * α0) :=
      exponentSumGrowthExponent_le_iff.1 (hβ αnn hαnnle)
    have hpb := hsumbound N₂ T₂ F₂ a₂ b₂ hN₂1 hT₂1 hT₂unbounded hpa
      (hF.comp_strictMono hφ2) (fun j ↦ hab _)
    refine ⟨ψ, hψ, ?_⟩
    exact exponentPairBound_of_isPowerBounded hN₂1 hT₂1 hpa hpb
  · -- If the lengths do not tend to infinity, extract bounded lengths and use the trivial bound.
    have hfreq : ∃ M : ℝ, ∃ᶠ j in atTop, N (φ j) < M := by
      simpa only [Filter.tendsto_atTop, not_forall, Filter.not_eventually, not_le,
        Function.comp_apply] using hNtop
    obtain ⟨M, hM⟩ := hfreq
    obtain ⟨ψ, hψ, hψbound⟩ := extraction_of_frequently_atTop hM
    refine ⟨ψ, hψ, VariableObject.fixed k, VariableObject.fixed l,
      IsEqUpToInfinitesimal.refl _, IsEqUpToInfinitesimal.refl _, ?_⟩
    let Q : ℝ := max M 1
    have hQ : 1 ≤ Q := le_max_right _ _
    refine Asymptotics.IsBigO.of_bound (3 * Q) ?_
    filter_upwards [] with j
    have hNpos : 0 < N (φ (ψ j)) := zero_lt_one.trans_le (hN _)
    have hTpos : 0 < T (φ (ψ j)) := hNpos.trans_le (hNT _)
    have hdiv : 1 ≤ T (φ (ψ j)) / N (φ (ψ j)) :=
      (le_div_iff₀ hNpos).mpr (by simpa using hNT (φ (ψ j)))
    have hkpow := Real.one_le_rpow hdiv hk
    have hlpow := Real.one_le_rpow (hN (φ (ψ j))) (hk.trans hkl)
    have hmajor : 1 ≤ (T (φ (ψ j)) / N (φ (ψ j))) ^ k * N (φ (ψ j)) ^ l := by
      nlinarith
    have hS := norm_exponentialSum_le_three_mul F T N a b hN hab (φ (ψ j))
    have hNQ : N (φ (ψ j)) ≤ Q := (hψbound j).le.trans (le_max_left _ _)
    simp only [VariableObject.fixed, Real.norm_eq_abs,
      abs_of_pos (mul_pos (Real.rpow_pos_of_pos (div_pos hTpos hNpos) _)
        (Real.rpow_pos_of_pos hNpos _))]
    nlinarith

/-- **Duality between exponent pairs and `β`** (blueprint Lemma `beta-duality`). A point
`(k, ℓ)` of the triangle `{0 ≤ k ≤ 1/2 ≤ ℓ ≤ 1, k + ℓ ≤ 1}` is an exponent pair if and only if
`β(α) ≤ k + (ℓ - k) α` for all `0 ≤ α ≤ 1`. -/
theorem isExponentPair_iff_forall_exponentSumGrowthExponent_le
    {k l : ℝ} (htriangle : 0 ≤ k ∧ k ≤ 1 / 2 ∧ 1 / 2 ≤ l ∧ l ≤ 1 ∧ k + l ≤ 1) :
    IsExponentPair k l ↔
      ∀ α : ℝ≥0, α ≤ 1 → exponentSumGrowthExponent α ≤ k + (l - k) * α := by
  obtain ⟨hk, hk2, hl2, hl1, hkl⟩ := htriangle
  refine ⟨fun h _ hα ↦ h.exponentSumGrowthExponent_le hα, fun hβ ↦
    ⟨⟨hk, hk2, hl2, hl1, hkl⟩,
      isExponentPairBound_of_forall_exponentSumGrowthExponent_le hk
        (hk2.trans hl2) hβ⟩⟩

/-! ## The set of exponent pairs -/

/-- The set of exponent pairs, as a subset of the plane. -/
def exponentPairs : Set (ℝ × ℝ) := {p | IsExponentPair p.1 p.2}

/-- Membership in the set of exponent pairs. -/
@[simp] theorem mem_exponentPairs {p : ℝ × ℝ} :
    p ∈ exponentPairs ↔ IsExponentPair p.1 p.2 := Iff.rfl

/-- The set of exponent pairs is the intersection of the triangle with the half-planes coming
from the linear bounds on `β`. -/
theorem exponentPairs_eq_inter :
    exponentPairs =
      ({p : ℝ × ℝ | 0 ≤ p.1} ∩ {p : ℝ × ℝ | p.1 ≤ 1 / 2} ∩ {p : ℝ × ℝ | 1 / 2 ≤ p.2} ∩
        {p : ℝ × ℝ | p.2 ≤ 1} ∩ {p : ℝ × ℝ | p.1 + p.2 ≤ 1}) ∩
        ⋂ α : {α : ℝ≥0 // α ≤ 1},
          {p : ℝ × ℝ | exponentSumGrowthExponent α ≤ p.1 + (p.2 - p.1) * ((α : ℝ≥0) : ℝ)} := by
  ext p
  simp only [mem_exponentPairs, Set.mem_inter_iff, Set.mem_setOf_eq, Set.mem_iInter]
  constructor
  · intro hp
    exact ⟨⟨⟨⟨⟨hp.k_nonneg, hp.k_le_half⟩, hp.half_le_l⟩, hp.l_le_one⟩, hp.add_le_one⟩,
      fun α ↦ hp.exponentSumGrowthExponent_le α.2⟩
  · rintro ⟨⟨⟨⟨⟨h1, h2⟩, h3⟩, h4⟩, h5⟩, h6⟩
    exact (isExponentPair_iff_forall_exponentSumGrowthExponent_le
      ⟨h1, h2, h3, h4, h5⟩).2 fun α hα ↦ h6 ⟨α, hα⟩

/-- The set of exponent pairs is closed (blueprint Corollary `exp-pair-closed`). -/
theorem isClosed_exponentPairs : IsClosed exponentPairs := by
  rw [exponentPairs_eq_inter]
  refine IsClosed.inter ?_ (isClosed_iInter fun α ↦ ?_)
  · exact ((((isClosed_le continuous_const continuous_fst).inter
      (isClosed_le continuous_fst continuous_const)).inter
      (isClosed_le continuous_const continuous_snd)).inter
      (isClosed_le continuous_snd continuous_const)).inter
      (isClosed_le (continuous_fst.add continuous_snd) continuous_const)
  · exact isClosed_le continuous_const
      (continuous_fst.add ((continuous_snd.sub continuous_fst).mul continuous_const))

/-- The set of exponent pairs is convex (blueprint Corollary `exp-pair-closed`). -/
theorem convex_exponentPairs : Convex ℝ exponentPairs := by
  have hfst : IsLinearMap ℝ (fun p : ℝ × ℝ ↦ p.1) := ⟨fun _ _ ↦ rfl, fun _ _ ↦ rfl⟩
  have hsnd : IsLinearMap ℝ (fun p : ℝ × ℝ ↦ p.2) := ⟨fun _ _ ↦ rfl, fun _ _ ↦ rfl⟩
  have hadd : IsLinearMap ℝ (fun p : ℝ × ℝ ↦ p.1 + p.2) :=
    ⟨fun p q ↦ by simp [Prod.fst_add, Prod.snd_add]; ring,
      fun c p ↦ by simp [Prod.smul_def]; ring⟩
  have hlin : ∀ t : ℝ, IsLinearMap ℝ (fun p : ℝ × ℝ ↦ p.1 + (p.2 - p.1) * t) := fun t ↦
    ⟨fun p q ↦ by simp [Prod.fst_add, Prod.snd_add]; ring,
      fun c p ↦ by simp [Prod.smul_def]; ring⟩
  rw [exponentPairs_eq_inter]
  refine Convex.inter ?_ (convex_iInter fun α ↦ ?_)
  · exact ((((convex_halfSpace_ge hfst 0).inter
      (convex_halfSpace_le hfst (1 / 2))).inter
      (convex_halfSpace_ge hsnd (1 / 2))).inter
      (convex_halfSpace_le hsnd 1)).inter
      (convex_halfSpace_le hadd 1)
  · exact convex_halfSpace_ge (hlin ((α : ℝ≥0) : ℝ)) (exponentSumGrowthExponent α)

/-! ## Trivial exponent pairs -/

/-- `(0, 1)` is an exponent pair (blueprint Proposition `exp-pair-trivial`). -/
theorem isExponentPair_zero_one : IsExponentPair 0 1 := by
  refine (isExponentPair_iff_forall_exponentSumGrowthExponent_le (by norm_num)).2 ?_
  intro α _
  simpa using exponentSumGrowthExponent_le_self α

/-- The pair `(1/2, 1/2)` is an exponent pair as soon as `β(α) ≤ 1/2` holds on the upper half
`[1/2, 1]` of the scale range; on the lower half the trivial bound `β(α) ≤ α` suffices. -/
theorem isExponentPair_half_half_of_exponentSumGrowthExponent_le_half
    (hβ : ∀ α : ℝ≥0, 1 / 2 ≤ (α : ℝ) → α ≤ 1 → exponentSumGrowthExponent α ≤ 1 / 2) :
    IsExponentPair (1 / 2) (1 / 2) := by
  refine (isExponentPair_iff_forall_exponentSumGrowthExponent_le (by norm_num)).2 ?_
  intro α hα
  simp only [sub_self, zero_mul, add_zero]
  rcases le_or_gt ((α : ℝ)) (1 / 2) with h | h
  · exact (exponentSumGrowthExponent_le_self α).trans h
  · exact hβ α h.le hα

/-- `(1/2, 1/2)` is an exponent pair (blueprint Proposition `exp-pair-trivial`). -/
theorem isExponentPair_half_half : IsExponentPair (1 / 2) (1 / 2) := by
  apply isExponentPair_half_half_of_exponentSumGrowthExponent_le_half
  intro α hhalf hα
  by_cases hαone : α = 1
  · subst α
    simp
  have hαlt : α < 1 := lt_of_le_of_ne hα hαone
  let γ : ℝ≥0 := 1 - α
  have hγpos : 0 < γ := tsub_pos_of_lt hαlt
  have hαpos : 0 < α := by
    exact_mod_cast lt_of_lt_of_le (by norm_num : (0 : ℝ) < 1 / 2) hhalf
  have hγlt : γ < 1 := tsub_lt_self (by norm_num) hαpos
  have href := exponentSumGrowthExponent_reflection γ hγpos hγlt
  have hγα : 1 - γ = α := tsub_tsub_cancel_of_le hα
  rw [hγα] at href
  rw [href]
  linarith [exponentSumGrowthExponent_le_self γ]

/-! ## The exponent-pair conjecture -/

/-- The exponent-pair conjecture (blueprint Conjecture `exp-pair-conj`). -/
def ExponentPairConjecture : Prop :=
  IsExponentPair 0 (1 / 2)

/-- The exponent-pair conjecture is equivalent to every point of the exponent-pair triangle
being an exponent pair (blueprint Conjecture `exp-pair-conj`). -/
theorem exponentPairConjecture_iff_triangle :
    ExponentPairConjecture ↔
      ∀ k l : ℝ, (0 ≤ k ∧ k ≤ 1 / 2 ∧ 1 / 2 ≤ l ∧ l ≤ 1 ∧ k + l ≤ 1) →
        IsExponentPair k l := by
  constructor
  · intro h k l htriangle
    obtain ⟨hk, hkhalf, hlhalf, hl, hkl⟩ := htriangle
    -- The vertex weights are `2(1-k-ℓ)`, `2ℓ-1`, and `2k`.
    let w : Fin 3 → ℝ := ![2 * (1 - k - l), 2 * l - 1, 2 * k]
    let v : Fin 3 → ℝ × ℝ := ![(0, 1 / 2), (0, 1), (1 / 2, 1 / 2)]
    have hw : ∀ i ∈ Finset.univ, 0 ≤ w i := by
      intro i _
      fin_cases i <;> dsimp [w] <;> linarith
    have hsum : ∑ i, w i = 1 := by
      rw [Fin.sum_univ_three]
      change 2 * (1 - k - l) + (2 * l - 1) + 2 * k = 1
      ring
    have hv : ∀ i ∈ Finset.univ, v i ∈ exponentPairs := by
      intro i _
      fin_cases i
      · exact h
      · exact isExponentPair_zero_one
      · exact isExponentPair_half_half
    have hmem := convex_exponentPairs.sum_mem hw hsum hv
    have heq : (∑ i, w i • v i) = (k, l) := by
      simp only [Fin.sum_univ_three]
      ext <;> dsimp [w, v] <;> ring
    change (k, l) ∈ exponentPairs
    simpa only [heq] using hmem
  · intro h
    exact h 0 (1 / 2) (by norm_num)

/-- The exponent-pair conjecture is equivalent to `β(α) = α / 2` on `[0, 1]`
(blueprint Lemma `exp-pair-conj-beta`). -/
theorem exponentPairConjecture_iff_beta :
    ExponentPairConjecture ↔
      ∀ α : ℝ≥0, α ≤ 1 → exponentSumGrowthExponent α = (α : ℝ) / 2 := by
  rw [ExponentPairConjecture,
    isExponentPair_iff_forall_exponentSumGrowthExponent_le
      (k := 0) (l := 1 / 2) (by norm_num)]
  constructor
  · intro h α hα
    refine le_antisymm ?_ (exponentSumGrowthExponent_mem_Icc hα).1
    have := h α hα
    linarith
  · intro h α hα
    have := h α hα
    linarith

/-! ## The van der Corput A-transform -/

/-- The van der Corput `A`-transform `A(k, ℓ) = (k / (2k + 2), ℓ / (2k + 2) + 1/2)`. -/
def vanDerCorputA (p : ℝ × ℝ) : ℝ × ℝ :=
  (p.1 / (2 * p.1 + 2), p.2 / (2 * p.1 + 2) + 1 / 2)

/-- The van der Corput `A`-transform maps the exponent pair triangle into itself. -/
theorem vanDerCorputA_mem_triangle {k l : ℝ}
    (htriangle : 0 ≤ k ∧ k ≤ 1 / 2 ∧ 1 / 2 ≤ l ∧ l ≤ 1 ∧ k + l ≤ 1) :
    0 ≤ (vanDerCorputA (k, l)).1 ∧ (vanDerCorputA (k, l)).1 ≤ 1 / 2 ∧
      1 / 2 ≤ (vanDerCorputA (k, l)).2 ∧ (vanDerCorputA (k, l)).2 ≤ 1 ∧
      (vanDerCorputA (k, l)).1 + (vanDerCorputA (k, l)).2 ≤ 1 := by
  obtain ⟨hk, hk2, hl2, hl1, hkl⟩ := htriangle
  have hpos : (0 : ℝ) < 2 * k + 2 := by linarith
  have hl0 : (0 : ℝ) ≤ l := by linarith
  have hhalf : l / (2 * k + 2) ≤ 1 / 2 := by
    rw [div_le_iff₀ hpos]; linarith
  have hsum : k / (2 * k + 2) + l / (2 * k + 2) = (k + l) / (2 * k + 2) := by ring
  have hsum' : (k + l) / (2 * k + 2) ≤ 1 / 2 := by
    rw [div_le_iff₀ hpos]; linarith
  dsimp [vanDerCorputA]
  have hl0' : (0 : ℝ) ≤ l / (2 * k + 2) := by positivity
  refine ⟨by positivity, ?_, by linarith, by linarith, by linarith⟩
  rw [div_le_iff₀ hpos]; linarith

/-- The van der Corput `A`-process preserves exponent pairs (blueprint Proposition `vdc-a`). -/
theorem IsExponentPair.vanDerCorputA {k l : ℝ} (hp : IsExponentPair k l) :
    IsExponentPair (vanDerCorputA (k, l)).1 (vanDerCorputA (k, l)).2 := by
  apply (isExponentPair_iff_forall_exponentSumGrowthExponent_le
    (vanDerCorputA_mem_triangle hp.1)).mpr
  intro α hα
  let L : ℝ := (k + (l + 1) * (α : ℝ)) / (2 * k + 2)
  have hd : 0 < 2 * k + 2 := by linarith [hp.k_nonneg]
  have hL : (2 * k + 2) * L = k + (l + 1) * (α : ℝ) := by
    dsimp [L]
    exact mul_div_cancel₀ _ hd.ne'
  have htarget : (Expdb.vanDerCorputA (k, l)).1 +
      ((Expdb.vanDerCorputA (k, l)).2 - (Expdb.vanDerCorputA (k, l)).1) * (α : ℝ) = L := by
    dsimp [Expdb.vanDerCorputA, L]
    ring_nf
    field_simp [show 1 + k ≠ 0 by linarith [hp.k_nonneg]]
    ring
  rw [htarget]
  by_cases htrivial : (α : ℝ) ≤ L
  · exact (exponentSumGrowthExponent_le_self α).trans htrivial
  by_cases hhalf : 1 / 2 ≤ L
  · have hb := isExponentPair_half_half.exponentSumGrowthExponent_le hα
    simp only [sub_self, zero_mul, add_zero] at hb
    exact hb.trans hhalf
  have hLα : L < (α : ℝ) := lt_of_not_ge htrivial
  have hLhalf : L < 1 / 2 := lt_of_not_ge hhalf
  have hLnonneg : 0 ≤ L := by
    dsimp [L]
    exact div_nonneg (by nlinarith [hp.k_nonneg, hp.half_le_l, α.coe_nonneg]) hd.le
  have hαpos : 0 < (α : ℝ) := hLnonneg.trans_lt hLα
  have hsmall : (l + 1) * (α : ℝ) < 1 := by
    nlinarith [mul_pos hd (sub_pos.mpr hLhalf)]
  have hαtwo : (α : ℝ) ≤ 2 / 3 := by
    nlinarith [hp.half_le_l]
  have hLlower : 3 * (α : ℝ) / 4 ≤ L := by
    dsimp [L]
    rw [le_div_iff₀ hd]
    nlinarith [mul_nonneg hp.k_nonneg (show 0 ≤ 1 - 3 * (α : ℝ) / 2 by linarith),
      mul_nonneg (sub_nonneg.mpr hp.half_le_l) α.coe_nonneg]
  let h : ℝ := 2 * (α : ℝ) - 2 * L
  have hh : 0 ≤ h := by dsimp [h]; linarith
  have hinterval : 2 * (α : ℝ) - 1 ≤ h := by dsimp [h]; linarith
  have hdiag : 2 * (α : ℝ) - h = 2 * L := by dsimp [h]; ring
  have hshort : 2 * h ≤ 2 * L := by dsimp [h]; linarith
  let g : ℝ → ℝ := fun t ↦ (t + 1 - α) *
    exponentSumGrowthExponentReal ((α : ℝ) / (t + 1 - α)) + t
  have hsup : sSup (g '' Set.Icc (2 * (α : ℝ) - 1) h) ≤
      k + (l - 2 * k) * (α : ℝ) + (k + 1) * h := by
    apply csSup_le ((show (Set.Icc (2 * (α : ℝ) - 1) h).Nonempty from
      ⟨h, hinterval, le_rfl⟩).image g)
    rintro _ ⟨t, ht, rfl⟩
    have hD : 0 < t + 1 - (α : ℝ) := by linarith [ht.1]
    have hx : 0 ≤ (α : ℝ) / (t + 1 - α) := div_nonneg α.coe_nonneg hD.le
    have hxone : ((α : ℝ) / (t + 1 - α)).toNNReal ≤ 1 := by
      rw [← NNReal.coe_le_coe, Real.coe_toNNReal _ hx, NNReal.coe_one]
      exact (div_le_one hD).mpr (by linarith [ht.1])
    have hb := hp.exponentSumGrowthExponent_le hxone
    rw [Real.coe_toNNReal _ hx] at hb
    have hmul := mul_le_mul_of_nonneg_left hb hD.le
    rw [mul_add, mul_comm (l - k), ← mul_assoc, mul_div_cancel₀ _ hD.ne'] at hmul
    dsimp [g, exponentSumGrowthExponentReal]
    nlinarith [mul_nonneg (show 0 ≤ k + 1 by linarith [hp.k_nonneg])
      (sub_nonneg.mpr ht.2)]
  have hthird : (α : ℝ) - h +
      (k + (l - 2 * k) * (α : ℝ) + (k + 1) * h) = 2 * L := by
    dsimp [h]
    nlinarith [hL]
  have hb := vanDerCorputA_exponentSumGrowthExponent α h
    (by exact_mod_cast hαtwo) hh
  change 2 * exponentSumGrowthExponent α ≤
    max (max (2 * (α : ℝ) - h) (2 * h))
      ((α : ℝ) - h + sSup (g '' Set.Icc (2 * (α : ℝ) - 1) h)) at hb
  have hmax : max (max (2 * (α : ℝ) - h) (2 * h))
      ((α : ℝ) - h + sSup (g '' Set.Icc (2 * (α : ℝ) - 1) h)) ≤ 2 * L := by
    apply max_le (max_le hdiag.le hshort)
    linarith [hsup, hthird]
  linarith

/-! ## The van der Corput B-transform -/

/-- The van der Corput `B`-transform `B(k, ℓ) = (ℓ - 1/2, k + 1/2)`. -/
def vanDerCorputB (p : ℝ × ℝ) : ℝ × ℝ :=
  (p.2 - 1 / 2, p.1 + 1 / 2)

/-- The van der Corput `B`-transform maps the exponent pair triangle into itself. -/
theorem vanDerCorputB_mem_triangle {k l : ℝ}
    (htriangle : 0 ≤ k ∧ k ≤ 1 / 2 ∧ 1 / 2 ≤ l ∧ l ≤ 1 ∧ k + l ≤ 1) :
    0 ≤ (vanDerCorputB (k, l)).1 ∧ (vanDerCorputB (k, l)).1 ≤ 1 / 2 ∧
      1 / 2 ≤ (vanDerCorputB (k, l)).2 ∧ (vanDerCorputB (k, l)).2 ≤ 1 ∧
      (vanDerCorputB (k, l)).1 + (vanDerCorputB (k, l)).2 ≤ 1 := by
  obtain ⟨hk, hk2, hl2, hl1, hkl⟩ := htriangle
  dsimp [vanDerCorputB]
  refine ⟨by linarith, by linarith, by linarith, by linarith, by linarith⟩

/-- The van der Corput `B`-process preserves exponent pairs (blueprint Proposition `vdc-b`). -/
theorem IsExponentPair.vanDerCorputB {k l : ℝ} (hp : IsExponentPair k l) :
    IsExponentPair (vanDerCorputB (k, l)).1 (vanDerCorputB (k, l)).2 := by
  apply (isExponentPair_iff_forall_exponentSumGrowthExponent_le
    (vanDerCorputB_mem_triangle hp.1)).mpr
  intro α hα
  dsimp [Expdb.vanDerCorputB]
  by_cases hzero : α = 0
  · subst α
    simp only [exponentSumGrowthExponent_zero, NNReal.coe_zero, mul_zero, add_zero]
    linarith [hp.half_le_l]
  by_cases hone : α = 1
  · subst α
    simp only [exponentSumGrowthExponent_one, NNReal.coe_one, mul_one]
    linarith [hp.k_nonneg]
  have hαpos : 0 < α := pos_iff_ne_zero.mpr hzero
  have hαlt : α < 1 := lt_of_le_of_ne hα hone
  have hb := hp.exponentSumGrowthExponent_le (α := 1 - α) (tsub_le_self)
  rw [NNReal.coe_sub hα, NNReal.coe_one] at hb
  have href := exponentSumGrowthExponent_reflection α hαpos hαlt
  nlinarith

end Expdb
