/-
  ANTEDB Blueprint -- Chapter 2: Basic Notation
  ===============================================
-/

import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.Topology.Algebra.Order.LiminfLimsup
import Mathlib.Analysis.Normed.Field.Basic
import Mathlib.Analysis.Asymptotics.Asymptotics

open Filter Topology Asymptotics

-- ===========================================================
-- SECTION 1: Phase Function e(θ) = exp(2πiθ)
-- ===========================================================

/-- Definition: e(θ) := e^(2πiθ) as defined in Blueprint p. 4 -/
noncomputable def e (θ : ℝ) : ℂ :=
  Complex.exp (2 * Real.pi * θ * Complex.I)

/-- Base case: e(0) = 1 -/
lemma e_zero : e 0 = 1 := by
  simp [e]

/-- Absolute value / Norm: |e(θ)| = 1 for all θ -/
lemma norm_e (θ : ℝ) : Complex.abs (e θ) = 1 := by
  have h : (2 * Real.pi * θ * Complex.I).re = 0 := by simp
  rw [Complex.abs_exp, h, Real.exp_zero]

/-- Homomorphism property: e(θ₁ + θ₂) = e(θ₁) * e(θ₂) -/
lemma e_add (θ₁ θ₂ : ℝ) : e (θ₁ + θ₂) = e θ₁ * e θ₂ := by
  simp [e, ← Complex.exp_add]

/-- Periodicity over integers: e(n) = 1 for any n : ℤ -/
lemma e_int (n : ℤ) : e n = 1 := by
  simp [e]
  have h : (2 * Real.pi * (n : ℝ) * Complex.I) = (n : ℂ) * (2 * Real.pi * Complex.I) := by push_cast; ring
  rw [h, Complex.exp_int_mul_two_pi_mul_I]

-- ===========================================================
-- SECTION 2 & 3: Asymptotic Notation
-- ===========================================================

/-- Bounded Sequence: defined as O(1) atTop -/
def IsBoundedSeq (X : ℕ → ℝ) : Prop :=
  IsBigO atTop X (fun _ => (1 : ℝ))

/-- Infinitesimal Sequence: defined as o(1) atTop -/
def IsInfinitesimal (X : ℕ → ℝ) : Prop :=
  IsLittleO atTop X (fun _ => (1 : ℝ))

-- Custom notation matching the Blueprint's style
notation X " =O= " Y => IsBigO atTop X Y
notation X " =o= " Y => IsLittleO atTop X Y

-- ===========================================================
-- SECTION 4: Separated Sequences
-- ===========================================================

/-- 1-Separated Sets: distance between distinct elements is at least 1 -/
def IsSeparated (W : Finset ℝ) : Prop :=
  ∀ t ∈ W, ∀ t' ∈ W, t ≠ t' → 1 ≤ |t - t'|

/-- λ-Separated Sets: distance between distinct elements is at least λ -/
def IsLambdaSeparated (λ : ℝ) (W : Finset ℝ) : Prop :=
  ∀ t ∈ W, ∀ t' ∈ W, t ≠ t' → λ ≤ |t - t'|

/-- 1-separated is structurally identical to λ-separated with λ = 1 -/
lemma isSeparated_iff_isLambdaSeparated_one (W : Finset ℝ) :
    IsSeparated W ↔ IsLambdaSeparated 1 W := by
  rfl

-- ===========================================================
-- SECTION 5: 1-Bounded Sequences
-- ===========================================================

/-- 1-Bounded complex sequence: |aₙ| ≤ 1 for all n -/
def IsOneBounded (a : ℕ → ℂ) : Prop :=
  ∀ n, Complex.abs (a n) ≤ 1

/-- The phase sequence e(θₙ) is always 1-bounded -/
lemma e_is_one_bounded (θ : ℕ → ℝ) : IsOneBounded (fun n => e (θ n)) := by
  intro n
  rw [norm_e]

-- ===========================================================
-- Cheap nonstandard notation
-- ===========================================================

/-- An infinitesimal sequence: a sequence that converges to 0 -/
def IsInfinitesimal (X : ℕ → ℝ) : Prop :=
  Tendsto X atTop (nhds 0)

/-- X ≤ Y + o(1) in a strict sense:
    There exists an infinitesimal sequence ε_i such that x_i ≤ y_i + ε_i eventually. -/
def EventuallyLeUpToInfinitesimal (X Y : ℕ → ℝ) : Prop :=
  ∃ ε : ℕ → ℝ, IsInfinitesimal ε ∧ 
               (∀ᶠ i in atTop, X i ≤ Y i + ε i)

-- Notation shorthand
notation X " ≤o " Y => EventuallyLeUpToInfinitesimal X Y

-- ===========================================================
-- Underspill Principle (Blueprint p. 6)
-- ===========================================================

/-- Underspill Principle:
    X ≤ Y + o(1)  ↔  For every constant ε > 0, X ≤ Y + ε + o(1) -/
theorem underspill (X Y : ℕ → ℝ) :
    (X ≤o Y) ↔ 
    (∀ ε : ℝ, ε > 0 → X ≤o (fun i => Y i + ε)) := by
  constructor

  -- =======================
  -- Forward Direction (→)
  -- =======================
  · intro ⟨εseq, hεseq_inf, hεseq_bound⟩ ε hε
    -- We choose the exact same sequence εseq
    -- x_i ≤ y_i + εseq_i ≤ y_i + ε + εseq_i
    use εseq
    constructor
    · exact hεseq_inf
    · filter_upwards [hεseq_bound] with i hi
      -- hi : x_i ≤ y_i + εseq_i
      -- Goal: x_i ≤ (y_i + ε) + εseq_i
      linarith

  -- =======================
  -- Backward Direction (←)
  -- =======================
  · intro h
    -- We want to construct an infinitesimal sequence d_i such that x_i ≤ y_i + d_i
    -- Strategy: For each c > 0, by hypothesis with ε = c/2:
    --   x_i ≤ y_i + c/2 + d_i where d_i → 0
    --   For sufficiently large i: d_i < c/2
    --   Therefore: x_i ≤ y_i + c
    -- This implies: x_i - y_i ≤ c for all c > 0
    -- We build the sequence explicitly.
    
    -- For each n : ℕ, use ε = 1/(n+1)
    -- From the hypothesis, we obtain d^n_i such that x_i ≤ y_i + 1/(n+1) + d^n_i
    -- We define ε_i = inf_{n} (1/(n+1) + d^n_i)
    -- However, this is complex, so we use a more direct approach:
    
    -- We define z_i = max(x_i - y_i, 0)
    -- and prove that z_i → 0
    
    -- First, we prove: ∀ c > 0, ∀ᶠ i, x_i - y_i < c
    have key : ∀ c : ℝ, c > 0 → ∀ᶠ i in atTop, X i - Y i < c := by
      intro c hc
      -- Use the hypothesis with ε = c/2
      have hc2 : c / 2 > 0 := by linarith
      obtain ⟨dseq, hdseq_inf, hdseq_bound⟩ := h (c / 2) hc2
      -- dseq → 0, so ∀ᶠ i, |dseq i| < c/2
      rw [Metric.tendsto_atTop] at hdseq_inf
      have hdseq_small := hdseq_inf (c / 2) hc2
      -- For sufficiently large i: x_i ≤ y_i + c/2 + dseq_i and |dseq_i| < c/2
      filter_upwards [hdseq_bound, hdseq_small] with i hi_bound hi_small
      -- hi_bound : x_i ≤ y_i + c/2 + dseq_i
      -- hi_small : dist (dseq i) 0 < c/2, i.e., |dseq_i| < c/2
      rw [Real.dist_eq] at hi_small
      simp at hi_small
      -- dseq_i < c/2 follows from |dseq_i| < c/2
      have hdseq_lt : dseq i < c / 2 := by
        exact lt_of_abs_lt hi_small
      linarith
    
    -- Now we construct the infinitesimal sequence
    -- We use the sequence z_i = max(x_i - y_i, 0)
    -- and prove that it converges to 0
    
    -- Alternatively, more simply: we use x_i - y_i directly
    -- and prove that (x_i - y_i)⁺ → 0, then conclude
    
    -- For simplicity, we show the existence of ε_i = max(x_i - y_i, 1/i) approximately
    -- But the simplest proof uses the squeeze theorem
    
    -- We define ε_i explicitly via: for any i, take 1/(i+1) as an approximation
    -- If x_i ≤ y_i + 1/(i+1) + d_i where d_i → 0
    
    -- Direct proof: we show x_i - y_i → 0 by definition
    use fun i => max (X i - Y i) 0
    constructor
    · -- We prove that max(x_i - y_i, 0) → 0
      rw [Metric.tendsto_atTop]
      intro δ hδ
      -- From key with c = δ
      have h_ev := key δ hδ
      -- We also need x_i - y_i > -δ, but this is not guaranteed
      -- In fact, max(z, 0) ≤ |z|, so it suffices that |x_i - y_i| < δ
      -- But key only provides x_i - y_i < δ
      -- We use key with c = δ
      filter_upwards [h_ev] with i hi
      -- hi : x_i - y_i < δ
      rw [Real.dist_eq]
      simp
      constructor
      · -- max(x_i - y_i, 0) < δ
        constructor
        · linarith
        · linarith
      · -- -δ < max(x_i - y_i, 0)
        have : max (X i - Y i) 0 ≥ 0 := le_max_right _ _
        linarith
    · -- We prove x_i ≤ y_i + max(x_i - y_i, 0)
      apply Filter.Eventually.of_forall
      intro i
      have : max (X i - Y i) 0 ≥ X i - Y i := le_max_left _ _
      linarith
