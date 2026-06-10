/-
  ANTEDB Blueprint -- Chapter 2: Basic Notation
  ===============================================
  Source: "Database of known results on analytic number theory exponents"
  Tao, Trudgian, Yang (2026)
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
-- SECTION 6: Underspill Principle (Blueprint p. 6)
-- ===========================================================

/-- Underspill Principle: If the difference is o(1), 
    then X_i ≤ Y_i + ε holds eventually for all ε > 0. -/
lemma underspill (X Y : ℕ → ℝ) (h : IsInfinitesimal (fun i => X i - Y i)) :
    ∀ ε > 0, ∀ᶠ i in atTop, X i ≤ Y i + ε := by
  intro ε hε
  have h_tendsto : Tendsto (fun i => X i - Y i) atTop (nhds 0) := h.tendsto_nhds
  have h_mem : {x : ℝ | x < ε} ∈ nhds (0 : ℝ) := Iio_mem_nhds hε
  have h_ev := h_tendsto h_mem
  filter_upwards [h_ev] with i hi
  simp at hi
  linarith

-- ===========================================================
-- SECTION 7: Automatic Uniformity
-- ===========================================================

/-- Formulation of Proposition 2.1 (i): If a sequence is pointwise bounded for each x,
    it is uniformly bounded by a constant C independent of both x and i. -/
def AutomaticUniformityProperty (f : ℕ → ℕ → ℝ) : Prop :=
  (∀ x, IsBoundedSeq (fun i => f i x)) → ∃ C : ℝ, ∀ x i, |f i x| ≤ C

end
