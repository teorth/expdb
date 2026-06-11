/-
  ANTEDB Blueprint -- Chapter 2: Basic Notation
  ===============================================
-/

import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.Topology.Algebra.Order.LiminfLimsup
import Mathlib.Analysis.Normed.Field.Basic
import Mathlib.Analysis.Asymptotics.Asymptotics
import Mathlib.Data.EReal.Basic
import Mathlib.Order.Filter.Basic
import Mathlib.Topology.MetricSpace.Sequences

open Filter Topology Asymptotics Real

-- ===========================================================
--  Function e(θ) = exp(2πiθ)
-- ===========================================================

/-- Definition: e(θ) := e^(2πiθ) as defined in Blueprint p. 4 -/
noncomputable def e (θ : ℝ) : ℂ :=
  Complex.exp (2 * Real.pi * θ * Complex.I)

/-- Base case: e(0) = 1
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
-- Empty Supremum / Infimum Conventions 
-- ===========================================================

/-- Convention: empty supremum = -∞ (⊥ in EReal) -/
lemma blueprintSup_empty_convention :
    Sup (∅ : Set EReal) = ⊥ :=
  EReal.sSup_empty

/-- Convention: empty infimum = +∞ (⊤ in EReal) -/
lemma blueprintInf_empty_convention :
    Inf (∅ : Set EReal) = ⊤ :=
  EReal.sInf_empty

/-- sup_{σ₀ ≤ σ < σ₁} f(σ) = -∞ when σ₁ < σ₀ -/
noncomputable def blueprintSup (σ₀ σ₁ : ℝ) (f : ℝ → EReal) : EReal :=
  Sup (f '' {σ : ℝ | σ₀ ≤ σ ∧ σ < σ₁})

lemma blueprintSup_empty {σ₀ σ₁ : ℝ} (f : ℝ → EReal) (h : σ₁ < σ₀) :
    blueprintSup σ₀ σ₁ f = ⊥ := by
  have h_empty : {σ : ℝ | σ₀ ≤ σ ∧ σ < σ₁} = ∅ := by
    ext σ
    simp only [Set.mem_setOf_eq, Set.mem_empty_iff_false, iff_false]
    intro ⟨h1, h2⟩
    linarith
  simp [blueprintSup, h_empty]
  exact EReal.sSup_empty

-- ===========================================================
-- N^(-∞) = 0 Convention 
-- ===========================================================

/-- Extended real power: handles N^r for r : EReal.
    The only special case stated in the blueprint is N^(-∞) = 0 when N > 1.
    For r = +∞ we leave it as the natural limit (not specified by blueprint). -/
noncomputable def blueprintPower (N : ℝ) (r : EReal) : ℝ :=
  if r = ⊥ then 0
  else N ^ r.toReal

/-- Blueprint convention: N^(-∞) = 0 for N > 1 -/
lemma blueprintPower_bot {N : ℝ} (hN : N > 1) :
    blueprintPower N ⊥ = 0 := by
  simp [blueprintPower]

/-- For finite exponents, blueprintPower agrees with real power -/
lemma blueprintPower_coe {N : ℝ} (r : ℝ) :
    blueprintPower N (r : EReal) = N ^ r := by
  simp [blueprintPower, EReal.coe_ne_bot, EReal.toReal_coe]

-- ===========================================================
-- Indicator Function 
-- ===========================================================

/-- 1_I(n) = 1 if n ∈ I, else 0 -/
def indicatorFunction {α : Type*} [DecidableEq α] (I : Set α) 
    [DecidablePred (· ∈ I)] (n : α) : ℝ :=
  if n ∈ I then 1 else 0

/-- Indicator is 0 or 1 -/
lemma indicatorFunction_values {α : Type*} [DecidableEq α] 
    (I : Set α) [DecidablePred (· ∈ I)] (n : α) :
    indicatorFunction I n = 0 ∨ indicatorFunction I n = 1 := by
  simp [indicatorFunction]
  split_ifs <;> simp

-- ===========================================================
-- Cardinality |W| for Finsets 
-- ===========================================================

/-- We use Finset.card for cardinality, written |W| in the blueprint.
    In Lean we use W.card or Finset.card W to avoid notation conflicts. -/
example (W : Finset ℝ) : W.card = Finset.card W := rfl

-- ===========================================================
--  Separated Sequences
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
-- 1-Bounded Sequences
-- ===========================================================

/-- 1-Bounded complex sequence: |aₙ| ≤ 1 for all n -/
def IsOneBounded (a : ℕ → ℂ) : Prop :=
  ∀ n, Complex.abs (a n) ≤ 1

/-- The phase sequence e(θₙ) is always 1-bounded -/
lemma e_is_one_bounded (θ : ℕ → ℝ) : IsOneBounded (fun n => e (θ n)) := by
  intro n
  rw [norm_e]

-- ===========================================================
-- Asymptotic Notation 
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

-- ============================================================
--  Auxiliary lemmas for subsequence extraction
-- ============================================================

/-- From a property holding arbitrarily late, extract a strictly
    increasing sequence on which it holds. -/
private lemma extract_strictMono_subseq {P : ℕ → Prop}
    (h : ∀ j, ∃ i ≥ j, P i) :
    ∃ φ : ℕ → ℕ, StrictMono φ ∧ ∀ n, P (φ n) := by
  have step : ∀ prev, ∃ i > prev, P i :=
    fun prev => let ⟨i, hi, hP⟩ := h (prev + 1); ⟨i, by omega, hP⟩
  let φ : ℕ → ℕ :=
    fun n => n.rec (h 0).choose (fun _ prev => (step prev).choose)
  exact ⟨strictMono_nat_of_lt_succ fun n => (step (φ n)).choose_spec.1,
         fun n => n.casesOn (h 0).choose_spec.2
                             (fun _ => (step (φ _)).choose_spec.2)⟩

/-- Extract a strictly increasing φ and bad elements x(n) ∈ E(φ n)
    with |f(φ n)(x n)| > n.  Used in the proof of Proposition 2.1(i). -/
private lemma extract_bad_seq_i
    (E : ℕ → Finset ℝ) (f : ∀ i, (E i : Set ℝ) → ℂ)
    (bad : ∀ j, ∃ i ≥ j, ∃ x : (E i : Set ℝ), (j : ℝ) < Complex.abs (f i x)) :
    ∃ φ : ℕ → ℕ, StrictMono φ ∧
    ∃ x : ∀ n, (E (φ n) : Set ℝ), ∀ n, (n : ℝ) < Complex.abs (f (φ n) (x n)) := by
  have step : ∀ n prev, ∃ i > prev, ∃ x : (E i : Set ℝ), (n : ℝ) < Complex.abs (f i x) :=
    fun n prev => by
      obtain ⟨i, hi, x, hx⟩ := bad (max (prev + 1) n)
      exact ⟨i, by linarith [le_max_left (prev+1) n, hi],
             x, by exact_mod_cast lt_of_le_of_lt (le_max_right _ _) hx⟩
  let data : ℕ → Σ i, (E i : Set ℝ) :=
    fun n => n.rec ⟨(step 0 0).choose, (step 0 0).choose_spec.2.choose⟩
      fun n p => ⟨(step (n+1) p.1).choose, (step (n+1) p.1).choose_spec.2.choose⟩
  exact ⟨strictMono_nat_of_lt_succ fun n => (step (n+1) (data n).1).choose_spec.1,
         fun n => (data n).2,
         fun n => n.casesOn (step 0 0).choose_spec.2.choose_spec
                             (fun _ => (step _ _).choose_spec.2.choose_spec)⟩

/-- Extract a strictly increasing φ and bad elements x(n) ∈ E(φ n)
    with |f(φ n)(x n)| ≥ 1/(n+1).  Used in the proof of Proposition 2.1(ii). -/
private lemma extract_bad_seq_ii
    (E : ℕ → Finset ℝ) (f : ∀ i, (E i : Set ℝ) → ℂ)
    (bad : ∀ j, ∃ i ≥ j, ∃ x : (E i : Set ℝ), (1:ℝ)/(j+1) ≤ Complex.abs (f i x)) :
    ∃ φ : ℕ → ℕ, StrictMono φ ∧
    ∃ x : ∀ n, (E (φ n) : Set ℝ), ∀ n, (1:ℝ)/(n+1) ≤ Complex.abs (f (φ n) (x n)) := by
  have step : ∀ n prev, ∃ i > prev, ∃ x : (E i : Set ℝ), (1:ℝ)/(n+1) ≤ Complex.abs (f i x) :=
    fun n prev => by
      obtain ⟨i, hi, x, hx⟩ := bad (max (prev+1) n)
      refine ⟨i, by linarith [le_max_left (prev+1) n, hi], x, ?_⟩
      calc (1:ℝ)/(n+1) ≤ 1/(max (prev+1) n + 1) := by
              gcongr; exact_mod_cast Nat.succ_le_succ (le_max_right _ _)
           _ ≤ Complex.abs (f i x) := hx
  let data : ℕ → Σ i, (E i : Set ℝ) :=
    fun n => n.rec ⟨(step 0 0).choose, (step 0 0).choose_spec.2.choose⟩
      fun n p => ⟨(step (n+1) p.1).choose, (step (n+1) p.1).choose_spec.2.choose⟩
  exact ⟨strictMono_nat_of_lt_succ fun n => (step (n+1) (data n).1).choose_spec.1,
         fun n => (data n).2,
         fun n => n.casesOn (step 0 0).choose_spec.2.choose_spec
                             (fun _ => (step _ _).choose_spec.2.choose_spec)⟩

/-- Build a strictly increasing threshold sequence φ such that
    |f(φ n)(x)| ≤ 1/(n+1) for all x ∈ E(φ n). -/
private lemma build_increasing_thresholds
    (E : ℕ → Finset ℝ) (f : ∀ i, (E i : Set ℝ) → ℂ)
    (scale : ∀ n, 0 < n → ∃ i_n, ∀ i ≥ i_n, ∀ x : (E i : Set ℝ),
             Complex.abs (f i x) ≤ 1/n) :
    ∃ φ : ℕ → ℕ, StrictMono φ ∧
    ∀ n, ∀ x : (E (φ n) : Set ℝ), Complex.abs (f (φ n) x) ≤ 1/(n+1) := by
  choose i_seq hi_seq using fun n => scale (n+1) (Nat.succ_pos n)
  let φ : ℕ → ℕ :=
    fun n => n.rec (i_seq 0) (fun k p => max (i_seq (k+1)) (p+1))
  refine ⟨strictMono_nat_of_lt_succ fun n =>
           Nat.lt_of_lt_of_le (Nat.lt_succ_self _) (le_max_right _ _),
         fun n x => hi_seq n (φ n) ?_ x⟩
  cases n with
  | zero => exact le_refl _
  | succ n => exact le_max_left _ _

-- ===========================================================
-- Underspill Principle 
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

-- ============================================================
-- Asymptotic relations  X = O(Y),  X ≪ Y,  X ≍ Y
-- ============================================================

/-- X = O(Y): there exists a fixed C with |X| ≤ C·Y eventually. -/
def IsBigOSeq (X Y : ℕ → ℝ) : Prop :=
  ∃ C : ℝ, 0 < C ∧ ∀ᶠ i in atTop, |X i| ≤ C * Y i

/-- X ≪ Y  (X is much less than Y):
    same as X = O(Y) in the blueprint's variable-quantity sense. -/
def IsVeryLT (X Y : ℕ → ℝ) : Prop := IsBigOSeq X Y

notation X " ≪ " Y => IsVeryLT X Y
notation Y " ≫ " X => IsVeryLT X Y

/-- X = o(Y): there exists an infinitesimal c with |X| ≤ c·Y eventually. -/
def IsLittleOSeq (X Y : ℕ → ℝ) : Prop :=
  ∃ c : ℕ → ℝ, IsInfinitesimal c ∧ ∀ᶠ i in atTop, |X i| ≤ c i * Y i

/-- X ≍ Y  (X and Y are comparable):
    X ≪ Y and Y ≪ X,  i.e. X = O(Y) and Y = O(X). -/
def IsAsymptoticallyEquiv (X Y : ℕ → ℝ) : Prop :=
  (X ≪ Y) ∧ (Y ≪ X)

notation X " ≍ " Y => IsAsymptoticallyEquiv X Y

-- ============================================================
-- Pointwise-bounded and pointwise-infinitesimal functions
-- ============================================================

/-- f is pointwise O(1): for every variable sequence (x_i) ∈ E_i,
    the values (f_i(x_i)) are eventually bounded. -/
def IsPointwiseBounded (E : ℕ → Set ℝ) (f : ∀ i, E i → ℂ) : Prop :=
  ∀ x : ∀ i, E i, ∃ C : ℝ, ∀ᶠ i in atTop, Complex.abs (f i (x i)) ≤ C

/-- f is pointwise o(1): for every variable sequence (x_i) ∈ E_i,
    the values (f_i(x_i)) tend to 0. -/
def IsPointwiseInfinitesimal (E : ℕ → Set ℝ) (f : ∀ i, E i → ℂ) : Prop :=
  ∀ x : ∀ i, E i, IsInfinitesimal (fun i => Complex.abs (f i (x i)))

-- ============================================================
-- Proposition 2.1 — Automatic uniformity   
-- ============================================================

-- Helper: rewrite |f(φ m)(y(φ m))| as |f(φ m)(x_bad m)| avoiding cast issues.
private lemma abs_y_eq_abs_x_bad
    {E : ℕ → Finset ℝ} {f : ∀ i, (E i : Set ℝ) → ℂ}
    {φ : ℕ → ℕ} (hφ : StrictMono φ)
    {x_bad : ∀ n, (E (φ n) : Set ℝ)}
    {default_elem : ∀ i, (E i : Set ℝ)}
    (m : ℕ) :
    let y : ∀ j, (E j : Set ℝ) := fun j =>
      if h : ∃ n, φ n = j then
        (show (E (φ h.choose) : Set ℝ) = E j by rw [h.choose_spec]) ▸ x_bad h.choose
      else default_elem j
    Complex.abs (f (φ m) (y (φ m))) = Complex.abs (f (φ m) (x_bad m)) := by
  simp only []
  split_ifs with h
  · generalize_proofs hp
    have : h.choose = m := hφ.injective hp
    rw [this]
  · exact absurd ⟨m, rfl⟩ h

/-- **Proposition 2.1(i) — Automatic uniform bound.**
    If f(x) = O(1) for every variable x ∈ E, then after passing to a
    subsequence there exists a *fixed* C with |f(x)| ≤ C for all x ∈ E. -/
theorem automatic_uniformity_i
    (E : ℕ → Finset ℝ) (hE : ∀ i, (E i).Nonempty)
    (f : ∀ i, (E i : Set ℝ) → ℂ)
    (hf : IsPointwiseBounded (fun i => (E i : Set ℝ)) f) :
    ∃ φ : ℕ → ℕ, StrictMono φ ∧
    ∃ C : ℝ, ∀ᶠ i in atTop, ∀ x : (E (φ i) : Set ℝ),
    Complex.abs (f (φ i) x) ≤ C := by
  by_contra h_fail
  push_neg at h_fail
  -- Build a bad sequence: for each j find i ≥ j and x with |f i x| > j
  have bad : ∀ j, ∃ i ≥ j, ∃ x : (E i : Set ℝ), (j:ℝ) < Complex.abs (f i x) := fun j => by
    obtain ⟨i, hi, x, hx⟩ := (h_fail id strictMono_id j).exists; exact ⟨i, hi, x, hx⟩
  obtain ⟨φ, hφ, x_bad, hx_bad⟩ := extract_bad_seq_i E f bad
  -- Extend x_bad to a full variable sequence y
  let default_elem : ∀ i, (E i : Set ℝ) :=
    fun i => ⟨(hE i).choose, (hE i).choose_mem⟩
  classical
  let y : ∀ j, (E j : Set ℝ) := fun j =>
    if h : ∃ n, φ n = j then
      (show (E (φ h.choose) : Set ℝ) = E j by rw [h.choose_spec]) ▸ x_bad h.choose
    else default_elem j
  -- Apply pointwise bound to y
  obtain ⟨C_y, hC_y⟩ := hf y
  rw [Filter.eventually_atTop] at hC_y
  obtain ⟨j₀, hj₀⟩ := hC_y
  obtain ⟨n₁, hn₁⟩ := exists_nat_gt C_y
  -- Find m with φ(m) ≥ j₀ and m > n₁
  obtain ⟨m, hm_ge, hm_large⟩ : ∃ m, φ m ≥ j₀ ∧ m > n₁ := by
    obtain ⟨m, hm⟩ := (hφ.tendsto_atTop).eventually (eventually_ge_atTop j₀) |>.exists
    exact ⟨max m (n₁+1), le_trans hm (hφ.monotone (le_max_left _ _)), by omega⟩
  -- Derive contradiction
  have heq := abs_y_eq_abs_x_bad hφ m (default_elem := default_elem)
  linarith [hj₀ (φ m) hm_ge (y (φ m)),
            hx_bad m,
            show C_y < (m:ℝ) from by exact_mod_cast hm_large,
            heq ▸ hj₀ (φ m) hm_ge (y (φ m))]

/-- **Proposition 2.1(ii) — Automatic uniform infinitesimal.**
    If f(x) = o(1) for every variable x ∈ E, then after passing to a
    subsequence there exists an *infinitesimal* c with |f(x)| ≤ c for all x ∈ E. -/
theorem automatic_uniformity_ii
    (E : ℕ → Finset ℝ) (hE : ∀ i, (E i).Nonempty)
    (f : ∀ i, (E i : Set ℝ) → ℂ)
    (hf : IsPointwiseInfinitesimal (fun i => (E i : Set ℝ)) f) :
    ∃ φ : ℕ → ℕ, StrictMono φ ∧
    ∃ c : ℕ → ℝ, IsInfinitesimal c ∧
    ∀ᶠ i in atTop, ∀ x : (E (φ i) : Set ℝ),
    Complex.abs (f (φ i) x) ≤ c i := by
  -- Step 1: for each n ≥ 1, the bound 1/n eventually holds uniformly
  have scale : ∀ n, 0 < n → ∃ i_n, ∀ i ≥ i_n, ∀ x : (E i : Set ℝ),
      Complex.abs (f i x) ≤ 1/n := by
    intro n hn
    by_contra h_fail; push_neg at h_fail
    have bad : ∀ j, ∃ i ≥ j, ∃ x : (E i : Set ℝ), (1:ℝ)/n < Complex.abs (f i x) :=
      fun j => by obtain ⟨i, hi, x, hx⟩ := h_fail j; exact ⟨i, hi, x, hx⟩
    obtain ⟨φ, hφ, x_bad, hx_bad⟩ :=
      extract_bad_seq_ii E f (fun j => let ⟨i,hi,x,hx⟩ := bad j; ⟨i,hi,x,le_of_lt hx⟩)
    let default_elem : ∀ i, (E i : Set ℝ) :=
      fun i => ⟨(hE i).choose, (hE i).choose_mem⟩
    classical
    let y : ∀ j, (E j : Set ℝ) := fun j =>
      if h : ∃ m, φ m = j then
        (show (E (φ h.choose) : Set ℝ) = E j by rw [h.choose_spec]) ▸ x_bad h.choose
      else default_elem j
    have hfy := hf y
    rw [IsInfinitesimal, Metric.tendsto_atTop] at hfy
    obtain ⟨N₀, hN₀⟩ := hfy (1/(2*n)) (by positivity)
    obtain ⟨m, hm⟩ := (hφ.tendsto_atTop).eventually (eventually_ge_atTop N₀) |>.exists
    have heq := abs_y_eq_abs_x_bad hφ m (default_elem := default_elem)
    have h1 : Complex.abs (f (φ m) (y (φ m))) < 1/(2*n) := by
      have := hN₀ (φ m) hm
      rwa [Real.dist_eq, abs_of_nonneg (Complex.abs.nonneg _), sub_zero] at this
    linarith [hx_bad m, heq ▸ h1,
              show (1:ℝ)/(2*n) < 1/n by
                apply div_lt_div_of_pos_left one_pos (by positivity); linarith]
  -- Step 2: build strictly increasing thresholds and conclude
  obtain ⟨φ, hφ, hφ_bd⟩ := build_increasing_thresholds E f scale
  refine ⟨φ, hφ, fun n => 1/(↑n+1), ?_, Filter.eventually_atTop.mpr ⟨0, fun n _ x => hφ_bd n x⟩⟩
  rw [IsInfinitesimal]
  exact (tendsto_const_nhds.div_atTop
    (Filter.Tendsto.atTop_add tendsto_natCast_atTop_atTop tendsto_const_nhds)).congr
    (fun n => by ring_nf)
