module

public import Expdb.Basic.Asymptotics

/-!
# Automatic uniformity

This module formalizes the automatic uniformity results of Chapter 2 of the ANTEDB blueprint:
pointwise boundedness or infinitesimality along every variable sequence can be made uniform
after passing to a subsequence.
-/

@[expose] public section

open Filter Topology

namespace Expdb

variable {α : Type*} [SeminormedAddCommGroup α]

/-! ### Auxiliary lemmas for subsequence extraction -/

/-- From a property holding arbitrarily late, extract a strictly
    increasing sequence on which it holds. -/
private lemma extract_strictMono_subseq {P : ℕ → Prop}
    (h : ∀ j, ∃ i ≥ j, P i) :
    ∃ φ : ℕ → ℕ, StrictMono φ ∧ ∀ n, P (φ n) := by
  have step : ∀ prev, ∃ i > prev, P i :=
    fun prev => let ⟨i, hi, hP⟩ := h (prev + 1); ⟨i, by omega, hP⟩
  let φ : ℕ → ℕ :=
    fun n => n.rec (h 0).choose (fun _ prev => (step prev).choose)
  exact ⟨φ, strictMono_nat_of_lt_succ fun n => (step (φ n)).choose_spec.1,
         fun n => n.casesOn (h 0).choose_spec.2
                             (fun _ => (step (φ _)).choose_spec.2)⟩

/-- Extract a strictly increasing φ and bad elements x(n) ∈ E(φ n)
    with |f(φ n)(x n)| > n.  Used in the proof of Proposition 2.1(i). -/
private lemma extract_bad_seq_i
    (E : VariableObject (Set ℝ)) (f : VariableFunction (fun i ↦ E i) α)
    (bad : ∀ j, ∃ i ≥ j, ∃ x : E i, (j : ℝ) < ‖f i x‖) :
    ∃ φ : ℕ → ℕ, StrictMono φ ∧
    ∃ x : ∀ n, E (φ n), ∀ n, n < ‖f (φ n) (x n)‖ := by
  have step : ∀ (n prev : ℕ), ∃ i > prev, ∃ x : E i, (n : ℝ) < ‖f i x‖ :=
    fun n prev => by
      obtain ⟨i, hi, x, hx⟩ := bad (max (prev + 1) n)
      exact ⟨i, by have := le_max_left (prev+1) n; omega,
             x, lt_of_le_of_lt (by exact_mod_cast le_max_right (prev+1) n) hx⟩
  let data : ℕ → Σ i, E i :=
    fun n => n.rec ⟨(step 0 0).choose, (step 0 0).choose_spec.2.choose⟩
      fun n p => ⟨(step (n+1) p.1).choose, (step (n+1) p.1).choose_spec.2.choose⟩
  exact ⟨fun n => (data n).1,
         strictMono_nat_of_lt_succ fun n => (step (n+1) (data n).1).choose_spec.1,
         fun n => (data n).2,
         fun n => n.casesOn (step 0 0).choose_spec.2.choose_spec
                             (fun _ => (step _ _).choose_spec.2.choose_spec)⟩

/-- For a fixed threshold ε, extract a strictly increasing φ and bad elements
    x(n) ∈ E(φ n) with |f(φ n)(x n)| > ε. Used in Proposition 2.1(ii). -/
private lemma extract_bad_seq_ii
    (E : VariableObject (Set ℝ)) (f : VariableFunction (fun i ↦ E i) α)
    {ε : ℝ}
    (bad : ∀ j, ∃ i ≥ j, ∃ x : E i, ε < ‖f i x‖) :
    ∃ φ : ℕ → ℕ, StrictMono φ ∧
    ∃ x : ∀ n, E (φ n), ∀ n, ε < ‖f (φ n) (x n)‖ := by
  obtain ⟨φ, hφ, hbad⟩ := extract_strictMono_subseq
    (P := fun i => ∃ x : E i, ε < ‖f i x‖) bad
  choose x hx using hbad
  exact ⟨φ, hφ, x, hx⟩

/-- Build a strictly increasing threshold sequence φ such that
    |f(φ n)(x)| ≤ 1/(n+1) for all x ∈ E(φ n). -/
private lemma build_increasing_thresholds
    (E : VariableObject (Set ℝ)) (f : VariableFunction (fun i ↦ E i) α)
    (scale : ∀ n : ℕ, 0 < n → ∃ i_n, ∀ i ≥ i_n, ∀ x : E i,
             ‖f i x‖ ≤ 1/n) :
    ∃ φ : ℕ → ℕ, StrictMono φ ∧
    ∀ n, ∀ x : E (φ n), ‖f (φ n) x‖ ≤ 1/(n+1) := by
  choose i_seq hi_seq using fun n => scale (n+1) (Nat.succ_pos n)
  let φ : ℕ → ℕ :=
    fun n => n.rec (i_seq 0) (fun k p => max (i_seq (k+1)) (p+1))
  refine ⟨φ, strictMono_nat_of_lt_succ fun n =>
           Nat.lt_of_lt_of_le (Nat.lt_succ_self _) (le_max_right _ _),
         fun n x => ?_⟩
  have hge : φ n ≥ i_seq n := by
    cases n with
    | zero => exact le_refl _
    | succ n => exact le_max_left _ _
  have hb := hi_seq n (φ n) hge x
  rwa [Nat.cast_add, Nat.cast_one] at hb

/-! ### Proposition 2.1: automatic uniformity -/

open Classical in
private noncomputable def extend_subsequence
    (E : VariableObject (Set ℝ)) (hE : ∀ i, (E i).Nonempty)
    (φ : ℕ → ℕ) (x : ∀ n, E (φ n)) : ∀ j, E j :=
  fun j =>
    if h : ∃ n, φ n = j then
      (show E (φ h.choose) = E j by rw [h.choose_spec]) ▸ x h.choose
    else ⟨(hE j).choose, (hE j).choose_spec⟩

private lemma norm_extend_subsequence_apply
    {E : VariableObject (Set ℝ)} (hE : ∀ i, (E i).Nonempty)
    {f : VariableFunction (fun i ↦ E i) α} {φ : ℕ → ℕ} (hφ : StrictMono φ)
    (x : ∀ n, E (φ n)) (m : ℕ) :
    ‖f (φ m) (extend_subsequence E hE φ x (φ m))‖ = ‖f (φ m) (x m)‖ := by
  classical
  simp only [extend_subsequence]
  split_ifs with h
  · have hm : h.choose = m := hφ.injective h.choose_spec
    have hx : x h.choose ≍ x m := by rw [hm]
    apply congrArg
    apply congrArg
    apply eq_of_heq
    exact rec_heq_of_heq _ hx
  · exact absurd ⟨m, rfl⟩ h

/-- **Proposition 2.1(i) — Automatic uniform bound.**
    If f(x) = O(1) for every variable x ∈ E, then after passing to a
    subsequence there exists a *fixed* C with |f(x)| ≤ C for all x ∈ E. -/
theorem automatic_uniformity_of_pointwise_bounded
    (E : VariableObject (Set ℝ)) (hE : ∀ i, (E i).Nonempty)
    (f : VariableFunction (fun i ↦ E i) α)
    (hf : f.IsPointwiseBounded) :
    ∃ φ : ℕ → ℕ, StrictMono φ ∧
    ∃ C : ℝ, ∀ i, ∀ x : E (φ i),
    ‖f (φ i) x‖ ≤ C := by
  have h_eventual :
      ∃ φ : ℕ → ℕ, StrictMono φ ∧
      ∃ C : ℝ, ∀ᶠ i in atTop, ∀ x : E (φ i), ‖f (φ i) x‖ ≤ C := by
    by_contra h_fail
    push Not at h_fail
    -- Build a bad sequence: for each j find i ≥ j and x with |f i x| > j
    have bad : ∀ j, ∃ i ≥ j, ∃ x : E i, (j:ℝ) < ‖f i x‖ := fun j => by
      rcases Filter.frequently_atTop.mp (h_fail id strictMono_id j) j with ⟨i, hi, x, hx⟩
      exact ⟨i, hi, x, hx⟩
    obtain ⟨φ, hφ, x_bad, hx_bad⟩ := extract_bad_seq_i E f bad
    let y : ∀ j, E j := extend_subsequence E hE φ x_bad
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
    have heq := norm_extend_subsequence_apply hE (f := f) hφ x_bad m
    linarith [hj₀ (φ m) hm_ge,
              hx_bad m,
              show C_y < (m:ℝ) from
              lt_trans hn₁ (by exact_mod_cast hm_large),
              heq ▸ hj₀ (φ m) hm_ge]
  obtain ⟨φ, hφ, C, hC⟩ := h_eventual
  rw [Filter.eventually_atTop] at hC
  obtain ⟨N, hN⟩ := hC
  refine ⟨fun n => φ (n + N), ?_, C, ?_⟩
  · intro m n hmn
    exact hφ (Nat.add_lt_add_right hmn N)
  · intro n x
    exact hN (n + N) (by omega) x

/-- **Proposition 2.1(ii) — Automatic uniform infinitesimal.**
    If f(x) = o(1) for every variable x ∈ E, then after passing to a
    subsequence there exists an *infinitesimal* c with |f(x)| ≤ c for all x ∈ E. -/
theorem automatic_uniformity_of_pointwise_infinitesimal
    (E : VariableObject (Set ℝ)) (hE : ∀ i, (E i).Nonempty)
    (f : VariableFunction (fun i ↦ E i) α)
    (hf : f.IsPointwiseInfinitesimal) :
    ∃ φ : ℕ → ℕ, StrictMono φ ∧
    ∃ c : VariableObject ℝ, c.IsInfinitesimal ∧
    ∀ i, ∀ x : E (φ i),
    ‖f (φ i) x‖ ≤ c i := by
  -- Step 1: for each n ≥ 1, the bound 1/n eventually holds uniformly
  have scale : ∀ n : ℕ, 0 < n → ∃ i_n, ∀ i ≥ i_n, ∀ x : E i,
      ‖f i x‖ ≤ 1/n := by
    intro n hn
    by_contra h_fail; push Not at h_fail
    have bad : ∀ j, ∃ i ≥ j, ∃ x : E i, (1:ℝ)/n < ‖f i x‖ :=
      fun j => by obtain ⟨i, hi, x, hx⟩ := h_fail j; exact ⟨i, hi, x, hx⟩
    obtain ⟨φ, hφ, x_bad, hx_bad⟩ :=
      extract_bad_seq_ii E f bad
    let y : ∀ j, E j := extend_subsequence E hE φ x_bad
    have hfy := hf y
    rw [VariableObject.IsInfinitesimal] at hfy
    rw [Metric.tendsto_atTop] at hfy
    obtain ⟨N₀, hN₀⟩ := hfy (1/(2*n)) (by positivity)
    obtain ⟨m, hm⟩ := (hφ.tendsto_atTop).eventually (eventually_ge_atTop N₀) |>.exists
    have heq := norm_extend_subsequence_apply hE (f := f) hφ x_bad m
    have h1 : ‖f (φ m) (y (φ m))‖ < 1/(2*n) := by
      have := hN₀ (φ m) hm
      rwa [Real.dist_eq, sub_zero, abs_of_nonneg (norm_nonneg _)] at this
    linarith [hx_bad m, heq ▸ h1,
              show (1:ℝ)/(2*n) < 1/n by
                have hn : (0 : ℝ) < n := by positivity
                exact one_div_lt_one_div_of_lt hn (by nlinarith)]
  -- Step 2: build strictly increasing thresholds and conclude
  obtain ⟨φ, hφ, hφ_bd⟩ := build_increasing_thresholds E f scale
  refine ⟨φ, hφ, fun n => 1/(↑n+1), ?_, hφ_bd⟩
  rw [VariableObject.IsInfinitesimal]
  have hc : Tendsto (fun n : ℕ ↦ 1 / ((n : ℝ) + 1)) atTop (nhds 0) :=
    tendsto_one_div_add_atTop_nhds_zero_nat
  simpa using hc.norm

/-! ### Full-tail uniformity -/

/-- Boundedness along every variable choice is equivalent to an eventual bound uniform over
the original variable sets. This strengthens `automatic_uniformity_of_pointwise_bounded`. -/
theorem VariableFunction.isPointwiseBounded_iff_eventually_uniform
    (E : VariableObject (Set ℝ)) (hE : ∀ i, (E i).Nonempty)
    (f : VariableFunction (fun i ↦ E i) α) :
    f.IsPointwiseBounded ↔
      ∃ C : ℝ, ∀ᶠ i in atTop, ∀ x : E i, ‖f i x‖ ≤ C := by
  constructor
  · intro hf
    by_contra h
    push Not at h
    have bad : ∀ j, ∃ i ≥ j, ∃ x : E i, (j : ℝ) < ‖f i x‖ := fun j => by
      rcases Filter.frequently_atTop.mp (h j) j with ⟨i, hi, x, hx⟩
      exact ⟨i, hi, x, hx⟩
    obtain ⟨φ, hφ, x, hx⟩ := extract_bad_seq_i E f bad
    let y : ∀ j, E j := extend_subsequence E hE φ x
    obtain ⟨C, hC⟩ := hf y
    rw [Filter.eventually_atTop] at hC
    obtain ⟨j₀, hj₀⟩ := hC
    obtain ⟨n₁, hn₁⟩ := exists_nat_gt C
    obtain ⟨m, hm_ge, hm_large⟩ : ∃ m, φ m ≥ j₀ ∧ m > n₁ := by
      obtain ⟨m, hm⟩ := (hφ.tendsto_atTop).eventually (eventually_ge_atTop j₀) |>.exists
      exact ⟨max m (n₁ + 1), le_trans hm (hφ.monotone (le_max_left _ _)), by omega⟩
    have heq := norm_extend_subsequence_apply hE (f := f) hφ x m
    linarith [hj₀ (φ m) hm_ge, hx m,
      show C < (m : ℝ) from
        lt_trans hn₁ (by exact_mod_cast hm_large),
      heq ▸ hj₀ (φ m) hm_ge]
  · rintro ⟨C, hC⟩ x
    exact ⟨C, hC.mono fun i hi ↦ hi (x i)⟩

/-- Infinitesimality along every variable choice is equivalent to convergence that is eventually
uniform over the original variable sets. This strengthens
`automatic_uniformity_of_pointwise_infinitesimal`. -/
theorem VariableFunction.isPointwiseInfinitesimal_iff_forall_pos_uniform
    (E : VariableObject (Set ℝ)) (hE : ∀ i, (E i).Nonempty)
    (f : VariableFunction (fun i ↦ E i) α) :
    f.IsPointwiseInfinitesimal ↔
      ∀ ε : ℝ, 0 < ε →
        ∀ᶠ i in atTop, ∀ x : E i, ‖f i x‖ < ε := by
  constructor
  · intro hf ε hε
    classical
    let x₀ : ∀ i, E i := fun i ↦ ⟨(hE i).choose, (hE i).choose_spec⟩
    let x : ∀ i, E i := fun i ↦
      if h : ∃ y : E i, ε ≤ ‖f i y‖ then h.choose else x₀ i
    have hx (i : ℕ) (h : ∃ y : E i, ε ≤ ‖f i y‖) :
        ε ≤ ‖f i (x i)‖ := by
      simp only [x, dif_pos h]
      exact h.choose_spec
    have hsmall :=
      (VariableObject.isInfinitesimal_iff_forall_pos
        (fun i ↦ f i (x i))).1 (hf x) ε hε
    by_contra h
    have hfrequent :
        ∃ᶠ i in atTop, ε ≤ ‖f i (x i)‖ := (not_eventually.mp h).mono fun i hi ↦ by
      apply hx i
      simpa only [not_forall, not_lt] using hi
    obtain ⟨i, hilarge, hismall⟩ := (hfrequent.and_eventually hsmall).exists
    exact (not_lt_of_ge hilarge) hismall
  · intro h x
    apply (VariableObject.isInfinitesimal_iff_forall_pos
      (fun i ↦ f i (x i))).2
    intro ε hε
    filter_upwards [h ε hε] with i hi
    exact hi (x i)

end Expdb
