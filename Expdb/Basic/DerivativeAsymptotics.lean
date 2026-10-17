module

public import Expdb.Basic.Asymptotics
import Mathlib.Analysis.Asymptotics.Lemmas
public import Mathlib.Analysis.Calculus.ContDiff.Basic
public import Mathlib.Analysis.Calculus.IteratedDeriv.Lemmas
import Mathlib.Analysis.Calculus.MeanValue
import Mathlib.Analysis.SpecificLimits.Basic
public import Mathlib.Order.Filter.AtTopBot.Basic
public import Mathlib.Topology.Algebra.Order.Field

/-!
# Derivative bounds and convergence along sequences

`DerivativesClose` records equality up to an infinitesimal and boundedness of the comparison
derivatives through a specified order at varying points. It uses the existing `=o` relation
and `VariableObject.IsBounded`; `PositiveDerivativesClose` omits the function values.

The product, composition, square-root and inverse rules propagate these bounds and convergence.
They are used by the dual phase construction and the stationary-phase amplitude estimates.
-/

@[expose] public section

namespace Expdb

open Filter Topology
open scoped Expdb ContDiff

noncomputable section

/-! ## Auxiliary bounds -/

private theorem const_eq_and_bounded (c : ℝ) :
    VariableObject.fixed c =o VariableObject.fixed c ∧
      VariableObject.IsBounded (VariableObject.fixed c) :=
  ⟨IsEqUpToInfinitesimal.refl _, VariableObject.IsBounded.const c⟩

private theorem mul_eq_and_bounded {u₁ v₁ u₂ v₂ : VariableObject ℝ}
    (h₁ : u₁ =o v₁ ∧ v₁.IsBounded) (h₂ : u₂ =o v₂ ∧ v₂.IsBounded) :
    (fun i ↦ u₁ i * u₂ i) =o (fun i ↦ v₁ i * v₂ i) ∧
      VariableObject.IsBounded (fun i ↦ v₁ i * v₂ i) := by
  exact ⟨by simpa only [IsEqUpToInfinitesimal, VariableObject.IsInfinitesimal, Pi.sub_apply,
      Pi.mul_apply] using h₁.1.mul h₂.1 h₁.2 h₂.2,
    by simpa only [VariableObject.IsBounded, Pi.mul_apply] using h₁.2.mul h₂.2⟩

private theorem sum_eq_and_bounded {ι : Type*} (s : Finset ι)
    {U V : ι → VariableObject ℝ} (h : ∀ k ∈ s, U k =o V k ∧ (V k).IsBounded) :
    (fun i ↦ ∑ k ∈ s, U k i) =o (fun i ↦ ∑ k ∈ s, V k i) ∧
      VariableObject.IsBounded (fun i ↦ ∑ k ∈ s, V k i) :=
  ⟨IsEqUpToInfinitesimal.sum s fun k hk ↦ (h k hk).1,
    VariableObject.IsBounded.sum s fun k hk ↦ (h k hk).2⟩

private theorem congr_eq_and_bounded {u v u' v' : VariableObject ℝ}
    (h : u =o v ∧ v.IsBounded) (hu : u =ᶠ[atTop] u') (hv : v =ᶠ[atTop] v') :
    u' =o v' ∧ v'.IsBounded := by
  refine ⟨h.1.congr hu hv, ?_⟩
  obtain ⟨C, hC⟩ := h.2
  refine ⟨C, ?_⟩
  filter_upwards [hC, hv] with i hi hvi
  rwa [← hvi]

/-! ## Derivative bounds and convergence -/

/-- Derivative differences tend to zero through order `n` at two sequences of points,
and the comparison derivatives are eventually bounded. -/
def DerivativesClose (f g : ℕ → ℝ → ℝ) (x y : ℕ → ℝ) (n : ℕ) : Prop :=
  ∀ k ≤ n,
    (fun i ↦ iteratedDeriv k (f i) (x i)) =o (fun i ↦ iteratedDeriv k (g i) (y i)) ∧
      VariableObject.IsBounded (fun i ↦ iteratedDeriv k (g i) (y i))

/-- Derivative differences of orders `1` through `n` tend to zero at two sequences of points,
and the comparison derivatives are eventually bounded. -/
def PositiveDerivativesClose (f g : ℕ → ℝ → ℝ) (x y : ℕ → ℝ) (n : ℕ) : Prop :=
  ∀ k, 1 ≤ k → k ≤ n →
    (fun i ↦ iteratedDeriv k (f i) (x i)) =o (fun i ↦ iteratedDeriv k (g i) (y i)) ∧
      VariableObject.IsBounded (fun i ↦ iteratedDeriv k (g i) (y i))

/-- Uniform derivative convergence with bounded limits gives `DerivativesClose`
at asymptotically equal points. -/
theorem DerivativesClose.of_uniform {f : ℕ → ℝ → ℝ} {f₀ : ℝ → ℝ} (hf₀ : ContDiff ℝ ∞ f₀)
    (hconv : ∀ k, ∀ ε > 0, ∀ᶠ i in atTop, ∀ x,
      ‖iteratedDeriv k (f i) x - iteratedDeriv k f₀ x‖ ≤ ε)
    (hbdd : ∀ k, ∃ K, ∀ x, ‖iteratedDeriv k f₀ x‖ ≤ K)
    {p p' : ℕ → ℝ} (hpp : Tendsto (fun i ↦ p i - p' i) atTop (𝓝 0)) (n : ℕ) :
    DerivativesClose f (fun _ ↦ f₀) p p' n := by
  intro k _
  obtain ⟨K, hK⟩ := hbdd k
  obtain ⟨K', hK'⟩ := hbdd (k + 1)
  refine ⟨?_, ⟨K, ?_⟩⟩
  · rw [isEqUpToInfinitesimal_iff_tendsto_sub]
    have h1 : Tendsto (fun i ↦ iteratedDeriv k (f i) (p i) - iteratedDeriv k f₀ (p i))
        atTop (𝓝 0) := by
      rw [Metric.tendsto_nhds]
      intro ε hε
      filter_upwards [hconv k (ε / 2) (half_pos hε)] with i hi
      rw [dist_zero_right]
      exact lt_of_le_of_lt (hi _) (half_lt_self hε)
    have hlip : ∀ x y, ‖iteratedDeriv k f₀ x - iteratedDeriv k f₀ y‖ ≤ K' * ‖x - y‖ := by
      intro x y
      exact Convex.norm_image_sub_le_of_norm_deriv_le (s := Set.univ)
        (fun z _ ↦ (hf₀.differentiable_iteratedDeriv k (by exact_mod_cast WithTop.coe_lt_top _))
          z)
        (fun z _ ↦ by rw [← iteratedDeriv_succ]; exact hK' z) convex_univ trivial trivial
    have h2 : Tendsto (fun i ↦ iteratedDeriv k f₀ (p i) - iteratedDeriv k f₀ (p' i))
        atTop (𝓝 0) := by
      apply squeeze_zero_norm (fun i ↦ hlip _ _)
      simpa using hpp.norm.const_mul K'
    simpa using h1.add h2
  · exact Eventually.of_forall fun i ↦ hK _

/-- Uniform derivative bounds give `DerivativesClose` for a family compared with itself. -/
theorem DerivativesClose.of_uniform_bound {f : ℕ → ℝ → ℝ} {z : ℕ → ℝ} {n : ℕ}
    (h : ∀ k ≤ n, ∃ K, ∀ᶠ i in atTop, ∀ x, ‖iteratedDeriv k (f i) x‖ ≤ K) :
    DerivativesClose f f z z n := by
  intro k hk
  obtain ⟨K, hK⟩ := h k hk
  exact ⟨IsEqUpToInfinitesimal.refl _, K, hK.mono fun i hi ↦ hi (z i)⟩

theorem DerivativesClose.mono {f g : ℕ → ℝ → ℝ} {x y : ℕ → ℝ} {m n : ℕ}
    (h : DerivativesClose f g x y n) (hmn : m ≤ n) : DerivativesClose f g x y m :=
  fun k hk ↦ h k (hk.trans hmn)

theorem DerivativesClose.positiveDerivativesClose {f g : ℕ → ℝ → ℝ} {x y : ℕ → ℝ} {n : ℕ}
    (h : DerivativesClose f g x y n) : PositiveDerivativesClose f g x y n :=
  fun k _ hk ↦ h k hk

theorem PositiveDerivativesClose.mono {f g : ℕ → ℝ → ℝ} {x y : ℕ → ℝ} {m n : ℕ}
    (h : PositiveDerivativesClose f g x y n) (hmn : m ≤ n) : PositiveDerivativesClose f g x y m :=
  fun k hk1 hk ↦ h k hk1 (hk.trans hmn)

/-- Differentiation shifts the orders in `DerivativesClose`. -/
theorem DerivativesClose.shift {f g : ℕ → ℝ → ℝ} {x y : ℕ → ℝ} {n : ℕ}
    (h : DerivativesClose f g x y (n + 1)) :
    DerivativesClose (fun i ↦ deriv (f i)) (fun i ↦ deriv (g i)) x y n := by
  intro k hk
  simpa only [← iteratedDeriv_succ'] using h (k + 1) (by omega)

/-- Combine convergence and boundedness of the values with those of the derivatives. -/
theorem DerivativesClose.unshift {f g : ℕ → ℝ → ℝ} {x y : ℕ → ℝ} {n : ℕ}
    (h0 : (fun i ↦ f i (x i)) =o (fun i ↦ g i (y i)))
    (hb0 : VariableObject.IsBounded (fun i ↦ g i (y i)))
    (h : DerivativesClose (fun i ↦ deriv (f i)) (fun i ↦ deriv (g i)) x y n) :
    DerivativesClose f g x y (n + 1) := by
  intro k hk
  rcases k with _ | k
  · exact ⟨h0, hb0⟩
  · simpa only [iteratedDeriv_succ'] using h k (by omega)

/-- Products preserve derivative bounds and convergence. -/
theorem DerivativesClose.mul {f g f' g' : ℕ → ℝ → ℝ} {x y : ℕ → ℝ} {n : ℕ}
    (hf : ∀ i, ContDiff ℝ ∞ (f i)) (hg : ∀ i, ContDiff ℝ ∞ (g i))
    (hf' : ∀ i, ContDiff ℝ ∞ (f' i)) (hg' : ∀ i, ContDiff ℝ ∞ (g' i))
    (h₁ : DerivativesClose f f' x y n) (h₂ : DerivativesClose g g' x y n) :
    DerivativesClose (fun i z ↦ f i z * g i z) (fun i z ↦ f' i z * g' i z) x y n := by
  intro k hk
  simp only [iteratedDeriv_fun_mul ((hf _).contDiffAt.of_le (by exact_mod_cast le_top))
    ((hg _).contDiffAt.of_le (by exact_mod_cast le_top)), iteratedDeriv_fun_mul
    ((hf' _).contDiffAt.of_le (by exact_mod_cast le_top))
    ((hg' _).contDiffAt.of_le (by exact_mod_cast le_top))]
  apply sum_eq_and_bounded
  intro m hm
  have hm' := Finset.mem_range.1 hm
  exact mul_eq_and_bounded (mul_eq_and_bounded (const_eq_and_bounded _)
    (h₁ m (by omega))) (h₂ (k - m) (by omega))

/-- The derivative of a composition, differentiated `j` more times. -/
private theorem iteratedDeriv_succ_comp {g f : ℝ → ℝ} (hg : ContDiff ℝ ∞ g) (hf : ContDiff ℝ ∞ f)
    (j : ℕ) (x : ℝ) :
    iteratedDeriv (j + 1) (g ∘ f) x =
      ∑ m ∈ Finset.range (j + 1), (j.choose m : ℝ) *
        iteratedDeriv m (deriv g ∘ f) x * iteratedDeriv (j - m + 1) f x := by
  have hd : deriv (g ∘ f) = fun y ↦ (deriv g ∘ f) y * deriv f y := by
    funext y
    exact deriv_comp y ((hg.differentiable (by simp)) _) ((hf.differentiable (by simp)) _)
  rw [iteratedDeriv_succ', hd,
    iteratedDeriv_fun_mul
      (((hg.deriv' (n := ∞)).comp hf).contDiffAt.of_le (by exact_mod_cast le_top))
      ((hf.deriv' (n := ∞)).contDiffAt.of_le (by exact_mod_cast le_top))]
  refine Finset.sum_congr rfl fun m _ ↦ ?_
  rw [← iteratedDeriv_succ']

/-- The chain rule for derivative bounds and convergence. -/
theorem DerivativesClose.comp {n : ℕ} :
    ∀ {g g' f f' : ℕ → ℝ → ℝ} {x x' : ℕ → ℝ},
    (∀ i, ContDiff ℝ ∞ (g i)) → (∀ i, ContDiff ℝ ∞ (g' i)) →
    (∀ i, ContDiff ℝ ∞ (f i)) → (∀ i, ContDiff ℝ ∞ (f' i)) →
    DerivativesClose g g' (fun i ↦ f i (x i)) (fun i ↦ f' i (x' i)) n →
    PositiveDerivativesClose f f' x x' n →
    DerivativesClose (fun i ↦ g i ∘ f i) (fun i ↦ g' i ∘ f' i) x x' n := by
  induction n with
  | zero =>
    intro g g' f f' x x' _ _ _ _ hg _ k hk
    obtain rfl : k = 0 := by omega
    simpa using hg 0 le_rfl
  | succ n ih =>
    intro g g' f f' x x' hgs hg's hfs hf's hg hf k hk
    rcases k with _ | j
    · simpa using hg 0 (by omega)
    · have hIH := ih (g := fun i ↦ deriv (g i)) (g' := fun i ↦ deriv (g' i))
        (fun i ↦ (hgs i).deriv' (n := ∞)) (fun i ↦ (hg's i).deriv' (n := ∞))
        hfs hf's hg.shift (hf.mono (by omega))
      simp only [iteratedDeriv_succ_comp (hgs _) (hfs _),
        iteratedDeriv_succ_comp (hg's _) (hf's _)]
      apply sum_eq_and_bounded
      intro m hm
      have hm' := Finset.mem_range.1 hm
      exact mul_eq_and_bounded (mul_eq_and_bounded (const_eq_and_bounded _)
        (hIH m (by omega))) (hf (j - m + 1) (by omega) (by omega))

/-- If `u * a = b`, derivative bounds and convergence for `a` and `b` imply them for `u`,
provided `a` stays away from zero. -/
theorem DerivativesClose.of_mul_eq {u u' a a' b b' : ℕ → ℝ → ℝ} {x x' : ℕ → ℝ} {n : ℕ}
    (hu : ∀ i, ContDiff ℝ ∞ (u i)) (hu' : ∀ i, ContDiff ℝ ∞ (u' i))
    (ha : ∀ i, ContDiff ℝ ∞ (a i)) (ha' : ∀ i, ContDiff ℝ ∞ (a' i))
    (heq : ∀ i y, u i y * a i y = b i y) (heq' : ∀ i y, u' i y * a' i y = b' i y)
    (haj : DerivativesClose a a' x x' n) (hbj : DerivativesClose b b' x x' n)
    {c : ℝ} (hc : 0 < c) (hlow : ∀ᶠ i in atTop, c ≤ ‖a' i (x' i)‖) :
    DerivativesClose u u' x x' n := by
  have hclose : ((fun i ↦ a i (x i)) =o (fun i ↦ a' i (x' i)) ∧
      VariableObject.IsBounded (fun i ↦ a' i (x' i))) := by
    simpa using haj 0 (Nat.zero_le _)
  have hainv : ((fun i ↦ (a i (x i))⁻¹) =o (fun i ↦ (a' i (x' i))⁻¹) ∧
      VariableObject.IsBounded (fun i ↦ (a' i (x' i))⁻¹)) :=
    ⟨hclose.1.inv hc hlow, VariableObject.isBounded_inv_of_norm_ge hc hlow⟩
  have ha0 : ∀ᶠ i in atTop, a i (x i) ≠ 0 ∧ a' i (x' i) ≠ 0 := by
    have hsmall : ∀ᶠ i in atTop, ‖a i (x i) - a' i (x' i)‖ ≤ c / 2 :=
      ((isEqUpToInfinitesimal_iff_tendsto_sub _ _).1 hclose.1).norm.eventually
        (Iic_mem_nhds (by simpa using half_pos hc))
    filter_upwards [hsmall, hlow] with i h1 h2
    constructor
    · intro h0
      rw [h0, zero_sub, norm_neg] at h1
      linarith
    · intro h0
      rw [h0, norm_zero] at h2
      linarith
  -- Solve the Leibniz identity for the top derivative of `u`.
  have hformula : ∀ (v A B : ℝ → ℝ), ContDiff ℝ ∞ v → ContDiff ℝ ∞ A →
      (∀ y, v y * A y = B y) → ∀ k z, A z ≠ 0 →
      iteratedDeriv k v z = (iteratedDeriv k B z -
        ∑ m ∈ Finset.range k, (k.choose m : ℝ) * iteratedDeriv m v z *
          iteratedDeriv (k - m) A z) * (A z)⁻¹ := by
    intro v A B hv hA hvAB k z hz
    have hB : B = fun y ↦ v y * A y := funext fun y ↦ (hvAB y).symm
    rw [hB, iteratedDeriv_fun_mul ((hv).contDiffAt.of_le (by exact_mod_cast le_top))
      ((hA).contDiffAt.of_le (by exact_mod_cast le_top)), Finset.sum_range_succ]
    simp only [Nat.choose_self, Nat.cast_one, one_mul, Nat.sub_self, iteratedDeriv_zero]
    field_simp
    ring
  have key : DerivativesClose u u' x x' n := by
    intro k
    induction k using Nat.strong_induction_on with
    | _ k ih =>
      intro hk
      have hsum := sum_eq_and_bounded (Finset.range k) (fun m hm ↦ by
        have hm' := Finset.mem_range.1 hm
        exact mul_eq_and_bounded (mul_eq_and_bounded (const_eq_and_bounded (k.choose m : ℝ))
          (ih m hm' (by omega))) (haj (k - m) (by omega)))
      have hexpr := mul_eq_and_bounded
        ⟨(hbj k hk).1.sub hsum.1, (hbj k hk).2.sub hsum.2⟩ hainv
      refine congr_eq_and_bounded hexpr ?_ ?_
      · filter_upwards [ha0] with i hi
        exact (hformula (u i) (a i) (b i) (hu i) (ha i) (heq i) k (x i) hi.1).symm
      · filter_upwards [ha0] with i hi
        exact (hformula (u' i) (a' i) (b' i) (hu' i) (ha' i) (heq' i) k (x' i) hi.2).symm
  exact key

/-- Nonnegative square roots preserve derivative bounds and convergence when they stay
away from zero. -/
theorem DerivativesClose.of_mul_self_eq {s s' h h' : ℕ → ℝ → ℝ} {x x' : ℕ → ℝ} {n : ℕ}
    (hs : ∀ i, ContDiff ℝ ∞ (s i)) (hs' : ∀ i, ContDiff ℝ ∞ (s' i))
    (hpos : ∀ i y, 0 ≤ s i y) (hpos' : ∀ i y, 0 ≤ s' i y)
    (heq : ∀ i y, s i y * s i y = h i y) (heq' : ∀ i y, s' i y * s' i y = h' i y)
    (hhj : DerivativesClose h h' x x' n)
    {c : ℝ} (hc : 0 < c) (hlow : ∀ᶠ i in atTop, c ≤ s' i (x' i)) :
    DerivativesClose s s' x x' n := by
  have hh0 : ((fun i ↦ h i (x i)) =o (fun i ↦ h' i (x' i)) ∧
      VariableObject.IsBounded (fun i ↦ h' i (x' i))) := by
    simpa using hhj 0 (Nat.zero_le _)
  have h0 : ((fun i ↦ s i (x i)) =o (fun i ↦ s' i (x' i)) ∧
      VariableObject.IsBounded (fun i ↦ s' i (x' i))) := by
    have hbd : VariableObject.IsBounded (fun i ↦ s' i (x' i)) := by
      obtain ⟨C, hC⟩ := hh0.2
      refine ⟨C + 1, ?_⟩
      filter_upwards [hC] with i hi
      rw [← heq'] at hi
      rw [Real.norm_of_nonneg (hpos' _ _)]
      rw [Real.norm_of_nonneg (mul_self_nonneg _)] at hi
      nlinarith [hpos' i (x' i)]
    refine ⟨?_, hbd⟩
    rw [isEqUpToInfinitesimal_iff_tendsto_sub]
    have hq : IsBoundedUnder (· ≤ ·) atTop
        (fun i ↦ ‖(s i (x i) + s' i (x' i))⁻¹‖) := by
      refine ⟨c⁻¹, ?_⟩
      simp only [eventually_map]
      filter_upwards [hlow] with i hi
      rw [norm_inv, Real.norm_of_nonneg (add_nonneg (hpos _ _) (hpos' _ _))]
      exact inv_anti₀ hc (by linarith [hpos i (x i)])
    have := ((isEqUpToInfinitesimal_iff_tendsto_sub _ _).1 hh0.1).zero_mul_isBoundedUnder_le hq
    refine this.congr' ?_
    filter_upwards [hlow] with i hi
    have hne : s i (x i) + s' i (x' i) ≠ 0 := by linarith [hpos i (x i)]
    rw [← heq, ← heq']
    field_simp
    ring
  have hs0 : ((fun i ↦ (2 * s i (x i))⁻¹) =o (fun i ↦ (2 * s' i (x' i))⁻¹) ∧
      VariableObject.IsBounded (fun i ↦ (2 * s' i (x' i))⁻¹)) := by
    have heq : (fun i ↦ 2 * s i (x i)) =o (fun i ↦ 2 * s' i (x' i)) := by
      simpa only [IsEqUpToInfinitesimal, VariableObject.IsInfinitesimal, Pi.sub_apply,
        Pi.smul_apply, smul_eq_mul] using h0.1.const_smul (2 : ℝ)
    have hlow2 : ∀ᶠ i in atTop, 2 * c ≤ ‖2 * s' i (x' i)‖ := by
      filter_upwards [hlow] with i hi
      rw [Real.norm_of_nonneg (by linarith [hpos' i (x' i)])]
      linarith
    exact ⟨heq.inv (by positivity) hlow2,
      VariableObject.isBounded_inv_of_norm_ge (by positivity) hlow2⟩
  have hlowl : ∀ᶠ i in atTop, s i (x i) ≠ 0 := by
    have hsmall : ∀ᶠ i in atTop, ‖s i (x i) - s' i (x' i)‖ ≤ c / 2 :=
      ((isEqUpToInfinitesimal_iff_tendsto_sub _ _).1 h0.1).norm.eventually
        (Iic_mem_nhds (by simpa using half_pos hc))
    filter_upwards [hsmall, hlow] with i h1 h2 h3
    rw [h3, zero_sub, norm_neg, Real.norm_of_nonneg (hpos' _ _)] at h1
    linarith
  have hlowr : ∀ᶠ i in atTop, s' i (x' i) ≠ 0 := by
    filter_upwards [hlow] with i hi h0
    rw [h0] at hi
    linarith
  have hformula : ∀ (v H : ℝ → ℝ), ContDiff ℝ ∞ v → (∀ y, v y * v y = H y) →
      ∀ k z, v z ≠ 0 →
      iteratedDeriv (k + 1) v z = (iteratedDeriv (k + 1) H z -
        ∑ m ∈ Finset.range k, ((k + 1).choose (m + 1) : ℝ) * iteratedDeriv (m + 1) v z *
          iteratedDeriv (k + 1 - (m + 1)) v z) * (2 * v z)⁻¹ := by
    intro v H hv hvH k z hz
    have hH : H = fun y ↦ v y * v y := funext fun y ↦ (hvH y).symm
    rw [hH, iteratedDeriv_fun_mul ((hv).contDiffAt.of_le (by exact_mod_cast le_top))
      ((hv).contDiffAt.of_le (by exact_mod_cast le_top)), Finset.sum_range_succ',
      Finset.sum_range_succ]
    simp only [Nat.choose_self, Nat.cast_one, one_mul, Nat.sub_self, iteratedDeriv_zero,
      Nat.choose_zero_right, Nat.sub_zero]
    field_simp
    ring
  have key : DerivativesClose s s' x x' n := by
    intro k
    induction k using Nat.strong_induction_on with
    | _ k ih =>
      intro hk
      rcases k with _ | k
      · simpa using h0
      · have hsum := sum_eq_and_bounded (Finset.range k) (fun m hm ↦ by
          have hm' := Finset.mem_range.1 hm
          exact mul_eq_and_bounded (mul_eq_and_bounded
            (const_eq_and_bounded ((k + 1).choose (m + 1) : ℝ))
            (ih (m + 1) (by omega) (by omega)))
            (ih (k + 1 - (m + 1)) (by omega) (by omega)))
        have hexpr := mul_eq_and_bounded
          ⟨(hhj (k + 1) hk).1.sub hsum.1, (hhj (k + 1) hk).2.sub hsum.2⟩ hs0
        refine congr_eq_and_bounded hexpr ?_ ?_
        · filter_upwards [hlowl] with i hi
          exact (hformula (s i) (h i) (hs i) (heq i) k (x i) hi).symm
        · filter_upwards [hlowr] with i hi
          exact (hformula (s' i) (h' i) (hs' i) (heq' i) k (x' i) hi).symm
  exact key

/-- If `Y ∘ X = id`, bounds and convergence of the derivatives of `X` through order `n + 1`
follow from those of its values and of the derivatives of `Y'` through order `n` at the images. -/
theorem DerivativesClose.inverse {X X' Y Y' : ℕ → ℝ → ℝ} {x x' : ℕ → ℝ} {n : ℕ}
    (hX : ∀ i, ContDiff ℝ ∞ (X i)) (hX' : ∀ i, ContDiff ℝ ∞ (X' i))
    (hY : ∀ i, ContDiff ℝ ∞ (Y i)) (hY' : ∀ i, ContDiff ℝ ∞ (Y' i))
    (hinv : ∀ i y, Y i (X i y) = y) (hinv' : ∀ i y, Y' i (X' i y) = y)
    (hYj : DerivativesClose (fun i ↦ deriv (Y i)) (fun i ↦ deriv (Y' i))
      (fun i ↦ X i (x i)) (fun i ↦ X' i (x' i)) n)
    (h0 : (fun i ↦ X i (x i)) =o (fun i ↦ X' i (x' i)))
    (hb0 : VariableObject.IsBounded (fun i ↦ X' i (x' i)))
    {c : ℝ} (hc : 0 < c) (hlow : ∀ᶠ i in atTop, c ≤ ‖deriv (Y' i) (X' i (x' i))‖) :
    DerivativesClose X X' x x' (n + 1) := by
  have hprod : ∀ (X₀ Y₀ : ℝ → ℝ), ContDiff ℝ ∞ X₀ → ContDiff ℝ ∞ Y₀ →
      (∀ y, Y₀ (X₀ y) = y) → ∀ y, deriv X₀ y * (deriv Y₀ ∘ X₀) y = 1 := by
    intro X₀ Y₀ hX₀ hY₀ hinv₀ y
    have hcomp : deriv (Y₀ ∘ X₀) y = 1 := by
      have : Y₀ ∘ X₀ = id := funext hinv₀
      rw [this, deriv_id]
    rw [deriv_comp y ((hY₀.differentiable (by simp)) _) ((hX₀.differentiable (by simp)) _)]
      at hcomp
    simpa [mul_comm] using hcomp
  have step : ∀ m ≤ n, DerivativesClose X X' x x' m → DerivativesClose X X' x x' (m + 1) := by
    intro m hm hXm
    have hcomp := DerivativesClose.comp (g := fun i ↦ deriv (Y i)) (g' := fun i ↦ deriv (Y' i))
      (fun i ↦ (hY i).deriv' (n := ∞)) (fun i ↦ (hY' i).deriv' (n := ∞))
      hX hX' (hYj.mono hm) hXm.positiveDerivativesClose
    have hu := DerivativesClose.of_mul_eq (u := fun i ↦ deriv (X i)) (u' := fun i ↦ deriv (X' i))
      (b := fun _ _ ↦ 1) (b' := fun _ _ ↦ 1)
      (fun i ↦ (hX i).deriv' (n := ∞)) (fun i ↦ (hX' i).deriv' (n := ∞))
      (fun i ↦ ((hY i).deriv' (n := ∞)).comp (hX i))
      (fun i ↦ ((hY' i).deriv' (n := ∞)).comp (hX' i))
      (fun i y ↦ hprod (X i) (Y i) (hX i) (hY i) (hinv i) y)
      (fun i y ↦ hprod (X' i) (Y' i) (hX' i) (hY' i) (hinv' i) y)
      hcomp (by
        intro k _
        rcases k with _ | k
        · exact const_eq_and_bounded 1
        · simpa [iteratedDeriv_succ'] using const_eq_and_bounded 0) hc hlow
    exact DerivativesClose.unshift h0 hb0 hu
  have : ∀ m ≤ n + 1, DerivativesClose X X' x x' m := by
    intro m
    induction m with
    | zero =>
      intro _ k hk
      obtain rfl : k = 0 := by omega
      exact ⟨h0, hb0⟩
    | succ m ih => intro hm; exact step m (by omega) (ih (by omega))
  exact this (n + 1) le_rfl

end

end Expdb
