module

public import Expdb.Basic.Definitions
public import Mathlib.Analysis.SpecificLimits.Basic
public import Mathlib.Order.Filter.Basic
public import Mathlib.Topology.MetricSpace.Sequences

/-!
# Project-specific asymptotic notation

This module formalizes the infinitesimal equality and comparison relations and the underspill
principle from Chapter 2 of the ANTEDB blueprint. The notations `X =o Y` and `X ≤o Y` are
available after `open scoped Expdb`.

The blueprint also uses non-standard objects indexed by an ambient natural-number parameter.
These are represented by `VariableObject`. Their asymptotic properties are encoded in the
`VariableObject` namespace.

Variable functions whose domain type may depend on the ambient parameter are represented by
`VariableFunction`. Their pointwise asymptotic properties are encoded in the
`VariableFunction` namespace.
-/

@[expose] public section

open Filter Topology

namespace Expdb

universe u v

/-- A family indexed by the ambient natural-number parameter. -/
abbrev VariableObject (α : Type u) := ℕ → α

/-- A family of functions whose domain may depend on the ambient natural-number parameter. -/
abbrev VariableFunction (domain : VariableObject (Type u)) (codomain : Type v) :=
  ∀ i, domain i → codomain

namespace VariableObject

/-- The constant variable object with value `x`. -/
abbrev fixed {α : Type u} (x : α) : VariableObject α := fun _ ↦ x

variable {α : Type u} [SeminormedAddCommGroup α]

/-- A variable object whose norm is eventually bounded by a fixed constant. -/
def IsBounded (X : VariableObject α) : Prop :=
  ∃ C : ℝ, ∀ᶠ i in atTop, ‖X i‖ ≤ C

/-- A variable object whose norm tends to infinity. -/
def IsUnbounded (X : VariableObject α) : Prop :=
  Tendsto (fun i ↦ ‖X i‖) atTop atTop

/-- A variable object is unbounded iff its norm eventually exceeds every fixed threshold. -/
theorem isUnbounded_iff_forall_eventually_norm_ge (X : VariableObject α) :
    X.IsUnbounded ↔ ∀ C : ℝ, ∀ᶠ i in atTop, C ≤ ‖X i‖ := by
  rw [IsUnbounded, tendsto_atTop]

/-- A variable object whose norm tends to zero. -/
def IsInfinitesimal (X : VariableObject α) : Prop :=
  Tendsto (fun i ↦ ‖X i‖) atTop (nhds 0)

/-- A variable object is infinitesimal iff its norm is eventually below every positive
tolerance. -/
theorem isInfinitesimal_iff_forall_pos (X : VariableObject α) :
    X.IsInfinitesimal ↔
    ∀ δ : ℝ, 0 < δ → ∀ᶠ i in atTop, ‖X i‖ < δ := by
  rw [IsInfinitesimal, Metric.tendsto_nhds]
  simp only [dist_zero_right, norm_norm]

end VariableObject

namespace VariableFunction

variable {domain : VariableObject (Type u)}
variable {α : Type v} [SeminormedAddCommGroup α]

/-- A variable function that is bounded along every variable choice.
The choice may vary with the ambient index; it is not a fixed point of the domain. -/
def IsPointwiseBounded (f : VariableFunction domain α) : Prop :=
  ∀ x : ∀ i, domain i, VariableObject.IsBounded (fun i ↦ f i (x i))

/-- A variable function that is infinitesimal along every variable choice.
The choice may vary with the ambient index; it is not a fixed point of the domain. -/
def IsPointwiseInfinitesimal (f : VariableFunction domain α) : Prop :=
  ∀ x : ∀ i, domain i, VariableObject.IsInfinitesimal (fun i ↦ f i (x i))

end VariableFunction

section InfinitesimalRelations

variable {E : Type u} [SeminormedAddCommGroup E]

/-- The relation `X = Y + o(1)`, meaning that `X - Y` is infinitesimal. -/
def IsEqUpToInfinitesimal (X Y : VariableObject E) : Prop :=
  (X - Y).IsInfinitesimal

/-- The relation `X ≤ Y + o(1)`, meaning that the positive part of `X - Y` is
infinitesimal. -/
def IsLEUpToInfinitesimal (X Y : VariableObject ℝ) : Prop :=
  VariableObject.IsInfinitesimal (fun i ↦ max (X i - Y i) 0)

/-- `X =o Y` denotes `X = Y + o(1)`, not the little-o relation `X = o(Y)`. -/
scoped[Expdb] infix:50 " =o " => Expdb.IsEqUpToInfinitesimal

/-- `X ≤o Y` denotes `X ≤ Y + o(1)`. -/
scoped[Expdb] infix:50 " ≤o " => Expdb.IsLEUpToInfinitesimal

open scoped Expdb

private theorem isEqUpToInfinitesimal_iff_tendsto_sub
    (X Y : VariableObject E) :
    (X =o Y) ↔ Tendsto (fun i ↦ X i - Y i) atTop (nhds 0) := by
  rw [IsEqUpToInfinitesimal, VariableObject.IsInfinitesimal]
  change Tendsto (fun i ↦ ‖X i - Y i‖) atTop (nhds 0) ↔
    Tendsto (fun i ↦ X i - Y i) atTop (nhds 0)
  exact tendsto_zero_iff_norm_tendsto_zero.symm

/-- The relation `=o` is reflexive. -/
theorem IsEqUpToInfinitesimal.refl (X : VariableObject E) :
    X =o X := by
  rw [isEqUpToInfinitesimal_iff_tendsto_sub]
  simpa only [Pi.sub_apply, sub_self] using
    (tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ (0 : E)) atTop (nhds 0))

/-- The relation `=o` is symmetric. -/
theorem IsEqUpToInfinitesimal.symm
    {X Y : VariableObject E} (h : X =o Y) :
    Y =o X := by
  rw [isEqUpToInfinitesimal_iff_tendsto_sub] at h ⊢
  simpa only [Pi.sub_apply, neg_sub, neg_zero] using h.neg

/-- The relation `=o` is transitive. -/
theorem IsEqUpToInfinitesimal.trans
    {X Y Z : VariableObject E} (hXY : X =o Y) (hYZ : Y =o Z) :
    X =o Z := by
  rw [isEqUpToInfinitesimal_iff_tendsto_sub] at hXY hYZ ⊢
  simpa only [sub_add_sub_cancel, zero_add] using hXY.add hYZ

/-- The relation `=o` is preserved by negation. -/
theorem IsEqUpToInfinitesimal.neg
    {X Y : VariableObject E} (h : X =o Y) :
    -X =o -Y := by
  rw [isEqUpToInfinitesimal_iff_tendsto_sub] at h ⊢
  simpa [sub_eq_add_neg, add_comm] using h.neg

/-- The relation `=o` is preserved by addition. -/
theorem IsEqUpToInfinitesimal.add
    {X₁ X₂ Y₁ Y₂ : VariableObject E}
    (h₁ : X₁ =o Y₁) (h₂ : X₂ =o Y₂) :
    X₁ + X₂ =o Y₁ + Y₂ := by
  rw [isEqUpToInfinitesimal_iff_tendsto_sub] at h₁ h₂ ⊢
  simpa only [Pi.add_apply, add_sub_add_comm, zero_add] using h₁.add h₂

/-- The relation `=o` is preserved by subtraction. -/
theorem IsEqUpToInfinitesimal.sub
    {X₁ X₂ Y₁ Y₂ : VariableObject E}
    (h₁ : X₁ =o Y₁) (h₂ : X₂ =o Y₂) :
    X₁ - X₂ =o Y₁ - Y₂ := by
  simpa only [sub_eq_add_neg] using h₁.add h₂.neg

/-- The relation `=o` is preserved by multiplication by a fixed scalar. -/
theorem IsEqUpToInfinitesimal.const_smul
    {𝕜 : Type u} {F : Type v}
    [NormedField 𝕜] [SeminormedAddCommGroup F] [NormedSpace 𝕜 F]
    {X Y : VariableObject F} (h : X =o Y) (c : 𝕜) :
    c • X =o c • Y := by
  rw [isEqUpToInfinitesimal_iff_tendsto_sub] at h ⊢
  simpa only [Pi.sub_apply, Pi.smul_apply, smul_sub, smul_zero] using h.const_smul c

/-- `X =o Y` iff `X` and `Y` are eventually within every positive tolerance. -/
theorem isEqUpToInfinitesimal_iff_forall_pos
    (X Y : VariableObject E) :
    (X =o Y) ↔
    ∀ δ : ℝ, 0 < δ → ∀ᶠ i in atTop, ‖X i - Y i‖ < δ := by
  rw [IsEqUpToInfinitesimal, VariableObject.isInfinitesimal_iff_forall_pos]
  simp only [Pi.sub_apply]

end InfinitesimalRelations

/-- `X ≤o Y` iff `X i ≤ Y i + δ` eventually for every `δ > 0`. -/
theorem isLEUpToInfinitesimal_iff_forall_pos (X Y : VariableObject ℝ) :
    (X ≤o Y) ↔
    ∀ δ : ℝ, 0 < δ → ∀ᶠ i in atTop, X i ≤ Y i + δ := by
  rw [IsLEUpToInfinitesimal, VariableObject.isInfinitesimal_iff_forall_pos]
  constructor
  · intro h δ hδ
    filter_upwards [h δ hδ] with i hi
    rw [Real.norm_eq_abs, abs_of_nonneg (le_max_right _ _)] at hi
    linarith [le_max_left (X i - Y i) 0]
  · intro h δ hδ
    have hδ2 : 0 < δ / 2 := by linarith
    filter_upwards [h (δ / 2) hδ2] with i hi
    rw [Real.norm_eq_abs, abs_of_nonneg (le_max_right _ _)]
    exact max_lt (by linarith) hδ

/-- For real variable objects, `X =o Y` iff `X ≤o Y` and `Y ≤o X`. -/
theorem isEqUpToInfinitesimal_iff_leUpToInfinitesimal
    (X Y : VariableObject ℝ) :
    (X =o Y) ↔ (X ≤o Y ∧ Y ≤o X) := by
  rw [isEqUpToInfinitesimal_iff_forall_pos,
    isLEUpToInfinitesimal_iff_forall_pos,
    isLEUpToInfinitesimal_iff_forall_pos]
  constructor
  · intro h
    constructor
    · intro δ hδ
      filter_upwards [h δ hδ] with i hi
      rw [Real.norm_eq_abs] at hi
      linarith [lt_of_le_of_lt (le_abs_self (X i - Y i)) hi]
    · intro δ hδ
      filter_upwards [h δ hδ] with i hi
      rw [Real.norm_eq_abs] at hi
      linarith [lt_of_le_of_lt (neg_le_abs (X i - Y i)) hi]
  · rintro ⟨hXY, hYX⟩ δ hδ
    have hδ2 : 0 < δ / 2 := by linarith
    filter_upwards [hXY (δ / 2) hδ2, hYX (δ / 2) hδ2] with i hiXY hiYX
    rw [Real.norm_eq_abs]
    rw [abs_lt]
    constructor <;> linarith

/-- **Underspill principle.** The relation `X ≤ Y + o(1)` holds if and only if
`X ≤ Y + ε + o(1)` for every fixed `ε > 0`. -/
theorem underspill (X Y : VariableObject ℝ) :
    (X ≤o Y) ↔
    (∀ ε : ℝ, ε > 0 → X ≤o (fun i => Y i + ε)) := by
  constructor
  · intro h ε hε
    apply (isLEUpToInfinitesimal_iff_forall_pos X (fun i => Y i + ε)).2
    intro δ hδ
    filter_upwards [(isLEUpToInfinitesimal_iff_forall_pos X Y).1 h δ hδ] with i hi
    linarith
  · intro h
    apply (isLEUpToInfinitesimal_iff_forall_pos X Y).2
    intro ε hε
    have hε2 : 0 < ε / 2 := by linarith
    have hbound := (isLEUpToInfinitesimal_iff_forall_pos
      X (fun i => Y i + ε / 2)).1 (h (ε / 2) hε2) (ε / 2) hε2
    filter_upwards [hbound] with i hi
    linarith

end Expdb
