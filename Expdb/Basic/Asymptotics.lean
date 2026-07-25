module

public import Expdb.Basic.Definitions
public import Mathlib.Analysis.SpecificLimits.Basic
public import Mathlib.Order.Filter.Basic
public import Mathlib.Topology.MetricSpace.Sequences

/-!
# Project-specific asymptotic notation

This module formalizes the infinitesimal comparison relation and underspill principle from
Chapter 2 of the ANTEDB blueprint. The notation `X ≤o Y` is available after
`open scoped Expdb`.

The blueprint also uses non-standard objects indexed by an ambient natural-number parameter.
These are represented by `VariableObject`. Some of their asymptotic properties are encoded as:
`VariableObject.IsBounded`, `VariableObject.IsUnbounded`,
and `VariableObject.IsInfinitesimal`.

Variable functions whose domain type may depend on the ambient parameter are represented by
`VariableFunction`. Their pointwise asymptotic properties are encoded in the
`VariableFunction` namespace.
-/

@[expose] public section

open Filter Topology

namespace Expdb

universe u v

/-- A variable object is a family indexed by the ambient natural-number parameter. -/
abbrev VariableObject (α : Type u) := ℕ → α

/-- A variable function is a family of functions whose domain type may depend on the ambient
natural-number parameter. -/
abbrev VariableFunction (domain : VariableObject (Type u)) (codomain : Type v) :=
  ∀ i, domain i → codomain

namespace VariableObject

/-- Regard a fixed object as a variable object. -/
def fixed {α : Type u} (x : α) : VariableObject α := fun _ ↦ x

/-- A fixed object has the same value at every ambient index. -/
@[simp]
theorem fixed_apply {α : Type u} (x : α) (i : ℕ) : fixed x i = x := rfl

variable {α : Type u} [SeminormedAddCommGroup α]

/-- A variable object is bounded when its norm is eventually bounded by a fixed constant. -/
def IsBounded (X : VariableObject α) : Prop :=
  ∃ C : ℝ, ∀ᶠ i in atTop, ‖X i‖ ≤ C

/-- A variable object is unbounded when its norm tends to infinity. -/
def IsUnbounded (X : VariableObject α) : Prop :=
  Tendsto (fun i ↦ ‖X i‖) atTop atTop

/-- A variable object is infinitesimal when its norm tends to zero. -/
def IsInfinitesimal (X : VariableObject α) : Prop :=
  Tendsto (fun i ↦ ‖X i‖) atTop (nhds 0)

end VariableObject

namespace VariableFunction

variable {domain : VariableObject (Type u)}
variable {α : Type v} [SeminormedAddCommGroup α]

/-- A variable function is pointwise `O(1)` when evaluation along every variable choice is
bounded. -/
def IsPointwiseBounded (f : VariableFunction domain α) : Prop :=
  ∀ x : ∀ i, domain i, VariableObject.IsBounded (fun i ↦ f i (x i))

/-- A variable function is pointwise `o(1)` when evaluation along every variable choice is
infinitesimal. -/
def IsPointwiseInfinitesimal (f : VariableFunction domain α) : Prop :=
  ∀ x : ∀ i, domain i, VariableObject.IsInfinitesimal (fun i ↦ f i (x i))

end VariableFunction

/-- The one-sided asymptotic relation `X ≤ Y + o(1)` from the blueprint.

It holds when there is a real error sequence tending to zero such that
`X i ≤ Y i + err i` eventually. Equivalently, for every fixed `δ > 0`, one eventually has
`X i ≤ Y i + δ`.

This does not assert that `X - Y` tends to zero; it only requires the positive part of `X - Y`
to tend to zero. -/
def IsLEUpToInfinitesimal (X Y : VariableObject ℝ) : Prop :=
  ∃ err : VariableObject ℝ, Tendsto err atTop (nhds 0) ∧
    ∀ᶠ i in atTop, X i ≤ Y i + err i

/-- `X ≤o Y` denotes the complete blueprint expression `X ≤ Y + o(1)`; it is not the
little-o relation `X = o(Y)`. -/
scoped[Expdb] notation X " ≤o " Y => IsLEUpToInfinitesimal X Y

open scoped Expdb

/-- The relation `X ≤ Y + o(1)` is equivalent to `X i ≤ Y i + δ` eventually for every fixed
positive `δ`. -/
theorem isLEUpToInfinitesimal_iff_forall_pos (X Y : VariableObject ℝ) :
    (X ≤o Y) ↔
    ∀ δ : ℝ, 0 < δ → ∀ᶠ i in atTop, X i ≤ Y i + δ := by
  constructor
  · rintro ⟨err, herr, hXY⟩ δ hδ
    rw [Metric.tendsto_nhds] at herr
    have herr_small := herr δ hδ
    filter_upwards [hXY, herr_small] with i hi hierr
    rw [Real.dist_eq, sub_zero] at hierr
    have herr_lt : err i < δ := lt_of_le_of_lt (le_abs_self _) hierr
    linarith
  · intro h
    refine ⟨fun i => max (X i - Y i) 0, ?_, Filter.Eventually.of_forall fun i => ?_⟩
    · rw [Metric.tendsto_nhds]
      intro δ hδ
      have hδ2 : 0 < δ / 2 := by linarith
      filter_upwards [h (δ / 2) hδ2] with i hi
      rw [Real.dist_eq, sub_zero, abs_of_nonneg (le_max_right _ _)]
      exact max_lt (by linarith) hδ
    · have hi : X i - Y i ≤ max (X i - Y i) 0 := le_max_left _ _
      linarith

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
