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
These are represented by `Variable`. Their asymptotic properties can be expressed using Mathlib's
filter API:
- a bounded variable `X` satisfies
 `∃ C : ℝ, ∀ᶠ i in atTop, ‖X i‖ ≤ C`;
- an unbounded variable `X` satisfies
 `Tendsto (fun i => ‖X i‖) atTop atTop`;
- an infinitesimal variable `X` satisfies
 `Tendsto X atTop (nhds 0)`.

If these conditions recur sufficiently often in later chapters, they may be given the
project-specific names `IsBoundedVariable`, `IsUnboundedVariable`, and
`IsInfinitesimalVariable`.
-/

@[expose] public section

open Filter Topology

namespace Expdb

universe u

/-- A variable object is a family indexed by the ambient natural-number parameter. -/
abbrev Variable (α : Type u) := ℕ → α

namespace Variable

/-- Regard a fixed object as a constant variable object. -/
def constant {α : Type u} (x : α) : Variable α := fun _ ↦ x

/-- A constant variable object has the same value at every ambient index. -/
@[simp]
theorem constant_apply {α : Type u} (x : α) (i : ℕ) : constant x i = x := rfl

end Variable

/-- The one-sided asymptotic relation `X ≤ Y + o(1)` from the blueprint.

It holds when there is a real error sequence tending to zero such that
`X i ≤ Y i + err i` eventually. Equivalently, for every fixed `δ > 0`, one eventually has
`X i ≤ Y i + δ`.

This does not assert that `X - Y` tends to zero; it only requires the positive part of `X - Y`
to tend to zero. -/
def IsLEUpToInfinitesimal (X Y : Variable ℝ) : Prop :=
  ∃ err : Variable ℝ, Tendsto err atTop (nhds 0) ∧
    ∀ᶠ i in atTop, X i ≤ Y i + err i

/-- `X ≤o Y` denotes the complete blueprint expression `X ≤ Y + o(1)`; it is not the
little-o relation `X = o(Y)`. -/
scoped[Expdb] notation X " ≤o " Y => IsLEUpToInfinitesimal X Y

open scoped Expdb

/-- The relation `X ≤ Y + o(1)` is equivalent to `X i ≤ Y i + δ` eventually for every fixed
positive `δ`. -/
theorem isLEUpToInfinitesimal_iff_forall_pos (X Y : Variable ℝ) :
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
theorem underspill (X Y : Variable ℝ) :
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
