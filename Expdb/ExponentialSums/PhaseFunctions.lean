module

public import Expdb.Basic.Asymptotics
public import Mathlib.Analysis.Calculus.IteratedDeriv.Defs
public import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Expdb.Basic.AutomaticUniformity
import Expdb.Mathlib.IteratedDeriv
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv

/-!
# Phase functions

This module formalizes the phase function definitions from the blueprint's Exponential sum growth
exponents chapter (`beta-chapter`).
A phase function is smooth on `[1, 2]`. A model phase function is a variable family
of phase functions whose successive derivatives asymptotically agree with those of
`u ↦ u ^ (-σ)` for some fixed positive exponent `σ`.
-/

@[expose] public section

open Filter Topology
open scoped ContDiff

noncomputable section

namespace Expdb

/-- The interval on which phase functions are considered. -/
def phaseInterval : Set ℝ := Set.Icc 1 2

/-- A phase function is a variable real-valued function that is smooth on `[1, 2]` at every
ambient index. -/
def IsPhaseFunction (F : VariableFunction (VariableObject.fixed ℝ) ℝ) : Prop :=
  ∀ i : ℕ, ContDiffOn ℝ ∞ (F i) phaseInterval

/-- The reference phase `u ↦ u ^ (-σ)` whose derivatives a model phase follows. -/
def modelPhase (σ : ℝ) : ℝ → ℝ := fun u ↦ u ^ (-σ)

/-- The fixed-data model-phase error at derivative order `p`. -/
def modelPhaseErrorAt (F : ℝ → ℝ) (σ : ℝ) (p : ℕ) (u : ℝ) : ℝ :=
  iteratedDerivWithin (p + 1) F phaseInterval u -
    iteratedDerivWithin p (modelPhase σ) phaseInterval u

/-- The variable error function in the model phase condition at derivative order `p`. -/
def modelPhaseError
    (F : VariableFunction (VariableObject.fixed ℝ) ℝ) (σ : ℝ) (p : ℕ) :
    VariableFunction (VariableObject.fixed phaseInterval) ℝ :=
  fun i u ↦ modelPhaseErrorAt (F i) σ p u

@[simp] theorem modelPhaseError_apply
    (F : VariableFunction (VariableObject.fixed ℝ) ℝ) (σ : ℝ) (p i : ℕ)
    (u : phaseInterval) :
    modelPhaseError F σ p i u = modelPhaseErrorAt (F i) σ p u :=
  rfl

/-- The phase interval has unique derivatives. -/
theorem uniqueDiffOn_phaseInterval : UniqueDiffOn ℝ phaseInterval :=
  uniqueDiffOn_Icc (by norm_num [phaseInterval])

/-- Every reference model phase is smooth on the phase interval. -/
theorem contDiffOn_modelPhase (σ : ℝ) :
    ContDiffOn ℝ ∞ (modelPhase σ) phaseInterval := by
  intro u hu
  exact (Real.contDiffAt_rpow_const_of_ne
    (ne_of_gt (zero_lt_one.trans_le hu.1))).contDiffWithinAt

/-- The derivatives of the reference model phase have the expected falling-factorial
formula on the phase interval. -/
theorem iteratedDerivWithin_modelPhase (σ : ℝ) (p : ℕ) {u : ℝ}
    (hu : u ∈ phaseInterval) :
    iteratedDerivWithin p (modelPhase σ) phaseInterval u =
      (descPochhammer ℝ p).eval (-σ) * u ^ (-σ - p) := by
  have hu0 : u ≠ 0 := ne_of_gt (zero_lt_one.trans_le hu.1)
  have hsmooth : ContDiffAt ℝ p (modelPhase σ) u := by
    change ContDiffAt ℝ p (fun x : ℝ ↦ x ^ (-σ)) u
    exact Real.contDiffAt_rpow_const_of_ne hu0
  rw [iteratedDerivWithin_eq_iteratedDeriv uniqueDiffOn_phaseInterval hsmooth hu,
    iteratedDeriv_eq_iterate]
  exact Real.iter_deriv_rpow_const (-σ) u p

/-- On `[1, 2]`, the falling-factorial coefficient bounds the corresponding derivative
of a model phase with nonnegative exponent. -/
theorem norm_iteratedDerivWithin_modelPhase_le {σ : ℝ} (hσ : 0 ≤ σ) (p : ℕ)
    {u : ℝ} (hu : u ∈ phaseInterval) :
    ‖iteratedDerivWithin p (modelPhase σ) phaseInterval u‖ ≤
      ‖(descPochhammer ℝ p).eval (-σ)‖ := by
  rw [iteratedDerivWithin_modelPhase σ p hu, norm_mul,
    Real.norm_of_nonneg (Real.rpow_nonneg (zero_le_one.trans hu.1) _)]
  simpa only [mul_one] using mul_le_mul_of_nonneg_left
    (Real.rpow_le_one_of_one_le_of_nonpos hu.1 (by
      have hp : 0 ≤ (p : ℝ) := Nat.cast_nonneg p
      linarith))
    (norm_nonneg ((descPochhammer ℝ p).eval (-σ)))

/-- The iterated derivatives of a reference model phase vary continuously on the phase
interval. -/
theorem continuousOn_iteratedDerivWithin_modelPhase (σ : ℝ) (p : ℕ) :
    ContinuousOn (iteratedDerivWithin p (modelPhase σ) phaseInterval) phaseInterval :=
  (contDiffOn_modelPhase σ).continuousOn_iteratedDerivWithin
    (ENat.natCast_le_of_coe_top_le_withTop le_rfl p) uniqueDiffOn_phaseInterval

/-- A fixed phase function whose model-phase errors through order `P` are at most `δ`. -/
def IsApproximateModelPhaseFunction
    (F : ℝ → ℝ) (σ : ℝ) (P : ℕ) (δ : ℝ) : Prop :=
  ContDiffOn ℝ ∞ F phaseInterval ∧
    ∀ p : ℕ, p ≤ P → ∀ u : phaseInterval,
      ‖modelPhaseErrorAt F σ p u‖ ≤ δ

/-- The derivatives of orders `1, …, n` of `F` are bounded by `K` on the phase interval. -/
def HasPhaseDerivBound (F : ℝ → ℝ) (n : ℕ) (K : ℝ) : Prop :=
  ∀ k : ℕ, 1 ≤ k → k ≤ n → ∀ u ∈ phaseInterval,
    ‖iteratedDerivWithin k F phaseInterval u‖ ≤ K

/-- The first derivative of `F` is at least `c` on the phase interval. -/
def HasPhaseFirstDerivLowerBound (F : ℝ → ℝ) (c : ℝ) : Prop :=
  ∀ u ∈ phaseInterval, c ≤ iteratedDerivWithin 1 F phaseInterval u

/-- A derivative bound may be restricted to a smaller order. -/
theorem HasPhaseDerivBound.of_le {F : ℝ → ℝ} {m n : ℕ} {K : ℝ}
    (hF : HasPhaseDerivBound F n K) (hmn : m ≤ n) :
    HasPhaseDerivBound F m K :=
  fun k hk hkm u hu ↦ hF k hk (hkm.trans hmn) u hu

/-- The data defining a model phase when its reference exponent `σ` has already been fixed. -/
structure IsModelPhaseFunctionWith
    (F : VariableFunction (VariableObject.fixed ℝ) ℝ) (σ : ℝ) : Prop where
  /-- Every member of the family is smooth on the phase interval. -/
  isPhaseFunction : IsPhaseFunction F
  /-- Every fixed-order model-phase error is choicewise infinitesimal. -/
  error_isChoicewiseInfinitesimal :
    ∀ p : ℕ, (modelPhaseError F σ p).IsChoicewiseInfinitesimal

/-- A variable phase function is a model phase function when, for some fixed `σ > 0`,
the error between its `(p + 1)`st derivative and the `p`th derivative of `u ↦ u ^ (-σ)` is
choicewise infinitesimal for every fixed derivative order `p`. -/
def IsModelPhaseFunction (F : VariableFunction (VariableObject.fixed ℝ) ℝ) : Prop :=
  IsPhaseFunction F ∧
    ∃ σ : ℝ, 0 < σ ∧
      ∀ p : ℕ, (modelPhaseError F σ p).IsChoicewiseInfinitesimal

/-- Fixed model-phase data with a positive reference exponent gives a model phase function. -/
theorem IsModelPhaseFunctionWith.isModelPhaseFunction
    {F : VariableFunction (VariableObject.fixed ℝ) ℝ} {σ : ℝ}
    (hF : IsModelPhaseFunctionWith F σ) (hσ : 0 < σ) :
    IsModelPhaseFunction F :=
  ⟨hF.isPhaseFunction, σ, hσ, hF.error_isChoicewiseInfinitesimal⟩

/-- A model phase function admits fixed-exponent model-phase data. -/
theorem IsModelPhaseFunction.exists_with
    {F : VariableFunction (VariableObject.fixed ℝ) ℝ}
    (hF : IsModelPhaseFunction F) :
    ∃ σ : ℝ, 0 < σ ∧ IsModelPhaseFunctionWith F σ := by
  rcases hF with ⟨hphase, σ, hσ, herror⟩
  exact ⟨σ, hσ, ⟨hphase, herror⟩⟩

/-- Adding a linear function changes only the zeroth model-phase error. -/
theorem modelPhaseErrorAt_add_linear
    {F : ℝ → ℝ} {σ c u : ℝ} {p : ℕ}
    (hF : ContDiffOn ℝ ∞ F phaseInterval) (hu : u ∈ phaseInterval) :
    modelPhaseErrorAt (fun x ↦ F x + c * x) σ p u =
      modelPhaseErrorAt F σ p u + if p = 0 then c else 0 := by
  change modelPhaseErrorAt (F + fun x ↦ c * x) σ p u = _
  rw [modelPhaseErrorAt, modelPhaseErrorAt,
    iteratedDerivWithin_add hu uniqueDiffOn_phaseInterval
      ((hF u hu).of_le
        (ENat.natCast_le_of_coe_top_le_withTop le_rfl (p + 1)))
      ((show ContDiffOn ℝ ∞ (fun x : ℝ ↦ c * x) phaseInterval by fun_prop)
        u hu |>.of_le
          (ENat.natCast_le_of_coe_top_le_withTop le_rfl (p + 1)))]
  rw [iteratedDerivWithin_const_mul_field,
    iteratedDerivWithin_fun_id hu uniqueDiffOn_phaseInterval]
  by_cases hp : p = 0
  · subst p
    simp
    ring
  · simp only [if_neg (by omega : p + 1 ≠ 0), if_neg (by omega : p + 1 ≠ 1),
      mul_zero, add_zero, hp]
    simp

/-- Adding a controlled linear term increases the model-phase error by at most the norm of its
coefficient. -/
theorem IsApproximateModelPhaseFunction.add_linear
    {F : ℝ → ℝ} {σ δ η c : ℝ} {P : ℕ}
    (hF : IsApproximateModelPhaseFunction F σ P δ) (hc : ‖c‖ ≤ η) :
    IsApproximateModelPhaseFunction (fun u ↦ F u + c * u) σ P (δ + η) := by
  refine ⟨hF.1.add (by fun_prop), ?_⟩
  intro p hp u
  have hη : 0 ≤ η := (norm_nonneg c).trans hc
  rw [modelPhaseErrorAt_add_linear hF.1 u.property]
  calc
    ‖modelPhaseErrorAt F σ p u + if p = 0 then c else 0‖ ≤
        ‖modelPhaseErrorAt F σ p u‖ + ‖if p = 0 then c else 0‖ :=
      norm_add_le _ _
    _ ≤ δ + η := add_le_add (hF.2 p hp u) (by
      by_cases h : p = 0
      · simpa [h] using hc
      · simp [h, hη])

/-- Model-phase data is stable under addition of an infinitesimal linear term. -/
theorem IsModelPhaseFunctionWith.add_infinitesimal_linear
    {F : VariableFunction (VariableObject.fixed ℝ) ℝ} {σ : ℝ}
    (hF : IsModelPhaseFunctionWith F σ)
    {c : VariableObject ℝ} (hc : c.IsInfinitesimal) :
    IsModelPhaseFunctionWith (fun i u ↦ F i u + c i * u) σ := by
  refine ⟨fun i ↦ (hF.isPhaseFunction i).add (by fun_prop), ?_⟩
  intro p u
  rw [VariableObject.IsInfinitesimal]
  have herr := hF.error_isChoicewiseInfinitesimal p u
  rw [VariableObject.IsInfinitesimal] at herr
  have hc' : Tendsto (fun i ↦ ‖c i‖) atTop (𝓝 0) := hc
  apply squeeze_zero' (g := fun i ↦
    ‖modelPhaseErrorAt (F i) σ p (u i)‖ + ‖c i‖)
  · exact Filter.Eventually.of_forall fun _ ↦ norm_nonneg _
  · filter_upwards [] with i
    rw [modelPhaseError_apply,
      modelPhaseErrorAt_add_linear (hF.isPhaseFunction i) (u i).property]
    split_ifs
    · exact norm_add_le _ _
    · simp
  · simpa [modelPhaseError] using herr.add hc'

/-- Model-phase data is stable under affine reparametrizations that tend to the identity and map
the phase interval to itself. -/
theorem IsModelPhaseFunctionWith.comp_affine_tendsto_id
    {F : VariableFunction (VariableObject.fixed ℝ) ℝ} {σ : ℝ}
    (hF : IsModelPhaseFunctionWith F σ) {c d : VariableObject ℝ}
    (hc : Tendsto c atTop (𝓝 1)) (hd : Tendsto d atTop (𝓝 0))
    (hmap : ∀ i, Set.MapsTo (fun u ↦ c i * u + d i) phaseInterval phaseInterval) :
    IsModelPhaseFunctionWith (fun i u ↦ F i (c i * u + d i)) σ := by
  refine ⟨fun i ↦ (hF.isPhaseFunction i).comp (by fun_prop) (hmap i), ?_⟩
  intro p u
  rw [VariableObject.IsInfinitesimal]
  let v : ∀ i, phaseInterval :=
    fun i ↦ ⟨c i * (u i : ℝ) + d i, hmap i (u i).property⟩
  have herr := hF.error_isChoicewiseInfinitesimal p v
  rw [VariableObject.IsInfinitesimal] at herr
  have hcsub : Tendsto (fun i ↦ c i - 1) atTop (𝓝 0) := by
    simpa using hc.sub
      (tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ (1 : ℝ)) atTop (𝓝 1))
  have huBounded : IsBoundedUnder (· ≤ ·) atTop
      ((‖·‖) ∘ fun i ↦ (u i : ℝ)) :=
    Filter.isBoundedUnder_of_eventually_le <| Filter.Eventually.of_forall fun i ↦ by
      rw [Function.comp_apply, Real.norm_eq_abs, abs_of_nonneg
        (zero_le_one.trans (u i).property.1)]
      exact (u i).property.2
  have hmul : Tendsto (fun i ↦ (c i - 1) * (u i : ℝ)) atTop (𝓝 0) :=
    hcsub.zero_mul_isBoundedUnder_le huBounded
  have hvsub : Tendsto (fun i ↦ (v i : ℝ) - (u i : ℝ)) atTop (𝓝 0) := by
    convert hmul.add hd using 1
    · ext i
      dsimp [v]
      ring_nf
    · ring_nf
  have hcpow : Tendsto (fun i ↦ c i ^ (p + 1)) atTop (𝓝 1) := by
    simpa using hc.pow (p + 1)
  have href := continuousOn_iteratedDerivWithin_modelPhase σ p
  have hrefUniform := isCompact_Icc.uniformContinuousOn_of_continuous href
  have hrefclose : Tendsto (fun i ↦
      iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) -
        iteratedDerivWithin p (modelPhase σ) phaseInterval (u i)) atTop (𝓝 0) := by
    rw [Metric.uniformContinuousOn_iff] at hrefUniform
    rw [Metric.tendsto_nhds]
    intro ε hε
    obtain ⟨δ, hδ, hclose⟩ := hrefUniform ε hε
    have hevent : ∀ᶠ i in atTop, ‖(v i : ℝ) - (u i : ℝ)‖ < δ := by
      simpa [Real.dist_eq] using (Metric.tendsto_nhds.1 hvsub) δ hδ
    filter_upwards [hevent] with i hi
    simpa [dist_eq_norm] using
      hclose (v i) (v i).property (u i) (u i).property
        (by simpa [Real.dist_eq] using hi)
  have hzero : Tendsto (fun i ↦ c i ^ (p + 1) - 1) atTop (𝓝 0) := by
    simpa using hcpow.sub
      (tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ (1 : ℝ)) atTop (𝓝 1))
  obtain ⟨D, hD⟩ := isCompact_Icc.exists_bound_of_continuousOn href
  have hrefBounded : IsBoundedUnder (· ≤ ·) atTop
      ((‖·‖) ∘ fun i ↦
        iteratedDerivWithin p (modelPhase σ) phaseInterval (v i)) :=
    Filter.isBoundedUnder_of_eventually_le <| Filter.Eventually.of_forall fun i ↦
      hD (v i) (v i).property
  have hscale : Tendsto (fun i ↦
      (c i ^ (p + 1) - 1) *
        iteratedDerivWithin p (modelPhase σ) phaseInterval (v i)) atTop (𝓝 0) :=
    hzero.zero_mul_isBoundedUnder_le hrefBounded
  have hmajor : Tendsto (fun i ↦
      ‖c i ^ (p + 1)‖ * ‖modelPhaseErrorAt (F i) σ p (v i)‖ +
        ‖(c i ^ (p + 1) - 1) *
            iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) +
          (iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) -
            iteratedDerivWithin p (modelPhase σ) phaseInterval (u i))‖)
      atTop (𝓝 0) := by
    simpa using (hcpow.norm.mul herr).add (hscale.add hrefclose).norm
  apply squeeze_zero' (g := fun i ↦
    ‖c i ^ (p + 1)‖ * ‖modelPhaseErrorAt (F i) σ p (v i)‖ +
      ‖(c i ^ (p + 1) - 1) *
          iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) +
        (iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) -
          iteratedDerivWithin p (modelPhase σ) phaseInterval (u i))‖)
  · exact Filter.Eventually.of_forall fun _ ↦ norm_nonneg _
  · filter_upwards [] with i
    change ‖modelPhaseErrorAt (fun x ↦ F i (c i * x + d i)) σ p (u i)‖ ≤ _
    rw [modelPhaseErrorAt,
      iteratedDerivWithin_comp_affine_of_mapsTo (s := phaseInterval) (t := phaseInterval)
        ((hF.isPhaseFunction i).of_le
          (ENat.natCast_le_of_coe_top_le_withTop le_rfl (p + 1)))
        uniqueDiffOn_phaseInterval uniqueDiffOn_phaseInterval (u i).property (hmap i)]
    change ‖c i ^ (p + 1) *
        iteratedDerivWithin (p + 1) (F i) phaseInterval (v i) - _‖ ≤ _
    rw [show iteratedDerivWithin (p + 1) (F i) phaseInterval (v i) =
        modelPhaseErrorAt (F i) σ p (v i) +
          iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) by
      simp [modelPhaseErrorAt]]
    calc
      _ = ‖c i ^ (p + 1) * modelPhaseErrorAt (F i) σ p (v i) +
          ((c i ^ (p + 1) - 1) *
              iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) +
            (iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) -
              iteratedDerivWithin p (modelPhase σ) phaseInterval (u i)))‖ := by
        congr 1
        ring
      _ ≤ ‖c i ^ (p + 1)‖ * ‖modelPhaseErrorAt (F i) σ p (v i)‖ +
          ‖(c i ^ (p + 1) - 1) *
              iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) +
            (iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) -
              iteratedDerivWithin p (modelPhase σ) phaseInterval (u i))‖ := by
        simpa only [norm_mul] using norm_add_le
          (c i ^ (p + 1) * modelPhaseErrorAt (F i) σ p (v i))
          ((c i ^ (p + 1) - 1) *
              iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) +
            (iteratedDerivWithin p (modelPhase σ) phaseInterval (v i) -
              iteratedDerivWithin p (modelPhase σ) phaseInterval (u i)))
      _ = _ := by rfl
  · exact hmajor

/-- At every fixed finite order, fixed-exponent model-phase data is eventually a uniformly
approximate model phase. -/
theorem IsModelPhaseFunctionWith.eventually_isApproximate
    {F : VariableFunction (VariableObject.fixed ℝ) ℝ} {σ : ℝ}
    (hF : IsModelPhaseFunctionWith F σ) (P : ℕ) {δ : ℝ} (hδ : 0 < δ) :
    ∀ᶠ i in atTop, IsApproximateModelPhaseFunction (F i) σ P δ := by
  have hp : ∀ p ∈ Finset.range (P + 1), ∀ᶠ i in atTop,
      ∀ u : phaseInterval, ‖modelPhaseError F σ p i u‖ < δ :=
    fun p _ ↦
      (VariableFunction.isChoicewiseInfinitesimal_iff_forall_pos_uniform
        (VariableObject.fixed phaseInterval)
        (fun _ ↦ ⟨1, by simp [phaseInterval]⟩)
        (modelPhaseError F σ p)).1
        (hF.error_isChoicewiseInfinitesimal p) δ hδ
  have hall :
      ∀ᶠ i in atTop, ∀ p ∈ Finset.range (P + 1),
        ∀ u : phaseInterval, ‖modelPhaseError F σ p i u‖ < δ :=
    (Finset.range (P + 1)).eventually_all.mpr hp
  filter_upwards [hall] with i hi
  refine ⟨hF.isPhaseFunction i, ?_⟩
  intro p hp' u
  exact (hi p (Finset.mem_range.mpr (by omega)) u).le

private lemma modelPhaseError_sum_isChoicewiseInfinitesimal
    {F : VariableFunction (VariableObject.fixed ℝ) ℝ} {σ : ℝ}
    (herror : ∀ p : ℕ, (modelPhaseError F σ p).IsChoicewiseInfinitesimal)
    (P : ℕ) :
    VariableFunction.IsChoicewiseInfinitesimal
      ((fun i (u : phaseInterval) ↦
        ∑ p ∈ Finset.range (P + 1), ‖modelPhaseError F σ p i u‖) :
        VariableFunction (VariableObject.fixed phaseInterval) ℝ) := by
  intro u
  rw [VariableObject.IsInfinitesimal]
  have hsum : Tendsto
      (fun i ↦ ∑ p ∈ Finset.range (P + 1), ‖modelPhaseError F σ p i (u i)‖)
      atTop (nhds 0) := by
    simpa using tendsto_finsetSum (Finset.range (P + 1)) fun p _ ↦ herror p u
  simpa using hsum.norm

/-- For any fixed derivative cutoff, the errors of a model phase function have a common
infinitesimal bound, uniform over the phase interval and all derivative orders up to the cutoff,
after passing to a subsequence. -/
theorem IsModelPhaseFunction.exists_subsequence_uniform_error
    {F : VariableFunction (VariableObject.fixed ℝ) ℝ}
    (hF : IsModelPhaseFunction F) (P : ℕ) :
    ∃ σ : ℝ, 0 < σ ∧
      ∃ φ : ℕ → ℕ, StrictMono φ ∧
        ∃ c : VariableObject ℝ, c.IsInfinitesimal ∧
          ∀ i p, p ≤ P →
            ∀ u : phaseInterval,
              ‖modelPhaseError F σ p (φ i) u‖ ≤ c i := by
  rcases hF with ⟨_, σ, hσ, herror⟩
  let g : VariableFunction (VariableObject.fixed phaseInterval) ℝ :=
    fun i u ↦ ∑ p ∈ Finset.range (P + 1), ‖modelPhaseError F σ p i u‖
  have hg : g.IsChoicewiseInfinitesimal :=
    modelPhaseError_sum_isChoicewiseInfinitesimal herror P
  obtain ⟨φ, hφ, c, hc, hbound⟩ :=
    automatic_uniformity_of_choicewise_infinitesimal
      (E := VariableObject.fixed phaseInterval)
      (fun _ ↦ ⟨1, by simp [phaseInterval]⟩) g hg
  refine ⟨σ, hσ, φ, hφ, c, hc, ?_⟩
  intro i p hp u
  have hmem : p ∈ Finset.range (P + 1) := Finset.mem_range.mpr (by omega)
  calc
    ‖modelPhaseError F σ p (φ i) u‖ ≤
        ∑ q ∈ Finset.range (P + 1), ‖modelPhaseError F σ q (φ i) u‖ :=
      Finset.single_le_sum
        (fun q _ ↦ norm_nonneg (modelPhaseError F σ q (φ i) u)) hmem
    _ = ‖g (φ i) u‖ := by
      rw [Real.norm_eq_abs, abs_of_nonneg (Finset.sum_nonneg fun _ _ ↦ norm_nonneg _)]
    _ ≤ c i := hbound i u

end Expdb
