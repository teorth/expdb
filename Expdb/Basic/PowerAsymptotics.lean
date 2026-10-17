module

public import Expdb.Basic.Asymptotics
public import Mathlib.Analysis.Asymptotics.Defs
public import Mathlib.Analysis.Asymptotics.Lemmas
public import Mathlib.Analysis.SpecialFunctions.Log.Base
public import Mathlib.Analysis.SpecialFunctions.Pow.Asymptotics

/-!
# Power asymptotics

This module collects the project-wide notation and elementary calculus for quantities of the
form `T ^ (α + o(1))`. It is independent of exponential sums and can therefore be reused by
all later asymptotic arguments.
-/

@[expose] public section

open Filter Topology
open scoped Expdb

noncomputable section

namespace Expdb

/-! ## Definitions -/

/-- `N = T ^ (α + o(1))`, expressed using a variable exponent that is equal to `α` up to an
infinitesimal. -/
def IsPowerAsymptotic (N T : VariableObject ℝ) (α : ℝ) : Prop :=
  ∃ exponent : VariableObject ℝ,
    exponent =o VariableObject.fixed α ∧
      ∀ᶠ i in atTop, N i = T i ^ exponent i

/-- `X ≪ T ^ (β + o(1))`, expressed using a variable exponent that is equal to `β` up to
an infinitesimal. -/
def IsPowerBounded {E : Type*} [SeminormedAddCommGroup E]
    (X : VariableObject E) (T : VariableObject ℝ) (β : ℝ) : Prop :=
  ∃ exponent : VariableObject ℝ,
    exponent =o VariableObject.fixed β ∧
      Asymptotics.IsBigO atTop X (fun i ↦ T i ^ exponent i)

/-! ## Power asymptotics -/

/-- Restrict a power asymptotic along a sequence tending to infinity. -/
theorem IsPowerAsymptotic.comp_tendsto
    {N T : VariableObject ℝ} {α : ℝ} (h : IsPowerAsymptotic N T α)
    {φ : ℕ → ℕ} (hφ : Tendsto φ atTop atTop) :
    IsPowerAsymptotic (N ∘ φ) (T ∘ φ) α := by
  obtain ⟨e, he, heq⟩ := h
  refine ⟨e ∘ φ, ?_, hφ.eventually heq⟩
  rw [IsEqUpToInfinitesimal, VariableObject.IsInfinitesimal] at he ⊢
  exact he.comp hφ

private theorem isEqUpToInfinitesimal_fixed_iff_tendsto {e : VariableObject ℝ} {α : ℝ} :
    e =o VariableObject.fixed α ↔ Tendsto e atTop (𝓝 α) := by
  rw [IsEqUpToInfinitesimal, VariableObject.IsInfinitesimal]
  constructor
  · intro he
    have hzero : Tendsto (fun i ↦ e i - α) atTop (𝓝 0) :=
      tendsto_zero_iff_norm_tendsto_zero.mpr he
    simpa using hzero.add (tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ α) atTop (𝓝 α))
  · intro he
    simpa using (he.sub (tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ α) atTop (𝓝 α))).norm

/-- Changing to a height of positive growth divides the scale exponent. -/
theorem IsPowerAsymptotic.change_base {N S T : VariableObject ℝ} {α κ : ℝ}
    (hN : IsPowerAsymptotic N T α) (hS : IsPowerAsymptotic S T κ)
    (hκ : 0 < κ) (hT : ∀ i, 1 ≤ T i) :
    IsPowerAsymptotic N S (α / κ) := by
  obtain ⟨a, ha, hNa⟩ := hN
  obtain ⟨b, hb, hSb⟩ := hS
  have hatop := isEqUpToInfinitesimal_fixed_iff_tendsto.mp ha
  have hbtop := isEqUpToInfinitesimal_fixed_iff_tendsto.mp hb
  refine ⟨a / b, isEqUpToInfinitesimal_fixed_iff_tendsto.mpr (hatop.div hbtop hκ.ne'), ?_⟩
  have hbpos := hbtop.eventually (lt_mem_nhds hκ)
  filter_upwards [hNa, hSb, hbpos] with i hiN hiS hiB
  rw [hiN, hiS, ← Real.rpow_mul (zero_le_one.trans (hT i))]
  congr 1
  dsimp
  field_simp

/-- If `N = T ^ (α + o(1))`, then `N` eventually lies between the powers with exponents
`α - ε` and `α + ε`. -/
theorem IsPowerAsymptotic.eventually_between
    {N T : VariableObject ℝ} {α ε : ℝ}
    (hNT : IsPowerAsymptotic N T α)
    (hT : ∀ᶠ i in atTop, 1 ≤ T i)
    (hε : 0 < ε) :
    ∀ᶠ i in atTop, T i ^ (α - ε) ≤ N i ∧ N i ≤ T i ^ (α + ε) := by
  obtain ⟨exponent, hexponent, hN⟩ := hNT
  obtain ⟨hexponent_le, halpha_le⟩ :=
    (isEqUpToInfinitesimal_iff_leUpToInfinitesimal
      exponent (VariableObject.fixed α)).1 hexponent
  have hupper := (isLEUpToInfinitesimal_iff_forall_pos
    exponent (VariableObject.fixed α)).1 hexponent_le ε hε
  have hlower := (isLEUpToInfinitesimal_iff_forall_pos
    (VariableObject.fixed α) exponent).1 halpha_le ε hε
  filter_upwards [hT, hN, hupper, hlower] with i hiT hiN hiupper hilower
  rw [hiN]
  exact
    ⟨Real.rpow_le_rpow_of_exponent_le hiT (by linarith),
      Real.rpow_le_rpow_of_exponent_le hiT (by linarith)⟩

/-- Convergence of the logarithmic exponent gives a power asymptotic. -/
theorem isPowerAsymptotic_of_logb_tendsto
    {N T : VariableObject ℝ} {α : ℝ}
    (hT : ∀ᶠ i in atTop, 1 < T i) (hN : ∀ᶠ i in atTop, 0 < N i)
    (hexponent : Tendsto (fun i ↦ Real.logb (T i) (N i)) atTop (nhds α)) :
    IsPowerAsymptotic N T α := by
  let exponent : VariableObject ℝ := fun i ↦ Real.logb (T i) (N i)
  refine ⟨exponent, ?_, ?_⟩
  · rw [IsEqUpToInfinitesimal, VariableObject.IsInfinitesimal]
    have hsub : Tendsto (fun i ↦ Real.logb (T i) (N i) - α) atTop (nhds 0) := by
      simpa using hexponent.sub
        (tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ α) atTop (nhds α))
    convert (continuous_norm.tendsto 0).comp hsub using 1 <;>
      simp [Function.comp_def, exponent, VariableObject.fixed]
  · filter_upwards [hT, hN] with i hiT hiN
    exact (Real.rpow_logb (zero_lt_one.trans hiT) hiT.ne' hiN).symm

/-- A variable power sandwich with an infinitesimal error determines a power asymptotic. No
sign condition on the error is needed. -/
theorem isPowerAsymptotic_of_between_of_infinitesimal
    {N T δ : VariableObject ℝ} {α : ℝ}
    (hT : ∀ i, 1 < T i) (hN : ∀ i, 0 < N i)
    (hδ : δ.IsInfinitesimal)
    (hbetween : ∀ i,
      T i ^ (α - δ i) ≤ N i ∧ N i ≤ T i ^ (α + δ i)) :
    IsPowerAsymptotic N T α := by
  let exponent : VariableObject ℝ := fun i ↦ Real.logb (T i) (N i)
  refine ⟨exponent, ?_, Filter.Eventually.of_forall fun i ↦ ?_⟩
  · apply (isEqUpToInfinitesimal_iff_forall_pos
    exponent (VariableObject.fixed α)).2
    intro ε hε
    have hδsmall :=
      (VariableObject.isInfinitesimal_iff_forall_pos δ).1 hδ ε hε
    filter_upwards [hδsmall] with i hi
    have hlower : α - δ i ≤ exponent i :=
      (Real.le_logb_iff_rpow_le (hT i) (hN i)).2 (hbetween i).1
    have hupper : exponent i ≤ α + δ i :=
      (Real.logb_le_iff_le_rpow (hT i) (hN i)).2 (hbetween i).2
    rw [Real.norm_eq_abs, abs_lt]
    rw [Real.norm_eq_abs] at hi
    have hδupper : δ i < ε := (le_abs_self (δ i)).trans_lt hi
    constructor <;> dsimp [exponent, VariableObject.fixed] at * <;> linarith
  · exact (Real.rpow_logb (lt_trans zero_lt_one (hT i))
      (hT i).ne' (hN i)).symm

/-- A variable power sandwich with a nonnegative infinitesimal error determines a power
asymptotic. This is the compatibility form of
`isPowerAsymptotic_of_between_of_infinitesimal`. -/
theorem isPowerAsymptotic_of_between
    {N T δ : VariableObject ℝ} {α : ℝ}
    (hT : ∀ i, 1 < T i) (hN : ∀ i, 0 < N i)
    (_hδnonneg : ∀ i, 0 ≤ δ i) (hδ : δ.IsInfinitesimal)
    (hbetween : ∀ i,
      T i ^ (α - δ i) ≤ N i ∧ N i ≤ T i ^ (α + δ i)) :
    IsPowerAsymptotic N T α :=
  isPowerAsymptotic_of_between_of_infinitesimal hT hN hδ hbetween

/-- If `T` lies between fixed positive multiples of `N ^ α⁻¹`, then the logarithm of `N` to
base `T` tends to `α`. -/
private theorem tendsto_logb_of_between_const_rpow
    {N T : VariableObject ℝ} {α A B : ℝ}
    (hα : 0 < α) (hA : 0 < A) (hB : 0 < B)
    (hNtop : Tendsto N atTop atTop)
    (hT : ∀ᶠ i in atTop, A * N i ^ α⁻¹ ≤ T i ∧ T i ≤ B * N i ^ α⁻¹) :
    Tendsto (fun i ↦ Real.logb (T i) (N i)) atTop (nhds α) := by
  have hlogNtop : Tendsto (fun i ↦ Real.log (N i)) atTop atTop :=
    Real.tendsto_log_atTop.comp hNtop
  have hNgt : ∀ᶠ i in atTop, 1 < N i := hNtop.eventually (eventually_gt_atTop 1)
  have hlogNpos : ∀ᶠ i in atTop, 0 < Real.log (N i) := by
    filter_upwards [hNgt] with i hi
    exact Real.log_pos hi
  have hpowtop : Tendsto (fun i ↦ N i ^ α⁻¹) atTop atTop :=
    (tendsto_rpow_atTop (inv_pos.mpr hα)).comp hNtop
  have hpowlarge : ∀ᶠ i in atTop, 1 / A < N i ^ α⁻¹ :=
    hpowtop.eventually (eventually_gt_atTop (1 / A))
  have hTgt : ∀ᶠ i in atTop, 1 < T i := by
    filter_upwards [hpowlarge, hT] with i hi hTi
    have : 1 < A * N i ^ α⁻¹ := by
      simpa [mul_comm] using (div_lt_iff₀ hA).1 hi
    exact this.trans_le hTi.1
  have hratio : Tendsto (fun i ↦ Real.log (T i) / Real.log (N i)) atTop (nhds α⁻¹) := by
    have hlowerLimit : Tendsto
        (fun i ↦ Real.log A / Real.log (N i) + α⁻¹) atTop (nhds α⁻¹) := by
      simpa using (tendsto_const_nhds.div_atTop hlogNtop).add
        (tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ α⁻¹) atTop (nhds α⁻¹))
    have hupperLimit : Tendsto
        (fun i ↦ Real.log B / Real.log (N i) + α⁻¹) atTop (nhds α⁻¹) := by
      simpa using (tendsto_const_nhds.div_atTop hlogNtop).add
        (tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ α⁻¹) atTop (nhds α⁻¹))
    apply tendsto_of_tendsto_of_tendsto_of_le_of_le' hlowerLimit hupperLimit
    · filter_upwards [hNgt, hlogNpos, hT] with i hNi hlogNi hTi
      have hbase : 0 < N i := zero_lt_one.trans hNi
      have hlower := Real.log_le_log (mul_pos hA (Real.rpow_pos_of_pos hbase _)) hTi.1
      rw [Real.log_mul hA.ne' (Real.rpow_pos_of_pos hbase _).ne',
          Real.log_rpow hbase] at hlower
      apply (le_div_iff₀ hlogNi).2
      rw [add_mul, div_mul_cancel₀ _ hlogNi.ne']
      simpa [mul_comm] using hlower
    · filter_upwards [hNgt, hlogNpos, hT] with i hNi hlogNi hTi
      have hbase : 0 < N i := zero_lt_one.trans hNi
      have hupper := Real.log_le_log (lt_of_lt_of_le
        (mul_pos hA (Real.rpow_pos_of_pos hbase _)) hTi.1) hTi.2
      rw [Real.log_mul hB.ne' (Real.rpow_pos_of_pos hbase _).ne',
          Real.log_rpow hbase] at hupper
      apply (div_le_iff₀ hlogNi).2
      rw [add_mul, div_mul_cancel₀ _ hlogNi.ne']
      simpa [mul_comm] using hupper
  have hinv := hratio.inv₀ (inv_ne_zero hα.ne')
  have hinv' : Tendsto (fun i ↦ Real.logb (T i) (N i)) atTop (nhds α⁻¹⁻¹) :=
    hinv.congr' <| by
      filter_upwards [hlogNpos, hTgt] with i hlogNi hTi
      dsimp [Real.logb]
      field_simp [hlogNi.ne', (Real.log_pos hTi).ne']
  simpa using hinv'

/-- Fixed positive multiplicative uncertainty in `T = N ^ α⁻¹` does not affect the
power asymptotic `N = T ^ (α + o(1))`. -/
theorem isPowerAsymptotic_of_between_const_rpow
    {N T : VariableObject ℝ} {α A B : ℝ}
    (hα : 0 < α) (hA : 0 < A) (hB : 0 < B)
    (hNtop : Tendsto N atTop atTop)
    (hT : ∀ᶠ i in atTop, A * N i ^ α⁻¹ ≤ T i ∧ T i ≤ B * N i ^ α⁻¹) :
    IsPowerAsymptotic N T α := by
  have hNpos : ∀ᶠ i in atTop, 0 < N i :=
    hNtop.eventually (eventually_gt_atTop 0)
  have hpowtop : Tendsto (fun i ↦ N i ^ α⁻¹) atTop atTop :=
    (tendsto_rpow_atTop (inv_pos.mpr hα)).comp hNtop
  have hpowlarge : ∀ᶠ i in atTop, 1 / A < N i ^ α⁻¹ :=
    hpowtop.eventually (eventually_gt_atTop (1 / A))
  have hTgt : ∀ᶠ i in atTop, 1 < T i := by
    filter_upwards [hpowlarge, hT] with i hi hTi
    have : 1 < A * N i ^ α⁻¹ := by
      simpa [mul_comm] using (div_lt_iff₀ hA).1 hi
    exact this.trans_le hTi.1
  exact isPowerAsymptotic_of_logb_tendsto hTgt
    hNpos (tendsto_logb_of_between_const_rpow hα hA hB hNtop hT)

/-! ## Algebra and limits -/

/-- A base is power-asymptotic to itself with exponent one. -/
theorem isPowerAsymptotic_self (T : VariableObject ℝ) :
    IsPowerAsymptotic T T 1 :=
  ⟨VariableObject.fixed 1, IsEqUpToInfinitesimal.refl _,
    Filter.Eventually.of_forall fun i ↦ (Real.rpow_one (T i)).symm⟩

/-- Products of powers add their asymptotic exponents. -/
theorem IsPowerAsymptotic.mul
    {A B T : VariableObject ℝ} {α β : ℝ}
    (hA : IsPowerAsymptotic A T α) (hB : IsPowerAsymptotic B T β)
    (hT : ∀ i, 1 ≤ T i) :
    IsPowerAsymptotic (A * B) T (α + β) := by
  obtain ⟨a, ha, hAeq⟩ := hA
  obtain ⟨b, hb, hBeq⟩ := hB
  refine ⟨a + b, ?_, ?_⟩
  · have hab := ha.add hb
    convert hab using 1
    · ext i
      rfl
  · filter_upwards [hAeq, hBeq] with i hiA hiB
    simp only [Pi.mul_apply, hiA, hiB]
    simpa only [Pi.add_apply] using
      (Real.rpow_add (zero_lt_one.trans_le (hT i)) (a i) (b i)).symm

/-- Quotients of powers subtract their asymptotic exponents. -/
theorem IsPowerAsymptotic.div
    {A B T : VariableObject ℝ} {α β : ℝ}
    (hA : IsPowerAsymptotic A T α) (hB : IsPowerAsymptotic B T β)
    (hT : ∀ i, 1 ≤ T i) :
    IsPowerAsymptotic (A / B) T (α - β) := by
  obtain ⟨a, ha, hAeq⟩ := hA
  obtain ⟨b, hb, hBeq⟩ := hB
  refine ⟨a - b, ?_, ?_⟩
  · have hab := ha.sub hb
    convert hab using 1
    · ext i
      rfl
  · filter_upwards [hAeq, hBeq] with i hiA hiB
    simp only [Pi.div_apply, hiA, hiB]
    simpa only [Pi.sub_apply] using
      (Real.rpow_sub (zero_lt_one.trans_le (hT i)) (a i) (b i)).symm

/-- A negative power-asymptotic exponent forces convergence to zero. -/
theorem IsPowerAsymptotic.tendsto_zero_of_neg
    {X T : VariableObject ℝ} {α : ℝ} (hX : IsPowerAsymptotic X T α)
    (hα : α < 0) (hT : ∀ i, 1 ≤ T i) (hTunbounded : T.IsUnbounded) :
    Tendsto X atTop (𝓝 0) := by
  have hTtop : Tendsto T atTop atTop :=
    (VariableObject.isUnbounded_iff_tendsto_atTop
      (fun i ↦ zero_le_one.trans (hT i))).mp hTunbounded
  have hbetween := hX.eventually_between (Filter.Eventually.of_forall hT)
    (show 0 < -α / 2 by linarith)
  have hpowzero : Tendsto (fun i ↦ T i ^ (α + -α / 2)) atTop (𝓝 0) := by
    have : 0 < -(α + -α / 2) := by linarith
    simpa [Function.comp_def, neg_neg] using
      (tendsto_rpow_neg_atTop this).comp hTtop
  apply squeeze_zero'
  · filter_upwards [hbetween] with i hi
    exact (Real.rpow_nonneg (zero_le_one.trans (hT i)) _).trans hi.1
  · exact hbetween.mono fun _ hi ↦ hi.2
  · exact hpowzero

/-- A positive power-asymptotic exponent forces convergence to infinity. -/
theorem IsPowerAsymptotic.tendsto_atTop_of_pos
    {X T : VariableObject ℝ} {α : ℝ} (hX : IsPowerAsymptotic X T α)
    (hα : 0 < α) (hT : ∀ i, 1 ≤ T i) (hTunbounded : T.IsUnbounded) :
    Tendsto X atTop atTop := by
  have hTtop : Tendsto T atTop atTop :=
    (VariableObject.isUnbounded_iff_tendsto_atTop
      (fun i ↦ zero_le_one.trans (hT i))).mp hTunbounded
  have hbetween := hX.eventually_between (Filter.Eventually.of_forall hT)
    (show 0 < α / 2 by linarith)
  apply tendsto_atTop_mono' atTop (hbetween.mono fun _ hi ↦ hi.1)
  exact (tendsto_rpow_atTop (show 0 < α - α / 2 by linarith)).comp hTtop

/-- The logarithmic exponent of a power asymptotic converges to its specified exponent. -/
theorem IsPowerAsymptotic.tendsto_logb {N T : VariableObject ℝ} {α : ℝ}
    (hNT : IsPowerAsymptotic N T α) (hT : ∀ᶠ i in atTop, 1 < T i) :
    Tendsto (fun i ↦ Real.logb (T i) (N i)) atTop (𝓝 α) := by
  obtain ⟨e, he, heq⟩ := hNT
  apply (isEqUpToInfinitesimal_fixed_iff_tendsto.mp he).congr'
  filter_upwards [heq, hT] with i hi hiT
  rw [hi, Real.logb_rpow (zero_lt_one.trans hiT) hiT.ne']

/-- A quantity tending to a positive constant has power-asymptotic exponent zero. -/
theorem isPowerAsymptotic_zero_of_tendsto_pos
    {X T : VariableObject ℝ} {c : ℝ} (hc : 0 < c) (hX : Tendsto X atTop (𝓝 c))
    (hT : ∀ i, 1 ≤ T i) (hTunbounded : T.IsUnbounded) :
    IsPowerAsymptotic X T 0 := by
  have hTtop : Tendsto T atTop atTop :=
    (VariableObject.isUnbounded_iff_tendsto_atTop
      (fun i ↦ zero_le_one.trans (hT i))).mp hTunbounded
  have hTgt : ∀ᶠ i in atTop, 1 < T i := hTtop.eventually (eventually_gt_atTop 1)
  have hXpos : ∀ᶠ i in atTop, 0 < X i := hX.eventually (eventually_gt_nhds hc)
  apply isPowerAsymptotic_of_logb_tendsto hTgt hXpos
  have hlogX : Tendsto (fun i ↦ Real.log (X i)) atTop (𝓝 (Real.log c)) :=
    (Real.continuousAt_log hc.ne').tendsto.comp hX
  have hlogTtop := Real.tendsto_log_atTop.comp hTtop
  simpa [Real.logb] using hlogX.div_atTop hlogTtop

/-- A quantity tending to one has power-asymptotic exponent zero relative to any unbounded
base that is at least one. -/
theorem isPowerAsymptotic_zero_of_tendsto_one
    {X T : VariableObject ℝ} (hX : Tendsto X atTop (𝓝 1))
    (hT : ∀ i, 1 ≤ T i) (hTunbounded : T.IsUnbounded) :
    IsPowerAsymptotic X T 0 :=
  isPowerAsymptotic_zero_of_tendsto_pos zero_lt_one hX hT hTunbounded

/-! ## Natural-valued power scales -/

/-- The positive natural scale obtained by taking the floor of `T ^ κ`. -/
def floorRpow (T : VariableObject ℝ) (κ : ℝ) : VariableObject ℕ :=
  fun i ↦ max 1 ⌊T i ^ κ⌋₊

/-- The positive natural power scale is nonzero. -/
theorem floorRpow_pos (T : VariableObject ℝ) (κ : ℝ) (i : ℕ) :
    0 < floorRpow T κ i := by
  simp [floorRpow]

/-- Flooring a divergent positive power introduces only a multiplicative `1 + o(1)` error. -/
theorem tendsto_floorRpow_cast_div_rpow_one
    {T : VariableObject ℝ} {κ : ℝ} (hκ : 0 < κ)
    (hT : Tendsto T atTop atTop) :
    Tendsto (fun i ↦ (floorRpow T κ i : ℝ) / T i ^ κ) atTop (𝓝 1) := by
  have hpow : Tendsto (fun i ↦ T i ^ κ) atTop atTop :=
    (tendsto_rpow_atTop hκ).comp hT
  have hlarge : ∀ᶠ i in atTop, 1 ≤ T i ^ κ := hpow.eventually (eventually_ge_atTop 1)
  have hfloor : ∀ᶠ i in atTop,
      floorRpow T κ i = ⌊T i ^ κ⌋₊ := by
    filter_upwards [hlarge] with i hi
    exact max_eq_right ((Nat.one_le_floor_iff _).mpr hi)
  apply (tendsto_nat_floor_div_atTop.comp hpow).congr'
  filter_upwards [hfloor] with i hi
  simp only [Function.comp_apply]
  rw [hi]

/-- The cast of `floorRpow T κ` is power-asymptotic to exponent `κ`. -/
theorem isPowerAsymptotic_floorRpow
    {T : VariableObject ℝ} {κ : ℝ} (hκ : 0 < κ)
    (hT : ∀ i, 1 ≤ T i) (hTunbounded : T.IsUnbounded) :
    IsPowerAsymptotic (fun i ↦ (floorRpow T κ i : ℝ)) T κ := by
  have hTtop : Tendsto T atTop atTop :=
    (VariableObject.isUnbounded_iff_tendsto_atTop
      (fun i ↦ zero_le_one.trans (hT i))).mp hTunbounded
  have hratio := tendsto_floorRpow_cast_div_rpow_one hκ hTtop
  have hTgt : ∀ᶠ i in atTop, 1 < T i := hTtop.eventually (eventually_gt_atTop 1)
  have hqpos : ∀ᶠ i in atTop, 0 < (floorRpow T κ i : ℝ) :=
    Filter.Eventually.of_forall fun i ↦ by exact_mod_cast floorRpow_pos T κ i
  apply isPowerAsymptotic_of_logb_tendsto hTgt hqpos
  have hlogratio : Tendsto (fun i ↦
      Real.log ((floorRpow T κ i : ℝ) / T i ^ κ)) atTop (𝓝 0) := by
    simpa [Function.comp_def] using
      (Real.continuousAt_log one_ne_zero).tendsto.comp hratio
  have hlogTtop := Real.tendsto_log_atTop.comp hTtop
  have hsmall := hlogratio.div_atTop hlogTtop
  have hlimit : Tendsto (fun i ↦ κ +
      Real.log ((floorRpow T κ i : ℝ) / T i ^ κ) / (Real.log ∘ T) i)
      atTop (𝓝 κ) := by
    simpa using (tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ κ) atTop (𝓝 κ)).add hsmall
  apply hlimit.congr'
  filter_upwards [hTgt] with i hi
  rw [Real.logb]
  have hTi : 0 < T i := zero_lt_one.trans hi
  have hqi : 0 < (floorRpow T κ i : ℝ) := by exact_mod_cast floorRpow_pos T κ i
  rw [Real.log_div hqi.ne' (Real.rpow_pos_of_pos hTi κ).ne', Real.log_rpow hTi]
  dsimp only [Function.comp_apply]
  field_simp [(Real.log_pos hi).ne']
  ring

/-! ## Power bounds -/

/-- A constant has power exponent zero. -/
theorem isPowerBounded_const (T : VariableObject ℝ) (K : ℝ) :
    IsPowerBounded (fun _ ↦ K) T 0 := by
  refine ⟨VariableObject.fixed 0, IsEqUpToInfinitesimal.refl _, ?_⟩
  simp only [VariableObject.fixed, Real.rpow_zero]
  exact Asymptotics.isBigO_const_const K one_ne_zero atTop

/-- Restrict a power bound along a sequence tending to infinity. -/
theorem IsPowerBounded.comp_tendsto {E : Type*} [SeminormedAddCommGroup E]
    {X : VariableObject E} {T : VariableObject ℝ} {β : ℝ}
    (h : IsPowerBounded X T β) {φ : ℕ → ℕ} (hφ : Tendsto φ atTop atTop) :
    IsPowerBounded (X ∘ φ) (T ∘ φ) β := by
  obtain ⟨e, he, hbound⟩ := h
  refine ⟨e ∘ φ, ?_, hbound.comp_tendsto hφ⟩
  rw [IsEqUpToInfinitesimal, VariableObject.IsInfinitesimal] at he ⊢
  exact he.comp hφ

/-- A power bound in a height of growth `κ` becomes a bound of exponent `κ * β`
in the original height. -/
theorem IsPowerBounded.change_base {E : Type*} [SeminormedAddCommGroup E]
    {X : VariableObject E} {S T : VariableObject ℝ} {β κ : ℝ}
    (hX : IsPowerBounded X S β) (hS : IsPowerAsymptotic S T κ)
    (hT : ∀ i, 1 ≤ T i) : IsPowerBounded X T (κ * β) := by
  obtain ⟨a, ha, hXa⟩ := hX
  obtain ⟨b, hb, hSb⟩ := hS
  have hatop := isEqUpToInfinitesimal_fixed_iff_tendsto.mp ha
  have hbtop := isEqUpToInfinitesimal_fixed_iff_tendsto.mp hb
  refine ⟨b * a, isEqUpToInfinitesimal_fixed_iff_tendsto.mpr (hbtop.mul hatop), ?_⟩
  apply hXa.congr' Filter.EventuallyEq.rfl
  filter_upwards [hSb] with i hiS
  rw [hiS, ← Real.rpow_mul (zero_le_one.trans (hT i))]
  rfl

/-- A power asymptotic also gives the corresponding upper bound. -/
theorem IsPowerAsymptotic.isPowerBounded {N T : VariableObject ℝ} {α : ℝ}
    (h : IsPowerAsymptotic N T α) : IsPowerBounded N T α := by
  obtain ⟨e, he, heq⟩ := h
  have heq' : N =ᶠ[atTop] (fun i ↦ T i ^ e i) := heq
  exact ⟨e, he, (Asymptotics.isBigO_refl _ atTop).congr' heq'.symm Filter.EventuallyEq.rfl⟩

/-- Taking norms preserves power bounds. -/
theorem IsPowerBounded.norm {E : Type*} [SeminormedAddCommGroup E]
    {X : VariableObject E} {T : VariableObject ℝ} {β : ℝ}
    (h : IsPowerBounded X T β) : IsPowerBounded (fun i ↦ ‖X i‖) T β := by
  obtain ⟨e, he, hbound⟩ := h
  exact ⟨e, he, hbound.norm_left⟩

/-- Multiplying real power bounds adds their exponents. -/
theorem IsPowerBounded.mul {X Y T : VariableObject ℝ} {β γ : ℝ}
    (hX : IsPowerBounded X T β) (hY : IsPowerBounded Y T γ)
    (hT : ∀ i, 1 ≤ T i) : IsPowerBounded (X * Y) T (β + γ) := by
  obtain ⟨e, he, hXe⟩ := hX
  obtain ⟨f, hf, hYf⟩ := hY
  refine ⟨e + f, he.add hf, ?_⟩
  apply (hXe.mul hYf).congr' Filter.EventuallyEq.rfl
  filter_upwards [] with i
  simpa using (Real.rpow_add (zero_lt_one.trans_le (hT i)) (e i) (f i)).symm

/-- Scaling by a constant preserves a power bound. -/
theorem IsPowerBounded.const_mul {X T : VariableObject ℝ} {β : ℝ} (hX : IsPowerBounded X T β)
    (K : ℝ) : IsPowerBounded (fun i ↦ K * X i) T β := by
  obtain ⟨e, he, hXe⟩ := hX
  exact ⟨e, he, hXe.const_mul_left K⟩

/-- Strictly increasing a fixed exponent makes a real power little-o along a base tending to
infinity. -/
theorem isLittleO_rpow_comp_tendsto_of_lt
    {T : VariableObject ℝ} {β γ : ℝ}
    (hT : Tendsto T atTop atTop) (hβγ : β < γ) :
    (fun i ↦ T i ^ β) =o[atTop] (fun i ↦ T i ^ γ) := by
  have hzero : ∀ᶠ i in atTop, T i ^ γ = 0 → T i ^ β = 0 := by
    filter_upwards [hT.eventually (eventually_gt_atTop 0)] with i hi
    exact fun hzero ↦ ((Real.rpow_pos_of_pos hi γ).ne' hzero).elim
  rw [Asymptotics.isLittleO_iff_tendsto' hzero]
  have hlimit := (tendsto_rpow_neg_atTop (sub_pos.mpr hβγ)).comp hT
  apply hlimit.congr'
  filter_upwards [hT.eventually (eventually_gt_atTop 0)] with i hi
  change T i ^ (-(γ - β)) = T i ^ β / T i ^ γ
  rw [← Real.rpow_sub hi]
  congr 1
  ring

/-- If `T` is at least one and unbounded, then `X ≪ T ^ (β + o(1))` iff
`X = O(T ^ (β + ε))` for every fixed `ε > 0`. -/
theorem isPowerBounded_iff_forall_pos
    {E : Type*} [SeminormedAddCommGroup E]
    (X : VariableObject E) (T : VariableObject ℝ) (β : ℝ)
    (hT : ∀ i, 1 ≤ T i) (hTunbounded : T.IsUnbounded) :
    IsPowerBounded X T β ↔
      ∀ ε : ℝ, 0 < ε →
        Asymptotics.IsBigO atTop X (fun i ↦ T i ^ (β + ε)) := by
  have hTtop : Tendsto T atTop atTop :=
    (VariableObject.isUnbounded_iff_tendsto_atTop
      fun i ↦ zero_le_one.trans (hT i)).1 hTunbounded
  have hTgt : ∀ᶠ i in atTop, 1 < T i :=
    (tendsto_atTop.1 hTtop 2).mono fun _ hi ↦ by linarith
  constructor
  · rintro ⟨exponent, hexponent, hX⟩ ε hε
    refine hX.trans (Asymptotics.IsBigO.of_bound 1 ?_)
    have hexponent_upper :=
      (isEqUpToInfinitesimal_iff_forall_pos
        exponent (VariableObject.fixed β)).1 hexponent ε hε
    filter_upwards [hexponent_upper] with i hi
    rw [one_mul, Real.norm_eq_abs,
      abs_of_nonneg (Real.rpow_nonneg (le_trans zero_le_one (hT i)) _),
      Real.norm_eq_abs,
      abs_of_nonneg (Real.rpow_nonneg (le_trans zero_le_one (hT i)) _)]
    apply Real.rpow_le_rpow_of_exponent_le (hT i)
    rw [Real.norm_eq_abs] at hi
    linarith [le_abs_self (exponent i - β)]
  · intro hX
    let exponent : VariableObject ℝ := fun i ↦
      if ‖X i‖ = 0 then β else
        max β (Real.log ‖X i‖ / Real.log (T i))
    refine ⟨exponent, ?_, Asymptotics.IsBigO.of_bound 1 ?_⟩
    · apply (isEqUpToInfinitesimal_iff_forall_pos
        exponent (VariableObject.fixed β)).2
      intro δ hδ
      have hδ3 : 0 < δ / 3 := by linarith
      obtain ⟨C, hCpos, hC⟩ := (hX (δ / 3) hδ3).exists_pos
      have hC' : ∀ᶠ i in atTop,
          ‖X i‖ ≤ C * T i ^ (β + δ / 3) := hC.bound.mono fun i hi ↦ by
        simpa only [Real.norm_eq_abs,
          abs_of_nonneg (Real.rpow_nonneg (zero_le_one.trans (hT i)) _)] using hi
      have hpowtop : Tendsto (fun i ↦ T i ^ (δ / 3)) atTop atTop :=
        (tendsto_rpow_atTop hδ3).comp hTtop
      have hCpower : ∀ᶠ i in atTop, C ≤ T i ^ (δ / 3) :=
        tendsto_atTop.1 hpowtop C
      filter_upwards [hTgt, hC', hCpower] with i hiT hiX hiC
      by_cases hXi : ‖X i‖ = 0
      · simp [exponent, hXi, hδ]
      · have hXpos : 0 < ‖X i‖ := lt_of_le_of_ne (norm_nonneg _) (Ne.symm hXi)
        have hiTpos : 0 < T i := lt_trans zero_lt_one hiT
        have hXpower : ‖X i‖ ≤ T i ^ (β + 2 * δ / 3) := by
          calc
            ‖X i‖ ≤ C * T i ^ (β + δ / 3) := hiX
            _ ≤ T i ^ (δ / 3) * T i ^ (β + δ / 3) := by gcongr
            _ = T i ^ (β + 2 * δ / 3) := by
              rw [← Real.rpow_add hiTpos]
              congr 1
              ring
        have hratio :
            Real.log ‖X i‖ / Real.log (T i) ≤ β + 2 * δ / 3 := by
          apply (div_le_iff₀ (Real.log_pos hiT)).2
          exact (Real.le_rpow_iff_log_le hXpos hiTpos).1 hXpower
        rw [Real.norm_eq_abs, abs_lt]
        constructor
        · simp only [exponent, hXi, ↓reduceIte, VariableObject.fixed]
          linarith [le_max_left β
            (Real.log ‖X i‖ / Real.log (T i))]
        · simp only [exponent, hXi, ↓reduceIte, VariableObject.fixed]
          have := max_lt (show β < β + δ by linarith)
            (lt_of_le_of_lt hratio (show β + 2 * δ / 3 < β + δ by linarith))
          linarith
    · filter_upwards [hTgt] with i hiT
      rw [one_mul, Real.norm_eq_abs,
        abs_of_nonneg (Real.rpow_nonneg (le_trans zero_le_one (hT i)) _)]
      by_cases hXi : ‖X i‖ = 0
      · rw [hXi]
        exact Real.rpow_nonneg (le_trans zero_le_one (hT i)) _
      · have hXpos : 0 < ‖X i‖ := lt_of_le_of_ne (norm_nonneg _) (Ne.symm hXi)
        have hiTpos : 0 < T i := lt_trans zero_lt_one hiT
        apply (Real.le_rpow_iff_log_le hXpos hiTpos).2
        apply (div_le_iff₀ (Real.log_pos hiT)).1
        simp only [exponent, hXi, ↓reduceIte]
        exact le_max_right β (Real.log ‖X i‖ / Real.log (T i))

/-- Increasing the exponent preserves a power bound when the base is eventually at least one. -/
theorem IsPowerBounded.mono
    {E : Type*} [SeminormedAddCommGroup E]
    {X : VariableObject E} {T : VariableObject ℝ} {β γ : ℝ}
    (hX : IsPowerBounded X T β) (hT : ∀ᶠ i in atTop, 1 ≤ T i)
    (hβγ : β ≤ γ) :
    IsPowerBounded X T γ := by
  obtain ⟨exponent, hexponent, hX⟩ := hX
  let shiftedExponent : VariableObject ℝ :=
    fun i ↦ exponent i + (γ - β)
  refine ⟨shiftedExponent, ?_, hX.trans (Asymptotics.IsBigO.of_bound 1 ?_)⟩
  · have hshift := hexponent.add
      (IsEqUpToInfinitesimal.refl (VariableObject.fixed (γ - β)))
    convert hshift using 1
    · ext i
      dsimp [shiftedExponent]
    · ext i
      dsimp [shiftedExponent]
      ring
  filter_upwards [hT] with i hiT
  have hT0 : 0 ≤ T i := zero_le_one.trans hiT
  rw [one_mul, Real.norm_eq_abs, abs_of_nonneg (Real.rpow_nonneg hT0 _),
    Real.norm_eq_abs, abs_of_nonneg (Real.rpow_nonneg hT0 _)]
  exact Real.rpow_le_rpow_of_exponent_le hiT (by
    dsimp [shiftedExponent]
    linarith)

/-- A power lower bound for an unbounded object forces its exponent to be no larger than any
admissible power-bound exponent. -/
theorem exponent_le_of_isPowerBounded_of_eventually_norm_ge_rpow
    {E : Type*} [SeminormedAddCommGroup E]
    {X : VariableObject E} {T : VariableObject ℝ} {β γ c : ℝ}
    (hT : ∀ i, 1 ≤ T i) (hTunbounded : T.IsUnbounded)
    (hbound : IsPowerBounded X T β) (hc : 0 < c)
    (hlower : ∀ᶠ i in atTop, c * T i ^ γ ≤ ‖X i‖) :
    γ ≤ β := by
  by_contra hγβ
  have hβγ : β < γ := lt_of_not_ge hγβ
  let ε := (γ - β) / 2
  have hε : 0 < ε := by dsimp [ε]; linarith
  have hO := (isPowerBounded_iff_forall_pos X T β hT hTunbounded).1 hbound ε hε
  have hTtop : Tendsto T atTop atTop :=
    (VariableObject.isUnbounded_iff_tendsto_atTop
      fun i ↦ zero_le_one.trans (hT i)).1 hTunbounded
  have hpowO : (fun i ↦ T i ^ γ) =O[atTop] X := by
    apply Asymptotics.IsBigO.of_bound (1 / c)
    filter_upwards [hlower] with i hi
    rw [Real.norm_eq_abs, abs_of_nonneg
      (Real.rpow_nonneg (zero_le_one.trans (hT i)) _)]
    calc
      T i ^ γ ≤ ‖X i‖ / c := (le_div_iff₀ hc).2 (by simpa [mul_comm] using hi)
      _ = 1 / c * ‖X i‖ := by ring
  have hsmall : (fun i ↦ T i ^ (β + ε)) =o[atTop] (fun i ↦ T i ^ γ) :=
    isLittleO_rpow_comp_tendsto_of_lt hTtop (by dsimp [ε]; linarith)
  have hnonzero : ∃ᶠ i in atTop, T i ^ (β + ε) ≠ 0 :=
    (Filter.Eventually.of_forall fun i ↦
      (Real.rpow_pos_of_pos (zero_lt_one.trans_le (hT i)) _).ne').frequently
  exact hsmall.not_isBigO hnonzero (hpowO.trans hO)

/-- Underspill for power bounds: bounds at every strictly larger exponent give a bound
at the limiting exponent. -/
theorem IsPowerBounded.of_forall_gt
    {E : Type*} [SeminormedAddCommGroup E]
    {X : VariableObject E} {T : VariableObject ℝ} {β : ℝ}
    (hT : ∀ i, 1 ≤ T i) (hTunbounded : T.IsUnbounded)
    (h : ∀ γ : ℝ, β < γ → IsPowerBounded X T γ) :
    IsPowerBounded X T β := by
  apply (isPowerBounded_iff_forall_pos X T β hT hTunbounded).2
  intro ε hε
  have hb := (isPowerBounded_iff_forall_pos X T (β + ε / 2) hT hTunbounded).1
    (h (β + ε / 2) (by linarith)) (ε / 2) (by linarith)
  simpa only [show β + ε / 2 + ε / 2 = β + ε by ring] using hb

/-- A sum of real power bounds has exponent the maximum of the two exponents. -/
theorem IsPowerBounded.add {X Y T : VariableObject ℝ} {β γ : ℝ}
    (hX : IsPowerBounded X T β) (hY : IsPowerBounded Y T γ)
    (hT : ∀ i, 1 ≤ T i) : IsPowerBounded (X + Y) T (max β γ) := by
  obtain ⟨e, he, hXe⟩ := hX
  obtain ⟨f, hf, hYf⟩ := hY
  let m : VariableObject ℝ := fun i ↦ max (e i) (f i)
  have hm : m =o VariableObject.fixed (max β γ) :=
    isEqUpToInfinitesimal_fixed_iff_tendsto.mpr
      ((isEqUpToInfinitesimal_fixed_iff_tendsto.mp he).max
        (isEqUpToInfinitesimal_fixed_iff_tendsto.mp hf))
  have hupper (g : VariableObject ℝ) (hg : ∀ i, g i ≤ m i) :
      Asymptotics.IsBigO atTop (fun i ↦ T i ^ g i) (fun i ↦ T i ^ m i) := by
    apply Asymptotics.IsBigO.of_bound 1
    filter_upwards [] with i
    simp only [one_mul, Real.norm_of_nonneg (Real.rpow_nonneg (zero_le_one.trans (hT i)) _)]
    exact Real.rpow_le_rpow_of_exponent_le (hT i) (hg i)
  exact ⟨m, hm, (hXe.trans (hupper e fun i ↦ le_max_left _ _)).add
    (hYf.trans (hupper f fun i ↦ le_max_right _ _))⟩

/-- Finite sums of power-bounded quantities are power bounded with the same exponent. -/
theorem isPowerBounded_finset_sum {ι : Type*} (s : Finset ι) {X : ι → VariableObject ℝ}
    {T : VariableObject ℝ} {β : ℝ} (hT : ∀ i, 1 ≤ T i)
    (hX : ∀ j ∈ s, IsPowerBounded (X j) T β) :
    IsPowerBounded (fun i ↦ ∑ j ∈ s, X j i) T β := by
  classical
  induction s using Finset.induction_on with
  | empty =>
    simp only [Finset.sum_empty]
    exact ⟨VariableObject.fixed β, IsEqUpToInfinitesimal.refl _,
      Asymptotics.isBigO_zero _ _⟩
  | insert j s hj ih =>
    simp only [Finset.sum_insert hj]
    have h := (hX j (Finset.mem_insert_self j s)).add
      (ih fun k hk ↦ hX k (Finset.mem_insert_of_mem hk)) hT
    rw [max_self] at h
    exact h
/-- The logarithm of a quantity between one and a fixed power of the height has
power exponent zero. -/
theorem isPowerBounded_log_of_le_const_rpow {X T : VariableObject ℝ} {C p : ℝ}
    (hC : 1 ≤ C) (hp : 0 ≤ p) (hT : ∀ i, 1 ≤ T i) (hTunb : T.IsUnbounded)
    (hX1 : ∀ᶠ i in atTop, 1 ≤ X i) (hXle : ∀ᶠ i in atTop, X i ≤ C * T i ^ p) :
    IsPowerBounded (fun i ↦ Real.log (X i)) T 0 := by
  apply (isPowerBounded_iff_forall_pos _ T 0 hT hTunb).2
  intro ε hε
  apply Asymptotics.IsBigO.of_bound (Real.log C + p / ε)
  filter_upwards [hX1, hXle] with i h1 h2
  have hTpos : 0 < T i := zero_lt_one.trans_le (hT i)
  have hCpos : 0 < C := zero_lt_one.trans_le hC
  have hlog : Real.log (X i) ≤ Real.log C + p * Real.log (T i) := by
    calc Real.log (X i) ≤ Real.log (C * T i ^ p) :=
          Real.log_le_log (zero_lt_one.trans_le h1) h2
      _ = Real.log C + p * Real.log (T i) := by
          rw [Real.log_mul hCpos.ne' (Real.rpow_pos_of_pos hTpos _).ne', Real.log_rpow hTpos]
  have hlogT := Real.log_le_rpow_div hTpos.le hε
  have hTε := Real.one_le_rpow (hT i) hε.le
  have hlogC := Real.log_nonneg hC
  rw [Real.norm_of_nonneg (Real.log_nonneg h1), zero_add,
    Real.norm_of_nonneg (Real.rpow_nonneg hTpos.le _)]
  calc Real.log (X i) ≤ Real.log C + p * Real.log (T i) := hlog
    _ ≤ Real.log C * T i ^ ε + p * (T i ^ ε / ε) :=
      add_le_add (by nlinarith) (mul_le_mul_of_nonneg_left hlogT hp)
    _ = _ := by ring

/-- An eventual pointwise norm bound transfers a power bound. -/
theorem IsPowerBounded.of_norm_le {E : Type*} [SeminormedAddCommGroup E]
    {X : VariableObject E} {Y T : VariableObject ℝ} {β : ℝ}
    (hY : IsPowerBounded Y T β) (hXY : ∀ᶠ i in atTop, ‖X i‖ ≤ Y i) :
    IsPowerBounded X T β := by
  obtain ⟨e, he, hYe⟩ := hY
  refine ⟨e, he, (Asymptotics.IsBigO.of_bound 1 ?_).trans hYe⟩
  exact hXY.mono fun i hi ↦ by simpa using hi.trans (le_abs_self (Y i))

/-- Taking the square root of a bound for a squared norm halves its exponent. -/
theorem IsPowerBounded.of_norm_sq {E : Type*} [SeminormedAddCommGroup E]
    {X : VariableObject E} {T : VariableObject ℝ} {β : ℝ}
    (hX : IsPowerBounded (fun i ↦ ‖X i‖ ^ 2) T (2 * β))
    (hT : ∀ i, 1 ≤ T i) (hTunbounded : T.IsUnbounded) : IsPowerBounded X T β := by
  apply (isPowerBounded_iff_forall_pos X T β hT hTunbounded).2
  intro ε hε
  have hO := (isPowerBounded_iff_forall_pos _ T (2 * β) hT hTunbounded).1 hX
    (2 * ε) (by linarith)
  have hpow (i : ℕ) : (T i ^ (β + ε)) ^ 2 = T i ^ (2 * β + 2 * ε) := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul (zero_le_one.trans (hT i))]
    congr 1
    norm_num
    ring
  have hsq : ((fun i ↦ ‖X i‖) ^ 2) =O[atTop] ((fun i ↦ T i ^ (β + ε)) ^ 2) := by
    change (fun i ↦ ‖X i‖ ^ 2) =O[atTop] (fun i ↦ (T i ^ (β + ε)) ^ 2)
    simpa only [hpow] using hO
  exact Asymptotics.isBigO_norm_left.mp (hsq.of_pow (by decide : (2 : ℕ) ≠ 0))

/-- To prove a power bound, it suffices to find it on a further subsequence of every
subsequence. Failure of a bound would give a subsequence with a larger power lower bound. -/
theorem isPowerBounded_of_subsequences {E : Type*} [SeminormedAddCommGroup E]
    {X : VariableObject E} {T : VariableObject ℝ} {β : ℝ}
    (hT : ∀ i, 1 ≤ T i) (hTunbounded : T.IsUnbounded)
    (hsub : ∀ φ : ℕ → ℕ, StrictMono φ →
      ∃ ψ : ℕ → ℕ, StrictMono ψ ∧
        IsPowerBounded (X ∘ φ ∘ ψ) (T ∘ φ ∘ ψ) β) :
    IsPowerBounded X T β := by
  apply (isPowerBounded_iff_forall_pos X T β hT hTunbounded).2
  intro ε hε
  by_contra hO
  have hbad : ∀ k, ∃ i ≥ k, T i ^ (β + ε) < ‖X i‖ := by
    intro k
    by_contra hk
    push Not at hk
    apply hO
    apply Asymptotics.IsBigO.of_bound 1
    filter_upwards [eventually_ge_atTop k] with i hi
    simpa only [one_mul, Real.norm_eq_abs,
      abs_of_nonneg (Real.rpow_nonneg (zero_le_one.trans (hT i)) _)] using hk i hi
  obtain ⟨φ, hφ, hφbad⟩ := extraction_of_frequently_atTop (frequently_atTop.2 hbad)
  obtain ⟨ψ, hψ, hbound⟩ := hsub φ hφ
  have hT' : ∀ i, 1 ≤ (T ∘ φ ∘ ψ) i := fun i ↦ hT _
  have hTtop := (VariableObject.isUnbounded_iff_tendsto_atTop
    (fun i ↦ zero_le_one.trans (hT i))).1 hTunbounded
  have hT'unbounded : VariableObject.IsUnbounded (T ∘ φ ∘ ψ) :=
    (VariableObject.isUnbounded_iff_tendsto_atTop
      (fun i ↦ zero_le_one.trans (hT' i))).2
      ((hTtop.comp hφ.tendsto_atTop).comp hψ.tendsto_atTop)
  have hle := exponent_le_of_isPowerBounded_of_eventually_norm_ge_rpow
    hT' hT'unbounded hbound zero_lt_one
    (Filter.Eventually.of_forall fun i ↦ by simpa using (hφbad (ψ i)).le)
  linarith

end Expdb
