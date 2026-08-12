module

public import Expdb.ExponentialSums.TrivialBounds

import Expdb.ExponentialSums.ExponentSumGrowthNonAsymptotic
import Expdb.ExponentialSums.ScaleTransfer
import Mathlib.Data.Finset.Max
import Mathlib.Topology.Semicontinuity.Basic

/-!
# Upper semicontinuity of the exponential-sum growth exponent

This file proves the upper semicontinuity assertion in the blueprint.  The main input is scale
transference: dilation and finite Fourier inversion increase the scale, while decomposition into
residue classes decreases it.  Both operations preserve the fixed model order of the phase.
-/

@[expose] public section

open Filter Topology
open scoped ContDiff Expdb FourierTransform NNReal

noncomputable section

namespace Expdb

/-! ## Scale transference -/

/-- On the subcritical interval, the exponential-sum growth exponent is monotone. -/
theorem exponentSumGrowthExponent_monoOn_Icc :
    MonotoneOn exponentSumGrowthExponent (Set.Icc 0 1) := by
  intro α hα γ hγ hαγ
  by_cases heq : α = γ
  · simp [heq]
  have hlt : (α : ℝ) < (γ : ℝ) := by exact_mod_cast lt_of_le_of_ne hαγ heq
  rw [exponentSumGrowthExponent_le_iff]
  intro N T F a b hN hT hTunbounded hNT hF hab
  apply (isPowerBounded_iff_forall_pos _ T _ hT hTunbounded).mpr
  intro ε hε
  obtain ⟨hphase, σ, hσ, herror⟩ := hF
  obtain ⟨δ, hδ, P, hP, C, hC, hfixed⟩ :=
    (exponentSumGrowthExponent_le_iff_nonAsymptotic.mp
      (le_rfl : exponentSumGrowthExponent γ ≤ exponentSumGrowthExponent γ))
      (ε / 2) (by linarith) σ hσ
  let η : ℝ := min (δ / 4) (((γ : ℝ) - (α : ℝ)) / 2)
  have hη : 0 < η := by dsimp [η]; positivity
  have hηδ : η < δ := lt_of_le_of_lt (min_le_left _ _) (by linarith)
  let κ : ℝ := (γ : ℝ) - (α : ℝ) - η
  have hκ : 0 < κ := by
    dsimp [κ, η]
    have := min_le_right (δ / 4) (((γ : ℝ) - (α : ℝ)) / 2)
    linarith
  let q : VariableObject ℕ := floorRpow T κ
  let Q : VariableObject ℝ := fun i ↦ (q i : ℝ)
  let M : VariableObject ℝ := Q * N
  have hTtop : Tendsto T atTop atTop :=
    (VariableObject.isUnbounded_iff_tendsto_atTop
      (fun i ↦ zero_le_one.trans (hT i))).mp hTunbounded
  have hQasym : IsPowerAsymptotic Q T κ :=
    isPowerAsymptotic_floorRpow hκ hT hTunbounded
  have hMasym0 : IsPowerAsymptotic M T (κ + (α : ℝ)) := hQasym.mul hNT hT
  have hMasym : IsPowerAsymptotic M T ((γ : ℝ) - η) := by
    convert hMasym0 using 1
    dsimp [κ]
    ring
  have hMdivAsym : IsPowerAsymptotic (M / T) T ((γ : ℝ) - η - 1) := by
    apply hMasym.div
    · refine ⟨VariableObject.fixed 1, IsEqUpToInfinitesimal.refl _, ?_⟩
      exact Filter.Eventually.of_forall fun i ↦ by simp [VariableObject.fixed]
    · exact hT
  have hMdivZero : Tendsto (M / T) atTop (𝓝 0) := by
    apply hMdivAsym.tendsto_zero_of_neg
    · have hγone : (γ : ℝ) ≤ 1 := by exact_mod_cast hγ.2
      linarith
    · exact hT
    · exact hTunbounded
  have hFapprox : ∀ᶠ i in atTop,
      IsApproximateModelPhaseFunction (F i) σ P (δ / 2) :=
    (IsModelPhaseFunctionWith.mk hphase herror).eventually_isApproximate P (by linarith)
  have hMdivSmall : ∀ᶠ i in atTop, M i / T i ≤ δ / 2 := by
    have hnear := (Metric.tendsto_nhds.1 hMdivZero) (δ / 2) (by linarith)
    filter_upwards [hnear] with i hi
    have hMnonneg : 0 ≤ M i := by
      change 0 ≤ (q i : ℝ) * N i
      exact mul_nonneg (Nat.cast_nonneg _) (zero_le_one.trans (hN i))
    rw [dist_zero_right] at hi
    change |M i / T i| < δ / 2 at hi
    rw [abs_of_nonneg (div_nonneg hMnonneg (zero_le_one.trans (hT i)))] at hi
    exact hi.le
  classical
  have hMbetween := hMasym.eventually_between
    (Filter.Eventually.of_forall hT) (sub_pos.mpr hηδ)
  have hthreshold : ∀ᶠ i in atTop, C ≤ T i :=
    hTtop.eventually (eventually_ge_atTop C)
  refine Asymptotics.IsBigO.of_bound C ?_
  filter_upwards [hFapprox, hMdivSmall, hMbetween, hthreshold]
    with i hiF hiMsmall hiM hiC
  have hqi : 0 < q i := floorRpow_pos T κ i
  letI : NeZero (q i) := ⟨hqi.ne'⟩
  have hMi : 1 ≤ M i := by
    dsimp [M, Q]
    have hqone : (1 : ℝ) ≤ q i := by exact_mod_cast (Nat.succ_le_iff.mp hqi)
    nlinarith [hN i]
  have hMbounds : T i ^ ((γ : ℝ) - δ) ≤ M i ∧ M i ≤ T i ^ ((γ : ℝ) + δ) := by
    constructor
    · convert hiM.1 using 1
      congr 1
      ring
    · exact hiM.2.trans (Real.rpow_le_rpow_of_exponent_le (hT i) (by linarith))
  have hinterval : M i ≤ ((q i * a i : ℕ) : ℝ) ∧
      ((q i * b i : ℕ) : ℝ) ≤ 2 * M i := by
    dsimp [M, Q]
    push_cast
    constructor <;> nlinarith [hqi.le, (hab i).1, (hab i).2]
  have hcoeffBound (k : ZMod (q i)) :
      ‖(k.val : ℝ) * N i / T i‖ ≤ M i / T i := by
    have hval : (k.val : ℝ) ≤ q i := by exact_mod_cast k.val_lt.le
    dsimp [M, Q]
    rw [abs_of_nonneg]
    · exact div_le_div_of_nonneg_right
        (mul_le_mul_of_nonneg_right hval (zero_le_one.trans (hN i)))
        (zero_le_one.trans (hT i))
    · exact div_nonneg
        (mul_nonneg (Nat.cast_nonneg _) (zero_le_one.trans (hN i)))
        (zero_le_one.trans (hT i))
  have hphaseBound (k : ZMod (q i)) :
      ‖exponentialSumAt
          (fun u ↦ F i u + (k.val : ℝ) * N i / T i * u)
          (T i) (M i) (q i * a i) (q i * b i)‖ ≤
        C * T i ^ (exponentSumGrowthExponent γ + ε / 2) := by
    have hiG : IsApproximateModelPhaseFunction
        (fun u ↦ F i u + (k.val : ℝ) * N i / T i * u) σ P δ := by
      have hk := hiF.add_linear ((hcoeffBound k).trans hiMsmall)
      convert hk using 1
      ring
    exact hfixed (T i) (M i) _ (q i * a i) (q i * b i)
      ⟨hiC, hMbounds.1, hMbounds.2, hiG, hinterval.1, hinterval.2⟩
  have hdilate := norm_exponentialSumAt_dilate_le (F i) (T i) (N i)
    (a i) (b i) (q i) hqi (ne_of_gt (zero_lt_one.trans_le (hT i)))
      (ne_of_gt (zero_lt_one.trans_le (hN i)))
  calc
    ‖exponentialSum F T N a b i‖ ≤
        ∑ _k : ZMod (q i),
          (C * T i ^ (exponentSumGrowthExponent γ + ε / 2)) / q i := by
      rw [exponentialSum_apply]
      refine hdilate.trans ?_
      apply Finset.sum_le_sum
      intro k _
      exact div_le_div_of_nonneg_right (hphaseBound k) (Nat.cast_nonneg _)
    _ = C * T i ^ (exponentSumGrowthExponent γ + ε / 2) := by
      rw [Finset.sum_const, show Finset.univ.card = q i from ZMod.card (q i),
        nsmul_eq_mul]
      field_simp [show (q i : ℝ) ≠ 0 by exact_mod_cast hqi.ne']
    _ ≤ C * ‖T i ^ (exponentSumGrowthExponent γ + ε)‖ := by
      rw [Real.norm_eq_abs,
        abs_of_nonneg (Real.rpow_nonneg (zero_le_one.trans (hT i)) _)]
      exact mul_le_mul_of_nonneg_left
        (Real.rpow_le_rpow_of_exponent_le (hT i) (by linarith))
        (zero_le_one.trans hC)

/-- Increasing the scale by `γ - α` costs at most that amount in the growth exponent. -/
theorem exponentSumGrowthExponent_le_add_sub
    {α γ : ℝ≥0} (hαγ : α ≤ γ) :
    exponentSumGrowthExponent γ ≤
      exponentSumGrowthExponent α + (γ : ℝ) - (α : ℝ) := by
  by_cases heq : α = γ
  · subst γ; simp
  by_cases hαzero : α = 0
  · subst α
    simpa using (exponentSumGrowthExponent_le_iff.mpr (isExponentSumBound_self γ))
  have hαpos : 0 < (α : ℝ) := by positivity
  have hgap : 0 < (γ : ℝ) - (α : ℝ) := by
    have hlt : (α : ℝ) < (γ : ℝ) := by
      exact_mod_cast lt_of_le_of_ne hαγ heq
    linarith
  rw [exponentSumGrowthExponent_le_iff]
  intro N T F a b hN hT hTunbounded hNT hF hab
  apply (isPowerBounded_iff_forall_pos _ T _ hT hTunbounded).mpr
  intro ε hε
  obtain ⟨hphase, σ, hσ, herror⟩ := hF
  obtain ⟨δ, hδ, P, hP, C, hC, hfixed⟩ :=
    (exponentSumGrowthExponent_le_iff_nonAsymptotic.mp
      (le_rfl : exponentSumGrowthExponent α ≤ exponentSumGrowthExponent α))
      (ε / 4) (by linarith) σ hσ
  let η : ℝ := min (δ / 4) (min (ε / 8) ((α : ℝ) / 4))
  have hη : 0 < η := by dsimp [η]; positivity
  have hηδ : η < δ := lt_of_le_of_lt (min_le_left _ _) (by linarith)
  have hηε : η ≤ ε / 8 :=
    (min_le_right (δ / 4) _).trans (min_le_left _ _)
  have hηα : η < (α : ℝ) := by
    have hηquarter : η ≤ (α : ℝ) / 4 :=
      (min_le_right (δ / 4) (min (ε / 8) ((α : ℝ) / 4))).trans
        (min_le_right (ε / 8) ((α : ℝ) / 4))
    linarith
  let κ : ℝ := (γ : ℝ) - (α : ℝ) + η
  have hκ : 0 < κ := by dsimp [κ]; linarith
  let q : VariableObject ℕ := floorRpow T κ
  let Q : VariableObject ℝ := fun i ↦ (q i : ℝ)
  have hTtop : Tendsto T atTop atTop :=
    (VariableObject.isUnbounded_iff_tendsto_atTop
      (fun i ↦ zero_le_one.trans (hT i))).mp hTunbounded
  have hQasym : IsPowerAsymptotic Q T κ :=
    isPowerAsymptotic_floorRpow hκ hT hTunbounded
  have hQdivNasym : IsPowerAsymptotic (Q / N) T (η - (α : ℝ)) := by
    have := hQasym.div hNT hT
    convert this using 1
    dsimp [κ]
    ring
  have hQdivNzero : Tendsto (Q / N) atTop (𝓝 0) := by
    apply hQdivNasym.tendsto_zero_of_neg
    · linarith
    · exact hT
    · exact hTunbounded
  have hqNevent : ∀ᶠ i in atTop, (q i : ℝ) ≤ N i := by
    have := (Metric.tendsto_nhds.1 hQdivNzero) 1 zero_lt_one
    filter_upwards [this] with i hi
    rw [dist_zero_right, Real.norm_eq_abs] at hi
    change |(q i : ℝ) / N i| < 1 at hi
    rw [
      abs_of_nonneg (div_nonneg (Nat.cast_nonneg _) (zero_le_one.trans (hN i)))] at hi
    have hNpos : 0 < N i := zero_lt_one.trans_le (hN i)
    exact (div_le_one hNpos).mp hi.le
  classical
  let A : ℕ → ℕ → ℕ := fun i r ↦ (a i - r) ⌈/⌉ q i
  let B : ℕ → ℕ → ℕ := fun i r ↦
    min ((b i - r) / q i) ⌊2 * residueScale (N i) (q i) r⌋₊
  have hmax (i : ℕ) : ∃ r ∈ Finset.range (q i), ∀ s ∈ Finset.range (q i),
      ‖exponentialSumAt (residuePhase (F i) (N i) s) (T i)
          (residueScale (N i) (q i) s) (A i s) (B i s)‖ ≤
      ‖exponentialSumAt (residuePhase (F i) (N i) r) (T i)
          (residueScale (N i) (q i) r) (A i r) (B i r)‖ :=
    Finset.exists_max_image (Finset.range (q i))
      (fun r ↦ ‖exponentialSumAt (residuePhase (F i) (N i) r) (T i)
        (residueScale (N i) (q i) r) (A i r) (B i r)‖)
      ⟨0, Finset.mem_range.mpr (floorRpow_pos T κ i)⟩
  let r : VariableObject ℕ := fun i ↦
    if (q i : ℝ) ≤ N i then (hmax i).choose else 0
  have hrq (i : ℕ) : r i < q i := by
    dsimp [r]
    split_ifs
    · exact Finset.mem_range.mp (hmax i).choose_spec.1
    · exact floorRpow_pos T κ i
  have hrN (i : ℕ) : (r i : ℝ) < N i := by
    dsimp [r]
    split_ifs with hi
    · exact lt_of_lt_of_le (by exact_mod_cast Finset.mem_range.mp (hmax i).choose_spec.1) hi
    · simpa only [Nat.cast_zero] using zero_lt_one.trans_le (hN i)
  let d : VariableObject ℝ := fun i ↦ (r i : ℝ) / N i
  let c : VariableObject ℝ := fun i ↦ 1 - d i
  have hdBound : ∀ i, 0 ≤ d i ∧ d i ≤ Q i / N i := by
    intro i
    constructor
    · dsimp [d]
      exact div_nonneg (Nat.cast_nonneg _) (zero_le_one.trans (hN i))
    · dsimp [d, Q]
      exact div_le_div_of_nonneg_right (by exact_mod_cast (hrq i).le)
        (zero_le_one.trans (hN i))
  have hdZero : Tendsto d atTop (𝓝 0) :=
    squeeze_zero' (Filter.Eventually.of_forall fun i ↦ (hdBound i).1)
      (Filter.Eventually.of_forall fun i ↦ (hdBound i).2) hQdivNzero
  have hcOne : Tendsto c atTop (𝓝 1) := by
    simpa [c] using tendsto_const_nhds.sub hdZero
  let M : VariableObject ℝ := fun i ↦ residueScale (N i) (q i) (r i)
  have hNasymDiv : IsPowerAsymptotic (N / Q) T ((γ : ℝ) - κ) := hNT.div hQasym hT
  have hcasym : IsPowerAsymptotic c T 0 :=
    isPowerAsymptotic_zero_of_tendsto_one hcOne hT hTunbounded
  have hMformula : M = (N / Q) * c := by
    funext i
    dsimp [M, residueScale, Q, c, d]
    field_simp [show (q i : ℝ) ≠ 0 by exact_mod_cast (floorRpow_pos T κ i).ne',
      show N i ≠ 0 by exact (ne_of_gt (zero_lt_one.trans_le (hN i)))]
  have hMasym0 : IsPowerAsymptotic M T (((γ : ℝ) - κ) + 0) := by
    rw [hMformula]
    exact hNasymDiv.mul hcasym hT
  have hMasym : IsPowerAsymptotic M T ((α : ℝ) - η) := by
    convert hMasym0 using 1
    dsimp [κ]
    ring
  let G : VariableFunction (VariableObject.fixed ℝ) ℝ := fun i ↦ residuePhase (F i) (N i) (r i)
  have hmap (i : ℕ) : Set.MapsTo (fun u ↦ c i * u + d i) phaseInterval phaseInterval := by
    intro u hu
    have hdI : d i ∈ Set.Ico (0 : ℝ) 1 := by
      exact ⟨(hdBound i).1,
        (div_lt_one (zero_lt_one.trans_le (hN i))).mpr (hrN i)⟩
    rcases hu with ⟨hu1, hu2⟩
    dsimp [c, phaseInterval] at *
    constructor <;> nlinarith [hdI.1, hdI.2]
  have hGdata : IsModelPhaseFunctionWith G σ := by
    have hdata :=
      (IsModelPhaseFunctionWith.mk hphase herror).comp_affine_tendsto_id hcOne hdZero hmap
    change IsModelPhaseFunctionWith
      (fun i u ↦ F i ((1 - (r i : ℝ) / N i) * u + (r i : ℝ) / N i)) σ
    simpa only [c, d] using hdata
  have hGapprox : ∀ᶠ i in atTop, IsApproximateModelPhaseFunction (G i) σ P δ :=
    hGdata.eventually_isApproximate P hδ
  have hMbetween := hMasym.eventually_between
    (Filter.Eventually.of_forall hT) (sub_pos.mpr hηδ)
  have hMtop : Tendsto M atTop atTop := by
    have hpowtop : Tendsto (fun i ↦ T i ^ (((α : ℝ) - η) / 2)) atTop atTop :=
      (tendsto_rpow_atTop (by linarith)).comp hTtop
    apply tendsto_atTop_mono' atTop _ hpowtop
    have hb := hMasym.eventually_between (Filter.Eventually.of_forall hT)
      (show 0 < ((α : ℝ) - η) / 2 by linarith)
    filter_upwards [hb] with i hi
    have hexp : ((α : ℝ) - η) / 2 =
        (α : ℝ) - η - ((α : ℝ) - η) / 2 := by ring
    rw [hexp]
    exact hi.1
  have hMone : ∀ᶠ i in atTop, 1 ≤ M i := hMtop.eventually (eventually_ge_atTop 1)
  have hthreshold : ∀ᶠ i in atTop, C ≤ T i := hTtop.eventually (eventually_ge_atTop C)
  have hqUpper : ∀ᶠ i in atTop, (q i : ℝ) ≤ T i ^ κ := by
    have hpowlarge : ∀ᶠ i in atTop, 1 ≤ T i ^ κ :=
      ((tendsto_rpow_atTop hκ).comp hTtop).eventually (eventually_ge_atTop 1)
    filter_upwards [hpowlarge] with i hi
    dsimp [q, floorRpow]
    rw [max_eq_right ((Nat.one_le_floor_iff _).mpr hi)]
    exact Nat.floor_le (Real.rpow_nonneg (zero_le_one.trans (hT i)) _)
  refine Asymptotics.IsBigO.of_bound (2 * C) ?_
  filter_upwards [hqNevent, hGapprox, hMbetween, hMone, hthreshold, hqUpper]
    with i hiqN hiG hiM hiMone hiC hiqUpper
  have hqi := floorRpow_pos T κ i
  have hMbounds : T i ^ ((α : ℝ) - δ) ≤ M i ∧ M i ≤ T i ^ ((α : ℝ) + δ) := by
    constructor
    · have hexp : (α : ℝ) - δ = (α : ℝ) - η - (δ - η) := by ring
      rw [hexp]
      exact hiM.1
    · exact hiM.2.trans (Real.rpow_le_rpow_of_exponent_le (hT i) (by linarith))
  have hA : M i ≤ (A i (r i) : ℝ) := by
    dsimp [M, residueScale, A]
    have hceil : a i - r i ≤ q i * ((a i - r i) ⌈/⌉ q i) :=
      le_smul_ceilDiv hqi
    have hra : r i ≤ a i := by
      exact_mod_cast (hrN i).le.trans (hab i).1
    have hceilR : ((a i - r i : ℕ) : ℝ) ≤
        (q i : ℝ) * (((a i - r i) ⌈/⌉ q i : ℕ) : ℝ) := by
      exact_mod_cast hceil
    rw [Nat.cast_sub hra] at hceilR
    rw [div_le_iff₀ (show (0 : ℝ) < q i by exact_mod_cast hqi)]
    nlinarith [(hab i).1]
  have hB : (B i (r i) : ℝ) ≤ 2 * M i := by
    dsimp [B, M]
    calc
      ((min ((b i - r i) / q i)
          ⌊2 * residueScale (N i) (q i) (r i)⌋₊ : ℕ) : ℝ) ≤
          (⌊2 * residueScale (N i) (q i) (r i)⌋₊ : ℕ) := by
        exact_mod_cast min_le_right ((b i - r i) / q i)
          ⌊2 * residueScale (N i) (q i) (r i)⌋₊
      _ ≤ 2 * residueScale (N i) (q i) (r i) :=
        Nat.floor_le (by linarith [hiMone])
  have hselected := hfixed (T i) (M i) (G i) (A i (r i)) (B i (r i))
    ⟨hiC, hMbounds.1, hMbounds.2, hiG, hA, hB⟩
  have hres := norm_exponentialSumAt_le_sum_residues (F i) (T i) (N i)
    (a i) (b i) (q i) hqi hiqN
      (ne_of_gt (zero_lt_one.trans_le (hN i))) (hab i).1 (hab i).2
  have hmaxEvent : ∀ s ∈ Finset.range (q i),
      ‖exponentialSumAt (residuePhase (F i) (N i) s) (T i)
          (residueScale (N i) (q i) s) (A i s) (B i s)‖ ≤
      ‖exponentialSumAt (G i) (T i) (M i) (A i (r i)) (B i (r i))‖ := by
    intro s hs
    have hri : r i = (hmax i).choose := by simp [r, hiqN]
    change ‖exponentialSumAt (residuePhase (F i) (N i) s) (T i)
        (residueScale (N i) (q i) s) (A i s) (B i s)‖ ≤
      ‖exponentialSumAt (residuePhase (F i) (N i) (r i)) (T i)
        (residueScale (N i) (q i) (r i)) (A i (r i)) (B i (r i))‖
    rw [hri]
    exact (hmax i).choose_spec.2 s hs
  have hsum : ‖exponentialSum F T N a b i‖ ≤
      (q i : ℝ) *
        (‖exponentialSumAt (G i) (T i) (M i) (A i (r i)) (B i (r i))‖ + 1) := by
    rw [exponentialSum_apply]
    refine hres.trans ?_
    calc
      _ ≤ ∑ _s ∈ Finset.range (q i),
          (‖exponentialSumAt (G i) (T i) (M i) (A i (r i)) (B i (r i))‖ + 1) := by
        gcongr with s hs
        simpa only [A, B] using hmaxEvent s hs
      _ = _ := by simp [mul_add, mul_comm]
  have hβnonneg : 0 ≤ exponentSumGrowthExponent α :=
    (isExponentSumBound_exponentSumGrowthExponent α).nonneg
  calc
    ‖exponentialSum F T N a b i‖ ≤
        (q i : ℝ) * (C * T i ^ (exponentSumGrowthExponent α + ε / 4) + 1) := by
      exact hsum.trans (mul_le_mul_of_nonneg_left
        (add_le_add hselected le_rfl) (Nat.cast_nonneg _))
    _ ≤ 2 * C * T i ^
        (exponentSumGrowthExponent α + (γ : ℝ) - (α : ℝ) + ε) := by
      have hone : 1 ≤ C * T i ^ (exponentSumGrowthExponent α + ε / 4) := by
        calc
          1 ≤ C := hC
          _ = C * 1 := by ring
          _ ≤ _ := by gcongr; exact Real.one_le_rpow (hT i) (by linarith)
      calc
        _ ≤ (q i : ℝ) * (2 * (C * T i ^ (exponentSumGrowthExponent α + ε / 4))) := by
          gcongr; linarith
        _ ≤ T i ^ κ * (2 * (C * T i ^ (exponentSumGrowthExponent α + ε / 4))) := by
          gcongr
        _ = 2 * C *
            (T i ^ κ * T i ^ (exponentSumGrowthExponent α + ε / 4)) := by ring
        _ = 2 * C * T i ^ (κ + exponentSumGrowthExponent α + ε / 4) := by
          congr 1
          rw [← Real.rpow_add (zero_lt_one.trans_le (hT i))]
          congr 1
          ring
        _ ≤ _ := by
          have hp : T i ^ (κ + exponentSumGrowthExponent α + ε / 4) ≤
              T i ^ (exponentSumGrowthExponent α + (γ : ℝ) - (α : ℝ) + ε) := by
            apply Real.rpow_le_rpow_of_exponent_le (hT i)
            dsimp [κ]
            linarith
          exact mul_le_mul_of_nonneg_left hp
            (mul_nonneg (by norm_num) (zero_le_one.trans hC))
    _ = 2 * C * ‖T i ^
        (exponentSumGrowthExponent α + (γ : ℝ) - (α : ℝ) + ε)‖ := by
      rw [Real.norm_eq_abs,
        abs_of_nonneg (Real.rpow_nonneg (zero_le_one.trans (hT i)) _)]

/-- The growth exponent is `1`-Lipschitz on the subcritical interval. -/
theorem lipschitzOnWith_exponentSumGrowthExponent :
    LipschitzOnWith 1 exponentSumGrowthExponent (Set.Icc 0 1) := by
  rw [lipschitzOnWith_iff_dist_le_mul]
  intro α hα γ hγ
  simp only [NNReal.coe_one, one_mul, Real.dist_eq]
  wlog hαγ : α ≤ γ generalizing α γ
  · rw [abs_sub_comm, dist_comm]
    exact this γ hγ α hα (le_of_not_ge hαγ)
  have hmono := exponentSumGrowthExponent_monoOn_Icc hα hγ hαγ
  have hslope := exponentSumGrowthExponent_le_add_sub hαγ
  rw [abs_of_nonpos (sub_nonpos.mpr hmono), NNReal.dist_eq,
    abs_of_nonpos (sub_nonpos.mpr (by exact_mod_cast hαγ))]
  linarith

/-- The exponential-sum growth exponent is continuous on the subcritical interval. -/
theorem continuousOn_exponentSumGrowthExponent :
    ContinuousOn exponentSumGrowthExponent (Set.Icc 0 1) :=
  lipschitzOnWith_exponentSumGrowthExponent.continuousOn

/-- The exponential-sum growth exponent is upper semicontinuous. -/
theorem upperSemicontinuous_exponentSumGrowthExponent :
    UpperSemicontinuous exponentSumGrowthExponent := by
  rw [upperSemicontinuous_iff]
  intro α y hy
  by_cases hαlt : α < 1
  · have hIcc : Set.Icc (0 : ℝ≥0) 1 ∈ nhds α :=
      Filter.mem_of_superset (Iio_mem_nhds hαlt) fun γ hγ ↦ ⟨zero_le, hγ.le⟩
    exact ((continuousOn_exponentSumGrowthExponent α ⟨zero_le, hαlt.le⟩).continuousAt
      hIcc).upperSemicontinuousAt y hy
  · have hαone : 1 ≤ α := le_of_not_gt hαlt
    by_cases hαeq : α = 1
    · subst α
      have hβnonneg := exponentSumGrowthExponent_nonneg 1
      let ρ : ℝ := min (y / 2) 1
      have hypos : 0 < y := lt_of_le_of_lt hβnonneg hy
      have hρ : 0 < ρ := by dsimp [ρ]; positivity
      rw [Metric.eventually_nhds_iff]
      refine ⟨ρ, hρ, ?_⟩
      intro γ hdist
      by_cases hγone : γ ≤ 1
      · exact lt_of_le_of_lt
          (exponentSumGrowthExponent_monoOn_Icc
            ⟨zero_le, hγone⟩ ⟨zero_le, le_rfl⟩ hγone) hy
      · have hγgt : 1 < γ := lt_of_not_ge hγone
        rw [exponentSumGrowthExponent_eq_sub_one hγgt]
        have hdist' : |(γ : ℝ) - 1| < ρ := by
          simpa [NNReal.dist_eq, abs_sub_comm] using hdist
        have hρy := min_le_left (y / 2) (1 : ℝ)
        have hdiffUpper := (abs_lt.mp hdist').2
        linarith
    · have hαgt : 1 < α := lt_of_le_of_ne hαone (Ne.symm hαeq)
      have heq : exponentSumGrowthExponent =ᶠ[nhds α]
          fun γ : ℝ≥0 ↦ (γ : ℝ) - 1 := by
        filter_upwards [Ioi_mem_nhds hαgt] with γ hγ
        exact exponentSumGrowthExponent_eq_sub_one hγ
      have hcont : ContinuousAt exponentSumGrowthExponent α :=
        (by fun_prop : ContinuousAt (fun γ : ℝ≥0 ↦ (γ : ℝ) - 1) α).congr_of_eventuallyEq
          heq
      exact hcont.upperSemicontinuousAt y hy

end Expdb
