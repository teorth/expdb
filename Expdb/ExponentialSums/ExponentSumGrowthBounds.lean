module

public import Expdb.ExponentialSums.UpperSemicontinuity
public import Expdb.ExponentialSums.AProcess
public import Expdb.ExponentialSums.BProcess
import Expdb.Mathlib.FourierCompletion
import Mathlib.Data.Nat.Basic

/-!
# Bounds for exponential sum growth exponents

This module records the classical transformations of the exponential sum growth exponent from
the blueprint's Exponential sum growth exponents chapter (`beta-chapter`).
-/

@[expose] public section

open Filter Topology
open scoped Expdb FourierTransform NNReal

noncomputable section

namespace Expdb

/-- Extend the exponential sum growth exponent to real inputs by truncating negative inputs to
zero. -/
def exponentSumGrowthExponentReal (α : ℝ) : ℝ :=
  exponentSumGrowthExponent α.toNNReal

/-- On nonnegative inputs, the real extension of `β` agrees with `β`. -/
@[simp] theorem exponentSumGrowthExponentReal_of_nonneg {α : ℝ} (hα : 0 ≤ α) :
    exponentSumGrowthExponentReal α = exponentSumGrowthExponent ⟨α, hα⟩ := by
  unfold exponentSumGrowthExponentReal
  apply congrArg exponentSumGrowthExponent
  apply NNReal.eq
  rw [Real.coe_toNNReal _ hα]
  rfl

/-! ## Classical transformations -/

/-- The van der Corput `A`-process for `β` (blueprint Lemma `vdca-beta`). -/
theorem vanDerCorputA_exponentSumGrowthExponent
    (α : ℝ≥0) (h : ℝ) (hα : α ≤ 2 / 3) (hh : 0 ≤ h) :
    2 * exponentSumGrowthExponent α ≤
      max (max (2 * (α : ℝ) - h) (2 * h))
        ((α : ℝ) - h + sSup
          ((fun h' : ℝ ↦
            (h' + 1 - α) *
                exponentSumGrowthExponentReal ((α : ℝ) / (h' + 1 - α)) + h') ''
            Set.Icc (2 * (α : ℝ) - 1) h)) := by
  let g : ℝ → ℝ := fun h' ↦ (h' + 1 - α) *
    exponentSumGrowthExponentReal ((α : ℝ) / (h' + 1 - α)) + h'
  let R : ℝ := max (max (2 * (α : ℝ) - h) (2 * h))
    ((α : ℝ) - h + sSup (g '' Set.Icc (2 * (α : ℝ) - 1) h))
  change 2 * exponentSumGrowthExponent α ≤ R
  have hself := exponentSumGrowthExponent_le_self α
  have hdiag : 2 * (α : ℝ) - h ≤ R :=
    (le_max_left _ _).trans (le_max_left _ _)
  have hhR : 2 * h ≤ R := (le_max_right _ _).trans (le_max_left _ _)
  by_cases hαh : (α : ℝ) ≤ h
  · linarith
  by_cases hhzero : h = 0
  · rw [hhzero] at hdiag
    linarith
  have hhpos : 0 < h := lt_of_le_of_ne hh (Ne.symm hhzero)
  have hhα : h < (α : ℝ) := lt_of_not_ge hαh
  have hαpos : 0 < (α : ℝ) := hhpos.trans hhα
  have hαR : (α : ℝ) ≤ 2 / 3 := by exact_mod_cast hα
  have hα1 : (α : ℝ) < 1 := by linarith
  have hbounded : BddAbove (g '' Set.Icc (2 * (α : ℝ) - 1) h) := by
    refine ⟨(α : ℝ) + h, ?_⟩
    rintro y ⟨h', hh', rfl⟩
    have hκ : 0 < h' + 1 - (α : ℝ) := by linarith [hh'.1]
    have hx : 0 ≤ (α : ℝ) / (h' + 1 - α) := (div_pos hαpos hκ).le
    have hβ := exponentSumGrowthExponent_le_self ((α : ℝ) / (h' + 1 - α)).toNNReal
    rw [Real.coe_toNNReal _ hx] at hβ
    have hproduct := mul_le_mul_of_nonneg_left hβ hκ.le
    rw [mul_div_cancel₀ _ hκ.ne'] at hproduct
    dsimp [g, exponentSumGrowthExponentReal]
    linarith [hh'.2]
  have hthird : (α : ℝ) - h + sSup (g '' Set.Icc (2 * (α : ℝ) - 1) h) ≤ R :=
    le_max_right _ _
  have hbound : IsExponentSumBound α (R / 2) := by
    apply isExponentSumBound_of_vanDerCorputA
      (by exact_mod_cast hαpos) (by exact_mod_cast hα1) hhpos hhα
    · linarith
    · intro h' hh'
      have hg : g h' ≤ sSup (g '' Set.Icc (2 * (α : ℝ) - 1) h) :=
        le_csSup hbounded ⟨h', hh', rfl⟩
      have heq : (1 + h' - α) *
          exponentSumGrowthExponent ((α : ℝ) / (1 + h' - α)).toNNReal + h' = g h' := by
        simp only [g, exponentSumGrowthExponentReal, add_comm 1 h']
      rw [heq]
      linarith
  have hβ := exponentSumGrowthExponent_le_iff.mpr hbound
  linarith

private theorem exponentSumGrowthExponent_le_half {α : ℝ≥0} (hαone : α ≤ 1) :
    exponentSumGrowthExponent α ≤ 1 / 2 := by
  have hinterior (α : ℝ≥0) (hα1 : α < 1) :
      exponentSumGrowthExponent α ≤ 1 / 2 := by
    by_cases hsmall : (α : ℝ) ≤ 1 / 2
    · exact (exponentSumGrowthExponent_le_self α).trans hsmall
    apply le_of_forall_pos_le_add
    intro ε hε
    let δ : ℝ := min ε ((2 * (α : ℝ) - 1) / 4)
    have hδ : 0 < δ := lt_min hε (by linarith)
    have hδε : δ ≤ ε := min_le_left _ _
    have hδsmall : δ ≤ (2 * (α : ℝ) - 1) / 4 := min_le_right _ _
    let h : ℝ := 2 * α - 1 - 2 * δ
    have hh : 0 < h := by dsimp [h]; linarith
    have hhα : h < (α : ℝ) := by
      dsimp [h]
      have : (α : ℝ) < 1 := by exact_mod_cast hα1
      linarith
    -- Below the small-shift cutoff the correlation hypothesis is empty.
    have hbound : IsExponentSumBound α (1 / 2 + δ) := by
      apply isExponentSumBound_of_vanDerCorputA
        (by change (0 : ℝ) < α; linarith) hα1 hh hhα
      · dsimp [h]; linarith
      · intro h' hh'
        dsimp [h] at hh'
        linarith [hh'.1, hh'.2]
    have hβ := exponentSumGrowthExponent_le_iff.mpr hbound
    linarith
  rcases lt_or_eq_of_le hαone with hlt | rfl
  · exact hinterior α hlt
  apply le_of_forall_pos_le_add
  intro ε hε
  let δ : ℝ := min ε (1 / 4)
  have hδ : 0 < δ := lt_min hε (by norm_num)
  have hδsmall : δ ≤ 1 / 4 := min_le_right _ _
  have hδε : δ ≤ ε := min_le_left _ _
  let x : ℝ≥0 := ⟨1 - δ, by linarith⟩
  have hxlt : x < 1 := by change 1 - δ < 1; linarith
  have hxhalf := hinterior x hxlt
  have htransfer := exponentSumGrowthExponent_le_add_sub (le_of_lt hxlt)
  change exponentSumGrowthExponent 1 ≤ exponentSumGrowthExponent x + 1 - (1 - δ)
    at htransfer
  linarith

/-- The first branch of the derivative bound holds throughout `0 ≤ α ≤ 1`. -/
theorem exponentSumGrowthExponent_le_vanDerCorput_of_le_one
    (α : ℝ≥0) (k : ℕ) (hk : 2 ≤ k) (hαone : α ≤ 1) :
    exponentSumGrowthExponent α ≤
      (α : ℝ) + (1 - k * (α : ℝ)) / ((2 : ℝ) ^ k - 2) := by
  induction k, hk using Nat.le_induction generalizing α with
  | base =>
    have h := exponentSumGrowthExponent_le_half hαone
    norm_num
    linarith
  | succ k hk ih =>
    have hkreal : (2 : ℝ) ≤ k := by exact_mod_cast hk
    have hq : (2 : ℝ) * k ≤ (2 : ℝ) ^ k := by
      exact_mod_cast Nat.mul_le_pow (by decide : 2 ≠ 1) k
    have hqpos : 0 < (2 : ℝ) ^ k - 2 := by linarith
    have hqone : 0 < (2 : ℝ) ^ k - 1 := by linarith
    have hαreal : (α : ℝ) ≤ 1 := by exact_mod_cast hαone
    rw [pow_succ]
    push_cast
    by_cases hsmall : ((k : ℝ) + 1) * α ≤ 1
    · exact (exponentSumGrowthExponent_le_self α).trans
        (le_add_of_nonneg_right (div_nonneg (by linarith) (by linarith)))
    by_cases heq : α = 1
    · subst α
      have hb := exponentSumGrowthExponent_le_half (le_refl (1 : ℝ≥0))
      have htarget : (1 / 2 : ℝ) ≤
          1 + (1 - ((k : ℝ) + 1) * 1) / ((2 : ℝ) ^ k * 2 - 2) := by
        have ht : (1 / 2 : ℝ) - 1 ≤
            (1 - ((k : ℝ) + 1) * 1) / ((2 : ℝ) ^ k * 2 - 2) := by
          apply (le_div_iff₀ (by linarith)).2
          nlinarith
        linarith
      exact hb.trans (by simpa using htarget)
    have hαlt : α < 1 := lt_of_le_of_ne hαone heq
    have hαpos : 0 < α := by
      change (0 : ℝ) < α
      nlinarith [α.coe_nonneg]
    -- Balance the diagonal bound with the correlation bound at the largest shift.
    let h : ℝ := (((k : ℝ) + 1) * α - 1) / ((2 : ℝ) ^ k - 1)
    have hh : 0 < h := div_pos (by linarith) hqone
    have hhα : h < (α : ℝ) := by
      apply (div_lt_iff₀ hqone).2
      have : (α : ℝ) < 1 := by exact_mod_cast hαlt
      nlinarith
    apply exponentSumGrowthExponent_le_iff.mpr
    apply isExponentSumBound_of_vanDerCorputA hαpos hαlt hh hhα
    · dsimp [h]
      field_simp
      nlinarith
    · intro s hs
      let w : ℝ := 1 + s - α
      have hwα : (α : ℝ) ≤ w := by dsimp [w]; linarith [hs.1]
      have hw : 0 < w := lt_of_lt_of_le (by exact_mod_cast hαpos) hwα
      have hx : 0 ≤ (α : ℝ) / w := div_nonneg α.coe_nonneg hw.le
      have hxone : ((α : ℝ) / w).toNNReal ≤ 1 := by
        rw [← NNReal.coe_le_coe, Real.coe_toNNReal _ hx, NNReal.coe_one]
        exact (div_le_one hw).2 hwα
      have hb := ih ((α : ℝ) / w).toNNReal hxone
      rw [Real.coe_toNNReal _ hx] at hb
      have hmul := add_le_add (mul_le_mul_of_nonneg_left hb hw.le) (le_refl s)
      change w * exponentSumGrowthExponent ((α : ℝ) / w).toNNReal + s ≤ _
      have halg : w * ((α : ℝ) / w + (1 - k * ((α : ℝ) / w)) / ((2 : ℝ) ^ k - 2)) + s =
          (α : ℝ) + (1 - ((k : ℝ) + 1) * α + s) / ((2 : ℝ) ^ k - 2) + s := by
        have hwne : 1 + s - (α : ℝ) ≠ 0 := hw.ne'
        dsimp [w]
        field_simp [hwne, hqpos.ne']
        ring
      rw [halg] at hmul
      have hmono : (α : ℝ) + (1 - ((k : ℝ) + 1) * α + s) /
          ((2 : ℝ) ^ k - 2) + s ≤
          (α : ℝ) + (1 - ((k : ℝ) + 1) * α + h) /
          ((2 : ℝ) ^ k - 2) + h := by
        apply add_le_add ?_ hs.2
        apply add_le_add le_rfl
        exact (div_le_div_iff_of_pos_right hqpos).2 (by linarith [hs.2])
      have hend : (α : ℝ) + (1 - ((k : ℝ) + 1) * α + h) /
          ((2 : ℝ) ^ k - 2) + h =
          2 * ((α : ℝ) + (1 - ((k : ℝ) + 1) * α) /
            ((2 : ℝ) ^ k * 2 - 2)) - α + h := by
        dsimp [h]
        field_simp
        ring
      exact hmul.trans (hmono.trans_eq hend)

/-- The van der Corput inequality for `β` (blueprint Proposition `beta-vdc`). -/
theorem exponentSumGrowthExponent_le_vanDerCorput
    (α : ℝ≥0) (k : ℕ) (hk : 2 ≤ k) (hα : 0 < α) :
    exponentSumGrowthExponent α ≤
      max ((α : ℝ) + (1 - k * (α : ℝ)) / ((2 : ℝ) ^ k - 2))
        ((1 - 2 ^ (2 - (k : ℤ))) * (α : ℝ) -
          (1 - k * (α : ℝ)) / ((2 : ℝ) ^ k - 2)) := by
  by_cases hαone : α ≤ 1
  · exact (exponentSumGrowthExponent_le_vanDerCorput_of_le_one α k hk hαone).trans
      (le_max_left _ _)
  have hαgt : 1 < α := lt_of_not_ge hαone
  rw [exponentSumGrowthExponent_eq_sub_one hαgt]
  apply le_max_of_le_right
  rcases eq_or_lt_of_le hk with rfl | hkgt
  · norm_num
    linarith
  rcases eq_or_lt_of_le (show 3 ≤ k by omega) with rfl | hkgt'
  · norm_num
    linarith
  have hkfour : (4 : ℝ) ≤ k := by exact_mod_cast (show 4 ≤ k by omega)
  have hq : (2 : ℝ) * k ≤ (2 : ℝ) ^ k := by
    exact_mod_cast Nat.mul_le_pow (by decide : 2 ≠ 1) k
  have hqpos : 0 < (2 : ℝ) ^ k := by positivity
  have hd : 0 < (2 : ℝ) ^ k - 2 := by linarith
  have hz : (2 : ℝ) ^ (2 - (k : ℤ)) = 4 / (2 : ℝ) ^ k := by
    rw [zpow_sub₀ (by norm_num)]
    norm_num
  rw [hz]
  have hr : 4 / (2 : ℝ) ^ k ≤ (k : ℝ) / ((2 : ℝ) ^ k - 2) := by
    apply (div_le_div_iff₀ hqpos hd).2
    nlinarith
  have hconst : 1 / ((2 : ℝ) ^ k - 2) ≤ 1 := by
    apply (div_le_one hd).2
    linarith
  have hm := mul_le_mul_of_nonneg_right hr (by exact_mod_cast hα.le : (0 : ℝ) ≤ α)
  have halg : (1 - 4 / (2 : ℝ) ^ k) * (α : ℝ) -
      (1 - k * (α : ℝ)) / ((2 : ℝ) ^ k - 2) =
      (α : ℝ) + ((k : ℝ) / ((2 : ℝ) ^ k - 2) - 4 / (2 : ℝ) ^ k) * α -
        1 / ((2 : ℝ) ^ k - 2) := by ring
  rw [halg]
  nlinarith

/-- The `k = 2` case displayed in blueprint Proposition `beta-vdc`. -/
theorem exponentSumGrowthExponent_le_vanDerCorput_two
    (α : ℝ≥0) (hα : 0 < α) :
    exponentSumGrowthExponent α ≤
      max (1 / 2) ((2 * (α : ℝ) - 1) / 2) := by
  have h := exponentSumGrowthExponent_le_vanDerCorput α 2 (by omega) hα
  norm_num only at h
  convert h using 1
  congr 1 <;> ring

/-- The `k = 3` case displayed in blueprint Proposition `beta-vdc`. -/
theorem exponentSumGrowthExponent_le_vanDerCorput_three
    (α : ℝ≥0) (hα : 0 < α) :
    exponentSumGrowthExponent α ≤
      max ((1 + 3 * (α : ℝ)) / 6) ((6 * (α : ℝ) - 1) / 6) := by
  have h := exponentSumGrowthExponent_le_vanDerCorput α 3 (by omega) hα
  norm_num only at h
  convert h using 1
  congr 1 <;> ring

/-- The `k = 4` case displayed in blueprint Proposition `beta-vdc`. -/
theorem exponentSumGrowthExponent_le_vanDerCorput_four
    (α : ℝ≥0) (hα : 0 < α) :
    exponentSumGrowthExponent α ≤
      max ((10 * (α : ℝ) + 1) / 14) ((29 * (α : ℝ) - 2) / 28) := by
  have h := exponentSumGrowthExponent_le_vanDerCorput α 4 (by omega) hα
  norm_num only at h
  convert h using 1
  congr 1 <;> ring

/-- The exponential sum growth exponent at one is `1 / 2`. -/
@[simp] theorem exponentSumGrowthExponent_one :
    exponentSumGrowthExponent 1 = 1 / 2 := by
  apply le_antisymm
  · have h := exponentSumGrowthExponent_le_vanDerCorput 1 2 (by norm_num) (by norm_num)
    norm_num at h ⊢
    exact h
  · exact (exponentSumGrowthExponent_mem_Icc (α := 1) (by norm_num)).1

/-- The optimized van der Corput inequality (blueprint Corollary `vdc-opt`). -/
theorem exponentSumGrowthExponent_le_vanDerCorput_sInf
    (α : ℝ≥0) (_hα : 0 < α) (hαone : α ≤ 1) :
    exponentSumGrowthExponent α ≤
      sInf {x : ℝ | ∃ k : ℕ, 2 ≤ k ∧
        x = (α : ℝ) + (1 - k * (α : ℝ)) / ((2 : ℝ) ^ k - 2)} := by
  apply le_csInf
  · exact ⟨_, 2, le_refl 2, rfl⟩
  · rintro x ⟨k, hk, rfl⟩
    exact exponentSumGrowthExponent_le_vanDerCorput_of_le_one α k hk hαone

/-- The three-term bound displayed in blueprint Corollary `vdc-opt`. -/
theorem exponentSumGrowthExponent_le_vanDerCorput_min
    (α : ℝ≥0) (_hα : 0 < α) (hαone : α ≤ 1) :
    exponentSumGrowthExponent α ≤
      min (min (1 / 2) ((1 + 3 * (α : ℝ)) / 6))
        ((10 * (α : ℝ) + 1) / 14) := by
  have htwo := exponentSumGrowthExponent_le_vanDerCorput_of_le_one α 2 (by omega) hαone
  have hthree := exponentSumGrowthExponent_le_vanDerCorput_of_le_one α 3 (by omega) hαone
  have hfour := exponentSumGrowthExponent_le_vanDerCorput_of_le_one α 4 (by omega) hαone
  norm_num only at htwo hthree hfour
  apply le_min (le_min ?_ ?_) ?_
  · convert htwo using 1
    · congr 1
    · ring
  · convert hthree using 1
    · congr 1
    · ring
  · convert hfour using 1
    · congr 1
    · ring

/-! ## Full dyadic intervals -/

private theorem completion_factor_bound {α : ℝ≥0} {N T : VariableObject ℝ}
    (hN : ∀ i, 1 ≤ N i) (hT : ∀ i, 1 ≤ T i) (hTunbounded : T.IsUnbounded)
    (hNT : IsPowerAsymptotic N T (α : ℝ)) :
    IsPowerBounded (fun i ↦ 3 + 2 * Real.log (⌊2 * N i⌋₊ + 1 : ℕ)) T 0 := by
  apply (isPowerBounded_iff_forall_pos _ T 0 hT hTunbounded).2
  intro ε hε
  let A : ℝ := 3 + 2 * Real.log 3
  have hA : 0 ≤ A := by dsimp [A]; positivity
  let C : ℝ := A + 2 * ((α : ℝ) + 1) / ε
  apply Asymptotics.IsBigO.of_bound C
  filter_upwards [hNT.eventually_between (Eventually.of_forall hT)
    (by norm_num : (0 : ℝ) < 1)] with i hi
  have hqpos : (0 : ℝ) < (⌊2 * N i⌋₊ + 1 : ℕ) := by positivity
  have hq : ((⌊2 * N i⌋₊ + 1 : ℕ) : ℝ) ≤ 3 * N i := by
    have hf := Nat.floor_le (by linarith [hN i] : 0 ≤ 2 * N i)
    push_cast
    linarith [hN i]
  have hscale : 3 * N i ≤ 3 * T i ^ ((α : ℝ) + 1) := by linarith [hi.2]
  have hlog : Real.log (⌊2 * N i⌋₊ + 1 : ℕ) ≤
      Real.log 3 + ((α : ℝ) + 1) * Real.log (T i) := by
    have htpos : 0 < T i := zero_lt_one.trans_le (hT i)
    have h := Real.log_le_log hqpos (hq.trans hscale)
    rw [Real.log_mul (by norm_num) (Real.rpow_pos_of_pos htpos _).ne', Real.log_rpow htpos] at h
    exact h
  have hpow : 1 ≤ T i ^ ε := Real.one_le_rpow (hT i) hε.le
  have hlogpow := Real.log_le_rpow_div (zero_le_one.trans (hT i)) hε
  have hL : 0 ≤ 3 + 2 * Real.log (⌊2 * N i⌋₊ + 1 : ℕ) := by
    positivity
  simp only [Real.norm_eq_abs, zero_add, abs_of_nonneg hL,
    abs_of_nonneg (Real.rpow_nonneg (zero_le_one.trans (hT i)) ε)]
  have hm := mul_le_mul_of_nonneg_left hlogpow (by positivity : 0 ≤ 2 * ((α : ℝ) + 1))
  have hApow := mul_le_mul_of_nonneg_left hpow hA
  dsimp [C, A] at hApow ⊢
  ring_nf at hm hApow hlog ⊢
  nlinarith

/-- Below scale one, finite Fourier completion turns full-interval bounds into
arbitrary-interval bounds. This also includes scale zero. -/
/- The blueprint uses smoothing and a Fourier integral. Here we use finite
Fourier completion instead: its logarithmic loss is T^{o(1)}, and its linear
twist preserves the model-phase conditions because N / T → 0 when α < 1. -/
private theorem isExponentSumBound_of_fullInterval_of_lt_one {α : ℝ≥0} {β : ℝ}
    (hbound : IsFullIntervalExponentSumBound α β) (hα : α < 1) :
    IsExponentSumBound α β := by
  intro N T F a b hN hT hTunbounded hNT hF hab
  let q (i : ℕ) : ℕ := ⌊2 * N i⌋₊ + 1
  have hq (i : ℕ) : 0 < q i := by dsimp [q]; omega
  have hcomplete (i : ℕ) := exists_norm_sum_Icc_le_log_mul_twist (hq i)
    (Nat.ceil_le.mpr (hab i).1) (Nat.le_floor (hab i).2)
    (by dsimp [q]; omega : ⌊2 * N i⌋₊ < q i)
    (fun n ↦ (𝐞 (T i * F i ((n : ℝ) / N i)) : ℂ))
  let r : VariableObject ℝ := fun i ↦ (hcomplete i).choose
  have hr (i : ℕ) : 0 ≤ r i ∧ r i ≤ 1 :=
    ⟨(hcomplete i).choose_spec.1, (hcomplete i).choose_spec.2.1⟩
  let c : VariableObject ℝ := fun i ↦ -(r i * (N i / T i))
  have hratio := (hNT.div (isPowerAsymptotic_self T) hT).tendsto_zero_of_neg
    (by
      have : (α : ℝ) < 1 := by exact_mod_cast hα
      linarith) hT hTunbounded
  have hc : c.IsInfinitesimal := by
    apply squeeze_zero' (g := fun i ↦ N i / T i)
    · exact Eventually.of_forall fun i ↦ norm_nonneg _
    · apply Eventually.of_forall
      intro i
      have hratio0 : 0 ≤ N i / T i :=
        div_nonneg (zero_le_one.trans (hN i)) (zero_le_one.trans (hT i))
      simp only [c, norm_neg, Real.norm_eq_abs, abs_mul, abs_of_nonneg (hr i).1,
        abs_of_nonneg hratio0]
      exact mul_le_of_le_one_left hratio0 (hr i).2
    · exact hratio
  obtain ⟨σ, hσ, hFσ⟩ := hF.exists_with
  let G : VariableFunction (VariableObject.fixed ℝ) ℝ := fun i u ↦ F i u + c i * u
  have hG : IsModelPhaseFunction G := (hFσ.add_infinitesimal_linear hc).isModelPhaseFunction hσ
  have hfull := hbound N T G hN hT hTunbounded hNT hG
  have hfactor := completion_factor_bound hN hT hTunbounded hNT
  have hproduct := hfactor.mul hfull.norm hT
  apply IsPowerBounded.of_norm_le (by simpa using hproduct)
  apply Eventually.of_forall
  intro i
  have he : ∑ n ∈ Finset.Icc ⌈N i⌉₊ ⌊2 * N i⌋₊,
      (𝐞 (T i * F i ((n : ℝ) / N i)) : ℂ) * (𝐞 (-(r i * n)) : ℂ) =
      exponentialSum G T N (fun i ↦ ⌈N i⌉₊) (fun i ↦ ⌊2 * N i⌋₊) i := by
    simp only [exponentialSum, exponentialSumAt, oscillatory, G, c]
    apply Finset.sum_congr rfl
    intro n _
    rw [← Circle.coe_mul, ← AddChar.map_add_eq_mul]
    congr 2
    have hNne : N i ≠ 0 := (zero_lt_one.trans_le (hN i)).ne'
    have hTne : T i ≠ 0 := (zero_lt_one.trans_le (hT i)).ne'
    field_simp [hNne, hTne]
  have hi := (hcomplete i).choose_spec.2.2
  rw [he] at hi
  simpa [exponentialSum, exponentialSumAt, oscillatory, q] using hi

/-- Arbitrary subintervals may be replaced by the full interval in the definition of `β`
(blueprint Lemma `interval-set`). -/
theorem isExponentSumBound_iff_fullInterval {α : ℝ≥0} {β : ℝ} :
    IsExponentSumBound α β ↔ IsFullIntervalExponentSumBound α β := by
  constructor
  · exact IsExponentSumBound.fullInterval
  intro hbound
  by_cases hαlt : α < 1
  · exact isExponentSumBound_of_fullInterval_of_lt_one hbound hαlt
  apply exponentSumGrowthExponent_le_iff.mp
  rcases eq_or_lt_of_le (le_of_not_gt hαlt) with rfl | hαgt
  · rw [exponentSumGrowthExponent_one]
    simpa using hbound.half_le (by norm_num) (le_refl (1 : ℝ≥0))
  · rw [exponentSumGrowthExponent_eq_sub_one hαgt]
    exact hbound.sub_one_le hαgt

/-- Characterization of `β` using only full dyadic intervals, combining blueprint Lemma
`interval-set` with the least-admissible-exponent characterization. -/
theorem exponentSumGrowthExponent_le_iff_fullInterval {α : ℝ≥0} {β : ℝ} :
    exponentSumGrowthExponent α ≤ β ↔ IsFullIntervalExponentSumBound α β := by
  rw [exponentSumGrowthExponent_le_iff, isExponentSumBound_iff_fullInterval]

/-! ## Reflection -/

/-- Reflection for the exponential sum growth exponent (blueprint Lemma `beta-reflect`). -/
theorem exponentSumGrowthExponent_reflection
    (α : ℝ≥0) (hα : 0 < α) (hα' : α < 1) :
    exponentSumGrowthExponent (1 - α) =
      1 / 2 - (α : ℝ) + exponentSumGrowthExponent α := by
  have hsub : ((1 - α : ℝ≥0) : ℝ) = 1 - α := by rw [NNReal.coe_sub hα'.le, NNReal.coe_one]
  apply le_antisymm
  · -- the B-process at scale `1 - α`, using the bound at scale `1 - (1 - α) = α`
    have h := isExponentSumBound_bProcess (α := 1 - α) (tsub_pos_of_lt hα')
      (tsub_lt_self one_pos hα) (B := exponentSumGrowthExponent α)
      (by rw [tsub_tsub_cancel_of_le hα'.le]
          exact isExponentSumBound_exponentSumGrowthExponent α)
    rw [hsub] at h
    refine (exponentSumGrowthExponent_le_iff.mpr h).trans (le_of_eq ?_)
    · ring
  · -- the B-process at scale `α`, using the bound at scale `1 - α`
    have h := isExponentSumBound_bProcess hα hα'
      (isExponentSumBound_exponentSumGrowthExponent (1 - α))
    have := exponentSumGrowthExponent_le_iff.mpr h
    linarith

end Expdb

end
