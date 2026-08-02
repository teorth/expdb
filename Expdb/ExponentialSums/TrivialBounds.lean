module

public import Expdb.ExponentialSums.ExponentSumGrowth

import Expdb.Basic.AutomaticUniformity
import Expdb.Fourier.L2Integral
import Expdb.Mathlib.EulerMaclaurin
import Mathlib.Analysis.Calculus.ContDiff.Bounds
import Mathlib.Analysis.Fourier.FourierTransformDeriv
import Mathlib.Analysis.SpecialFunctions.Pow.Asymptotics
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Base
import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.Analysis.Real.Pi.Bounds
import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts
import Mathlib.MeasureTheory.Integral.IntervalAverage

/-!
# Trivial bounds for exponential sum growth exponents

This module formalizes Lemma 4.4 of the ANTEDB blueprint.
-/

@[expose] public section

open Filter Topology
open scoped ContDiff Expdb FourierTransform NNReal

noncomputable section

namespace Expdb

private theorem contDiff_fourierChar : ContDiff ℝ ∞ (𝐞 · : ℝ → ℂ) := by
  rw [show (𝐞 · : ℝ → ℂ) = fun x : ℝ ↦
      Complex.exp (((2 * Real.pi * x : ℝ) : ℂ) * Complex.I) by
    funext x
    exact Real.fourierChar_apply x]
  have hreal : ContDiff ℝ ∞ (fun x : ℝ ↦ 2 * Real.pi * x) := by fun_prop
  exact ((Complex.ofRealCLM.contDiff.comp hreal).mul contDiff_const).cexp

private theorem iteratedDeriv_fourierChar (n : ℕ) (x : ℝ) :
    iteratedDeriv n (𝐞 · : ℝ → ℂ) x =
      ((2 * Real.pi : ℂ) * Complex.I) ^ n * 𝐞 x := by
  induction n with
  | zero => simp
  | succ n ih =>
      rw [iteratedDeriv_succ']
      have hderiv : deriv (𝐞 · : ℝ → ℂ) =
          fun y ↦ ((2 * Real.pi : ℂ) * Complex.I) * 𝐞 y := by
        funext y
        exact Real.deriv_fourierChar y
      rw [hderiv, iteratedDeriv_const_mul _
        (contDiff_fourierChar.contDiffAt.of_le
          (ENat.natCast_le_of_coe_top_le_withTop le_rfl n)), ih]
      ring

private lemma norm_iteratedDeriv_fourierChar_le
    {n k : ℕ} (hk : k ≤ n) (x : ℝ) :
    ‖iteratedDeriv k (𝐞 · : ℝ → ℂ) x‖ ≤ (2 * Real.pi + 1) ^ n := by
  rw [iteratedDeriv_fourierChar, norm_mul, norm_pow]
  have hcoeff : ‖(2 * Real.pi : ℂ) * Complex.I‖ = 2 * Real.pi := by
    simp [Complex.norm_real, Real.norm_eq_abs, Real.pi_pos.le]
  have hchar : ‖(𝐞 x : ℂ)‖ = 1 := by simp
  rw [hcoeff, hchar, mul_one]
  exact (pow_le_pow_left₀ (by positivity) (by linarith) k).trans
    (pow_le_pow_right₀ (by linarith [Real.pi_pos]) hk)

private theorem iteratedDerivWithin_comp_const_mul_of_mapsTo
    {f : ℝ → ℝ} {s t : Set ℝ} {x c : ℝ} {n : ℕ}
    (hf : ContDiffOn ℝ n f t) (hs : UniqueDiffOn ℝ s) (ht : UniqueDiffOn ℝ t)
    (hx : x ∈ s) (hst : Set.MapsTo (c * ·) s t) :
    iteratedDerivWithin n (fun y ↦ f (c * y)) s x =
      c ^ n * iteratedDerivWithin n f t (c * x) := by
  induction n generalizing x with
  | zero => simp
  | succ n ih =>
      have hcx : c * x ∈ t := hst hx
      have heq : Set.EqOn
          (iteratedDerivWithin n (fun y ↦ f (c * y)) s)
          (fun y ↦ c ^ n * iteratedDerivWithin n f t (c * y)) s :=
        fun y hy ↦ ih hf.of_succ hy
      have houter : DifferentiableWithinAt ℝ (iteratedDerivWithin n f t) t (c * x) :=
        hf.differentiableOn_iteratedDerivWithin (Nat.cast_lt.mpr n.lt_succ_self) ht _ hcx
      have hcomp : DifferentiableWithinAt ℝ
          (fun y ↦ iteratedDerivWithin n f t (c * y)) s x :=
        houter.comp x (by fun_prop) hst
      rw [iteratedDerivWithin_succ,
        derivWithin_congr heq (ih hf.of_succ hx),
        derivWithin_const_mul _ hcomp, iteratedDerivWithin_succ,
        ← Function.comp_def, derivWithin.scomp x houter (by fun_prop) hst,
        derivWithin_const_mul _ differentiableWithinAt_id,
        derivWithin_id' _ _ (hs.uniqueDiffWithinAt hx)]
      simp only [smul_eq_mul, mul_one, pow_succ]
      ring

private theorem norm_iteratedDerivWithin_oscillatory_le
    {F : ℝ → ℝ} {s : Set ℝ} {x T N K : ℝ} {n : ℕ}
    (hF : ContDiffOn ℝ n F phaseInterval) (hs : UniqueDiffOn ℝ s)
    (hx : x ∈ s) (hmap : Set.MapsTo (N⁻¹ * ·) s phaseInterval)
    (hT : 1 ≤ T) (hN : 1 ≤ N) (hK : 1 ≤ K)
    (hderiv : ∀ k : ℕ, 1 ≤ k → k ≤ n →
      ∀ u ∈ phaseInterval, ‖iteratedDerivWithin k F phaseInterval u‖ ≤ K) :
    ‖iteratedDerivWithin n (fun y ↦ (𝐞 (T * F (y / N)) : ℂ)) s x‖ ≤
      (n.factorial : ℝ) * (2 * Real.pi + 1) ^ n * (K * T / N) ^ n := by
  let q : ℝ → ℝ := fun y ↦ T * F (N⁻¹ * y)
  have hscale : ContDiffOn ℝ n (N⁻¹ * ·) s := by fun_prop
  have hFscale : ContDiffOn ℝ n (fun y ↦ F (N⁻¹ * y)) s :=
    by
      change ContDiffOn ℝ n (F ∘ fun y ↦ N⁻¹ * y) s
      exact hF.comp hscale hmap
  have hq : ContDiffOn ℝ n q s := by
    dsimp [q]
    exact contDiffOn_const.mul hFscale
  have houter : ContDiffOn ℝ n (𝐞 · : ℝ → ℂ) Set.univ :=
    (contDiff_fourierChar.of_le
      (ENat.natCast_le_of_coe_top_le_withTop le_rfl n)).contDiffOn
  have hC : ∀ k, k ≤ n →
      ‖iteratedFDerivWithin ℝ k (𝐞 · : ℝ → ℂ) Set.univ (q x)‖ ≤
        (2 * Real.pi + 1) ^ n := by
    intro k hk
    rw [norm_iteratedFDerivWithin_eq_norm_iteratedDerivWithin,
      iteratedDerivWithin_eq_iteratedDeriv uniqueDiffOn_univ
        (contDiff_fourierChar.contDiffAt.of_le
          (ENat.natCast_le_of_coe_top_le_withTop le_rfl k)) (Set.mem_univ _)]
    exact norm_iteratedDeriv_fourierChar_le hk _
  have hD : ∀ k, 1 ≤ k → k ≤ n →
      ‖iteratedFDerivWithin ℝ k q s x‖ ≤ (K * T / N) ^ k := by
    intro k hk hkn
    rw [norm_iteratedFDerivWithin_eq_norm_iteratedDerivWithin]
    have hrescale := iteratedDerivWithin_comp_const_mul_of_mapsTo
      (f := F) (c := N⁻¹) (n := k)
      (hF.of_le (Nat.cast_le.mpr hkn)) hs
      (uniqueDiffOn_Icc (by norm_num [phaseInterval])) hx hmap
    have hqderiv : iteratedDerivWithin k q s x =
        T * N⁻¹ ^ k * iteratedDerivWithin k F phaseInterval (N⁻¹ * x) := by
      rw [show q = fun y ↦ T * F (N⁻¹ * y) by rfl,
        iteratedDerivWithin_const_mul_field, hrescale]
      ring
    rw [hqderiv, norm_mul, norm_mul, Real.norm_eq_abs,
      abs_of_nonneg (le_trans zero_le_one hT), Real.norm_eq_abs,
      abs_of_nonneg (pow_nonneg (inv_nonneg.mpr (le_trans zero_le_one hN)) _)]
    have hFx : ‖iteratedDerivWithin k F phaseInterval (N⁻¹ * x)‖ ≤ K :=
      hderiv k hk hkn _ (hmap hx)
    calc
      T * N⁻¹ ^ k * ‖iteratedDerivWithin k F phaseInterval (N⁻¹ * x)‖ ≤
          T * N⁻¹ ^ k * K := by gcongr
      _ = (K * T) * N⁻¹ ^ k := by ring
      _ ≤ (K * T) ^ k * N⁻¹ ^ k := by
        gcongr
        simpa using pow_le_pow_right₀ (one_le_mul_of_one_le_of_one_le hK hT) hk
      _ = (K * T / N) ^ k := by simp only [div_eq_mul_inv, mul_pow]
  have hcomp := norm_iteratedFDerivWithin_comp_le
    houter hq (by rfl) uniqueDiffOn_univ hs (Set.mapsTo_univ q s) hx hC hD
  rw [norm_iteratedFDerivWithin_eq_norm_iteratedDerivWithin] at hcomp
  change ‖iteratedDerivWithin n (fun y ↦ (𝐞 (q y) : ℂ)) s x‖ ≤
    (n.factorial : ℝ) * (2 * Real.pi + 1) ^ n * (K * T / N) ^ n at hcomp
  simpa only [q, div_eq_mul_inv, mul_comm N⁻¹] using hcomp

private lemma map_Icc_natCast_int (a b : ℕ) :
    (Finset.Icc a b).map
      ⟨(fun n : ℕ ↦ (n : ℤ)), fun _ _ h ↦ Int.ofNat_inj.mp h⟩ =
      Finset.Icc (a : ℤ) (b : ℤ) := by
  ext z
  simp only [Finset.mem_map, Finset.mem_Icc]
  constructor
  · rintro ⟨n, hn, rfl⟩
    change (a : ℤ) ≤ (n : ℤ) ∧ (n : ℤ) ≤ (b : ℤ)
    exact_mod_cast hn
  · intro hz
    have hz0 : 0 ≤ z := le_trans (Int.natCast_nonneg a) hz.1
    have hzt : (z.toNat : ℤ) = z := Int.toNat_of_nonneg hz0
    refine ⟨z.toNat, ?_, hzt⟩
    constructor
    · exact_mod_cast (show (a : ℤ) ≤ (z.toNat : ℤ) by simpa [hzt] using hz.1)
    · exact_mod_cast (show (z.toNat : ℤ) ≤ (b : ℤ) by simpa [hzt] using hz.2)

private lemma sum_Icc_int_eq_sum_Icc_nat
    {E : Type*} [AddCommMonoid E] (f : ℤ → E) (a b : ℕ) :
    ∑ z ∈ Finset.Icc (a : ℤ) (b : ℤ), f z =
      ∑ n ∈ Finset.Icc a b, f n := by
  rw [← map_Icc_natCast_int]
  simp

private theorem oscillatory_sum_eq_eulerMaclaurin
    {F : ℝ → ℝ} {T N : ℝ} {a b s : ℕ} (hab : a ≤ b)
    (hF : ContDiffOn ℝ (s + 1) F phaseInterval)
    (hmap : Set.MapsTo (N⁻¹ * ·) (Set.Icc (a : ℝ) (b : ℝ)) phaseInterval) :
    ∑ n ∈ Finset.Icc a b, (𝐞 (T * F ((n : ℝ) / N)) : ℂ) =
      (∫ x in (a : ℝ)..b, (𝐞 (T * F (x / N)) : ℂ)) +
      (2⁻¹ : ℝ) •
        ((𝐞 (T * F ((a : ℝ) / N)) : ℂ) + (𝐞 (T * F ((b : ℝ) / N)) : ℂ)) +
      ∑ m ∈ Finset.range s, (-1 : ℝ) ^ m • EulerMaclaurin.saw (m + 2) 0 •
        (iteratedDerivWithin (m + 1) (fun x ↦ (𝐞 (T * F (x / N)) : ℂ))
            (Set.Icc (a : ℝ) (b : ℝ)) b -
          iteratedDerivWithin (m + 1) (fun x ↦ (𝐞 (T * F (x / N)) : ℂ))
            (Set.Icc (a : ℝ) (b : ℝ)) a) +
      (-1 : ℝ) ^ s •
        ∫ x in (a : ℝ)..b, EulerMaclaurin.saw (s + 1) x •
          iteratedDerivWithin (s + 1) (fun y ↦ (𝐞 (T * F (y / N)) : ℂ))
            (Set.Icc (a : ℝ) (b : ℝ)) x := by
  rcases hab.eq_or_lt with rfl | hab
  · simp only [Finset.Icc_self, Finset.sum_singleton, intervalIntegral.integral_same,
      sub_self, smul_zero, Finset.sum_const_zero]
    module
  let f : ℝ → ℂ := fun x ↦ 𝐞 (T * F (x / N))
  let interval : Set ℝ := Set.Icc (a : ℝ) (b : ℝ)
  have hscale : ContDiffOn ℝ (s + 1) (N⁻¹ * ·) interval := by fun_prop
  have hFscale : ContDiffOn ℝ (s + 1) (fun x ↦ F (N⁻¹ * x)) interval := by
    change ContDiffOn ℝ (s + 1) (F ∘ fun x ↦ N⁻¹ * x) interval
    exact hF.comp hscale hmap
  have hinner : ContDiffOn ℝ (s + 1) (fun x ↦ T * F (N⁻¹ * x)) interval :=
    contDiffOn_const.mul hFscale
  have hf : ContDiffOn ℝ (s + 1) f interval := by
    have heq : f = (𝐞 · : ℝ → ℂ) ∘ fun x ↦ T * F (N⁻¹ * x) := by
      funext x
      simp only [f, Function.comp_apply, div_eq_mul_inv]
      rw [mul_comm x N⁻¹]
    rw [heq]
    simpa only [Nat.cast_add, Nat.cast_one] using
      (contDiff_fourierChar.of_le
        (ENat.natCast_le_of_coe_top_le_withTop le_rfl (s + 1))).comp_contDiffOn hinner
  have hendNat : a + (b - a) = b := Nat.add_sub_of_le hab.le
  have hendInt : (a : ℤ) + (b - a : ℕ) = (b : ℤ) := by exact_mod_cast hendNat
  have hendReal : ((a : ℤ) : ℝ) + (b - a : ℕ) = ((b : ℤ) : ℝ) := by
    exact_mod_cast hendNat
  have hEM := EulerMaclaurin.sum_Icc_eq_integral_add
    (s := s) (a := (a : ℤ)) (n := b - a) (f := f) (t := interval)
    hf (uniqueDiffOn_Icc (by exact_mod_cast hab)) (by
      rw [hendReal]
      exact Set.Subset.rfl)
  rw [hendInt, hendReal, sum_Icc_int_eq_sum_Icc_nat] at hEM
  simpa only [f, interval, Int.cast_natCast, div_eq_mul_inv, mul_comm N⁻¹] using hEM

private theorem eventually_model_phase_deriv_bounds
    {F : VariableFunction (VariableObject.fixed ℝ) ℝ}
    (hF : IsModelPhaseFunction F) (P : ℕ) :
    ∃ c K : ℝ, 0 < c ∧ 1 ≤ K ∧
      ∀ᶠ i in atTop, ∀ u : phaseInterval,
        c ≤ iteratedDerivWithin 1 (F i) phaseInterval u ∧
          ∀ k : ℕ, 1 ≤ k → k ≤ P + 1 →
            ‖iteratedDerivWithin k (F i) phaseInterval u‖ ≤ K := by
  rcases hF with ⟨hphase, σ, hσ, herror⟩
  have href : ContDiffOn ℝ ∞ (fun u : ℝ ↦ u ^ (-σ)) phaseInterval := by
    intro u hu
    exact (Real.contDiffAt_rpow_const_of_ne (by
      exact ne_of_gt (lt_of_lt_of_le zero_lt_one hu.1))).contDiffWithinAt
  have hreference (p : ℕ) :
      ∃ B : ℝ, 0 ≤ B ∧ ∀ u : phaseInterval,
        ‖iteratedDerivWithin p (fun v : ℝ ↦ v ^ (-σ)) phaseInterval u‖ ≤ B := by
    have hcont : ContinuousOn
        (iteratedDerivWithin p (fun v : ℝ ↦ v ^ (-σ)) phaseInterval) phaseInterval :=
      href.continuousOn_iteratedDerivWithin
        (ENat.natCast_le_of_coe_top_le_withTop le_rfl p) (uniqueDiffOn_Icc (by
        norm_num [phaseInterval]))
    obtain ⟨B, hB⟩ := bddAbove_def.mp
      (isCompact_Icc.bddAbove_image hcont.norm)
    refine ⟨max B 0, le_max_right _ _, ?_⟩
    intro u
    exact (hB _ ⟨u, u.property, rfl⟩).trans (le_max_left _ _)
  choose B hBnonneg hB using hreference
  let K : ℝ := 1 + ∑ p ∈ Finset.range (P + 1), B p
  have hK : 1 ≤ K := by
    dsimp [K]
    exact le_add_of_nonneg_right (Finset.sum_nonneg fun p _ ↦ hBnonneg p)
  let c : ℝ := (2 : ℝ) ^ (-σ) / 2
  have hc : 0 < c := by dsimp [c]; positivity
  have herr (p : ℕ) (hp : p ≤ P) :
      ∀ᶠ i in atTop, ∀ u : phaseInterval,
        ‖modelPhaseError F σ p i u‖ < c :=
    (VariableFunction.isPointwiseInfinitesimal_iff_forall_pos_uniform
      (VariableObject.fixed phaseInterval)
      (fun _ ↦ ⟨1, by simp [phaseInterval]⟩)
      (modelPhaseError F σ p)).1 (herror p) c hc
  have hall : ∀ᶠ i in atTop, ∀ p ∈ Finset.range (P + 1),
      ∀ u : phaseInterval, ‖modelPhaseError F σ p i u‖ < c :=
    (Finset.range (P + 1)).eventually_all.mpr fun p hp ↦
      herr p (Nat.le_of_lt_succ (Finset.mem_range.mp hp))
  refine ⟨c, K, hc, hK, hall.mono ?_⟩
  intro i hi u
  have huPos : 0 < (u : ℝ) := zero_lt_one.trans_le u.property.1
  have huLower : (2 : ℝ) ^ (-σ) ≤ (u : ℝ) ^ (-σ) :=
    Real.rpow_le_rpow_of_nonpos huPos u.property.2 (neg_nonpos.mpr hσ.le)
  have herror0 := hi 0 (Finset.mem_range.mpr (Nat.zero_lt_succ P)) u
  have hfirst : c ≤ iteratedDerivWithin 1 (F i) phaseInterval u := by
    rw [Real.norm_eq_abs, abs_lt] at herror0
    simp only [modelPhaseError, zero_add, iteratedDerivWithin_zero] at herror0
    dsimp [c] at herror0 ⊢
    linarith
  refine ⟨hfirst, ?_⟩
  intro k hk hkP
  obtain ⟨p, rfl⟩ := Nat.exists_eq_add_of_le hk
  have hpP : p ≤ P := by omega
  have hpMem : p ∈ Finset.range (P + 1) := Finset.mem_range.mpr (by omega)
  have he := hi p hpMem u
  have htriangle :
      ‖iteratedDerivWithin (p + 1) (F i) phaseInterval u‖ ≤
        ‖modelPhaseError F σ p i u‖ +
          ‖iteratedDerivWithin p (fun v : ℝ ↦ v ^ (-σ)) phaseInterval u‖ := by
    rw [modelPhaseError] at he ⊢
    exact norm_le_norm_sub_add _ _
  calc
    ‖iteratedDerivWithin (1 + p) (F i) phaseInterval u‖ ≤ c + B p := by
      simpa [Nat.add_comm] using htriangle.trans (add_le_add he.le (hB p u))
    _ ≤ K := by
      dsimp [c, K]
      have hcOne : (2 : ℝ) ^ (-σ) / 2 ≤ 1 := by
        have hp : (2 : ℝ) ^ (-σ) ≤ 1 := by
          simpa using Real.rpow_le_one_of_one_le_of_nonpos (by norm_num : (1 : ℝ) ≤ 2)
            (neg_nonpos.mpr hσ.le)
        linarith [Real.rpow_nonneg (by norm_num : (0 : ℝ) ≤ 2) (-σ)]
      calc
        (2 : ℝ) ^ (-σ) / 2 + B p ≤ 1 + B p := by linarith
        _ ≤ 1 + ∑ q ∈ Finset.range (P + 1), B q := by
          gcongr
          exact Finset.single_le_sum (fun q _ ↦ hBnonneg q) hpMem

private theorem norm_integral_oscillatory_le
    {F : ℝ → ℝ} {A B T c K : ℝ}
    (hAB : A ≤ B) (hsub : Set.Icc A B ⊆ phaseInterval)
    (hF : ContDiffOn ℝ 2 F phaseInterval) (hT : 0 < T) (hc : 0 < c)
    (hfirst : ∀ x ∈ Set.Icc A B,
      c ≤ iteratedDerivWithin 1 F phaseInterval x)
    (hsecond : ∀ x ∈ Set.Icc A B,
      ‖iteratedDerivWithin 2 F phaseInterval x‖ ≤ K) :
    ‖∫ x in A..B, (𝐞 (T * F x) : ℂ)‖ ≤ (2 / c + K / c ^ 2) / T := by
  let d : ℝ → ℝ := iteratedDerivWithin 1 F phaseInterval
  let dd : ℝ → ℝ := iteratedDerivWithin 2 F phaseInterval
  let τ : ℂ := (2 * Real.pi : ℂ) * Complex.I
  let v : ℝ → ℂ := fun x ↦ 𝐞 (T * F x)
  let vp : ℝ → ℂ := fun x ↦ τ * T * d x * v x
  let z : ℝ → ℂ := fun x ↦ τ * T * d x
  let u : ℝ → ℂ := z⁻¹
  let up : ℝ → ℂ := fun x ↦ -(τ * T * dd x) / (z x) ^ 2
  have hunique : UniqueDiffOn ℝ phaseInterval := uniqueDiffOn_Icc (by
    norm_num [phaseInterval])
  have hτnorm : ‖τ‖ = 2 * Real.pi := by
    dsimp [τ]
    simp [Complex.norm_real, Real.norm_eq_abs, Real.pi_pos.le]
  have hτone : 1 ≤ ‖τ‖ := by rw [hτnorm]; linarith [Real.pi_gt_three]
  have hderivF (x : ℝ) (hx : x ∈ Set.Icc A B) :
      HasDerivWithinAt F (d x) (Set.Icc A B) x := by
    dsimp [d]
    rw [iteratedDerivWithin_one]
    exact ((hF.differentiableOn (by norm_num) x (hsub hx)).hasDerivWithinAt).mono hsub
  have hderivD (x : ℝ) (hx : x ∈ Set.Icc A B) :
      HasDerivWithinAt d (dd x) (Set.Icc A B) x := by
    simpa only [d, dd, show (2 : ℕ) = 1 + 1 by omega,
      iteratedDerivWithin_succ] using
        ((hF.differentiableOn_iteratedDerivWithin (m := 1) (by norm_num) hunique
          x (hsub hx)).hasDerivWithinAt).mono hsub
  have hv (x : ℝ) (hx : x ∈ Set.Icc A B) :
      HasDerivWithinAt v (vp x) (Set.Icc A B) x := by
    have hinner := (hderivF x hx).const_mul T
    have h := (Real.hasDerivAt_fourierChar (T * F x)).hasDerivWithinAt.scomp
      x hinner (Set.mapsTo_univ _ _)
    dsimp only [v]
    change HasDerivWithinAt ((𝐞 · : ℝ → ℂ) ∘ fun y ↦ T * F y)
      (vp x) (Set.Icc A B) x
    rw [show vp x = (T * d x : ℝ) •
        ((2 * Real.pi : ℂ) * Complex.I * 𝐞 (T * F x)) by
      dsimp [vp, τ]
      simp only [Complex.ofReal_mul]
      ring]
    exact h
  have hz (x : ℝ) (hx : x ∈ Set.Icc A B) :
      HasDerivWithinAt z (τ * T * dd x) (Set.Icc A B) x := by
    have h := (hderivD x hx).ofReal_comp.const_mul (τ * T)
    dsimp [z]
    simpa only [mul_assoc] using h
  have hz_ne (x : ℝ) (hx : x ∈ Set.Icc A B) : z x ≠ 0 := by
    have hdpos : 0 < d x := hc.trans_le (hfirst x hx)
    dsimp [z]
    exact mul_ne_zero (mul_ne_zero (norm_ne_zero_iff.mp
      (ne_of_gt (lt_of_lt_of_le zero_lt_one hτone)))
        (Complex.ofReal_ne_zero.mpr hT.ne'))
      (Complex.ofReal_ne_zero.mpr hdpos.ne')
  have hu (x : ℝ) (hx : x ∈ Set.Icc A B) :
      HasDerivWithinAt u (up x) (Set.Icc A B) x := by
    have h := (hz x hx).inv (hz_ne x hx)
    simpa only [u, up, Pi.inv_apply, div_eq_mul_inv] using h
  have hvCont : ContinuousOn v (Set.Icc A B) :=
    fun x hx ↦ (hv x hx).continuousWithinAt
  have huCont : ContinuousOn u (Set.Icc A B) :=
    fun x hx ↦ (hu x hx).continuousWithinAt
  have hvpCont : ContinuousOn vp (Set.Icc A B) := by
    have hdCont : ContinuousOn d (Set.Icc A B) := by
      dsimp [d]
      exact (hF.continuousOn_iteratedDerivWithin (by norm_num) hunique).mono hsub
    dsimp [vp]
    fun_prop
  have hupCont : ContinuousOn up (Set.Icc A B) := by
    have hdCont : ContinuousOn d (Set.Icc A B) := by
      dsimp [d]
      exact (hF.continuousOn_iteratedDerivWithin (by norm_num) hunique).mono hsub
    have hddCont : ContinuousOn dd (Set.Icc A B) := by
      dsimp [dd]
      exact (hF.continuousOn_iteratedDerivWithin (by norm_num) hunique).mono hsub
    have hdC : ContinuousOn (fun x ↦ (d x : ℂ)) (Set.Icc A B) :=
      Complex.continuous_ofReal.comp_continuousOn hdCont
    have hddC : ContinuousOn (fun x ↦ (dd x : ℂ)) (Set.Icc A B) :=
      Complex.continuous_ofReal.comp_continuousOn hddCont
    dsimp [up, z]
    exact ((continuousOn_const.mul continuousOn_const).mul hddC).neg.div
      (((continuousOn_const.mul continuousOn_const).mul hdC).pow 2)
      (fun x hx ↦ pow_ne_zero _ (hz_ne x hx))
  have hu' : ∀ x ∈ Set.uIcc A B, HasDerivWithinAt u (up x) (Set.uIcc A B) x := by
    simpa [Set.uIcc_of_le hAB] using hu
  have hv' : ∀ x ∈ Set.uIcc A B, HasDerivWithinAt v (vp x) (Set.uIcc A B) x := by
    simpa [Set.uIcc_of_le hAB] using hv
  have hupInt : IntervalIntegrable up MeasureTheory.volume A B := by
    have hc : ContinuousOn up (Set.uIcc A B) := by
      simpa [Set.uIcc_of_le hAB] using hupCont
    exact hc.intervalIntegrable
  have hvpInt : IntervalIntegrable vp MeasureTheory.volume A B := by
    have hc : ContinuousOn vp (Set.uIcc A B) := by
      simpa [Set.uIcc_of_le hAB] using hvpCont
    exact hc.intervalIntegrable
  have hparts := intervalIntegral.integral_mul_deriv_eq_deriv_mul_of_hasDerivWithinAt
    hu' hv' hupInt hvpInt
  have huv (x : ℝ) (hx : x ∈ Set.Icc A B) : u x * vp x = v x := by
    change (z x)⁻¹ * (z x * v x) = v x
    rw [← mul_assoc, inv_mul_cancel₀ (hz_ne x hx), one_mul]
  have hIntegral : (∫ x in A..B, v x) =
      u B * v B - u A * v A - ∫ x in A..B, up x * v x := by
    rw [← hparts]
    apply intervalIntegral.integral_congr
    intro x hx
    exact (huv x (by simpa [Set.uIcc_of_le hAB] using hx)).symm
  have huBound (x : ℝ) (hx : x ∈ Set.Icc A B) : ‖u x‖ ≤ 1 / (T * c) := by
    rw [show ‖u x‖ = (‖τ‖ * T * |d x|)⁻¹ by
      dsimp [u, z]
      rw [norm_inv, norm_mul, norm_mul]
      simp only [Complex.norm_real, Real.norm_eq_abs, abs_of_pos hT]]
    have hd : c ≤ |d x| := (hfirst x hx).trans (le_abs_self _)
    rw [one_div]
    apply (inv_le_inv₀
      (mul_pos (mul_pos (lt_of_lt_of_le zero_lt_one hτone) hT)
        (lt_of_lt_of_le hc hd)) (mul_pos hT hc)).2
    calc
      T * c ≤ ‖τ‖ * T * c := by
        calc
          T * c = 1 * (T * c) := by ring
          _ ≤ ‖τ‖ * (T * c) := mul_le_mul_of_nonneg_right hτone (by positivity)
          _ = ‖τ‖ * T * c := by ring
      _ ≤ ‖τ‖ * T * |d x| := by gcongr
  have hupBound (x : ℝ) (hx : x ∈ Set.Icc A B) :
      ‖up x‖ ≤ K / (T * c ^ 2) := by
    have hdpos : 0 < d x := hc.trans_le (hfirst x hx)
    have hzNorm : ‖z x‖ = ‖τ‖ * T * |d x| := by
      dsimp [z]
      rw [norm_mul, norm_mul]
      simp only [Complex.norm_real, Real.norm_eq_abs, abs_of_pos hT]
    have hnumNorm : ‖-(τ * T * dd x)‖ = ‖τ‖ * T * |dd x| := by
      rw [norm_neg, norm_mul, norm_mul]
      simp only [Complex.norm_real, Real.norm_eq_abs, abs_of_pos hT]
    rw [show ‖up x‖ = ‖-(τ * T * dd x)‖ / ‖z x‖ ^ 2 by
      simp only [up, norm_div, norm_pow], hnumNorm, hzNorm]
    rw [show ‖τ‖ * T * |dd x| / (‖τ‖ * T * |d x|) ^ 2 =
        |dd x| / (‖τ‖ * T * |d x| ^ 2) by
      field_simp [hT.ne', (lt_of_lt_of_le zero_lt_one hτone).ne', hdpos.ne']]
    have hden : T * c ^ 2 ≤ ‖τ‖ * T * |d x| ^ 2 := by
      have hd : c ≤ |d x| := (hfirst x hx).trans (le_abs_self _)
      calc
        T * c ^ 2 ≤ T * |d x| ^ 2 := by gcongr
        _ ≤ ‖τ‖ * T * |d x| ^ 2 := by
          have : 0 ≤ T * |d x| ^ 2 := by positivity
          calc
            T * |d x| ^ 2 = 1 * (T * |d x| ^ 2) := by ring
            _ ≤ ‖τ‖ * (T * |d x| ^ 2) := mul_le_mul_of_nonneg_right hτone this
            _ = ‖τ‖ * T * |d x| ^ 2 := by ring
    have hK0 : 0 ≤ K := (norm_nonneg (dd x)).trans (hsecond x hx)
    calc
      |dd x| / (‖τ‖ * T * |d x| ^ 2) ≤
          K / (‖τ‖ * T * |d x| ^ 2) :=
        div_le_div_of_nonneg_right (by simpa [Real.norm_eq_abs] using hsecond x hx) (by positivity)
      _ ≤ K / (T * c ^ 2) :=
        div_le_div_of_nonneg_left hK0 (by positivity) hden
  rw [hIntegral]
  calc
    ‖u B * v B - u A * v A - ∫ x in A..B, up x * v x‖ ≤
        ‖u B‖ + ‖u A‖ + ‖∫ x in A..B, up x * v x‖ := by
      calc
        _ ≤ ‖u B * v B - u A * v A‖ + ‖∫ x in A..B, up x * v x‖ :=
          norm_sub_le _ _
        _ ≤ (‖u B * v B‖ + ‖u A * v A‖) + ‖∫ x in A..B, up x * v x‖ := by
          gcongr
          exact norm_sub_le _ _
        _ = _ := by simp [v]
    _ ≤ 2 / (T * c) + (B - A) * (K / (T * c ^ 2)) := by
      have hA : A ∈ Set.Icc A B := ⟨le_rfl, hAB⟩
      have hB : B ∈ Set.Icc A B := ⟨hAB, le_rfl⟩
      have hint := intervalIntegral.norm_integral_le_of_norm_le_const
        (a := A) (b := B) (C := K / (T * c ^ 2)) (fun x hx ↦ by
          have hx' : x ∈ Set.Icc A B := by
            rw [Set.uIoc_of_le hAB] at hx
            exact ⟨hx.1.le, hx.2⟩
          calc
            ‖up x * v x‖ = ‖up x‖ := by simp [v]
            _ ≤ K / (T * c ^ 2) := hupBound x hx')
      rw [abs_of_nonneg (sub_nonneg.mpr hAB)] at hint
      calc
        ‖u B‖ + ‖u A‖ + ‖∫ x in A..B, up x * v x‖ ≤
            1 / (T * c) + 1 / (T * c) +
              K / (T * c ^ 2) * (B - A) := by
          exact add_le_add (add_le_add (huBound B hB) (huBound A hA)) hint
        _ = 2 / (T * c) + (B - A) * (K / (T * c ^ 2)) := by ring
    _ ≤ (2 / c + K / c ^ 2) / T := by
      have hlen : B - A ≤ 1 := by
        have hA := (hsub ⟨le_rfl, hAB⟩).1
        have hB := (hsub ⟨hAB, le_rfl⟩).2
        linarith
      have hK0 : 0 ≤ K := (norm_nonneg (dd A)).trans (hsecond A ⟨le_rfl, hAB⟩)
      field_simp
      nlinarith

private def oscillatorySumConstant (s : ℕ) (K c : ℝ) : ℝ :=
  1 +
    (∑ m ∈ Finset.range s, |EulerMaclaurin.saw (m + 2) 0| *
      (2 * (((m + 1).factorial : ℝ) * (2 * Real.pi + 1) ^ (m + 1)))) +
    EulerMaclaurin.sawBound (s + 1) *
      (((s + 1).factorial : ℝ) * (2 * Real.pi + 1) ^ (s + 1)) +
    (2 / c + K / c ^ 2)

private lemma oscillatorySumConstant_nonneg
    (s : ℕ) {K c : ℝ} (hK : 0 ≤ K) (hc : 0 < c) :
    0 ≤ oscillatorySumConstant s K c := by
  have hsum : 0 ≤ ∑ m ∈ Finset.range s, |EulerMaclaurin.saw (m + 2) 0| *
      (2 * (((m + 1).factorial : ℝ) * (2 * Real.pi + 1) ^ (m + 1))) :=
    Finset.sum_nonneg fun m _ ↦ by positivity
  have hrem : 0 ≤ EulerMaclaurin.sawBound (s + 1) *
      (((s + 1).factorial : ℝ) * (2 * Real.pi + 1) ^ (s + 1)) :=
    mul_nonneg EulerMaclaurin.sawBound_nonneg (by positivity)
  have hD : 0 ≤ 2 / c + K / c ^ 2 := by positivity
  dsimp only [oscillatorySumConstant]
  linarith

private theorem norm_oscillatory_sum_le
    (s : ℕ) {F : ℝ → ℝ} {T N K c : ℝ} {a b : ℕ}
    (hs : 1 ≤ s) (hab : a ≤ b) (hN : 1 ≤ N) (hT : 1 ≤ T) (hK : 1 ≤ K) (hc : 0 < c)
    (hF : ContDiffOn ℝ ∞ F phaseInterval)
    (hfirst : ∀ u ∈ phaseInterval,
      c ≤ iteratedDerivWithin 1 F phaseInterval u)
    (hderiv : ∀ k : ℕ, 1 ≤ k → k ≤ s + 1 → ∀ u ∈ phaseInterval,
      ‖iteratedDerivWithin k F phaseInterval u‖ ≤ K)
    (ha : N ≤ a) (hb : (b : ℝ) ≤ 2 * N)
    (hratio : K * T / N ≤ 1)
    (hremainder : N * (K * T / N) ^ (s + 1) ≤ 1) :
    ‖∑ n ∈ Finset.Icc a b, (𝐞 (T * F ((n : ℝ) / N)) : ℂ)‖ ≤
      oscillatorySumConstant s K c * (1 + N / T) := by
  have hNpos : 0 < N := zero_lt_one.trans_le hN
  have hTpos : 0 < T := zero_lt_one.trans_le hT
  have hKpos : 0 < K := zero_lt_one.trans_le hK
  have hrnonneg : 0 ≤ K * T / N := by positivity
  let interval : Set ℝ := Set.Icc (a : ℝ) (b : ℝ)
  have hmap : Set.MapsTo (N⁻¹ * ·) interval phaseInterval := by
    intro x hx
    simp only [interval, phaseInterval, Set.mem_Icc] at hx ⊢
    constructor
    · rw [inv_mul_eq_div, le_div_iff₀ hNpos]
      simpa using ha.trans hx.1
    · rw [inv_mul_eq_div, div_le_iff₀ hNpos]
      exact hx.2.trans hb
  by_cases heq : a = b
  · subst b
    simp only [Finset.Icc_self, Finset.sum_singleton]
    rw [show ‖(𝐞 (T * F ((a : ℝ) / N)) : ℂ)‖ = 1 by simp]
    have : 0 ≤ N / T := by positivity
    have hconst : 1 ≤ oscillatorySumConstant s K c := by
      have hsum : 0 ≤ ∑ m ∈ Finset.range s, |EulerMaclaurin.saw (m + 2) 0| *
          (2 * (((m + 1).factorial : ℝ) * (2 * Real.pi + 1) ^ (m + 1))) :=
        Finset.sum_nonneg fun m _ ↦ by positivity
      have hrem : 0 ≤ EulerMaclaurin.sawBound (s + 1) *
          (((s + 1).factorial : ℝ) * (2 * Real.pi + 1) ^ (s + 1)) := by
        exact mul_nonneg EulerMaclaurin.sawBound_nonneg (by positivity)
      have hD : 0 ≤ 2 / c + K / c ^ 2 := by positivity
      dsimp only [oscillatorySumConstant]
      linarith
    nlinarith
  have hablt : a < b := hab.lt_of_ne heq
  have hunique : UniqueDiffOn ℝ interval := uniqueDiffOn_Icc (by exact_mod_cast hablt)
  let A : ℕ → ℝ := fun k ↦ (k.factorial : ℝ) * (2 * Real.pi + 1) ^ k
  have hosc (k : ℕ) (hk : 1 ≤ k) (hkS : k ≤ s + 1) (x : ℝ) (hx : x ∈ interval) :
      ‖iteratedDerivWithin k (fun y ↦ (𝐞 (T * F (y / N)) : ℂ)) interval x‖ ≤
        A k * (K * T / N) ^ k := by
    simpa only [A] using norm_iteratedDerivWithin_oscillatory_le
      (hF.of_le (ENat.natCast_le_of_coe_top_le_withTop le_rfl k)) hunique hx hmap hT hN hK
      (fun q hq hqk u hu ↦ hderiv q hq (hqk.trans hkS) u hu)
  have hoscConst (k : ℕ) (hk : 1 ≤ k) (hkS : k ≤ s + 1)
      (x : ℝ) (hx : x ∈ interval) :
      ‖iteratedDerivWithin k (fun y ↦ (𝐞 (T * F (y / N)) : ℂ)) interval x‖ ≤ A k := by
    exact (hosc k hk hkS x hx).trans <| by
      have hA0 : 0 ≤ A k := by dsimp [A]; positivity
      simpa only [mul_one] using mul_le_mul_of_nonneg_left
        (pow_le_one₀ hrnonneg hratio) hA0
  have hIntegral :
      ‖∫ x in (a : ℝ)..b, (𝐞 (T * F (x / N)) : ℂ)‖ ≤
        N * ((2 / c + K / c ^ 2) / T) := by
    have hAB : (a : ℝ) / N ≤ (b : ℝ) / N := div_le_div_of_nonneg_right
      (by exact_mod_cast hab) hNpos.le
    have hsub : Set.Icc ((a : ℝ) / N) ((b : ℝ) / N) ⊆ phaseInterval := by
      intro x hx
      constructor
      · exact ((le_div_iff₀ hNpos).2 (by simpa using ha)).trans hx.1
      · exact hx.2.trans ((div_le_iff₀ hNpos).2 (by simpa using hb))
    have hy := norm_integral_oscillatory_le hAB hsub
      (hF.of_le (ENat.natCast_le_of_coe_top_le_withTop le_rfl 2)) hTpos hc
      (fun x hx ↦ hfirst x (hsub hx))
      (fun x hx ↦ hderiv 2 (by norm_num) (by omega) x (hsub hx))
    have hscale := intervalIntegral.integral_comp_div
      (fun y : ℝ ↦ (𝐞 (T * F y) : ℂ)) (a := (a : ℝ)) (b := (b : ℝ)) hNpos.ne'
    rw [hscale, norm_smul, Real.norm_eq_abs, abs_of_pos hNpos]
    exact mul_le_mul_of_nonneg_left hy hNpos.le
  have hEM := oscillatory_sum_eq_eulerMaclaurin (T := T) (N := N) hab
    (hF.of_le (ENat.natCast_le_of_coe_top_le_withTop le_rfl (s + 1))) hmap
  let Ccorr : ℝ := ∑ m ∈ Finset.range s,
    |EulerMaclaurin.saw (m + 2) 0| * (2 * A (m + 1))
  let Crem : ℝ := EulerMaclaurin.sawBound (s + 1) * A (s + 1)
  let D : ℝ := 2 / c + K / c ^ 2
  let C : ℝ := 1 + Ccorr + Crem + D
  have hA_nonneg (k : ℕ) : 0 ≤ A k := by
    dsimp [A]
    positivity
  have hCcorr : 0 ≤ Ccorr := Finset.sum_nonneg fun m _ ↦ by
    positivity
  have hCrem : 0 ≤ Crem := mul_nonneg EulerMaclaurin.sawBound_nonneg (hA_nonneg _)
  have hD : 0 < D := by
    dsimp [D]
    positivity
  have hC : 0 < C := by dsimp [C]; positivity
  have hendpoint :
      ‖(2⁻¹ : ℝ) •
        ((𝐞 (T * F ((a : ℝ) / N)) : ℂ) + (𝐞 (T * F ((b : ℝ) / N)) : ℂ))‖ ≤ 1 := by
    calc
      _ ≤ |(2⁻¹ : ℝ)| *
          (‖(𝐞 (T * F ((a : ℝ) / N)) : ℂ)‖ +
            ‖(𝐞 (T * F ((b : ℝ) / N)) : ℂ)‖) := by
        rw [norm_smul, Real.norm_eq_abs]
        gcongr
        exact norm_add_le _ _
      _ = 1 := by norm_num
  have hcorr :
      ‖∑ m ∈ Finset.range s, (-1 : ℝ) ^ m • EulerMaclaurin.saw (m + 2) 0 •
        (iteratedDerivWithin (m + 1) (fun x ↦ (𝐞 (T * F (x / N)) : ℂ)) interval b -
          iteratedDerivWithin (m + 1) (fun x ↦ (𝐞 (T * F (x / N)) : ℂ)) interval a)‖ ≤
        Ccorr := by
    calc
      _ ≤ ∑ m ∈ Finset.range s,
          ‖(-1 : ℝ) ^ m • EulerMaclaurin.saw (m + 2) 0 •
            (iteratedDerivWithin (m + 1) (fun x ↦ (𝐞 (T * F (x / N)) : ℂ)) interval b -
              iteratedDerivWithin (m + 1) (fun x ↦ (𝐞 (T * F (x / N)) : ℂ)) interval a)‖ :=
        norm_sum_le _ _
      _ ≤ Ccorr := by
        apply Finset.sum_le_sum
        intro m hm
        have hmS : m < s := Finset.mem_range.mp hm
        have haI : (a : ℝ) ∈ interval := ⟨le_rfl, by exact_mod_cast hab⟩
        have hbI : (b : ℝ) ∈ interval := ⟨by exact_mod_cast hab, le_rfl⟩
        rw [norm_smul, norm_smul, Real.norm_eq_abs, Real.norm_eq_abs,
          abs_pow, abs_neg, abs_one, one_pow, one_mul]
        calc
          |EulerMaclaurin.saw (m + 2) 0| * ‖_ - _‖ ≤
              |EulerMaclaurin.saw (m + 2) 0| * (A (m + 1) + A (m + 1)) := by
            gcongr
            exact (norm_sub_le _ _).trans (add_le_add
              (hoscConst (m + 1) (by omega) (by omega) b hbI)
              (hoscConst (m + 1) (by omega) (by omega) a haI))
          _ = |EulerMaclaurin.saw (m + 2) 0| * (2 * A (m + 1)) := by ring
  have hrem :
      ‖∫ x in (a : ℝ)..b, EulerMaclaurin.saw (s + 1) x •
        iteratedDerivWithin (s + 1) (fun y ↦ (𝐞 (T * F (y / N)) : ℂ)) interval x‖ ≤
        Crem := by
    have hend : ((a : ℤ) : ℝ) + (b - a : ℕ) = (b : ℝ) := by
      exact_mod_cast Nat.add_sub_of_le hab
    have hraw := EulerMaclaurin.norm_integral_saw_smul_iteratedDerivWithin_le
      (a := (a : ℤ)) (n := b - a) (s := s)
      (f := fun y ↦ (𝐞 (T * F (y / N)) : ℂ)) (t := interval)
      (C := A (s + 1) * (K * T / N) ^ (s + 1)) (fun x hx ↦ by
        apply hosc (s + 1) (by omega) le_rfl x
        rw [hend] at hx
        simpa only [interval, Int.cast_natCast] using hx)
    rw [hend] at hraw
    calc
      _ ≤ EulerMaclaurin.sawBound (s + 1) *
          (A (s + 1) * (K * T / N) ^ (s + 1)) * (b - a : ℕ) := hraw
      _ ≤ Crem := by
        have hlen : (b - a : ℝ) ≤ N := by
          linarith
        dsimp [Crem]
        have hprod : (b - a : ℝ) * (K * T / N) ^ (s + 1) ≤ 1 := by
          exact (mul_le_mul_of_nonneg_right hlen (pow_nonneg hrnonneg _)).trans hremainder
        calc
          EulerMaclaurin.sawBound (s + 1) *
              (A (s + 1) * (K * T / N) ^ (s + 1)) * (b - a : ℕ) =
              Crem * ((b - a : ℝ) * (K * T / N) ^ (s + 1)) := by
            dsimp [Crem]
            rw [Nat.cast_sub hab]
            ring
          _ ≤ Crem * 1 := mul_le_mul_of_nonneg_left hprod hCrem
          _ = Crem := mul_one _
  change _ ≤ C * (1 + N / T)
  rw [hEM]
  have hsum : ‖(∫ x in (a : ℝ)..b, (𝐞 (T * F (x / N)) : ℂ)) +
      (2⁻¹ : ℝ) • ((𝐞 (T * F ((a : ℝ) / N)) : ℂ) +
        (𝐞 (T * F ((b : ℝ) / N)) : ℂ)) +
      (∑ m ∈ Finset.range s, (-1 : ℝ) ^ m • EulerMaclaurin.saw (m + 2) 0 •
        (iteratedDerivWithin (m + 1) (fun x ↦ (𝐞 (T * F (x / N)) : ℂ)) interval b -
          iteratedDerivWithin (m + 1) (fun x ↦ (𝐞 (T * F (x / N)) : ℂ)) interval a)) +
      (-1 : ℝ) ^ s • (∫ x in (a : ℝ)..b, EulerMaclaurin.saw (s + 1) x •
        iteratedDerivWithin (s + 1) (fun y ↦ (𝐞 (T * F (y / N)) : ℂ)) interval x)‖ ≤
      D * (N / T) + 1 + Ccorr + Crem := by
    calc
      _ ≤ ‖∫ x in (a : ℝ)..b, (𝐞 (T * F (x / N)) : ℂ)‖ +
          ‖(2⁻¹ : ℝ) • ((𝐞 (T * F ((a : ℝ) / N)) : ℂ) +
            (𝐞 (T * F ((b : ℝ) / N)) : ℂ))‖ +
          ‖∑ m ∈ Finset.range s, (-1 : ℝ) ^ m • EulerMaclaurin.saw (m + 2) 0 •
            (iteratedDerivWithin (m + 1) (fun x ↦ (𝐞 (T * F (x / N)) : ℂ)) interval b -
              iteratedDerivWithin (m + 1) (fun x ↦ (𝐞 (T * F (x / N)) : ℂ)) interval a)‖ +
          ‖(-1 : ℝ) ^ s • (∫ x in (a : ℝ)..b, EulerMaclaurin.saw (s + 1) x •
            iteratedDerivWithin (s + 1) (fun y ↦ (𝐞 (T * F (y / N)) : ℂ)) interval x)‖ := by
        calc
          _ ≤ ‖_ + _ + _‖ + ‖_‖ := norm_add_le _ _
          _ ≤ (‖_ + _‖ + ‖_‖) + ‖_‖ := by
            gcongr
            exact norm_add_le _ _
          _ ≤ ((‖_‖ + ‖_‖) + ‖_‖) + ‖_‖ := by gcongr; exact norm_add_le _ _
      _ ≤ D * (N / T) + 1 + Ccorr + Crem := by
        have hI : ‖∫ x in (a : ℝ)..b, (𝐞 (T * F (x / N)) : ℂ)‖ ≤ D * (N / T) := by
          calc
            _ ≤ N * (D / T) := by simpa only [D] using hIntegral
            _ = D * (N / T) := by ring
        have hremSmul :
            ‖(-1 : ℝ) ^ s • (∫ x in (a : ℝ)..b, EulerMaclaurin.saw (s + 1) x •
              iteratedDerivWithin (s + 1) (fun y ↦ (𝐞 (T * F (y / N)) : ℂ)) interval x)‖ ≤
              Crem := by
          rw [norm_smul, Real.norm_eq_abs, abs_pow, abs_neg, abs_one, one_pow, one_mul]
          exact hrem
        linarith
  calc
    _ ≤ D * (N / T) + 1 + Ccorr + Crem := hsum
    _ ≤ C * (1 + N / T) := by
      have hNT : 0 ≤ N / T := by positivity
      dsimp [C]
      nlinarith

private theorem isExponentSumBound_sub_one
    {α : ℝ≥0} (hα : 1 < α) :
    IsExponentSumBound α ((α : ℝ) - 1) := by
  intro N T F a b hN hT hTunbounded hNT hF hab
  apply (isPowerBounded_iff_forall_pos
    (exponentialSum F T N a b) T ((α : ℝ) - 1) hT hTunbounded).2
  intro ε hε
  have hαR : (1 : ℝ) < α := by exact_mod_cast hα
  let η : ℝ := min (((α : ℝ) - 1) / 4) (ε / 4)
  have hη : 0 < η := lt_min (by linarith) (by linarith)
  have hηε : η < ε := lt_of_le_of_lt (min_le_right _ _) (by linarith)
  let δ : ℝ := (α : ℝ) - η - 1
  have hδ : 0 < δ := by
    dsimp [δ, η]
    have := min_le_left (((α : ℝ) - 1) / 4) (ε / 4)
    linarith
  obtain ⟨s, hs⟩ := exists_nat_ge ((1 + δ / 2) / δ)
  have hsone : 1 ≤ s := by
    by_contra hs0
    have : s = 0 := by omega
    subst s
    norm_num at hs
    have : 0 < (1 + δ / 2) / δ := by positivity
    linarith
  have hsδ : 1 + δ / 2 ≤ (s : ℝ) * δ := by
    exact (div_le_iff₀ hδ).1 (by simpa [mul_comm] using hs)
  let q : ℝ := δ / (2 * (s + 1 : ℕ))
  have hq : 0 < q := by dsimp [q]; positivity
  have hqδ : q ≤ δ := by
    dsimp [q]
    push_cast
    apply (div_le_iff₀ (by positivity : 0 < (2 : ℝ) * ((s : ℝ) + 1))).2
    have hs0 : (0 : ℝ) ≤ s := Nat.cast_nonneg s
    calc
      δ = δ * 1 := by ring
      _ ≤ δ * (2 * ((s : ℝ) + 1)) := by
        gcongr
        nlinarith
  obtain ⟨c, K, hc, hK, hderivEventually⟩ :=
    eventually_model_phase_deriv_bounds hF s
  have hTtop : Tendsto T atTop atTop := by
    rw [VariableObject.IsUnbounded] at hTunbounded
    exact hTunbounded.congr' <| Filter.Eventually.of_forall fun i ↦ by
      rw [Real.norm_eq_abs, abs_of_nonneg (le_trans zero_le_one (hT i))]
  have hKpower : ∀ᶠ i in atTop, K ≤ T i ^ q := by
    have h := tendsto_atTop.1 ((tendsto_rpow_atTop hq).comp hTtop) K
    exact h
  have hNTbounds := hNT.eventually_between
    (Filter.Eventually.of_forall hT) hη
  let C := oscillatorySumConstant s K c
  have hC : 0 ≤ C := oscillatorySumConstant_nonneg s (zero_le_one.trans hK) hc
  refine Asymptotics.IsBigO.of_bound (2 * C) ?_
  filter_upwards [hderivEventually, hKpower, hNTbounds] with i hiDeriv hiK hiNT
  have hTi : 0 < T i := zero_lt_one.trans_le (hT i)
  have hNi : 0 < N i := zero_lt_one.trans_le (hN i)
  have hqmul : q * (s + 1 : ℕ) = δ / 2 := by
    dsimp [q]
    field_simp
  have hKpow : K ^ (s + 1) ≤ T i ^ (δ / 2) := by
    calc
      K ^ (s + 1) ≤ (T i ^ q) ^ (s + 1) :=
        pow_le_pow_left₀ (zero_le_one.trans hK) hiK _
      _ = T i ^ (q * (s + 1 : ℕ)) :=
        (Real.rpow_mul_natCast hTi.le q (s + 1)).symm
      _ = T i ^ (δ / 2) := by rw [hqmul]
  have hratio : K * T i / N i ≤ 1 := by
    apply (div_le_one hNi).2
    calc
      K * T i ≤ T i ^ q * T i := mul_le_mul_of_nonneg_right hiK hTi.le
      _ = T i ^ (q + 1) := by
        simpa only [Real.rpow_one] using (Real.rpow_add hTi q 1).symm
      _ ≤ T i ^ (δ + 1) :=
        Real.rpow_le_rpow_of_exponent_le (hT i) (by linarith)
      _ = T i ^ ((α : ℝ) - η) := by
        congr 1
        dsimp [δ]
        ring
      _ ≤ N i := hiNT.1
  have hremainder : N i * (K * T i / N i) ^ (s + 1) ≤ 1 := by
    have hbase : T i ^ ((α : ℝ) - η) ≤ N i := hiNT.1
    have hNpow : T i ^ (((α : ℝ) - η) * s) ≤ N i ^ s := by
      calc
        T i ^ (((α : ℝ) - η) * s) = (T i ^ ((α : ℝ) - η)) ^ s :=
          Real.rpow_mul_natCast hTi.le ((α : ℝ) - η) s
        _ ≤ N i ^ s := pow_le_pow_left₀ (Real.rpow_nonneg hTi.le _) hbase _
    have hexponent : δ / 2 + (s + 1 : ℕ) ≤ ((α : ℝ) - η) * s := by
      dsimp [δ] at hsδ ⊢
      push_cast
      linarith
    have hnum : K ^ (s + 1) * T i ^ (s + 1) ≤ N i ^ s := by
      calc
        K ^ (s + 1) * T i ^ (s + 1) ≤
            T i ^ (δ / 2) * T i ^ (s + 1 : ℕ) := by
          gcongr
        _ = T i ^ (δ / 2 + (s + 1 : ℕ)) := by
          rw [← Real.rpow_natCast]
          exact (Real.rpow_add hTi (δ / 2) (s + 1 : ℕ)).symm
        _ ≤ T i ^ (((α : ℝ) - η) * s) :=
          Real.rpow_le_rpow_of_exponent_le (hT i) hexponent
        _ ≤ N i ^ s := hNpow
    rw [show N i * (K * T i / N i) ^ (s + 1) =
        (K ^ (s + 1) * T i ^ (s + 1)) / N i ^ s by
      rw [div_pow, mul_pow, pow_succ]
      field_simp [hNi.ne', pow_ne_zero _ hNi.ne']
      ring]
    exact (div_le_one (pow_pos hNi _)).2 hnum
  have hsum : ‖exponentialSum F T N a b i‖ ≤ C * (1 + N i / T i) := by
    change ‖∑ n ∈ Finset.Icc (a i) (b i),
      (𝐞 (T i * F i ((n : ℝ) / N i)) : ℂ)‖ ≤ C * (1 + N i / T i)
    by_cases habi : a i ≤ b i
    · exact norm_oscillatory_sum_le s hsone habi
        (hN i) (hT i) hK hc (hF.1 i)
        (fun u hu ↦ (hiDeriv ⟨u, hu⟩).1)
        (fun k hk hkS u hu ↦ (hiDeriv ⟨u, hu⟩).2 k hk hkS)
        (hab i).1 (hab i).2 hratio hremainder
    · rw [Finset.Icc_eq_empty habi, Finset.sum_empty, norm_zero]
      exact mul_nonneg hC (by positivity)
  have hTargetNonneg : 0 ≤ (α : ℝ) - 1 + ε := by linarith
  have hNTdiv : N i / T i ≤ T i ^ ((α : ℝ) - 1 + ε) := by
    calc
      N i / T i ≤ T i ^ ((α : ℝ) + η) / T i :=
        div_le_div_of_nonneg_right hiNT.2 hTi.le
      _ = T i ^ ((α : ℝ) + η - 1) := by
        simpa only [Real.rpow_one] using
          (Real.rpow_sub hTi ((α : ℝ) + η) 1).symm
      _ ≤ T i ^ ((α : ℝ) - 1 + ε) :=
        Real.rpow_le_rpow_of_exponent_le (hT i) (by linarith)
  have hone : 1 ≤ T i ^ ((α : ℝ) - 1 + ε) := by
    simpa using Real.one_le_rpow (hT i) hTargetNonneg
  change ‖exponentialSum F T N a b i‖ ≤
    2 * C * ‖T i ^ ((α : ℝ) - 1 + ε)‖
  rw [Real.norm_eq_abs, abs_of_nonneg (Real.rpow_nonneg hTi.le _)]
  calc
    ‖exponentialSum F T N a b i‖ ≤ C * (1 + N i / T i) := hsum
    _ ≤ C * (2 * T i ^ ((α : ℝ) - 1 + ε)) := by
      gcongr
      linarith
    _ = 2 * C * T i ^ ((α : ℝ) - 1 + ε) := by ring

private def oscillatoryErrorConstant (s : ℕ) : ℝ :=
  1 +
    (∑ m ∈ Finset.range s, |EulerMaclaurin.saw (m + 2) 0| *
      (2 * (((m + 1).factorial : ℝ) * (2 * Real.pi + 1) ^ (m + 1)))) +
    EulerMaclaurin.sawBound (s + 1) *
      (((s + 1).factorial : ℝ) * (2 * Real.pi + 1) ^ (s + 1))

private theorem norm_oscillatory_sum_sub_integral_le
    (s : ℕ) {F : ℝ → ℝ} {T N K : ℝ} {a b : ℕ}
    (hab : a < b) (hN : 1 ≤ N) (hT : 1 ≤ T) (hK : 1 ≤ K)
    (hF : ContDiffOn ℝ ∞ F phaseInterval)
    (hderiv : ∀ k : ℕ, 1 ≤ k → k ≤ s + 1 → ∀ u ∈ phaseInterval,
      ‖iteratedDerivWithin k F phaseInterval u‖ ≤ K)
    (ha : N ≤ a) (hb : (b : ℝ) ≤ 2 * N)
    (hratio : K * T / N ≤ 1)
    (hremainder : N * (K * T / N) ^ (s + 1) ≤ 1) :
    ‖(∑ n ∈ Finset.Icc a b, (𝐞 (T * F ((n : ℝ) / N)) : ℂ)) -
      ∫ x in (a : ℝ)..b, (𝐞 (T * F (x / N)) : ℂ)‖ ≤
        oscillatoryErrorConstant s := by
  have hNpos : 0 < N := zero_lt_one.trans_le hN
  have hrnonneg : 0 ≤ K * T / N := by positivity
  let interval : Set ℝ := Set.Icc (a : ℝ) (b : ℝ)
  have hmap : Set.MapsTo (N⁻¹ * ·) interval phaseInterval := by
    intro x hx
    simp only [interval, phaseInterval, Set.mem_Icc] at hx ⊢
    constructor
    · rw [inv_mul_eq_div, le_div_iff₀ hNpos]
      simpa using ha.trans hx.1
    · rw [inv_mul_eq_div, div_le_iff₀ hNpos]
      exact hx.2.trans hb
  have hunique : UniqueDiffOn ℝ interval := uniqueDiffOn_Icc (by exact_mod_cast hab)
  let A : ℕ → ℝ := fun k ↦ (k.factorial : ℝ) * (2 * Real.pi + 1) ^ k
  have hosc (k : ℕ) (hk : 1 ≤ k) (hkS : k ≤ s + 1)
      (x : ℝ) (hx : x ∈ interval) :
      ‖iteratedDerivWithin k (fun y ↦ (𝐞 (T * F (y / N)) : ℂ)) interval x‖ ≤ A k := by
    have h := norm_iteratedDerivWithin_oscillatory_le
      (hF.of_le (ENat.natCast_le_of_coe_top_le_withTop le_rfl k))
      hunique hx hmap hT hN hK
      (fun q hq hqk u hu ↦ hderiv q hq (hqk.trans hkS) u hu)
    exact h.trans <| by
      have hA0 : 0 ≤ A k := by dsimp [A]; positivity
      change A k * (K * T / N) ^ k ≤ A k
      simpa only [mul_one] using mul_le_mul_of_nonneg_left
        (pow_le_one₀ hrnonneg hratio) hA0
  have hEM := oscillatory_sum_eq_eulerMaclaurin
    (T := T) (N := N) hab.le
    (hF.of_le (ENat.natCast_le_of_coe_top_le_withTop le_rfl (s + 1))) hmap
  let Ccorr : ℝ := ∑ m ∈ Finset.range s,
    |EulerMaclaurin.saw (m + 2) 0| * (2 * A (m + 1))
  let Crem : ℝ := EulerMaclaurin.sawBound (s + 1) * A (s + 1)
  have hendpoint :
      ‖(2⁻¹ : ℝ) •
        ((𝐞 (T * F ((a : ℝ) / N)) : ℂ) + (𝐞 (T * F ((b : ℝ) / N)) : ℂ))‖ ≤ 1 := by
    calc
      _ ≤ |(2⁻¹ : ℝ)| *
          (‖(𝐞 (T * F ((a : ℝ) / N)) : ℂ)‖ +
            ‖(𝐞 (T * F ((b : ℝ) / N)) : ℂ)‖) := by
        rw [norm_smul, Real.norm_eq_abs]
        gcongr
        exact norm_add_le _ _
      _ = 1 := by norm_num
  have hcorr :
      ‖∑ m ∈ Finset.range s, (-1 : ℝ) ^ m • EulerMaclaurin.saw (m + 2) 0 •
        (iteratedDerivWithin (m + 1) (fun x ↦ (𝐞 (T * F (x / N)) : ℂ)) interval b -
          iteratedDerivWithin (m + 1) (fun x ↦ (𝐞 (T * F (x / N)) : ℂ)) interval a)‖ ≤
        Ccorr := by
    calc
      _ ≤ ∑ m ∈ Finset.range s,
          ‖(-1 : ℝ) ^ m • EulerMaclaurin.saw (m + 2) 0 •
            (iteratedDerivWithin (m + 1) (fun x ↦ (𝐞 (T * F (x / N)) : ℂ)) interval b -
              iteratedDerivWithin (m + 1) (fun x ↦ (𝐞 (T * F (x / N)) : ℂ)) interval a)‖ :=
        norm_sum_le _ _
      _ ≤ Ccorr := by
        apply Finset.sum_le_sum
        intro m hm
        have hmS := Finset.mem_range.mp hm
        have haI : (a : ℝ) ∈ interval := ⟨le_rfl, by exact_mod_cast hab.le⟩
        have hbI : (b : ℝ) ∈ interval := ⟨by exact_mod_cast hab.le, le_rfl⟩
        rw [norm_smul, norm_smul, Real.norm_eq_abs, Real.norm_eq_abs,
          abs_pow, abs_neg, abs_one, one_pow, one_mul]
        calc
          |EulerMaclaurin.saw (m + 2) 0| * ‖_ - _‖ ≤
              |EulerMaclaurin.saw (m + 2) 0| * (A (m + 1) + A (m + 1)) := by
            gcongr
            exact (norm_sub_le _ _).trans (add_le_add
              (hosc (m + 1) (by omega) (by omega) b hbI)
              (hosc (m + 1) (by omega) (by omega) a haI))
          _ = |EulerMaclaurin.saw (m + 2) 0| * (2 * A (m + 1)) := by ring
  have hrem :
      ‖∫ x in (a : ℝ)..b, EulerMaclaurin.saw (s + 1) x •
        iteratedDerivWithin (s + 1) (fun y ↦ (𝐞 (T * F (y / N)) : ℂ)) interval x‖ ≤
        Crem := by
    have hend : ((a : ℤ) : ℝ) + (b - a : ℕ) = (b : ℝ) := by
      exact_mod_cast Nat.add_sub_of_le hab.le
    have hraw := EulerMaclaurin.norm_integral_saw_smul_iteratedDerivWithin_le
      (a := (a : ℤ)) (n := b - a) (s := s)
      (f := fun y ↦ (𝐞 (T * F (y / N)) : ℂ)) (t := interval)
      (C := A (s + 1) * (K * T / N) ^ (s + 1)) (fun x hx ↦ by
        apply (norm_iteratedDerivWithin_oscillatory_le
          (hF.of_le (ENat.natCast_le_of_coe_top_le_withTop le_rfl (s + 1)))
          hunique _ hmap hT hN hK
          (fun q hq hqS u hu ↦ hderiv q hq hqS u hu))
        rw [hend] at hx
        simpa only [interval, Int.cast_natCast] using hx)
    rw [hend] at hraw
    calc
      _ ≤ EulerMaclaurin.sawBound (s + 1) *
          (A (s + 1) * (K * T / N) ^ (s + 1)) * (b - a : ℕ) := hraw
      _ ≤ Crem := by
        have hlen : (b - a : ℝ) ≤ N := by linarith
        have hprod : (b - a : ℝ) * (K * T / N) ^ (s + 1) ≤ 1 :=
          (mul_le_mul_of_nonneg_right hlen (pow_nonneg hrnonneg _)).trans hremainder
        dsimp [Crem]
        calc
          EulerMaclaurin.sawBound (s + 1) *
              (A (s + 1) * (K * T / N) ^ (s + 1)) * (b - a : ℕ) =
              Crem * ((b - a : ℝ) * (K * T / N) ^ (s + 1)) := by
            dsimp [Crem]
            rw [Nat.cast_sub hab.le]
            ring
          _ ≤ Crem := by
            have hCrem : 0 ≤ Crem := by
              dsimp [Crem, A]
              exact mul_nonneg EulerMaclaurin.sawBound_nonneg (by positivity)
            simpa only [mul_one] using mul_le_mul_of_nonneg_left hprod hCrem
  let I : ℂ := ∫ x in (a : ℝ)..b, (𝐞 (T * F (x / N)) : ℂ)
  let E : ℂ := (2⁻¹ : ℝ) •
    ((𝐞 (T * F ((a : ℝ) / N)) : ℂ) + (𝐞 (T * F ((b : ℝ) / N)) : ℂ))
  let Q : ℂ := ∑ m ∈ Finset.range s, (-1 : ℝ) ^ m • EulerMaclaurin.saw (m + 2) 0 •
    (iteratedDerivWithin (m + 1) (fun x ↦ (𝐞 (T * F (x / N)) : ℂ)) interval b -
      iteratedDerivWithin (m + 1) (fun x ↦ (𝐞 (T * F (x / N)) : ℂ)) interval a)
  let R : ℂ := (-1 : ℝ) ^ s •
    (∫ x in (a : ℝ)..b, EulerMaclaurin.saw (s + 1) x •
      iteratedDerivWithin (s + 1) (fun y ↦ (𝐞 (T * F (y / N)) : ℂ)) interval x)
  have hEM' : (∑ n ∈ Finset.Icc a b, (𝐞 (T * F ((n : ℝ) / N)) : ℂ)) =
      I + E + Q + R := by
    simpa only [I, E, Q, R, interval] using hEM
  rw [hEM']
  change ‖(I + E + Q + R) - I‖ ≤ oscillatoryErrorConstant s
  rw [show (I + E + Q + R) - I = E + Q + R by abel]
  have hE : ‖E‖ ≤ 1 := by simpa only [E] using hendpoint
  have hQ : ‖Q‖ ≤ Ccorr := by simpa only [Q] using hcorr
  have hremSmul :
      ‖(-1 : ℝ) ^ s • (∫ x in (a : ℝ)..b, EulerMaclaurin.saw (s + 1) x •
        iteratedDerivWithin (s + 1) (fun y ↦ (𝐞 (T * F (y / N)) : ℂ)) interval x)‖ ≤
        Crem := by
    rw [norm_smul, Real.norm_eq_abs, abs_pow, abs_neg, abs_one, one_pow, one_mul]
    exact hrem
  have hR : ‖R‖ ≤ Crem := by simpa only [R] using hremSmul
  calc
    ‖E + Q + R‖ ≤ ‖E + Q‖ + ‖R‖ := norm_add_le _ _
    _ ≤ ‖E‖ + ‖Q‖ + ‖R‖ := by
      gcongr
      exact norm_add_le _ _
    _ ≤ 1 + Ccorr + Crem := add_le_add (add_le_add hE hQ) hR
    _ = oscillatoryErrorConstant s := by rfl

private theorem norm_log_phase_integral_ge
    {N T : ℝ} (hN : 0 < N) (hT : 1 ≤ T) :
    (1 / (1 + 2 * Real.pi)) * (N / T) ≤
      ‖∫ x in N..2 * N, (𝐞 (T * Real.log (x / N)) : ℂ)‖ := by
  let τ : ℂ := (2 * Real.pi : ℂ) * Complex.I
  let v : ℝ → ℂ := fun x ↦ 𝐞 (T * Real.log (x / N))
  let denom : ℂ := 1 + τ * T
  let G : ℝ → ℂ := fun x ↦ x * v x / denom
  have hdenom : denom ≠ 0 := by
    intro h
    have hre := congrArg Complex.re h
    simp [denom, τ] at hre
  have hv (x : ℝ) (hx : 0 < x) :
      HasDerivAt v (τ * T / x * v x) x := by
    have hdiv : HasDerivAt (fun y : ℝ ↦ y / N) (1 / N) x := by
      simpa using (hasDerivAt_id x).div_const N
    have hlog := (Real.hasDerivAt_log (div_pos hx hN).ne').comp x hdiv
    have hinner : HasDerivAt (fun y ↦ T * Real.log (y / N)) (T / x) x := by
      convert hlog.const_mul T using 1 <;>
        first | rfl | field_simp [hN.ne', hx.ne']
    have hchar' := (Real.hasDerivAt_fourierChar (T * Real.log (x / N))).scomp x hinner
    dsimp only [v]
    change HasDerivAt ((𝐞 · : ℝ → ℂ) ∘ fun y ↦ T * Real.log (y / N))
      (τ * T / x * (𝐞 (T * Real.log (x / N)) : ℂ)) x
    rw [show τ * T / x * (𝐞 (T * Real.log (x / N)) : ℂ) =
        (T / x) • ((2 * Real.pi : ℂ) * Complex.I * 𝐞 (T * Real.log (x / N))) by
      dsimp [τ]
      simp only [Complex.ofReal_div]
      ring]
    exact hchar'
  have hG (x : ℝ) (hx : 0 < x) : HasDerivAt G (v x) x := by
    have hprod := (hasDerivAt_id x).ofReal_comp.mul (hv x hx)
    have hquot := hprod.div_const denom
    have hderivEq :
        ((1 : ℂ) * v x + x * (τ * T / x * v x)) / denom = v x := by
      apply (div_eq_iff hdenom).2
      dsimp [denom, τ]
      field_simp [hx.ne']
    rw [← hderivEq]
    simpa only [G, Pi.mul_apply, id_eq, Complex.ofReal_one] using hquot
  have hvCont : ContinuousOn v (Set.uIcc N (2 * N)) := by
    intro x hx
    have hxIcc : x ∈ Set.Icc N (2 * N) := by
      simpa [Set.uIcc_of_le (by linarith : N ≤ 2 * N)] using hx
    exact (hv x (hN.trans_le hxIcc.1)).continuousAt.continuousWithinAt
  have hFTC := intervalIntegral.integral_eq_sub_of_hasDerivAt
    (a := N) (b := 2 * N) (f := G) (f' := v)
    (fun x hx ↦ by
      have hxIcc : x ∈ Set.Icc N (2 * N) := by
        simpa [Set.uIcc_of_le (by linarith : N ≤ 2 * N)] using hx
      exact hG x (hN.trans_le hxIcc.1)) hvCont.intervalIntegrable
  have hformula : (∫ x in N..2 * N, v x) =
      (N : ℂ) * (2 * 𝐞 (T * Real.log 2) - 1) / denom := by
    rw [hFTC]
    dsimp [G, v]
    rw [show (2 * N) / N = 2 by field_simp [hN.ne'],
      show N / N = 1 by field_simp [hN.ne'], Real.log_one]
    simp only [mul_zero]
    rw [show (𝐞 (0 : ℝ) : ℂ) = 1 by norm_num]
    push_cast
    field_simp [hdenom]
  rw [show (∫ x in N..2 * N, (𝐞 (T * Real.log (x / N)) : ℂ)) =
      ∫ x in N..2 * N, v x by rfl, hformula, norm_div, norm_mul,
      Complex.norm_real, Real.norm_eq_abs, abs_of_pos hN]
  have hnum : 1 ≤ ‖(2 : ℂ) * 𝐞 (T * Real.log 2) - 1‖ := by
    calc
      1 = ‖(2 : ℂ) * 𝐞 (T * Real.log 2)‖ - ‖(1 : ℂ)‖ := by norm_num
      _ ≤ ‖(2 : ℂ) * 𝐞 (T * Real.log 2) - 1‖ := norm_sub_norm_le _ _
  have hden : ‖denom‖ ≤ (1 + 2 * Real.pi) * T := by
    calc
      ‖denom‖ ≤ ‖(1 : ℂ)‖ + ‖τ * T‖ := by
        dsimp [denom]
        exact norm_add_le _ _
      _ = 1 + 2 * Real.pi * T := by
        dsimp [τ]
        simp [Complex.norm_real, Real.norm_eq_abs, Real.pi_pos.le,
          abs_of_nonneg (le_trans zero_le_one hT)]
      _ ≤ (1 + 2 * Real.pi) * T := by
        nlinarith [Real.pi_pos]
  have hdenpos : 0 < ‖denom‖ := norm_pos_iff.mpr hdenom
  have hcoefpos : 0 < 1 + 2 * Real.pi := by positivity
  apply (le_div_iff₀ hdenpos).2
  calc
    (1 / (1 + 2 * Real.pi)) * (N / T) * ‖denom‖ ≤
        (1 / (1 + 2 * Real.pi)) * (N / T) * ((1 + 2 * Real.pi) * T) := by
      gcongr
    _ = N := by field_simp [hcoefpos.ne', (zero_lt_one.trans_le hT).ne']
    _ ≤ N * ‖(2 : ℂ) * 𝐞 (T * Real.log 2) - 1‖ := by
      simpa only [mul_one] using mul_le_mul_of_nonneg_left hnum hN.le

private theorem exponent_le_of_isPowerBounded_of_eventually_norm_ge_rpow
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
  have hδ : β + ε + ε = γ := by dsimp [ε]; ring
  have hO := (isPowerBounded_iff_forall_pos X T β hT hTunbounded).1 hbound ε hε
  obtain ⟨C, hC⟩ := hO.bound
  have hTtop : Tendsto T atTop atTop := by
    rw [VariableObject.IsUnbounded] at hTunbounded
    exact hTunbounded.congr' <| Filter.Eventually.of_forall fun i ↦ by
      rw [Real.norm_eq_abs, abs_of_nonneg (le_trans zero_le_one (hT i))]
  have hpowtop : Tendsto (fun i ↦ T i ^ ε) atTop atTop :=
    (tendsto_rpow_atTop hε).comp hTtop
  have hlarge : ∀ᶠ i in atTop, C < c * T i ^ ε := by
    have h := (tendsto_atTop.1 hpowtop (C / c + 1))
    filter_upwards [h] with i hi
    have hc' : C < c * (C / c + 1) := by
      field_simp
      linarith
    exact lt_of_lt_of_le hc' (mul_le_mul_of_nonneg_left hi hc.le)
  obtain ⟨i, hiLower, hiUpper, hiLarge⟩ :=
    (hlower.and (hC.and hlarge)).exists
  have hTi : 0 < T i := lt_of_lt_of_le zero_lt_one (hT i)
  have hApos : 0 < T i ^ (β + ε) := Real.rpow_pos_of_pos hTi _
  have hpower : T i ^ γ = T i ^ (β + ε) * T i ^ ε := by
    rw [← Real.rpow_add hTi, hδ]
  have hsmall : c * T i ^ ε ≤ C := by
    apply (mul_le_mul_iff_right₀ hApos).mp
    calc
      T i ^ (β + ε) * (c * T i ^ ε) = c * T i ^ γ := by rw [hpower]; ring
      _ ≤ ‖X i‖ := hiLower
      _ ≤ C * ‖T i ^ (β + ε)‖ := hiUpper
      _ = T i ^ (β + ε) * C := by
        rw [Real.norm_eq_abs, abs_of_nonneg (Real.rpow_nonneg hTi.le _)]
        ring
  exact (not_lt_of_ge hsmall) hiLarge

private theorem isPowerAsymptotic_of_logb_tendsto
    {N T : VariableObject ℝ} {α : ℝ}
    (hT : ∀ i, 1 < T i) (hN : ∀ i, 0 < N i)
    (hexponent : Tendsto (fun i ↦ Real.logb (T i) (N i)) atTop (nhds α)) :
    IsPowerAsymptotic N T α := by
  let exponent : VariableObject ℝ := fun i ↦ Real.logb (T i) (N i)
  refine ⟨exponent, ?_, Filter.Eventually.of_forall fun i ↦ ?_⟩
  · rw [IsEqUpToInfinitesimal, VariableObject.IsInfinitesimal]
    have hsub : Tendsto (fun i ↦ Real.logb (T i) (N i) - α) atTop (nhds 0) := by
      simpa using hexponent.sub
        (tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ α) atTop (nhds α))
    convert (continuous_norm.tendsto 0).comp hsub using 1 <;>
      simp [Function.comp_def, exponent, VariableObject.fixed]
  · exact (Real.rpow_logb (zero_lt_one.trans (hT i)) (hT i).ne' (hN i)).symm

private theorem tendsto_logb_of_between_const_rpow
    {N T : VariableObject ℝ} {α D : ℝ}
    (hα : 0 < α) (hD : 1 ≤ D)
    (hNtop : Tendsto N atTop atTop) (hN : ∀ i, 2 ≤ N i)
    (hT : ∀ i, D * N i ^ α⁻¹ ≤ T i ∧ T i ≤ 2 * (D * N i ^ α⁻¹)) :
    Tendsto (fun i ↦ Real.logb (T i) (N i)) atTop (nhds α) := by
  have hlogNtop : Tendsto (fun i ↦ Real.log (N i)) atTop atTop :=
    Real.tendsto_log_atTop.comp hNtop
  have hlogNpos (i : ℕ) : 0 < Real.log (N i) := Real.log_pos (lt_of_lt_of_le one_lt_two (hN i))
  have hDpos : 0 < D := zero_lt_one.trans_le hD
  have hratio : Tendsto (fun i ↦ Real.log (T i) / Real.log (N i)) atTop (nhds α⁻¹) := by
    have hlowerLimit : Tendsto
        (fun i ↦ Real.log D / Real.log (N i) + α⁻¹) atTop (nhds α⁻¹) := by
      simpa using (tendsto_const_nhds.div_atTop hlogNtop).add
        (tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ α⁻¹) atTop (nhds α⁻¹))
    have hupperLimit : Tendsto
        (fun i ↦ Real.log (2 * D) / Real.log (N i) + α⁻¹) atTop (nhds α⁻¹) := by
      simpa using (tendsto_const_nhds.div_atTop hlogNtop).add
        (tendsto_const_nhds : Tendsto (fun _ : ℕ ↦ α⁻¹) atTop (nhds α⁻¹))
    apply tendsto_of_tendsto_of_tendsto_of_le_of_le' hlowerLimit hupperLimit
    · exact Filter.Eventually.of_forall fun i ↦ by
        have hbase : 0 < N i := zero_lt_two.trans_le (hN i)
        have hlower := Real.log_le_log (mul_pos hDpos (Real.rpow_pos_of_pos hbase _)) (hT i).1
        rw [Real.log_mul hDpos.ne' (Real.rpow_pos_of_pos hbase _).ne',
          Real.log_rpow hbase] at hlower
        apply (le_div_iff₀ (hlogNpos i)).2
        rw [add_mul, div_mul_cancel₀ _ (hlogNpos i).ne']
        simpa [mul_comm] using hlower
    · exact Filter.Eventually.of_forall fun i ↦ by
        have hbase : 0 < N i := zero_lt_two.trans_le (hN i)
        have htwoD : 0 < 2 * D := mul_pos (by norm_num) hDpos
        have hupper := Real.log_le_log (lt_of_lt_of_le
          (mul_pos hDpos (Real.rpow_pos_of_pos hbase _)) (hT i).1) (hT i).2
        rw [show 2 * (D * N i ^ α⁻¹) = (2 * D) * N i ^ α⁻¹ by ring,
          Real.log_mul htwoD.ne' (Real.rpow_pos_of_pos hbase _).ne',
          Real.log_rpow hbase] at hupper
        apply (div_le_iff₀ (hlogNpos i)).2
        rw [add_mul, div_mul_cancel₀ _ (hlogNpos i).ne']
        simpa [mul_comm] using hupper
  have hinv := hratio.inv₀ (inv_ne_zero hα.ne')
  have hTgt (i : ℕ) : 1 < T i := by
    have hbase : 1 < N i := one_lt_two.trans_le (hN i)
    have hrpow : 1 < N i ^ α⁻¹ := Real.one_lt_rpow hbase (inv_pos.mpr hα)
    exact hrpow.trans_le <| (le_mul_of_one_le_left
      (Real.rpow_nonneg (zero_le_one.trans hbase.le) _) hD).trans (hT i).1
  have hinv' : Tendsto (fun i ↦ Real.logb (T i) (N i)) atTop (nhds α⁻¹⁻¹) :=
    hinv.congr' <| Filter.Eventually.of_forall fun i ↦ by
      dsimp [Real.logb]
      field_simp [(hlogNpos i).ne', (Real.log_pos (hTgt i)).ne']
  simpa using hinv'

private theorem sub_one_le_exponentSumGrowthExponent
    {α : ℝ≥0} (hα : 1 < α) :
    (α : ℝ) - 1 ≤ exponentSumGrowthExponent α := by
  have hαR : (1 : ℝ) < α := by exact_mod_cast hα
  let δ : ℝ := (α : ℝ) - 1
  have hδ : 0 < δ := by dsimp [δ]; linarith
  obtain ⟨s, hs⟩ := exists_nat_ge ((1 + δ / 2) / δ)
  have hsone : 1 ≤ s := by
    by_contra hs0
    have : s = 0 := by omega
    subst s
    norm_num at hs
    have : 0 < (1 + δ / 2) / δ := by positivity
    linarith
  have hsδ : 1 + δ / 2 ≤ (s : ℝ) * δ :=
    (div_le_iff₀ hδ).1 (by simpa [mul_comm] using hs)
  let q : ℝ := δ / (2 * (s + 1 : ℕ))
  have hq : 0 < q := by dsimp [q]; positivity
  have hqδ : q ≤ δ := by
    dsimp [q]
    push_cast
    apply (div_le_iff₀ (by positivity : 0 < (2 : ℝ) * ((s : ℝ) + 1))).2
    have hs0 : (0 : ℝ) ≤ s := Nat.cast_nonneg s
    calc
      δ = δ * 1 := by ring
      _ ≤ δ * (2 * ((s : ℝ) + 1)) := by
        gcongr
        nlinarith
  let N : VariableObject ℝ := fun i ↦ (i : ℝ) + 1
  let T : VariableObject ℝ := fun i ↦ N i ^ ((α : ℝ)⁻¹)
  let a : VariableObject ℕ := fun i ↦ i + 1
  let b : VariableObject ℕ := fun i ↦ 2 * (i + 1)
  have hN : ∀ i, 1 ≤ N i := by
    intro i
    dsimp [N]
    exact le_add_of_nonneg_left (Nat.cast_nonneg i)
  have hT : ∀ i, 1 ≤ T i := by
    intro i
    dsimp [T]
    exact Real.one_le_rpow (hN i) (inv_nonneg.mpr (zero_le_one.trans hαR.le))
  have hNtop : Tendsto N atTop atTop := by
    exact tendsto_atTop_add_const_right atTop 1 tendsto_natCast_atTop_atTop
  have hTtop : Tendsto T atTop atTop :=
    (tendsto_rpow_atTop (inv_pos.mpr (zero_lt_one.trans hαR))).comp hNtop
  have hTunbounded : T.IsUnbounded := by
    rw [VariableObject.IsUnbounded]
    convert hTtop using 1
    ext i
    rw [Real.norm_eq_abs, abs_of_nonneg (le_trans zero_le_one (hT i))]
  have hNTexact (i : ℕ) : N i = T i ^ (α : ℝ) := by
    dsimp [T]
    rw [← Real.rpow_mul (le_trans zero_le_one (hN i)),
      inv_mul_cancel₀ (zero_lt_one.trans hαR).ne',
      Real.rpow_one]
  have hNT : IsPowerAsymptotic N T (α : ℝ) := by
    refine ⟨VariableObject.fixed (α : ℝ), IsEqUpToInfinitesimal.refl _, ?_⟩
    exact Filter.Eventually.of_forall hNTexact
  have hab : ∀ i, N i ≤ (a i : ℝ) ∧ (b i : ℝ) ≤ 2 * N i := by
    intro i
    dsimp [N, a, b]
    push_cast
    exact ⟨le_rfl, le_rfl⟩
  have hbound := isPowerBounded_logPhase α hN hT hTunbounded hNT hab
  obtain ⟨c, K, hc, hK, hderivEventually⟩ :=
    eventually_model_phase_deriv_bounds isModelPhaseFunction_log s
  have hKpower : ∀ᶠ i in atTop, K ≤ T i ^ q :=
    tendsto_atTop.1 ((tendsto_rpow_atTop hq).comp hTtop) K
  let d : ℝ := 1 / (1 + 2 * Real.pi)
  have hd : 0 < d := by dsimp [d]; positivity
  let E : ℝ := oscillatoryErrorConstant s
  have hE : 0 ≤ E := by
    dsimp [E, oscillatoryErrorConstant]
    have hsum : 0 ≤ ∑ m ∈ Finset.range s, |EulerMaclaurin.saw (m + 2) 0| *
        (2 * (((m + 1).factorial : ℝ) * (2 * Real.pi + 1) ^ (m + 1))) :=
      Finset.sum_nonneg fun m _ ↦ by positivity
    have hrem : 0 ≤ EulerMaclaurin.sawBound (s + 1) *
        (((s + 1).factorial : ℝ) * (2 * Real.pi + 1) ^ (s + 1)) :=
      mul_nonneg EulerMaclaurin.sawBound_nonneg (by positivity)
    linarith
  have herrorSmall : ∀ᶠ i in atTop, 2 * E / d ≤ T i ^ ((α : ℝ) - 1) :=
    tendsto_atTop.1 ((tendsto_rpow_atTop (by linarith : 0 < (α : ℝ) - 1)).comp hTtop)
      (2 * E / d)
  have hlower : ∀ᶠ i in atTop,
      (d / 2) * T i ^ ((α : ℝ) - 1) ≤
        ‖exponentialSum logPhase T N a b i‖ := by
    filter_upwards [hderivEventually, hKpower, herrorSmall] with i hiDeriv hiK hiSmall
    have hTi : 0 < T i := zero_lt_one.trans_le (hT i)
    have hNi : 0 < N i := zero_lt_one.trans_le (hN i)
    have hqmul : q * (s + 1 : ℕ) = δ / 2 := by
      dsimp [q]
      field_simp
    have hKpow : K ^ (s + 1) ≤ T i ^ (δ / 2) := by
      calc
        K ^ (s + 1) ≤ (T i ^ q) ^ (s + 1) :=
          pow_le_pow_left₀ (zero_le_one.trans hK) hiK _
        _ = T i ^ (q * (s + 1 : ℕ)) :=
          (Real.rpow_mul_natCast hTi.le q (s + 1)).symm
        _ = T i ^ (δ / 2) := by rw [hqmul]
    have hratio : K * T i / N i ≤ 1 := by
      apply (div_le_one hNi).2
      rw [hNTexact]
      calc
        K * T i ≤ T i ^ q * T i := mul_le_mul_of_nonneg_right hiK hTi.le
        _ = T i ^ (q + 1) := by
          simpa only [Real.rpow_one] using (Real.rpow_add hTi q 1).symm
        _ ≤ T i ^ (δ + 1) :=
          Real.rpow_le_rpow_of_exponent_le (hT i) (by linarith)
        _ = T i ^ (α : ℝ) := by dsimp [δ]; congr 1; ring
    have hremainder : N i * (K * T i / N i) ^ (s + 1) ≤ 1 := by
      have hexponent : δ / 2 + (s + 1 : ℕ) ≤ (α : ℝ) * s := by
        dsimp [δ] at hsδ ⊢
        push_cast
        linarith
      have hnum : K ^ (s + 1) * T i ^ (s + 1) ≤ N i ^ s := by
        rw [hNTexact]
        calc
          K ^ (s + 1) * T i ^ (s + 1) ≤
              T i ^ (δ / 2) * T i ^ (s + 1 : ℕ) := by gcongr
          _ = T i ^ (δ / 2 + (s + 1 : ℕ)) := by
            rw [← Real.rpow_natCast]
            exact (Real.rpow_add hTi (δ / 2) (s + 1 : ℕ)).symm
          _ ≤ T i ^ ((α : ℝ) * s) :=
            Real.rpow_le_rpow_of_exponent_le (hT i) hexponent
          _ = (T i ^ (α : ℝ)) ^ s := Real.rpow_mul_natCast hTi.le (α : ℝ) s
      rw [show N i * (K * T i / N i) ^ (s + 1) =
          (K ^ (s + 1) * T i ^ (s + 1)) / N i ^ s by
        rw [div_pow, mul_pow, pow_succ]
        field_simp [hNi.ne', pow_ne_zero _ hNi.ne']
        ring]
      exact (div_le_one (pow_pos hNi _)).2 hnum
    have hablt : a i < b i := by dsimp [a, b]; omega
    have herr := norm_oscillatory_sum_sub_integral_le s hablt
      (hN i) (hT i) hK (isModelPhaseFunction_log.1 i)
      (fun k hk hkS u hu ↦ (hiDeriv ⟨u, hu⟩).2 k hk hkS)
      (hab i).1 (hab i).2 hratio hremainder
    have hInt := norm_log_phase_integral_ge hNi (hT i)
    have hpower : N i / T i = T i ^ ((α : ℝ) - 1) := by
      rw [hNTexact]
      simpa only [Real.rpow_one] using (Real.rpow_sub hTi (α : ℝ) 1).symm
    rw [hpower] at hInt
    have htriangle :
        ‖∫ x in (a i : ℝ)..b i,
          (𝐞 (T i * Real.log (x / N i)) : ℂ)‖ ≤
          ‖exponentialSum logPhase T N a b i‖ + E := by
      calc
        _ ≤ ‖exponentialSum logPhase T N a b i‖ +
            ‖exponentialSum logPhase T N a b i -
              (∫ x in (a i : ℝ)..b i,
                (𝐞 (T i * Real.log (x / N i)) : ℂ))‖ :=
          norm_le_norm_add_norm_sub _ _
        _ ≤ ‖exponentialSum logPhase T N a b i‖ + E := by
          gcongr
          simpa only [E, exponentialSum, logPhase] using herr
    have hIntEndpoints :
        (∫ x in (a i : ℝ)..b i,
          (𝐞 (T i * Real.log (x / N i)) : ℂ)) =
        ∫ x in N i..2 * N i,
          (𝐞 (T i * Real.log (x / N i)) : ℂ) := by
      congr 1 <;> dsimp [N, a, b] <;> push_cast <;> ring
    rw [hIntEndpoints] at htriangle
    have hEsmall : E ≤ d / 2 * T i ^ ((α : ℝ) - 1) := by
      calc
        E = d / 2 * (2 * E / d) := by field_simp [hd.ne']
        _ ≤ d / 2 * T i ^ ((α : ℝ) - 1) :=
          mul_le_mul_of_nonneg_left hiSmall (by positivity)
    linarith
  exact exponent_le_of_isPowerBounded_of_eventually_norm_ge_rpow
    hT hTunbounded hbound (half_pos hd) (by simpa [d] using hlower)

private lemma log_div_separated
    {N : ℝ} {a b : ℕ} (hN : 0 < N) (ha : N ≤ a) (hb : (b : ℝ) ≤ 2 * N) :
    IsSeparatedFamily (1 / (2 * N))
      (fun n : ↥(Finset.Icc a b) ↦ Real.log ((n : ℝ) / N)) := by
  intro m n hmn
  rw [Real.dist_eq]
  rcases lt_or_gt_of_ne (Subtype.coe_ne_coe.mpr hmn) with hmn' | hnm'
  · have hmnR : (m : ℝ) ≤ (n : ℝ) := by exact_mod_cast hmn'.le
    have hm_mem := Finset.mem_Icc.mp m.property
    have hmpos : 0 < (m : ℝ) := lt_of_lt_of_le hN (ha.trans (by exact_mod_cast hm_mem.1))
    have hnpos : 0 < (n : ℝ) := lt_trans hmpos (by exact_mod_cast hmn')
    rw [abs_of_nonpos (sub_nonpos.mpr (Real.log_le_log (by positivity)
      (div_le_div_of_nonneg_right hmnR hN.le))), neg_sub]
    have hlog : 1 - (m : ℝ) / n ≤ Real.log ((n : ℝ) / m) := by
      simpa [inv_div] using Real.one_sub_inv_le_log_of_pos (div_pos hnpos hmpos)
    have hdiff : 1 / (2 * N) ≤ 1 - (m : ℝ) / n := by
      have hgapNat : (m : ℕ) + 1 ≤ (n : ℕ) := Nat.add_one_le_iff.mpr hmn'
      have hgapCast : ((m : ℕ) : ℝ) + 1 ≤ ((n : ℕ) : ℝ) := by exact_mod_cast hgapNat
      have hgap : (1 : ℝ) ≤ (n : ℝ) - m := by linarith
      have hn_mem := Finset.mem_Icc.mp n.property
      have hnle : (n : ℝ) ≤ 2 * N := le_trans (by exact_mod_cast hn_mem.2) hb
      rw [one_sub_div hnpos.ne']
      calc
        1 / (2 * N) ≤ 1 / (n : ℝ) := one_div_le_one_div_of_le hnpos hnle
        _ ≤ ((n : ℝ) - m) / n := (div_le_div_iff_of_pos_right hnpos).2 hgap
    calc
      1 / (2 * N) ≤ Real.log ((n : ℝ) / m) := hdiff.trans hlog
      _ = Real.log ((n : ℝ) / N) - Real.log ((m : ℝ) / N) := by
        rw [Real.log_div (ne_of_gt hnpos) hN.ne', Real.log_div (ne_of_gt hmpos) hN.ne',
          Real.log_div (ne_of_gt hnpos) (ne_of_gt hmpos)]
        ring
  · have hnmR : (n : ℝ) ≤ (m : ℝ) := by exact_mod_cast hnm'.le
    have hn_mem := Finset.mem_Icc.mp n.property
    have hnpos : 0 < (n : ℝ) := lt_of_lt_of_le hN (ha.trans (by exact_mod_cast hn_mem.1))
    have hmpos : 0 < (m : ℝ) := lt_trans hnpos (by exact_mod_cast hnm')
    rw [abs_of_nonneg (sub_nonneg.mpr (Real.log_le_log (by positivity)
      (div_le_div_of_nonneg_right hnmR hN.le)))]
    have hlog : 1 - (n : ℝ) / m ≤ Real.log ((m : ℝ) / n) := by
      simpa [inv_div] using Real.one_sub_inv_le_log_of_pos (div_pos hmpos hnpos)
    have hdiff : 1 / (2 * N) ≤ 1 - (n : ℝ) / m := by
      have hgapNat : (n : ℕ) + 1 ≤ (m : ℕ) := Nat.add_one_le_iff.mpr hnm'
      have hgapCast : ((n : ℕ) : ℝ) + 1 ≤ ((m : ℕ) : ℝ) := by exact_mod_cast hgapNat
      have hgap : (1 : ℝ) ≤ (m : ℝ) - n := by linarith
      have hm_mem := Finset.mem_Icc.mp m.property
      have hmle : (m : ℝ) ≤ 2 * N := le_trans (by exact_mod_cast hm_mem.2) hb
      rw [one_sub_div hmpos.ne']
      calc
        1 / (2 * N) ≤ 1 / (m : ℝ) := one_div_le_one_div_of_le hmpos hmle
        _ ≤ ((m : ℝ) - n) / m := (div_le_div_iff_of_pos_right hmpos).2 hgap
    calc
      1 / (2 * N) ≤ Real.log ((m : ℝ) / n) := hdiff.trans hlog
      _ = Real.log ((m : ℝ) / N) - Real.log ((n : ℝ) / N) := by
        rw [Real.log_div (ne_of_gt hmpos) hN.ne', Real.log_div (ne_of_gt hnpos) hN.ne',
          Real.log_div (ne_of_gt hmpos) (ne_of_gt hnpos)]
        ring

private theorem exists_logPhase_sum_norm_sq_ge :
    ∃ C : ℝ, 0 < C ∧
      ∀ {R N : ℝ} {a b : ℕ},
        0 < R → 0 < N → N ≤ a → (b : ℝ) ≤ 2 * N → N ≤ R / (4 * C) →
        ∃ t ∈ Set.Icc R (2 * R),
          ((Finset.Icc a b).card : ℝ) / 2 ≤
            ‖∑ n ∈ Finset.Icc a b,
              (𝐞 (t * Real.log ((n : ℝ) / N)) : ℂ)‖ ^ 2 := by
  obtain ⟨C₀, hC₀, hestimate⟩ := l2_integral_estimate
  refine ⟨C₀, hC₀, ?_⟩
  intro R N a b hR hN ha hb hNR
  let ξ : ↥(Finset.Icc a b) → ℝ := fun n ↦ Real.log ((n : ℝ) / N)
  let coeff : ↥(Finset.Icc a b) → ℂ := fun _ ↦ 1
  have hsep : IsSeparatedFamily (1 / (2 * N)) ξ :=
    log_div_separated hN ha hb
  obtain ⟨θ, hθ, hint⟩ := hestimate coeff ξ (2 * N) (by positivity) hsep
    R (2 * R) R (by ring) (by linarith)
  let g : ℝ → ℝ := fun t ↦ ‖∑ n, coeff n * 𝐞 (ξ n * t)‖ ^ 2
  have hg : Continuous g := by
    dsimp [g, coeff, ξ]
    fun_prop
  have hcard : ∑ n, ‖coeff n‖ ^ 2 = ((Finset.Icc a b).card : ℝ) := by
    simp [coeff]
  have hfactor : R / 2 ≤ R + θ * (2 * N) := by
    have hθlower : -C₀ ≤ θ := (neg_le_of_abs_le hθ)
    have hCN : 2 * C₀ * N ≤ R / 2 := by
      apply (le_div_iff₀ (by positivity : 0 < (2 : ℝ))).2
      calc
        2 * C₀ * N * 2 = N * (4 * C₀) := by ring
        _ ≤ R := (le_div_iff₀ (by positivity : 0 < 4 * C₀)).1 hNR
    nlinarith
  have hint_lower : R / 2 * ((Finset.Icc a b).card : ℝ) ≤
      ∫ t in Set.Icc R (2 * R), g t := by
    rw [show (∫ t in Set.Icc R (2 * R), g t) =
        (R + θ * (2 * N)) * ((Finset.Icc a b).card : ℝ) by
      simpa only [g, hcard] using hint]
    exact mul_le_mul_of_nonneg_right hfactor (Nat.cast_nonneg _)
  obtain ⟨t, ht, htavg⟩ := exists_eq_interval_average
    (by linarith : R ≠ 2 * R) hg.continuousOn
  have ht' : t ∈ Set.Ioo R (2 * R) := by
    simpa [Set.uIoo_of_le (by linarith : R ≤ 2 * R)] using ht
  refine ⟨t, ⟨ht'.1.le, ht'.2.le⟩, ?_⟩
  have havg : ((Finset.Icc a b).card : ℝ) / 2 ≤ g t := by
    rw [htavg]
    rw [interval_average_eq_div]
    have hset_interval : (∫ x in R..2 * R, g x) =
        ∫ x in Set.Icc R (2 * R), g x := by
      rw [intervalIntegral.integral_of_le (by linarith),
        MeasureTheory.integral_Icc_eq_integral_Ioc]
    rw [hset_interval]
    rw [show 2 * R - R = R by ring]
    apply (le_div_iff₀ hR).2
    nlinarith [hint_lower]
  change ((Finset.Icc a b).card : ℝ) / 2 ≤
    ‖∑ n ∈ Finset.Icc a b, (𝐞 (t * Real.log ((n : ℝ) / N)) : ℂ)‖ ^ 2
  rw [Finset.sum_subtype (Finset.Icc a b) (fun _ ↦ Iff.rfl)]
  simpa only [g, coeff, ξ, one_mul, mul_comm] using havg

private theorem half_le_exponentSumGrowthExponent
    {α : ℝ≥0} (hα : 0 < α) (hαone : α ≤ 1) :
    (α : ℝ) / 2 ≤ exponentSumGrowthExponent α := by
  have hαR : 0 < (α : ℝ) := by exact_mod_cast hα
  have hαoneR : (α : ℝ) ≤ 1 := by exact_mod_cast hαone
  obtain ⟨C, hC, hlarge⟩ := exists_logPhase_sum_norm_sq_ge
  let D : ℝ := max (4 * C) 1
  have hD : 1 ≤ D := le_max_right _ _
  have hCD : 4 * C ≤ D := le_max_left _ _
  have hDpos : 0 < D := zero_lt_one.trans_le hD
  let n : VariableObject ℕ := fun i ↦ i + 2
  let N : VariableObject ℝ := fun i ↦ n i
  let a : VariableObject ℕ := n
  let b : VariableObject ℕ := fun i ↦ 2 * n i
  let R : VariableObject ℝ := fun i ↦ D * N i ^ (α : ℝ)⁻¹
  have hN (i : ℕ) : 2 ≤ N i := by
    dsimp [N, n]
    exact_mod_cast Nat.le_add_left 2 i
  have hNpos (i : ℕ) : 0 < N i := zero_lt_two.trans_le (hN i)
  have hNtop : Tendsto N atTop atTop := by
    dsimp [N, n]
    simpa only [Nat.cast_add, Nat.cast_ofNat] using
      tendsto_atTop_add_const_right atTop (2 : ℝ) tendsto_natCast_atTop_atTop
  have hinvα : 1 ≤ (α : ℝ)⁻¹ := (one_le_inv₀ hαR).2 hαoneR
  have hNR (i : ℕ) : N i ≤ R i / (4 * C) := by
    have hpow : N i ≤ N i ^ (α : ℝ)⁻¹ := by
      simpa only [Real.rpow_one] using
        Real.rpow_le_rpow_of_exponent_le (one_le_two.trans (hN i)) hinvα
    apply (le_div_iff₀ (by positivity : 0 < 4 * C)).2
    dsimp [R]
    calc
      N i * (4 * C) ≤ N i * D := mul_le_mul_of_nonneg_left hCD (hNpos i).le
      _ = D * N i := by ring
      _ ≤ D * N i ^ (α : ℝ)⁻¹ := mul_le_mul_of_nonneg_left hpow hDpos.le
  have hgood (i : ℕ) : ∃ t ∈ Set.Icc (R i) (2 * R i),
      ((Finset.Icc (a i) (b i)).card : ℝ) / 2 ≤
        ‖∑ m ∈ Finset.Icc (a i) (b i),
          (𝐞 (t * Real.log ((m : ℝ) / N i)) : ℂ)‖ ^ 2 := by
    apply hlarge (mul_pos hDpos (Real.rpow_pos_of_pos (hNpos i) _)) (hNpos i)
    · dsimp [N, a, n]
      norm_num
    · dsimp [N, b, n]
      push_cast
      norm_num
    · exact hNR i
  let T : VariableObject ℝ := fun i ↦ (hgood i).choose
  have hTmem (i : ℕ) : T i ∈ Set.Icc (R i) (2 * R i) := (hgood i).choose_spec.1
  have hTsquare (i : ℕ) : ((Finset.Icc (a i) (b i)).card : ℝ) / 2 ≤
      ‖∑ m ∈ Finset.Icc (a i) (b i),
        (𝐞 (T i * Real.log ((m : ℝ) / N i)) : ℂ)‖ ^ 2 :=
    (hgood i).choose_spec.2
  have hTbounds (i : ℕ) : D * N i ^ (α : ℝ)⁻¹ ≤ T i ∧
      T i ≤ 2 * (D * N i ^ (α : ℝ)⁻¹) := by
    simpa only [R] using Set.mem_Icc.mp (hTmem i)
  have hTgt (i : ℕ) : 1 < T i := by
    have hpow : 1 < N i ^ (α : ℝ)⁻¹ :=
      Real.one_lt_rpow (one_lt_two.trans_le (hN i)) (inv_pos.mpr hαR)
    exact hpow.trans_le <| (le_mul_of_one_le_left
      (Real.rpow_nonneg (hNpos i).le _) hD).trans (hTbounds i).1
  have hT : ∀ i, 1 ≤ T i := fun i ↦ (hTgt i).le
  have hRtop : Tendsto R atTop atTop := by
    exact (tendsto_rpow_atTop (inv_pos.mpr hαR)).comp hNtop |>.const_mul_atTop hDpos
  have hTtop : Tendsto T atTop atTop :=
    tendsto_atTop_mono' atTop (Filter.Eventually.of_forall fun i ↦ (hTmem i).1) hRtop
  have hTunbounded : T.IsUnbounded := by
    rw [VariableObject.IsUnbounded]
    convert hTtop using 1
    ext i
    rw [Real.norm_eq_abs, abs_of_nonneg (zero_le_one.trans (hT i))]
  have hlogb : Tendsto (fun i ↦ Real.logb (T i) (N i)) atTop (nhds (α : ℝ)) :=
    tendsto_logb_of_between_const_rpow hαR hD hNtop hN hTbounds
  have hNT : IsPowerAsymptotic N T (α : ℝ) :=
    isPowerAsymptotic_of_logb_tendsto hTgt hNpos hlogb
  have hab : ∀ i, N i ≤ (a i : ℝ) ∧ (b i : ℝ) ≤ 2 * N i := by
    intro i
    dsimp [N, a, b, n]
    push_cast
    exact ⟨le_rfl, le_rfl⟩
  have hbound := isPowerBounded_logPhase α (fun i ↦ (hN i).trans' one_le_two)
    hT hTunbounded hNT hab
  let L : ℝ := (2 * D) ^ ((α : ℝ) / 2)
  have hL : 0 < L := Real.rpow_pos_of_pos (mul_pos (by norm_num) hDpos) _
  have hlower : ∀ᶠ i in atTop,
      (1 / (2 * L)) * T i ^ ((α : ℝ) / 2) ≤
        ‖exponentialSum logPhase T N a b i‖ :=
    Filter.Eventually.of_forall fun i ↦ by
      have hcard : N i ≤ ((Finset.Icc (a i) (b i)).card : ℝ) := by
        have hcardNat : n i ≤ (Finset.Icc (a i) (b i)).card := by
          rw [Nat.card_Icc]
          dsimp [a, b]
          omega
        change ((n i : ℕ) : ℝ) ≤ ((Finset.Icc (a i) (b i)).card : ℝ)
        exact_mod_cast hcardNat
      have hsumSquare : N i / 2 ≤ ‖exponentialSum logPhase T N a b i‖ ^ 2 := by
        calc
          N i / 2 ≤ ((Finset.Icc (a i) (b i)).card : ℝ) / 2 := by gcongr
          _ ≤ ‖exponentialSum logPhase T N a b i‖ ^ 2 := by
            simpa only [exponentialSum, logPhase] using hTsquare i
      have hrootSquare : (N i ^ (1 / 2 : ℝ)) ^ 2 = N i := by
        rw [← Real.rpow_natCast]
        rw [← Real.rpow_mul (hNpos i).le]
        norm_num
      have hroot : N i ^ (1 / 2 : ℝ) / 2 ≤
          ‖exponentialSum logPhase T N a b i‖ := by
        have hrootNonneg := Real.rpow_nonneg (hNpos i).le (1 / 2 : ℝ)
        have hnormNonneg := norm_nonneg (exponentialSum logPhase T N a b i)
        nlinarith
      have hTpower : T i ^ ((α : ℝ) / 2) ≤ L * N i ^ (1 / 2 : ℝ) := by
        calc
          T i ^ ((α : ℝ) / 2) ≤
              (2 * (D * N i ^ (α : ℝ)⁻¹)) ^ ((α : ℝ) / 2) :=
            Real.rpow_le_rpow (zero_le_one.trans (hT i)) (hTbounds i).2 (by positivity)
          _ = L * N i ^ (1 / 2 : ℝ) := by
            rw [show 2 * (D * N i ^ (α : ℝ)⁻¹) = (2 * D) * N i ^ (α : ℝ)⁻¹ by ring,
              Real.mul_rpow (by positivity) (Real.rpow_nonneg (hNpos i).le _),
              ← Real.rpow_mul (hNpos i).le]
            dsimp [L]
            congr 2
            field_simp [hαR.ne']
      calc
        (1 / (2 * L)) * T i ^ ((α : ℝ) / 2) ≤
            (1 / (2 * L)) * (L * N i ^ (1 / 2 : ℝ)) := by gcongr
        _ = N i ^ (1 / 2 : ℝ) / 2 := by field_simp [hL.ne']
        _ ≤ ‖exponentialSum logPhase T N a b i‖ := hroot
  exact exponent_le_of_isPowerBounded_of_eventually_norm_ge_rpow
    hT hTunbounded hbound (by positivity) hlower

private theorem exponentSumGrowthExponent_nonneg (α : ℝ≥0) :
    0 ≤ exponentSumGrowthExponent α :=
  (isExponentSumBound_exponentSumGrowthExponent α).nonneg

private theorem exponentSumGrowthExponent_le_self (α : ℝ≥0) :
    exponentSumGrowthExponent α ≤ (α : ℝ) :=
  exponentSumGrowthExponent_le_iff.mpr (isExponentSumBound_self α)

/-- Above the transition point, the exponential sum growth exponent is `α - 1`. -/
theorem exponentSumGrowthExponent_eq_sub_one
    {α : ℝ≥0} (hα : 1 < α) :
    exponentSumGrowthExponent α = (α : ℝ) - 1 := by
  exact le_antisymm
    (exponentSumGrowthExponent_le_iff.mpr (isExponentSumBound_sub_one hα))
    (sub_one_le_exponentSumGrowthExponent hα)

/-- Up to the transition point, the exponential sum growth exponent lies in
the interval `[α / 2, α]`. -/
theorem exponentSumGrowthExponent_mem_Icc
    {α : ℝ≥0} (hα : α ≤ 1) :
    exponentSumGrowthExponent α ∈ Set.Icc ((α : ℝ) / 2) (α : ℝ) := by
  refine ⟨?_, exponentSumGrowthExponent_le_self α⟩
  by_cases hαzero : α = 0
  · subst α
    simpa using exponentSumGrowthExponent_nonneg 0
  · have hαpos : 0 < α := lt_of_le_of_ne (zero_le : (0 : ℝ≥0) ≤ α) (Ne.symm hαzero)
    exact half_le_exponentSumGrowthExponent hαpos hα

/-- The exponential sum growth exponent vanishes at scale zero. -/
@[simp] theorem exponentSumGrowthExponent_zero :
    exponentSumGrowthExponent 0 = 0 := by
  apply le_antisymm
  · simpa using exponentSumGrowthExponent_le_self 0
  · exact exponentSumGrowthExponent_nonneg 0

end Expdb
