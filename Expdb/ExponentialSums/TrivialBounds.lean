module

public import Expdb.ExponentialSums.ExponentSumGrowth

import Expdb.ExponentialSums.ExponentSumGrowthNonAsymptotic
import Expdb.ExponentialSums.OscillatoryBounds
import Expdb.Fourier.L2Integral
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv

/-!
# Trivial bounds for exponential sum growth exponents

This module formalizes the trivial bounds on `β` in the blueprint's Exponential sum growth
exponents chapter (`beta-chapter`).

For `α > 1`, a uniform Euler–Maclaurin estimate compares the exponential sum
of an arbitrary approximate model phase with its oscillatory integral. The
derivative bounds for model phases then give the upper bound `β(α) ≤ α - 1`.
For the reverse inequality, it suffices to use the logarithmic model phase,
since an admissible exponent must bound every model phase.
Its sum is essentially `∑ n ^ (iT)` and its integral has size comparable to `N / T`,
while Euler–Maclaurin contributes only a bounded error, giving `β(α) ≥ α - 1`.

For `0 ≤ α ≤ 1`, the triangle inequality gives `β(α) ≤ α`. The `L^2` integral estimate
lemma is applied to the separated values `log(n/N)` to find a parameter for
which the logarithmic-phase sum has size comparable to `N ^ (1 / 2)`, yielding
the lower bound `β(α) ≥ α / 2`.
-/

@[expose] public section

open Filter Topology
open scoped ContDiff Expdb FourierTransform NNReal

noncomputable section

namespace Expdb

/-! ## Choosing the Euler–Maclaurin order -/

private lemma exists_oscillatory_scale_parameters
    {δ : ℝ} (hδ : 0 < δ) :
    ∃ s : ℕ, ∃ q : ℝ,
      1 ≤ s ∧ 0 < q ∧ q ≤ δ ∧
      q * (s + 1 : ℕ) = δ / 2 ∧
      1 + δ / 2 ≤ (s : ℝ) * δ := by
  obtain ⟨s, hs⟩ := exists_nat_ge ((1 + δ / 2) / δ)
  have hsone : 1 ≤ s := by
    by_contra hs0
    have : s = 0 := by omega
    subst s
    norm_num at hs
    have : 0 < (1 + δ / 2) / δ := by positivity
    linarith
  let q : ℝ := δ / (2 * (s + 1 : ℕ))
  have hq : 0 < q := by dsimp [q]; positivity
  have hqδ : q ≤ δ := by
    dsimp [q]
    push_cast
    apply (div_le_iff₀ (by positivity : 0 < (2 : ℝ) * ((s : ℝ) + 1))).2
    have hs0 : (0 : ℝ) ≤ s := Nat.cast_nonneg s
    nlinarith
  have hqmul : q * (s + 1 : ℕ) = δ / 2 := by
    dsimp [q]
    field_simp
  have hsδ : 1 + δ / 2 ≤ (s : ℝ) * δ :=
    (div_le_iff₀ hδ).1 (by simpa [mul_comm] using hs)
  exact ⟨s, q, hsone, hq, hqδ, hqmul, hsδ⟩

private lemma oscillatory_scale_bounds
    {s : ℕ} {δ q K T N : ℝ}
    (hT : 1 ≤ T) (hN : 1 ≤ N) (hK : 1 ≤ K)
    (hqδ : q ≤ δ)
    (hqmul : q * (s + 1 : ℕ) = δ / 2)
    (hsδ : 1 + δ / 2 ≤ (s : ℝ) * δ)
    (hKpower : K ≤ T ^ q) (hNlower : T ^ (δ + 1) ≤ N) :
    K * T / N ≤ 1 ∧ N * (K * T / N) ^ (s + 1) ≤ 1 := by
  have hTpos : 0 < T := zero_lt_one.trans_le hT
  have hNpos : 0 < N := zero_lt_one.trans_le hN
  have hKpow : K ^ (s + 1) ≤ T ^ (δ / 2) := by
    calc
      K ^ (s + 1) ≤ (T ^ q) ^ (s + 1) :=
        pow_le_pow_left₀ (zero_le_one.trans hK) hKpower _
      _ = T ^ (q * (s + 1 : ℕ)) :=
        (Real.rpow_mul_natCast hTpos.le q (s + 1)).symm
      _ = T ^ (δ / 2) := by rw [hqmul]
  constructor
  · apply (div_le_one hNpos).2
    calc
      K * T ≤ T ^ q * T := mul_le_mul_of_nonneg_right hKpower hTpos.le
      _ = T ^ (q + 1) := by
        simpa only [Real.rpow_one] using (Real.rpow_add hTpos q 1).symm
      _ ≤ T ^ (δ + 1) := Real.rpow_le_rpow_of_exponent_le hT (by linarith)
      _ ≤ N := hNlower
  · have hNpow : T ^ ((δ + 1) * s) ≤ N ^ s := by
      calc
        T ^ ((δ + 1) * s) = (T ^ (δ + 1)) ^ s :=
          Real.rpow_mul_natCast hTpos.le (δ + 1) s
        _ ≤ N ^ s := pow_le_pow_left₀ (Real.rpow_nonneg hTpos.le _) hNlower _
    have hexponent : δ / 2 + (s + 1 : ℕ) ≤ (δ + 1) * s := by
      push_cast
      linarith
    have hnum : K ^ (s + 1) * T ^ (s + 1) ≤ N ^ s := by
      calc
        K ^ (s + 1) * T ^ (s + 1) ≤ T ^ (δ / 2) * T ^ (s + 1 : ℕ) := by gcongr
        _ = T ^ (δ / 2 + (s + 1 : ℕ)) := by
          rw [← Real.rpow_natCast]
          exact (Real.rpow_add hTpos (δ / 2) (s + 1 : ℕ)).symm
        _ ≤ T ^ ((δ + 1) * s) := Real.rpow_le_rpow_of_exponent_le hT hexponent
        _ ≤ N ^ s := hNpow
    rw [show N * (K * T / N) ^ (s + 1) =
        (K ^ (s + 1) * T ^ (s + 1)) / N ^ s by
      rw [div_pow, mul_pow, pow_succ]
      field_simp [hNpos.ne', pow_ne_zero _ hNpos.ne']
      ring]
    exact (div_le_one (pow_pos hNpos _)).2 hnum

/-! ## The upper bound for `α > 1` -/

private theorem exponentSumGrowthExponent_le_sub_one
    {α : ℝ≥0} (hα : 1 < α) :
    exponentSumGrowthExponent α ≤ (α : ℝ) - 1 := by
  rw [exponentSumGrowthExponent_le_iff_nonAsymptotic]
  intro ε hε σ hσ
  have hαR : (1 : ℝ) < α := by exact_mod_cast hα
  let c : ℝ := (2 : ℝ) ^ (-σ) / 2
  have hc : 0 < c := by dsimp [c]; positivity
  let η : ℝ := min (((α : ℝ) - 1) / 4) (min (ε / 4) (min c 1))
  have hη : 0 < η := by
    dsimp [η]
    exact lt_min (by linarith) (lt_min (by linarith) (lt_min hc zero_lt_one))
  have hηα : η ≤ ((α : ℝ) - 1) / 4 := min_le_left _ _
  have hηε : η ≤ ε / 4 := (min_le_right _ _).trans (min_le_left _ _)
  have hηc : η ≤ min c 1 := (min_le_right _ _).trans (min_le_right _ _)
  let δ : ℝ := (α : ℝ) - η - 1
  have hδ : 0 < δ := by dsimp [δ]; linarith
  obtain ⟨s, q, hs, hq, hqδ, hqmul, hsδ⟩ :=
    exists_oscillatory_scale_parameters hδ
  obtain ⟨K, hK, hphaseBounds⟩ := approximate_model_phase_deriv_bounds hσ s
  obtain ⟨M, hMone, hsumBound⟩ := exists_norm_oscillatory_sum_le s hs hK hc
  have hM : 0 ≤ M := zero_le_one.trans hMone
  let C : ℝ := max 1 (max (2 * M) (K ^ q⁻¹))
  have hC : 1 ≤ C := le_max_left _ _
  refine ⟨η, hη, s, hs, C, hC, ?_⟩
  intro T N F a b hTC hNlower hNupper hF ha hb
  have hTone : 1 ≤ T := hC.trans hTC
  have hTpos : 0 < T := zero_lt_one.trans_le hTone
  have hNone : 1 ≤ N := by
    calc
      1 ≤ T ^ ((α : ℝ) - η) := Real.one_le_rpow hTone (by linarith)
      _ ≤ N := hNlower
  have hKbase : K ^ q⁻¹ ≤ T := by
    calc
      K ^ q⁻¹ ≤ max (2 * M) (K ^ q⁻¹) := le_max_right _ _
      _ ≤ C := le_max_right _ _
      _ ≤ T := hTC
  have hKpower : K ≤ T ^ q := by
    calc
      K = (K ^ q⁻¹) ^ q := by
        rw [← Real.rpow_mul (zero_le_one.trans hK), inv_mul_cancel₀ hq.ne', Real.rpow_one]
      _ ≤ T ^ q := Real.rpow_le_rpow (Real.rpow_nonneg (zero_le_one.trans hK) _)
        hKbase hq.le
  have hNscale : T ^ (δ + 1) ≤ N := by
    simpa only [δ, sub_add_cancel] using hNlower
  obtain ⟨hratio, hremainder⟩ := oscillatory_scale_bounds
    hTone hNone hK hqδ hqmul hsδ hKpower hNscale
  obtain ⟨hfirst, hderiv⟩ := hphaseBounds hηc hF
  have htargetNonneg : 0 ≤ (α : ℝ) - 1 + ε := by linarith
  have hNTdiv : N / T ≤ T ^ ((α : ℝ) - 1 + ε) := by
    calc
      N / T ≤ T ^ ((α : ℝ) + η) / T :=
        div_le_div_of_nonneg_right hNupper hTpos.le
      _ = T ^ ((α : ℝ) + η - 1) := by
        simpa only [Real.rpow_one] using
          (Real.rpow_sub hTpos ((α : ℝ) + η) 1).symm
      _ ≤ T ^ ((α : ℝ) - 1 + ε) :=
        Real.rpow_le_rpow_of_exponent_le hTone (by linarith)
  have hone : 1 ≤ T ^ ((α : ℝ) - 1 + ε) := Real.one_le_rpow hTone htargetNonneg
  have hCM : 2 * M ≤ C := (le_max_left _ _).trans (le_max_right _ _)
  by_cases hab : a ≤ b
  · have hsum := hsumBound hab hNone hTone hF.1
      (fun u hu ↦ hfirst u hu) hderiv ha hb hratio hremainder
    change ‖exponentialSumAt F T N a b‖ ≤ C * T ^ ((α : ℝ) - 1 + ε)
    calc
      ‖exponentialSumAt F T N a b‖ ≤ M * (1 + N / T) := by
        simpa only [exponentialSumAt] using hsum
      _ ≤ M * (2 * T ^ ((α : ℝ) - 1 + ε)) := by
        gcongr
        linarith
      _ = (2 * M) * T ^ ((α : ℝ) - 1 + ε) := by ring
      _ ≤ C * T ^ ((α : ℝ) - 1 + ε) := by gcongr
  · rw [exponentialSumAt, Finset.Icc_eq_empty hab, Finset.sum_empty, norm_zero]
    exact mul_nonneg (zero_le_one.trans hC) (Real.rpow_nonneg hTpos.le _)

/-! ## The matching logarithmic lower bound for `α > 1` -/

/-- For each Euler–Maclaurin order, the logarithmic exponential sum differs from its
explicit main term by a uniformly bounded amount once the scale conditions hold. This is
the formal counterpart of the blueprint's estimate for `∑ n ^ (iT)`, with the factor `2π`
absorbed into Lean's Fourier character. -/
theorem exists_norm_logPhase_sum_sub_mainTerm_le (s : ℕ) :
    ∃ K C : ℝ, 1 ≤ K ∧ 1 ≤ C ∧
      ∀ {N T : ℝ} {a b : ℕ},
        a < b → 1 ≤ N → 1 ≤ T →
        (a : ℝ) = N → (b : ℝ) = 2 * N →
        K * T / N ≤ 1 → N * (K * T / N) ^ (s + 1) ≤ 1 →
        ‖(∑ n ∈ Finset.Icc a b, (𝐞 (T * Real.log ((n : ℝ) / N)) : ℂ)) -
          logPhaseMainTerm N T‖ ≤ C := by
  obtain ⟨K, hK, hphaseBounds⟩ := approximate_model_phase_deriv_bounds zero_lt_one s
  obtain ⟨_, hderiv⟩ := hphaseBounds (δ := 0) (F := Real.log) (by norm_num)
    (isApproximateModelPhaseFunction_log s)
  obtain ⟨C, hC, herror⟩ := exists_norm_oscillatory_sum_sub_integral_le s
  refine ⟨K, C, hK, hC, ?_⟩
  intro N T a b hab hN hT ha hb hratio hremainder
  have herr := herror hab hN hT hK
    (isModelPhaseFunction_log.1 0) hderiv (by exact_mod_cast ha.ge)
    (by exact_mod_cast hb.le) hratio hremainder
  have hint : (∫ x in (a : ℝ)..b,
      (𝐞 (T * Real.log (x / N)) : ℂ)) = logPhaseMainTerm N T := by
    rw [ha, hb]
    exact logPhase_integral_eq_mainTerm (zero_lt_one.trans_le hN)
  simp only [logPhase] at herr
  rw [hint] at herr
  exact herr

private theorem sub_one_le_exponentSumGrowthExponent
    {α : ℝ≥0} (hα : 1 < α) :
    (α : ℝ) - 1 ≤ exponentSumGrowthExponent α := by
  have hαR : (1 : ℝ) < α := by exact_mod_cast hα
  let δ : ℝ := (α : ℝ) - 1
  have hδ : 0 < δ := by dsimp [δ]; linarith
  obtain ⟨s, q, hsone, hq, hqδ, hqmul, hsδ⟩ :=
    exists_oscillatory_scale_parameters hδ
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
  obtain ⟨K, E, hK, hEone, hlogEstimate⟩ :=
    exists_norm_logPhase_sum_sub_mainTerm_le s
  have hKpower : ∀ᶠ i in atTop, K ≤ T i ^ q :=
    tendsto_atTop.1 ((tendsto_rpow_atTop hq).comp hTtop) K
  let d : ℝ := 1 / (1 + 2 * Real.pi)
  have hd : 0 < d := by dsimp [d]; positivity
  have hE : 0 ≤ E := zero_le_one.trans hEone
  have herrorSmall : ∀ᶠ i in atTop, 2 * E / d ≤ T i ^ ((α : ℝ) - 1) :=
    tendsto_atTop.1 ((tendsto_rpow_atTop (by linarith : 0 < (α : ℝ) - 1)).comp hTtop)
      (2 * E / d)
  have hlower : ∀ᶠ i in atTop,
      (d / 2) * T i ^ ((α : ℝ) - 1) ≤
        ‖exponentialSum logPhase T N a b i‖ := by
    filter_upwards [hKpower, herrorSmall] with i hiK hiSmall
    have hTi : 0 < T i := zero_lt_one.trans_le (hT i)
    have hNi : 0 < N i := zero_lt_one.trans_le (hN i)
    have hNscale : T i ^ (δ + 1) ≤ N i := by
      rw [hNTexact]
      apply le_of_eq
      congr 1
      dsimp [δ]
      ring
    obtain ⟨hratio, hremainder⟩ := oscillatory_scale_bounds
      (hT i) (hN i) hK hqδ hqmul hsδ hiK hNscale
    have hablt : a i < b i := by dsimp [a, b]; omega
    have haeq : (a i : ℝ) = N i := by dsimp [N, a]; push_cast; ring
    have hbeq : (b i : ℝ) = 2 * N i := by dsimp [N, b]; push_cast; ring
    have herr := hlogEstimate hablt (hN i) (hT i) haeq hbeq hratio hremainder
    have hInt := (norm_logPhaseMainTerm_bounds hNi (hT i)).1
    have hpower : N i / T i = T i ^ ((α : ℝ) - 1) := by
      rw [hNTexact]
      simpa only [Real.rpow_one] using (Real.rpow_sub hTi (α : ℝ) 1).symm
    rw [hpower] at hInt
    have htriangle :
        ‖logPhaseMainTerm (N i) (T i)‖ ≤
          ‖exponentialSum logPhase T N a b i‖ + E := by
      calc
        _ ≤ ‖exponentialSum logPhase T N a b i‖ +
            ‖exponentialSum logPhase T N a b i -
              logPhaseMainTerm (N i) (T i)‖ :=
          norm_le_norm_add_norm_sub _ _
        _ ≤ ‖exponentialSum logPhase T N a b i‖ + E := by
          gcongr
          simpa only [exponentialSum, logPhase] using herr
    have hEsmall : E ≤ d / 2 * T i ^ ((α : ℝ) - 1) := by
      calc
        E = d / 2 * (2 * E / d) := by field_simp [hd.ne']
        _ ≤ d / 2 * T i ^ ((α : ℝ) - 1) :=
          mul_le_mul_of_nonneg_left hiSmall (by positivity)
    linarith
  exact exponent_le_of_isPowerBounded_of_eventually_norm_ge_rpow
    hT hTunbounded hbound (half_pos hd) (by simpa [d] using hlower)

/-! ## The `L²` lower bound for `0 < α ≤ 1` -/

private theorem exists_logPhase_sum_norm_sq_ge :
    ∃ C : ℝ, 0 < C ∧
      ∀ {R N : ℝ} {a b : ℕ},
        0 < R → 0 < N → N ≤ a → (b : ℝ) ≤ 2 * N → N ≤ R / (4 * C) →
        ∃ t ∈ Set.Icc R (2 * R),
          ((Finset.Icc a b).card : ℝ) / 2 ≤
            ‖∑ n ∈ Finset.Icc a b,
              (𝐞 (t * Real.log ((n : ℝ) / N)) : ℂ)‖ ^ 2 := by
  obtain ⟨C₀, hC₀, hlarge⟩ := l2_integral_estimate_exists_norm_sq_ge
  refine ⟨C₀, hC₀, ?_⟩
  intro R N a b hR hN ha hb hNR
  let ξ : ↥(Finset.Icc a b) → ℝ := fun n ↦ Real.log ((n : ℝ) / N)
  let coeff : ↥(Finset.Icc a b) → ℂ := fun _ ↦ 1
  have hsep : IsSeparatedFamily (1 / (2 * N)) ξ :=
    log_div_separated hN ha hb
  have hlength : 2 * C₀ * (2 * N) ≤ 2 * R - R := by
    have hCN : N * (4 * C₀) ≤ R :=
      (le_div_iff₀ (by positivity : 0 < 4 * C₀)).1 hNR
    nlinarith
  obtain ⟨t, ht, htlarge⟩ := hlarge coeff ξ (2 * N) (by positivity) hsep
    R (2 * R) hlength
  have hmass : ∑ n, ‖coeff n‖ ^ 2 = ((Finset.Icc a b).card : ℝ) := by
    simp [coeff]
  rw [hmass] at htlarge
  refine ⟨t, ht, ?_⟩
  rw [Finset.sum_subtype (Finset.Icc a b) (fun _ ↦ Iff.rfl)]
  simpa only [coeff, ξ, one_mul, mul_comm] using htlarge

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
  have hNT : IsPowerAsymptotic N T (α : ℝ) :=
    isPowerAsymptotic_of_between_const_rpow (A := D) (B := 2 * D) hαR
      (zero_lt_one.trans_le hD) (mul_pos (by norm_num) (zero_lt_one.trans_le hD)) hNtop
      (Filter.Eventually.of_forall fun i ↦ by
        simpa only [mul_assoc] using hTbounds i)
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

/-! ## Trivial bounds -/

private theorem exponentSumGrowthExponent_nonneg (α : ℝ≥0) :
    0 ≤ exponentSumGrowthExponent α :=
  (isExponentSumBound_exponentSumGrowthExponent α).nonneg

private theorem exponentSumGrowthExponent_le_self (α : ℝ≥0) :
    exponentSumGrowthExponent α ≤ (α : ℝ) :=
  exponentSumGrowthExponent_le_iff.mpr (isExponentSumBound_self α)

/-- For `α > 1`, the exponential sum growth exponent is `α - 1`. -/
theorem exponentSumGrowthExponent_eq_sub_one
    {α : ℝ≥0} (hα : 1 < α) :
    exponentSumGrowthExponent α = (α : ℝ) - 1 := by
  exact le_antisymm
    (exponentSumGrowthExponent_le_sub_one hα)
    (sub_one_le_exponentSumGrowthExponent hα)

/-- For `0 ≤ α ≤ 1`, the exponential sum growth exponent lies in
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
