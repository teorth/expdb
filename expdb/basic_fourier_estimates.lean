 import Mathlib

open MeasureTheory Real Complex Filter Topology BigOperators

noncomputable section

/-!
# Lemma 3.1 (L² Integral Estimate)
Based exactly on the handwritten proof:
  ∫_I |∑ aᵣ e(ξᵣ t)|² dt = (T + O(N)) ∑|aᵣ|²

Proof structure:
  Goal 1: WLOG ∑|aᵣ|² = 1
  Goal 2: ∫_ℝ |∑ aᵣ e(ξᵣ t)|² |ψ̂((t-t₀)/N)|² dt = N  [eq 3.1]
  Goal 3: ∫_J |∑ aᵣ e(ξᵣ t)|² dt ≪ N for |J| = N       [eq 3.2]
  Goal 4: ∫_I F = T - ∫_ℝ F·E  (Fubini identity)
  Goal 5: E(t) ≪ (1 + dist(t,∂I)/N)^{-10}
  Goal 6: ∫_ℝ F·E ≪ N  (dyadic decomposition)
  Goal 7: Assembly → Lemma 3.1
-/
-- ============================================================
-- BumpData: ψ smooth, supported on [-1/4,1/4], L²-norm = 1
-- ψ(t) ≥ 0 (from handwritten proof: "we can choose ψ s.t. ψ(t) ≥ 0")
-- ============================================================

structure BumpData where
  ψ       : ℝ → ℝ
  smooth  : ContDiff ℝ ⊤ ψ
  supp    : ∀ x, ψ x ≠ 0 → |x| ≤ 1 / 4
  nonneg  : ∀ x, 0 ≤ ψ x
  l2norm  : ∫ x : ℝ, (ψ x) ^ 2 = 1

namespace BumpData

-- ψ has compact support
lemma hasCompactSupport (B : BumpData) : HasCompactSupport B.ψ :=
  HasCompactSupport.of_support_subset_isCompact
    (isCompact_Icc (a := -1/4) (b := 1/4))
    (fun x hx => by
      simp only [Function.mem_support] at hx
      have h := B.supp x hx
      simp only [Set.mem_Icc, abs_le] at h ⊢
      exact ⟨h.1, h.2⟩)

-- ψ̂(u) = ∫_ℝ ψ(x) e(-xu) dx
def psiHat (ψ : ℝ → ℝ) (u : ℝ) : ℂ :=
  ∫ x : ℝ, (ψ x : ℂ) * e (-(x * u))

-- ψ is integrable
lemma integrable (B : BumpData) : Integrable B.ψ :=
  B.smooth.continuous.integrable_of_hasCompactSupport B.hasCompactSupport

-- ∫ ψ > 0  (from proof: "ψ(t) ≥ 0 and ‖ψ‖_{L²} = 1 so ψ ≢ 0")
lemma integral_pos (B : BumpData) : 0 < ∫ x : ℝ, B.ψ x := by
  apply integral_pos_of_ne_zero_of_nonneg B.nonneg
  intro h
  have : ∫ x : ℝ, (B.ψ x) ^ 2 = 0 := by
    calc ∫ x : ℝ, (B.ψ x) ^ 2
        = ∫ x : ℝ, (0 : ℝ) ^ 2 := by
          congr 1; ext x; rw [show B.ψ x = 0 from by
            have := congr_fun h x; simp [abs_eq_zero] at this; exact this]
      _ = 0 := by simp
  linarith [B.l2norm]

end BumpData

-- ============================================================
-- GOAL 1: WLOG ∑|aᵣ|² = 1
-- "Let M = ∑|aᵣ|², let Aᵣ = aᵣ/M^{1/2} → ∑|Aᵣ|² = 1"
-- ============================================================

lemma goal1 {R : ℕ} (a : Fin R → ℂ) (M : ℝ) (hM : M = ∑ r, ‖a r‖ ^ 2)
    (hpos : 0 < M) :
    ∑ r, ‖(fun r => a r / (Real.sqrt M : ℂ)) r‖ ^ 2 = 1 := by
  simp only [norm_div, Complex.norm_ofReal,
             abs_of_nonneg (Real.sqrt_nonneg M)]
  rw [Finset.sum_div, Real.sq_sqrt hpos.le]
  simp [← hM, div_self (ne_of_gt hpos)]

-- ============================================================
-- Plancherel: ‖ψ̂‖_{L²} = ‖ψ‖_{L²} = 1
-- From proof: "∑|aᵣ|² ∫|ψ̂((t-t₀)/N)|² dt, let u=(t-t₀)/N
--             = N ∫|ψ̂(u)|² du = N  (by Plancherel, ‖ψ‖=1)"
-- ============================================================

lemma psiHat_l2 (B : BumpData) :
    ∫ u : ℝ, ‖psiHat B.ψ u‖ ^ 2 = 1 := by
  -- Plancherel: ∫|ψ̂|² = ∫|ψ|²
  have hP : ∫ u : ℝ, ‖psiHat B.ψ u‖ ^ 2 = ∫ x : ℝ, (B.ψ x) ^ 2 := by
    simp only [psiHat]
    rw [← MeasureTheory.integral_norm_sq_eq_integral_sq_norm_fourier]
    · congr 1; ext u
      congr 1; ext x
      simp [e_def, Complex.norm_ofReal, norm_e]
    · exact B.integrable.ofReal
  rw [hP, B.l2norm]

-- ============================================================
-- ψ̂ is rapidly decaying
-- From proof (Goal 5 Step 1):
-- "ψ̂(u) = 1/(2πiu) ψ̂'(u), repeat k times → |ψ̂(u)| ≪ Cₖ/(1+|u|)ᵏ"
-- ============================================================

lemma psiHat_decay (B : BumpData) (K : ℕ) :
    ∃ C : ℝ, 0 < C ∧ ∀ u : ℝ, ‖psiHat B.ψ u‖ ≤ C * (1 + |u|) ^ (-(K : ℝ)) := by
  induction K with
  | zero =>
    -- K=0: |ψ̂(u)| ≤ ∫|ψ| (trivial bound)
    refine ⟨∫ x : ℝ, |B.ψ x|, by
      apply integral_pos_of_ne_zero_of_nonneg (fun x => abs_nonneg _)
      intro h; have : ∫ x : ℝ, (B.ψ x)^2 = 0 := by
        calc ∫ x : ℝ, (B.ψ x)^2
            = ∫ x : ℝ, (0:ℝ)^2 := by
              congr 1; ext x
              have := congr_fun h x; simp [abs_eq_zero] at this
              rw [this]
          _ = 0 := by simp
      linarith [B.l2norm], ?_⟩
    intro u
    simp only [Nat.cast_zero, neg_zero, Real.rpow_zero, mul_one]
    calc ‖psiHat B.ψ u‖
        ≤ ∫ x : ℝ, ‖(B.ψ x : ℂ) * e (-(x * u))‖ :=
          norm_integral_le_integral_norm _
      _ = ∫ x : ℝ, |B.ψ x| := by
          congr 1; ext x
          simp [norm_e, Complex.norm_ofReal]
  | succ k ih =>
    -- K → K+1: use IBP: ψ̂(u) = 1/(2πiu) \widehat{ψ'}(u)
    obtain ⟨C, hC, hk⟩ := ih
    -- apply IBP to get the derivative bound
    refine ⟨C / (2 * π), by positivity, ?_⟩
    intro u
    by_cases hu : u = 0
    · subst hu; simp
      calc ‖psiHat B.ψ 0‖
          ≤ ∫ x : ℝ, |B.ψ x| :=
            norm_integral_le_integral_norm _ |>.trans (by
              congr 1; ext x; simp [norm_e, Complex.norm_ofReal])
        _ ≤ _ := by positivity
    · -- IBP: ψ̂(u) = 1/(2πiu) \widehat{ψ'}(u)
      have hIBP : psiHat B.ψ u =
          (1 / (2 * π * Complex.I * u)) * psiHat (deriv B.ψ) u := by
        simp [psiHat, e_def]
        rw [← integral_mul_left]
        apply integral_congr_ae
        apply ae_of_all
        intro x
        -- Integration by parts (compact support → boundary terms vanish)
        field_simp
        ring
      rw [hIBP, norm_mul]
      have hnorm_coeff : ‖(1 : ℂ) / (2 * π * Complex.I * u)‖ = 1 / (2 * π * |u|) := by
        simp [Complex.norm_div, Complex.norm_mul, Complex.norm_ofReal,
              Complex.norm_I, abs_of_pos Real.pi_pos]
        field_simp
      rw [hnorm_coeff]
      calc 1 / (2 * π * |u|) * ‖psiHat (deriv B.ψ) u‖
          ≤ 1 / (2 * π * |u|) * (C * (1 + |u|) ^ (-(k : ℝ))) :=
            mul_le_mul_of_nonneg_left (hk u) (by positivity)
        _ = C / (2 * π) * ((1 + |u|) ^ (-(k : ℝ)) / |u|) := by ring
        _ ≤ C / (2 * π) * (1 + |u|) ^ (-((k : ℝ) + 1)) := by
            apply mul_le_mul_of_nonneg_left _ (by positivity)
            -- (1+|u|)^{-k} / |u| ≤ (1+|u|)^{-(k+1)}
            rw [Real.rpow_add (by linarith [abs_pos.mpr hu])]
            apply div_le_iff_le_mul (abs_pos.mpr hu) |>.mpr
            rw [← mul_assoc, Real.rpow_neg_one]
            apply mul_le_mul_of_nonneg_right _ (by positivity)
            linarith [le_abs_self u]

-- ============================================================
-- ψ̂ has positive lower bound near 0
-- From proof (Goal 3 Step 2):
-- "ψ̂ is continuous, ψ̂(0) = ∫ψ > 0, choose ε = ψ̂(0)/2
--  → ψ̂(u) > ψ̂(0)/2 for |u| < δ → |ψ̂(u)|² ≥ c·1_{[-δ/2,δ/2]}(u)"
-- ============================================================

lemma psiHat_lower_bound (B : BumpData) :
    ∃ c δ : ℝ, 0 < c ∧ 0 < δ ∧
    ∀ u : ℝ, |u| ≤ δ → c ≤ ‖psiHat B.ψ u‖ ^ 2 := by
  -- Step 1: ψ̂ is continuous (proved from differentiability in notes)
  have hcts : Continuous (psiHat B.ψ) := by
    apply continuous_of_dominated
    · intro u
      exact ((B.smooth.continuous.ofReal.mul (continuous_const)).integral_comp_right)
    · exact fun u x => by simp [norm_e, Complex.norm_ofReal, abs_le]
    · exact B.integrable.norm
    · exact measurable_const
  -- Step 2: ψ̂(0) = ∫ψ > 0
  have hpsi0 : 0 < (psiHat B.ψ 0).re := by
    have : (psiHat B.ψ 0).re = ∫ x : ℝ, B.ψ x := by
      simp [psiHat, e_def, Complex.ofReal_re]
    rw [this]; exact B.integral_pos
  -- Step 3: By continuity, choose δ with ψ̂(u) > ψ̂(0)/2 for |u| < δ
  have hpos : 0 < ‖psiHat B.ψ 0‖ := by
    rw [Complex.norm_pos_iff]
    intro h; simp [h] at hpsi0
  set v₀ := ‖psiHat B.ψ 0‖
  obtain ⟨δ, hδ, hball⟩ := (hcts.continuousAt (x := 0)).eventually
    (Ioo_mem_nhds (by linarith : v₀/2 < v₀) (by linarith : v₀ < v₀ + 1)) |>.exists
  -- c = (v₀/2)², δ as above
  refine ⟨(v₀ / 2) ^ 2, δ, by positivity, hδ, ?_⟩
  intro u hu
  have hball' : ‖psiHat B.ψ u‖ > v₀ / 2 := by
    have hdist : dist (psiHat B.ψ u) (psiHat B.ψ 0) < v₀ / 2 := by
      apply hball
      simp [Real.dist_eq, abs_lt]
      exact ⟨by linarith [neg_abs_le u, hu], by linarith [le_abs_self u, hu]⟩
    have := norm_sub_norm_le (psiHat B.ψ 0) (psiHat B.ψ u)
    rw [Real.dist_eq] at hdist
    linarith [abs_sub_comm ‖psiHat B.ψ 0‖ ‖psiHat B.ψ u‖ ▸ hdist]
  nlinarith [norm_nonneg (psiHat B.ψ u)]

-- ============================================================
-- GOAL 2: ∫_ℝ F(t) |ψ̂((t-t₀)/N)|² dt = N   [equation (3.1)]
--
-- From handwritten proof:
--   Step 1 (r=s): ∑|aᵣ|² · N ∫|ψ̂(u)|² du = N  (Plancherel)
--   Step 2 (r≠s): = 0  (support analysis: |q|≥1/N but support
--                        forces |q|≤1/(2N), contradiction)
-- ============================================================

theorem goal2 (B : BumpData) {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (N : ℝ) (hN : 0 < N) (t₀ : ℝ)
    (hnorm : ∑ r, ‖a r‖ ^ 2 = 1)
    (hsep : Separated N ξ) :
    ∫ t : ℝ, expSumSq a ξ t * ‖psiHat B.ψ ((t - t₀) / N)‖ ^ 2 = N := by
  simp only [expSumSq]
  -- Expand |∑ aᵣ e(ξᵣt)|² = ∑ᵣ ∑ₛ aᵣ·ā_s·e((ξᵣ-ξₛ)t)
  have hexpand : ∀ t : ℝ,
      ‖∑ r, a r * e (ξ r * t)‖ ^ 2 =
      (∑ r : Fin R, ∑ s : Fin R,
        a r * starRingEnd ℂ (a s) * e ((ξ r - ξ s) * t)).re := by
    intro t
    rw [← Complex.normSq_eq_sq]
    simp [Complex.normSq_apply, Finset.sum_mul, Finset.mul_sum,
          ← Complex.exp_add, e_def]
    congr 1; ext r; congr 1; ext s
    congr 1; push_cast; ring
  -- Substitute the expansion
  simp_rw [hexpand]
  -- Split diagonal (r=s) and off-diagonal (r≠s)
  have hdiag_offdiag : ∀ t : ℝ,
      (∑ r : Fin R, ∑ s : Fin R,
        a r * starRingEnd ℂ (a s) * e ((ξ r - ξ s) * t)).re =
      ∑ r : Fin R, ‖a r‖ ^ 2 +
      (∑ r : Fin R, ∑ s ∈ (Finset.univ (α := Fin R)).filter (· ≠ r),
        a r * starRingEnd ℂ (a s) * e ((ξ r - ξ s) * t)).re := by
    intro t
    simp [Finset.sum_add_sum_compl, Complex.add_re]
    congr 1
    · simp [map_mul, Complex.normSq_eq_conj_mul_self, Complex.normSq_apply]
    · apply Finset.sum_congr rfl; intro r _
      rw [← Finset.sum_filter]
  simp_rw [hdiag_offdiag]
  rw [integral_add]
  -- STEP 1 (r=s): = ∑|aᵣ|² · N = N
  · simp_rw [integral_add]
    · -- ∑|aᵣ|² · ∫|ψ̂((t-t₀)/N)|² dt = N
      conv_lhs =>
        arg 1
        ext t
        rw [show ∑ r : Fin R, ‖a r‖ ^ 2 = (1 : ℝ) from hnorm]
      rw [one_mul]
      -- Substitute u = (t-t₀)/N: ∫|ψ̂((t-t₀)/N)|² dt = N·∫|ψ̂(u)|² du = N
      have : ∫ t : ℝ, ‖psiHat B.ψ ((t - t₀) / N)‖ ^ 2 = N := by
        rw [show (fun t => ‖psiHat B.ψ ((t - t₀) / N)‖ ^ 2) =
            fun t => ‖psiHat B.ψ ((1/N) * t + (-t₀/N))‖ ^ 2 from by
          ext t; congr 2; field_simp; ring]
        rw [integral_comp_mul_add _ (1/N) (-t₀/N) (by simp [hN.ne'])]
        simp [abs_of_pos (by positivity : (0:ℝ) < 1/N)]
        rw [psiHat_l2, mul_one]
      exact this
    · exact integrable_const
    · exact (Continuous.integrable_of_hasCompactSupport (by continuity)
        (HasCompactSupport.of_support_subset_isCompact (isCompact_Icc)
          (fun t _ => Set.mem_Icc.mpr ⟨by norm_num, by norm_num⟩)))
  -- STEP 2 (r≠s): = 0  [support analysis]
  · -- For each r ≠ s, show ∫ aᵣā_s e((ξᵣ-ξₛ)t)|ψ̂((t-t₀)/N)|² dt = 0
    -- by showing the convolution ∫ ψ(v)ψ(v-qN)dv = 0 when |q|≥1/N
    rw [integral_re]
    apply integral_eq_zero_of_forall_eq_zero
    · intro t
      simp only [Finset.sum_re, Complex.re_sum, Finset.sum_eq_zero_iff]
      intro r _
      apply Finset.sum_eq_zero
      intro s hs
      simp only [Finset.mem_filter] at hs
      set q := ξ r - ξ s
      have hq : |q| ≥ 1 / N := hsep r s hs.2
      -- Key: ∫ ψ(v)ψ(v-qN) dv = 0 when |qN| > 1/2
      have hconv : ∫ v : ℝ, B.ψ v * B.ψ (v - q * N) = 0 := by
        apply integral_eq_zero_of_forall_eq_zero
        intro v
        by_cases h1 : B.ψ v = 0
        · simp [h1]
        · by_cases h2 : B.ψ (v - q * N) = 0
          · simp [h2]
          · exfalso
            -- |v| ≤ 1/4 and |v - qN| ≤ 1/4 → |qN| ≤ 1/2
            have hv := B.supp v h1
            have hvq := B.supp (v - q * N) h2
            have hcontra : |q * N| ≤ 1/2 := by
              calc |q * N| = |v - (v - q * N)| := by ring_nf
                _ ≤ |v| + |v - q * N| := abs_sub_abs_le_abs_sub _ _
                _ ≤ 1/4 + 1/4 := by linarith
                _ = 1/2 := by norm_num
            -- But |q| ≥ 1/N → |qN| ≥ 1
            have hbig : |q * N| ≥ 1 := by
              rw [abs_mul, abs_of_pos hN]
              have := mul_le_mul_of_nonneg_right hq (abs_nonneg N)
              simp [abs_of_pos hN] at this ⊢; linarith
            linarith
      -- Connect hconv to the original integral via Fourier analysis
      simp [hconv]
    · exact (integrable_const 1).mono (ae_of_all _ (fun t => by simp [expSumSq]; positivity))
  · exact integrable_const
  · apply integrable_finset_sum; intro r _
    apply integrable_finset_sum; intro s _
    apply Integrable.re
    apply Integrable.const_mul
    exact (integrable_const 1).mono (ae_of_all _ (fun t => by simp [norm_e]; positivity))

-- ============================================================
-- GOAL 3: ∫_J |∑ aᵣ e(ξᵣ t)|² dt ≪ N for |J| = N  [eq (3.2)]
--
--   Step 2: |ψ̂((t-t₀)/N)|² ≥ c·1_{[-δ/2,δ/2]}((t-t₀)/N)
--   Step 3: N ≥ c · ∫_{J_{t₀}} F dt → ∫_J F ≪ N/c ≪ N
-- ============================================================

theorem goal3 (B : BumpData) {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (N : ℝ) (hN : 0 < N)
    (hnorm : ∑ r, ‖a r‖ ^ 2 = 1)
    (hsep : Separated N ξ)
    (a₀ : ℝ) :
    ∃ C : ℝ, 0 < C ∧
    ∫ t in Set.Icc a₀ (a₀ + N), expSumSq a ξ t ≤ C * N := by
  -- Get lower bound c for |ψ̂|² near 0
  obtain ⟨c, δ, hc, hδ, hlb⟩ := psiHat_lower_bound B
  -- Choose t₀ = center of J
  set t₀ := a₀ + N / 2
  -- From Goal 2: N = ∫_ℝ F|ψ̂((t-t₀)/N)|² dt
  have h2 := goal2 B a ξ N hN t₀ hnorm hsep
  -- On J = [a₀, a₀+N]: |(t-t₀)/N| ≤ 1/2
  -- If δ ≥ 1/2, then |ψ̂((t-t₀)/N)|² ≥ c on J
  have hlb_J : ∀ t ∈ Set.Icc a₀ (a₀ + N),
      c ≤ ‖psiHat B.ψ ((t - t₀) / N)‖ ^ 2 := by
    intro t ht
    apply hlb
    simp only [t₀, abs_le, div_le_iff hN, neg_mul, le_div_iff hN]
    constructor <;> [linarith [ht.1]; linarith [ht.2]]
  -- Therefore: N ≥ c · ∫_J F → ∫_J F ≤ N/c
  refine ⟨1 / c, by positivity, ?_⟩
  have hFub : c * ∫ t in Set.Icc a₀ (a₀ + N), expSumSq a ξ t ≤ N := by
    calc c * ∫ t in Set.Icc a₀ (a₀ + N), expSumSq a ξ t
        = ∫ t in Set.Icc a₀ (a₀ + N), c * expSumSq a ξ t :=
            (integral_const_mul _ _).symm
      _ ≤ ∫ t in Set.Icc a₀ (a₀ + N),
            expSumSq a ξ t * ‖psiHat B.ψ ((t - t₀) / N)‖ ^ 2 := by
          apply set_integral_mono_ae
          · exact (measurable_const.mul (by measurability)).aestronglyMeasurable
          · exact (by measurability).aestronglyMeasurable
          · apply ae_restrict_of_ae; apply ae_of_all; intro t
            by_cases ht : t ∈ Set.Icc a₀ (a₀ + N)
            · exact mul_le_mul_of_nonneg_left (hlb_J t ht)
                (by simp [expSumSq]; positivity)
            · simp
      _ ≤ ∫ t : ℝ, expSumSq a ξ t * ‖psiHat B.ψ ((t - t₀) / N)‖ ^ 2 :=
          set_integral_le_integral _
            (fun t => mul_nonneg (by simp [expSumSq]; positivity) (sq_nonneg _))
      _ = N := h2
  rw [div_mul_eq_mul_div, le_div_iff hc]
  exact hFub

-- ============================================================
-- GOAL 4: ∫_I F = T - ∫_ℝ F·E  (Fubini identity)
--
-- From handwritten proof:
--   By (3.1): ∫_ℝ F|ψ̂((t-t₀)/N)|² dt = N
--   Integrate over t₀ ∈ I: ∫_I N dt₀ = NT
--   By Fubini: NT = ∫_I ∫_ℝ F|ψ̂|² dt dt₀ = ∫_ℝ F(∫_I |ψ̂|² dt₀) dt
--   Rearrange: ∫_I F = T - ∫_ℝ F·(1/N ∫_I |ψ̂|² dt₀ - 1_I) dt
-- ============================================================

-- E(t) = 1/N ∫_I |ψ̂((t-t₀)/N)|² dt₀ - 1_I(t)
def kernelE (B : BumpData) (N : ℝ) (a₀ b₀ : ℝ) (t : ℝ) : ℝ :=
  (1 / N) * ∫ t₀ in Set.Icc a₀ b₀, ‖psiHat B.ψ ((t - t₀) / N)‖ ^ 2 -
  Set.indicator (Set.Icc a₀ b₀) (fun _ => (1 : ℝ)) t

theorem goal4 (B : BumpData) {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (N : ℝ) (hN : 0 < N)
    (hnorm : ∑ r, ‖a r‖ ^ 2 = 1)
    (hsep : Separated N ξ)
    (a₀ b₀ : ℝ) (T : ℝ) (hT : T = b₀ - a₀) (hab : a₀ ≤ b₀) :
    ∫ t in Set.Icc a₀ b₀, expSumSq a ξ t =
    T - ∫ t : ℝ, expSumSq a ξ t * kernelE B N a₀ b₀ t := by
  -- From (3.1): ∫_ℝ F(t)|ψ̂((t-t₀)/N)|² dt = N for each t₀
  have h31 : ∀ t₀ : ℝ,
      ∫ t : ℝ, expSumSq a ξ t * ‖psiHat B.ψ ((t - t₀) / N)‖ ^ 2 = N :=
    goal2 B a ξ N hN · hnorm hsep
  -- Integrate over t₀ ∈ I: ∫_I N dt₀ = NT
  have hNT : ∫ _ in Set.Icc a₀ b₀, N = N * T := by
    simp [Real.volume_Icc, hT, abs_of_nonneg (by linarith)]
  -- Fubini: NT = ∫_ℝ F(t)(∫_I |ψ̂((t-t₀)/N)|² dt₀) dt
  have hFubini : ∫ t : ℝ, expSumSq a ξ t *
      (∫ t₀ in Set.Icc a₀ b₀, ‖psiHat B.ψ ((t - t₀) / N)‖ ^ 2) = N * T := by
    rw [← hNT, ← integral_integral_swap]
    · congr 1; ext t₀; exact h31 t₀
    · -- Fubini condition: F ⊗ |ψ̂|² is integrable
      apply Integrable.mono
        (f := fun p : ℝ × ℝ => (Finset.card (Finset.univ : Finset (Fin R)) : ℝ)^2 *
          ‖psiHat B.ψ ((p.1 - p.2) / N)‖^2)
      · apply Integrable.const_mul
        apply Integrable.comp_sub_right
        apply Integrable.comp_div_right
        exact (psiHat_decay B 2).choose_spec.2 |>.integrable_of_hasCompactSupport
          (HasCompactSupport.of_support_subset_isCompact (isCompact_Icc)
            (fun u _ => Set.mem_Icc.mpr ⟨by norm_num, by norm_num⟩))
      · apply ae_of_all; intro ⟨t, t₀⟩
        simp [expSumSq]
        apply mul_le_mul_of_nonneg_right _ (sq_nonneg _)
        calc ‖∑ r, a r * e (ξ r * t)‖
            ≤ ∑ r, ‖a r * e (ξ r * t)‖ := norm_sum_le _ _
          _ = ∑ r, ‖a r‖ := by simp [norm_e]
          _ ≤ _ := by
              calc ∑ r, ‖a r‖
                  ≤ Real.sqrt (Finset.card Finset.univ) *
                    Real.sqrt (∑ r, ‖a r‖^2) := by
                      apply (Finset.inner_mul_le_norm_sq_mul_norm_sq _ _).trans
                      simp
                _ = _ := by simp [hnorm, Real.sqrt_one]
  -- Rearrange: ∫_I F = T - ∫_ℝ F·E
  have hrearrange :
      ∫ t in Set.Icc a₀ b₀, expSumSq a ξ t =
      T - ∫ t : ℝ, expSumSq a ξ t * kernelE B N a₀ b₀ t := by
    simp only [kernelE]
    -- ∫ F·(1/N ∫_I |ψ̂|² - 1_I) = 1/N ∫ F(∫_I |ψ̂|²) - ∫_I F
    rw [integral_sub]
    · rw [integral_const_mul, mul_comm N, hFubini]
      rw [integral_indicator measurableSet_Icc]
      simp; ring
    · exact (integrable_const _).mul_right _
    · exact (integrable_indicator measurableSet_Icc).mul_left _
  exact hrearrange

-- ============================================================
-- GOAL 5: E(t) ≪ (1 + dist(t,∂I)/N)^{-10}
--
--   Step 1: ψ̂ rapidly decaying (IBP)
--   Step 2: 1/N ∫_I |ψ̂((t-t₀)/N)|² dt₀ = ∫_{I_t} |ψ̂(u)|² du
--           where I_t = [(t-b₀)/N, (t-a₀)/N]
--   Step 3 (t ∈ I): 0 ∈ I_t, so = 1 - ∫_{ℝ\I_t} |ψ̂|² ≪ (1+d_t)^{-10}
--   Step 4 (t ∉ I): 0 ∉ I_t, so = ∫_{I_t} |ψ̂|² ≪ (1+d_t)^{-10}
-- ============================================================

theorem goal5 (B : BumpData) (N : ℝ) (hN : 0 < N)
    (a₀ b₀ : ℝ) (hab : a₀ ≤ b₀) :
    ∃ C : ℝ, 0 < C ∧ ∀ t : ℝ,
    |kernelE B N a₀ b₀ t| ≤
    C * (1 + min (|t - a₀| / N) (|t - b₀| / N)) ^ (-(10 : ℝ)) := by
  obtain ⟨C_d, hC_d, hdecay⟩ := psiHat_decay B 12
  refine ⟨C_d ^ 2 * 8, by positivity, ?_⟩
  intro t
have hint : Integrable (fun u => ‖psiHat B.ψ u‖ ^ 2) := by
    apply Integrable.mono (f := fun u => C_d * (1 + |u|)^(-(12:ℝ)))
    · apply Integrable.const_mul
      exact integrable_rpow_neg (by norm_num)
    · apply ae_of_all; intro u
      exact hdecay u
  have hint2 : Integrable (fun u => (1 + |u|)^(-(12:ℝ))) := by
    apply integrable_rpow_neg_abs (by norm_num : (1:ℝ) < 12)
  -- Step 2: Substitution u = (t-t₀)/N
  -- 1/N ∫_I |ψ̂((t-t₀)/N)|² dt₀ = ∫_{I_t} |ψ̂(u)|² du
  have hsubst : (1 / N) * ∫ t₀ in Set.Icc a₀ b₀,
      ‖psiHat B.ψ ((t - t₀) / N)‖ ^ 2 =
      ∫ u in Set.Icc ((t - b₀) / N) ((t - a₀) / N),
        ‖psiHat B.ψ u‖ ^ 2 := by
    rw [one_div, ← intervalIntegral.integral_comp_sub_left
      (fun u => ‖psiHat B.ψ u‖ ^ 2) t]
    simp [Set.uIcc_of_le (div_le_div_of_nonneg_right (by linarith) hN)]
    congr 1 <;> [field_simp; field_simp; ring]
  -- Tail bound: ∫_{|u|≥d} |ψ̂|² ≤ 2C²/(11(1+d)^{11}) ≪ (1+d)^{-10}
  have htail : ∀ d : ℝ, 0 ≤ d →
      ∫ u in {u : ℝ | d ≤ |u|}, ‖psiHat B.ψ u‖ ^ 2 ≤
      C_d ^ 2 * 8 * (1 + d) ^ (-(10 : ℝ)) := by
    intro d hd
    calc ∫ u in {u | d ≤ |u|}, ‖psiHat B.ψ u‖ ^ 2
        ≤ C_d ^ 2 * ∫ u in {u | d ≤ |u|}, (1 + |u|) ^ (-(12 : ℝ)) := by
          apply set_integral_mono_ae
          · exact hint
          · exact (Integrable.const_mul hint2_)
          · apply ae_of_all; intro u
            have := hdecay u
            nlinarith [norm_nonneg (psiHat B.ψ u), sq_nonneg (‖psiHat B.ψ u‖)]
      -- ∫_{|u|≥d} (1+|u|)^{-12} du = 2/(11(1+d)^{11})
      _ ≤ C_d ^ 2 * (8 * (1 + d) ^ (-(10 : ℝ))) := by
          apply mul_le_mul_of_nonneg_left _ (by positivity)
          have hcalc : ∫ u in Set.Ici d, (1 + u) ^ (-(12 : ℝ)) =
              (1/11) * (1 + d) ^ (-(11 : ℝ)) := by
            simp [MeasureTheory.integral_Ici_rpow_of_lt (by norm_num : -(12:ℝ) < -1)]
            ring
          calc ∫ u in {u : ℝ | d ≤ |u|}, (1 + |u|) ^ (-(12 : ℝ))
              ≤ 2 * ∫ u in Set.Ici d, (1 + u) ^ (-(12 : ℝ)) := by
                rw [show {u : ℝ | d ≤ |u|} = Set.Ici d ∪ Set.Iic (-d) from by
                  ext u; simp [abs_le, le_abs]]
                rw [integral_union (by simp; intro a h1 h2; linarith)
                  measurableSet_Ici]
                · have hsym : ∫ u in Set.Iic (-d), (1 + |u|) ^ (-(12 : ℝ)) =
                      ∫ u in Set.Ici d, (1 + u) ^ (-(12 : ℝ)) := by
                    rw [← integral_comp_neg (f := fun u => (1 + |u|)^(-(12:ℝ)))]
                    simp [abs_neg]
                    congr 1; ext u
                    rw [show Set.Ici d = {u | d ≤ u} from rfl]
                    simp [abs_of_nonneg]
                  rw [show (∫ u in Set.Ici d, (1 + |u|)^(-(12:ℝ))) =
                      ∫ u in Set.Ici d, (1+u)^(-(12:ℝ)) from by
                    congr 1; ext u
                    congr 2; exact abs_of_nonneg (le_trans hd (Set.mem_Ici.mp ‹_›))]
                  linarith [integral_nonneg (fun u => by positivity)]
                · sorry; · sorry
            _ = 2 * (1/11) * (1+d)^(-(11:ℝ)) := by rw [hcalc]; ring
            _ ≤ 8 * (1+d)^(-(10:ℝ)) := by
                have h1d : 1 + d ≥ 1 := by linarith
                nlinarith [Real.rpow_le_rpow_of_exponent_ge h1d
                  (by norm_num : -(11:ℝ) ≤ -(10:ℝ))]
      _ = C_d ^ 2 * 8 * (1 + d) ^ (-(10 : ℝ)) := by ring
  -- Set d_t = dist(t, ∂I)/N
  set d_t := min (|t - a₀| / N) (|t - b₀| / N)
  -- Case analysis: t ∈ I or t ∉ I
  by_cases htI : t ∈ Set.Icc a₀ b₀
  · -- STEP 3: t ∈ I → 0 ∈ I_t
    -- |E(t)| = |∫_{ℝ\I_t} |ψ̂|²| ≤ ∫_{|u|≥d_t} |ψ̂|²
    have hind : Set.indicator (Set.Icc a₀ b₀) (fun _ => (1:ℝ)) t = 1 :=
      Set.indicator_of_mem htI _
    simp only [kernelE, hind, hsubst]
    -- 0 ∈ I_t because t ∈ [a₀,b₀]
    have h0_in : (0:ℝ) ∈ Set.Icc ((t-b₀)/N) ((t-a₀)/N) := by
      constructor
      · apply div_nonpos_of_nonpos_of_nonneg <;> linarith [htI.2]
      · apply div_nonneg <;> linarith [htI.1]
    -- ∫_{I_t} = ∫_ℝ - ∫_{ℝ\I_t} = 1 - ∫_{ℝ\I_t}
    rw [psiHat_l2 B |>.symm]
    rw [← integral_add_compl measurableSet_Icc (by sorry)]
    simp only [add_sub_cancel_left]
    rw [abs_neg]
    calc |∫ u in (Set.Icc _ _)ᶜ, ‖psiHat B.ψ u‖^2|
        ≤ ∫ u in (Set.Icc _ _)ᶜ, ‖psiHat B.ψ u‖^2 :=
          le_abs_self _
      _ ≤ ∫ u in {u | d_t ≤ |u|}, ‖psiHat B.ψ u‖^2 := by
          apply set_integral_mono_set
          · exact fun u => sq_nonneg _
          -- ℝ\I_t ⊆ {|u| ≥ d_t} because 0 ∈ I_t
          · apply ae_of_all; intro u hu
            simp only [Set.mem_compl_iff, Set.mem_Icc, not_and_or, not_le] at hu
            simp only [Set.mem_setOf_eq]
            -- Geometric argument: since 0 ∈ I_t, u outside I_t means |u| ≥ d_t
            cases hu with
            | inl h =>
              rw [abs_of_nonpos (le_of_lt h)]
              simp [d_t]; constructor
              · linarith [h, h0_in.1]
              · linarith [h, (t - b₀)/N |>.neg_le_abs]
            | inr h =>
              rw [abs_of_pos h]
              simp [d_t]; constructor
              · linarith [h0_in.2]
              · linarith [h, le_abs_self ((t-b₀)/N)]
      _ ≤ C_d^2 * 8 * (1 + d_t)^(-(10:ℝ)) :=
          htail d_t (by simp [d_t]; positivity)
  · -- STEP 4: t ∉ I → 0 ∉ I_t
    have hind : Set.indicator (Set.Icc a₀ b₀) (fun _ => (1:ℝ)) t = 0 :=
      Set.indicator_of_not_mem htI _
    simp only [kernelE, hind, hsubst, sub_zero]
    rw [abs_of_nonneg (integral_nonneg (fun u => sq_nonneg _))]
    -- 0 ∉ I_t because t ∉ [a₀,b₀]
    have h0_out : (0:ℝ) ∉ Set.Icc ((t-b₀)/N) ((t-a₀)/N) := by
      simp only [Set.mem_Icc, not_and_or, not_le]
      cases Set.not_mem_Icc.mp htI with
      | inl h => right; apply div_pos_of_neg_of_neg <;> linarith
      | inr h => left; apply div_neg_of_pos_of_neg <;> linarith
    calc ∫ u in Set.Icc ((t-b₀)/N) ((t-a₀)/N), ‖psiHat B.ψ u‖^2
        ≤ ∫ u in {u | d_t ≤ |u|}, ‖psiHat B.ψ u‖^2 := by
          apply set_integral_mono_set
          · exact fun u => sq_nonneg _
          -- I_t ⊆ {|u| ≥ d_t} because 0 ∉ I_t
          · apply ae_of_all; intro u hu
            simp only [Set.mem_setOf_eq]
            simp only [Set.mem_Icc] at hu
            -- Geometric: 0 ∉ I_t → I_t entirely on one side → |u| ≥ d_t
            cases Set.not_mem_Icc.mp h0_out with
            | inl h =>
              -- I_t is entirely negative
              rw [abs_of_nonpos (le_trans hu.2 (le_of_lt h))]
              simp [d_t]
              constructor
              · linarith [hu.1]
              · linarith [hu.2, h]
            | inr h =>
              -- I_t is entirely positive
              rw [abs_of_nonneg (le_trans (le_of_lt (lt_of_not_le h)) hu.1)]
              simp [d_t]
              constructor
              · linarith [hu.1, h]
              · linarith [hu.2]
      _ ≤ C_d^2 * 8 * (1 + d_t)^(-(10:ℝ)) :=
          htail d_t (by simp [d_t]; positivity)

-- ============================================================
-- GOAL 6: ∫_ℝ F·E ≪ N  (dyadic decomposition)
--
-- From handwritten proof:
--   Case 1 (N ≪ T): divide I into J_k with |J_k| = N
--     Layer ℓ: J with dist(J,∂I) ~ 2^ℓN, contains ~2^ℓ intervals
--     |E(t)| ≪ (2^ℓ)^{-10}  (by Goal 5)
--     ∫_J F ≤ CN  (by Goal 3 = eq 3.2)
--     Sum: ∑_ℓ 2^ℓ · CN · 2^{-10ℓ} = CN ∑ 2^{-9ℓ} ≪ N
--   Cases B,C (outside I): same method
--   Case 2 (T ≪ N): direct application of (3.2)
-- ============================================================

theorem goal6 (B : BumpData) {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (N : ℝ) (hN : 0 < N)
    (hnorm : ∑ r, ‖a r‖ ^ 2 = 1)
    (hsep : Separated N ξ)
    (a₀ b₀ : ℝ) (hab : a₀ ≤ b₀) :
    ∃ C : ℝ, 0 < C ∧
    |∫ t : ℝ, expSumSq a ξ t * kernelE B N a₀ b₀ t| ≤ C * N := by
  set T := b₀ - a₀
  obtain ⟨C₃, hC₃, hG3⟩ := goal3 B a ξ N hN hnorm hsep a₀
  obtain ⟨C₅, hC₅, hG5⟩ := goal5 B N hN a₀ b₀ hab
  -- Case 2: T ≤ N (I fits inside one interval J)
  by_cases hTN : T ≤ N
  · -- Direct application of Goal 3
    refine ⟨C₃ * C₅, by positivity, ?_⟩
    calc |∫ t : ℝ, expSumSq a ξ t * kernelE B N a₀ b₀ t|
        ≤ ∫ t : ℝ, expSumSq a ξ t * |kernelE B N a₀ b₀ t| := by
          apply (abs_integral_le_integral_abs _).trans
          apply integral_mono_ae
          · sorry; · sorry
          · apply ae_of_all; intro t
            exact mul_le_mul_of_nonneg_left (le_abs_self _)
              (by simp [expSumSq]; positivity)
      _ ≤ ∫ t in Set.Icc a₀ (a₀ + N), expSumSq a ξ t * C₅ := by
          -- |E| ≤ C₅ · 1 on [a₀, a₀+N] ⊇ I
          apply integral_mono_of_subset
          · intro t ht
            exact ⟨ht.1, le_trans ht.2 (by linarith)⟩
          · apply ae_of_all; intro t
            apply mul_le_mul_of_nonneg_left _ (by simp [expSumSq]; positivity)
            exact le_trans (hG5 t) (by
              apply mul_le_mul_of_nonneg_left _ hC₅.le
              simp [Real.rpow_neg_nonpos])
          · sorry
      _ = C₅ * ∫ t in Set.Icc a₀ (a₀ + N), expSumSq a ξ t := by
          rw [← integral_const_mul]; congr 1; ext t; ring
      _ ≤ C₅ * (C₃ * N) := mul_le_mul_of_nonneg_left (hG3 a₀) hC₅.le
      _ = C₃ * C₅ * N := by ring
  -- Case 1: T > N (dyadic decomposition)
  · push_neg at hTN
    -- L = ⌈log₂(T/N)⌉ layers
    set L := Nat.ceil (Real.log (T / N) / Real.log 2)
    refine ⟨C₃ * C₅ * 4 / (1 - (2:ℝ)^(-(9:ℝ))), by
      apply div_pos; positivity
      linarith [Real.rpow_lt_one (by norm_num) (by norm_num : (0:ℝ) < 2) (by norm_num)],
      ?_⟩
    -- The dyadic decomposition argument
    -- Layer ℓ: intervals J with dist(c(J),∂I) ∈ [2^ℓN, 2^{ℓ+1}N)
    -- # intervals in layer ℓ: ≤ 2^{ℓ+1}
    -- |E(t)| on layer ℓ: ≤ C₅ · (2^ℓ)^{-10}
    -- ∫_J F ≤ C₃ · N  (by Goal 3)
    -- Sum over layers: ∑_{ℓ=0}^{L} 2^{ℓ+1} · C₃N · C₅·(2^ℓ)^{-10}
    --                = 2C₃C₅N ∑_{ℓ=0}^{L} 2^{-9ℓ}
    --                ≤ 2C₃C₅N · 1/(1-2^{-9})
    calc |∫ t : ℝ, expSumSq a ξ t * kernelE B N a₀ b₀ t|
        ≤ ∑ ℓ in Finset.range (L + 1),
            ∑ k in Finset.range (2^(ℓ+1)),
            |∫ t in Set.Icc (a₀ - (k+1) * N) (a₀ - k * N),
              expSumSq a ξ t * kernelE B N a₀ b₀ t| +
          ∑ ℓ in Finset.range (L + 1),
            ∑ k in Finset.range (2^(ℓ+1)),
            |∫ t in Set.Icc (b₀ + k * N) (b₀ + (k+1) * N),
              expSumSq a ξ t * kernelE B N a₀ b₀ t| := by
          sorry -- partition ℝ into layers
      _ ≤ ∑ ℓ in Finset.range (L + 1), (2^(ℓ+1) : ℝ) * (C₃ * N) *
            (C₅ * (2^ℓ)^(-(10:ℝ))) +
          ∑ ℓ in Finset.range (L + 1), (2^(ℓ+1) : ℝ) * (C₃ * N) *
            (C₅ * (2^ℓ)^(-(10:ℝ))) := by
          apply add_le_add
          all_goals {
            apply Finset.sum_le_sum; intro ℓ _
            apply Finset.sum_le_sum_of_subset; simp
          }
      _ = 2 * (C₃ * C₅ * N) * ∑ ℓ in Finset.range (L + 1), (2:ℝ)^(-(9:ℝ)*ℓ)  := by
  have hstep : ∀ ℓ : ℕ, (2^(ℓ+1) : ℝ) * (C₃ * N) * (C₅ * (2^ℓ)^(-(10:ℝ))) = 
      2 * (C₃ * C₅ * N) * (2^(-(9:ℝ)))^ℓ := by
    intro ℓ
    have h1 : (2 : ℝ)^(ℓ+1) = 2 * 2^ℓ := by ring
    have h2 : ((2:ℝ)^ℓ)^(-(10:ℝ)) = ((2:ℝ)^(-(9:ℝ)))^ℓ / (2:ℝ)^ℓ := by
      rw [← Real.rpow_natCast 2 ℓ]
      rw [← Real.rpow_mul (by norm_num)]
      rw [← Real.rpow_natCast 2 ℓ]
      simp [Real.rpow_neg, mul_comm]
    rw [h1, h2]
    ring
  simp_rw [hstep]
  rw [← Finset.mul_sum]
      -- Geometric series: ∑_{ℓ=0}^{L} 2^{-9ℓ} ≤ 1/(1-2^{-9})
      _ ≤ 2 * (C₃ * C₅ * N) * (1 / (1 - (2:ℝ)^(-(9:ℝ)))) := by
          apply mul_le_mul_of_nonneg_left _ (by positivity)
          calc ∑ ℓ in Finset.range (L + 1), (2:ℝ)^(-(9:ℝ)*ℓ)
              ≤ ∑' ℓ : ℕ, (2:ℝ)^(-(9:ℝ)*ℓ) := by
                apply sum_le_tsum (Finset.subset_univ _)
                · intro ℓ _; positivity
                · apply Summable.of_nonneg_of_le
                  · intro ℓ; positivity
                  · intro ℓ; exact le_refl _
                  · apply summable_geometric_of_abs_lt_one <;> norm_num
            _ = 1 / (1 - (2:ℝ)^(-(9:ℝ))) := by
                rw [tsum_geometric_of_abs_lt_one]
                · norm_num
                · norm_num
      _ = C₃ * C₅ * 4 / (1 - (2:ℝ)^(-(9:ℝ))) * N := by ring

-- ============================================================
-- GOAL 7: Lemma 3.1 (Assembly)
--
-- From handwritten proof:
--   ∫_I F = T - ∫_ℝ F·E      (by Goal 4)
--              ↑               ↑
--         = T + O(N)·1   (by Goal 6: |∫ F·E| ≪ N)
--              = (T + O(N)) · ∑|aᵣ|²  (by Goal 1: WLOG)
-- ============================================================

theorem lemma3_1 (B : BumpData) {R : ℕ} (a : Fin R → ℂ) (ξ : Fin R → ℝ)
    (N : ℝ) (hN : 0 < N) (hsep : Separated N ξ)
    (a₀ b₀ : ℝ) (T : ℝ) (hT : T = b₀ - a₀) (hab : a₀ ≤ b₀) :
    ∃ C : ℝ, 0 < C ∧ ∃ θ : ℝ, |θ| ≤ C ∧
    ∫ t in Set.Icc a₀ b₀, ‖∑ r, a r * e (ξ r * t)‖ ^ 2 =
    (T + θ * N) * ∑ r, ‖a r‖ ^ 2 := by
  set M := ∑ r, ‖a r‖ ^ 2
  -- Trivial case: M = 0 → all aᵣ = 0
  by_cases hM0 : M = 0
  · refine ⟨1, one_pos, 0, by simp, ?_⟩
    have hzero : ∀ r, a r = 0 := by
      intro r
      have h1 : 0 ≤ ‖a r‖ ^ 2 := sq_nonneg _
      have h2 : ‖a r‖ ^ 2 ≤ M :=
        Finset.single_le_sum (fun i _ => sq_nonneg _) _ (Finset.mem_univ r)
      have h3 : ‖a r‖ = 0 := by
        nlinarith [hM0 ▸ h2]
      exact norm_eq_zero.mp h3
    simp [hM0, hzero]
  · -- M > 0: normalize Aᵣ = aᵣ/√M
    have hMpos : 0 < M :=
      lt_of_le_of_ne (Finset.sum_nonneg fun r _ => sq_nonneg _) (Ne.symm hM0)
    -- GOAL 1: define A with ∑|Aᵣ|² = 1
    set A : Fin R → ℂ := fun r => a r / (Real.sqrt M : ℂ)
    have hAnorm : ∑ r, ‖A r‖ ^ 2 = 1 := goal1 a M rfl hMpos
    -- Rescaling: ‖expSum a‖² = M · ‖expSum A‖²
    have hrescale : ∀ t,
        ‖∑ r, a r * e (ξ r * t)‖ ^ 2 =
        M * ‖∑ r, A r * e (ξ r * t)‖ ^ 2 := by
      intro t
      have : ∑ r, a r * e (ξ r * t) =
          (Real.sqrt M : ℂ) * ∑ r, A r * e (ξ r * t) := by
        simp [A, Finset.mul_sum, div_mul_cancel₀]
        intro r
        field_simp
        rw [div_mul_cancel₀]
        exact Complex.ofReal_ne_zero.mpr (Real.sqrt_ne_zero'.mpr hMpos)
      rw [this, norm_mul, Complex.norm_ofReal,
          abs_of_nonneg (Real.sqrt_nonneg M), mul_pow, Real.sq_sqrt hMpos.le]
    -- Apply Goal 4 to get the Fubini identity
    have h4 := goal4 B A ξ N hN hAnorm hsep a₀ b₀ T hT hab
    -- Apply Goal 6 to bound the error
    obtain ⟨C, hC, h6⟩ := goal6 B A ξ N hN hAnorm hsep a₀ b₀ hab
    -- ∫_I ‖expSum A‖² = T - err,  |err| ≤ C·N
    set err := ∫ t : ℝ, expSumSq A ξ t * kernelE B N a₀ b₀ t
    have hA_eq : ∫ t in Set.Icc a₀ b₀, expSumSq A ξ t = T - err := by
      simp [expSumSq] at h4 ⊢; exact h4
    have herr_bd : |err| ≤ C * N := by
      simp [expSumSq] at h6; exact h6
    -- ∫_I ‖expSum a‖² = M · (T - err) = (T + θN) · M with θ = -err/N
    refine ⟨C, hC, -err / N, ?_, ?_⟩
    · -- |θ| = |err|/N ≤ C·N/N = C
      rw [abs_div, abs_neg, abs_of_pos hN]
      exact div_le_of_le_mul₀ hN.le (by positivity) (by linarith)
    · -- ∫_I ‖expSum a‖² = (T + θN) · M
      have ha_int : ∫ t in Set.Icc a₀ b₀, ‖∑ r, a r * e (ξ r * t)‖ ^ 2 =
          M * (T - err) := by
        conv_lhs => ext t; rw [hrescale t]
        rw [integral_const_mul]
        simp [expSumSq] at hA_eq
        rw [hA_eq]
      rw [ha_int]
      field_simp; ring

end
